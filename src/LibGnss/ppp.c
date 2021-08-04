/*------------------------------------------------------------------------------
* ppp.c : precise point positioning
*
*          Copyright (C) 2010-2018 by T.TAKASU, All rights reserved.
*
* options : -DIERS_MODEL  use IERS tide model
*           -DOUTSTAT_AMB output ambiguity parameters to solution status
*
* references :
*    [1] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*    [2] D.D.McCarthy and G.Petit, IERS Technical Note 32, IERS Conventions
*        2003, November 2003
*    [3] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*        Space Technology Library, 2004
*    [4] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*        May 2009
*    [5] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*        Code Biases, URA
*    [6] MacMillan et al., Atmospheric gradients and the VLBI terrestrial and
*        celestial reference frames, Geophys. Res. Let., 1997
*    [7] G.Petit and B.Luzum (eds), IERS Technical Note No. 36, IERS
*         Conventions (2010), 2010
*    [8] J.Kouba, A simplified yaw-attitude model for eclipsing GPS satellites,
*        GPS Solutions, 13:1-12, 2009
*    [9] F.Dilssner, GPS IIF-1 satellite antenna phase center and attitude
*        modeling, InsideGNSS, September, 2010
*    [10] F.Dilssner, The GLONASS-M satellite yaw-attitude model, Advances in
*        Space Research, 2010
*    [11] IGS MGEX (http://igs.org/mgex)
*
* version : $Revision:$ $Date:$
* history : 2010/07/20 1.0  new
*                           added api:
*                               tidedisp()
*           2010/12/11 1.1  enable exclusion of eclipsing satellite
*           2012/02/01 1.2  add gps-glonass h/w bias correction
*                           move windupcorr() to rtkcmn.c
*           2013/03/11 1.3  add otl and pole tides corrections
*                           involve iers model with -DIERS_MODEL
*                           change initial variances
*                           suppress acos domain error
*           2013/09/01 1.4  pole tide model by iers 2010
*                           add mode of ionosphere model off
*           2014/05/23 1.5  add output of trop gradient in solution status
*           2014/10/13 1.6  fix bug on P0(a[3]) computation in tide_oload()
*                           fix bug on m2 computation in tide_pole()
*           2015/03/19 1.7  fix bug on ionosphere correction for GLO and BDS
*           2015/05/10 1.8  add function to detect slip by MW-LC jump
*                           fix ppp solutin problem with large clock variance
*           2015/06/08 1.9  add precise satellite yaw-models
*                           cope with day-boundary problem of satellite clock
*           2015/07/31 1.10 fix bug on nan-solution without glonass nav-data
*                           pppoutsolsat() -> pppoutstat()
*           2015/11/13 1.11 add L5-receiver-dcb estimation
*                           merge post-residual validation by rnx2rtkp_test
*                           support support option opt->pppopt=-GAP_RESION=nnnn
*           2016/01/22 1.12 delete support for yaw-model bug
*                           add support for ura of ephemeris
*           2018/10/10 1.13 support api change of satexclude()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define SQR(x)      ((x)*(x))
//#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
//#define MAX(x,y)    ((x)>(y)?(x):(y))
//#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define MAX_ITER    8               /* max number of iterations */
#define MAX_STD_FIX 0.15            /* max std-dev (3d) to fix solution */
#define MIN_NSAT_SOL 4              /* min satellite number for solution */
#define THRES_REJECT 4.0             /* reject threshold of posfit-res (sigma) */

#define THRES_MW_JUMP 10.0

#define VAR_POS     SQR(60.0)       /* init variance receiver position (m^2) */
#define VAR_VEL     SQR(60.0)       /* init variance of receiver vel ((m/s)^2) */
#define VAR_ACC     SQR(60.0)       /* init variance of receiver acc ((m/ss)^2) */
#define VAR_CLK     SQR(60.0)       /* init variance receiver clock (m^2) */
#define VAR_DCB     SQR(60.0)       /* init variance dcb (m^2) */
#define VAR_IFCB    SQR(60.0)
#define VAR_ZWD     SQR(0.03)       /* init variance ztd (m^2) */
#define VAR_GRA     SQR(0.01)       /* init variance gradient (m^2) */
#define VAR_BIAS    SQR(60.0)       /* init variance phase-bias (m^2) */
#define VAR_IONO    SQR(60.0)       /* init variance iono-delay */
#define VAR_GLO_IFB SQR( 0.6)       /* variance of glonass ifb */

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.7             /* relative humidity for saastamoinen model */
#define GAP_RESION  120             /* default gap to reset ionos parameters (ep) */

#define EFACT_GPS_L5 10.0           /* error factor of GPS/QZS L5 */

#define MUDOT_GPS   (0.00836*D2R)   /* average angular velocity GPS (rad/s) */
#define MUDOT_GLO   (0.00888*D2R)   /* average angular velocity GLO (rad/s) */
#define EPS0_GPS    (13.5*D2R)      /* max shadow crossing angle GPS (rad) */
#define EPS0_GLO    (14.2*D2R)      /* max shadow crossing angle GLO (rad) */
#define T_POSTSHADOW 1800.0         /* post-shadow recovery time (s) */
#define QZS_EC_BETA 20.0            /* max beta angle for qzss Ec (deg) */

/* number and index of states */
#define NF(opt)      ((opt)->ionoopt==IONOOPT_IFLC?1:((opt)->ionoopt==IONOOPT_IF2?2:(opt)->nf))
#define NP(opt)      ((opt)->dynamics?9:3)
#define NC(opt)      ((opt)->sdopt?0:((opt)->bd3opt>=BD3OPT_BD2_3?NSYS+1:NSYS))
#define NRDCB(opt)   (((opt)->ionoopt==IONOOPT_UC_CONS&&(opt)->nf>=2)?((opt)->bd3opt>=BD3OPT_BD2_3?NSYS+1:NSYS):0)
#define NIFCB(opt)   ((opt)->nf>=3?((opt)->bd3opt>=BD3OPT_BD2_3?NSYS+1:NSYS):0)
#define NT(opt)      ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))
#define NI(opt)      ((opt)->ionoopt==IONOOPT_UC||(opt)->ionoopt==IONOOPT_UC_CONS?MAXSAT:0)
#define NR(opt)      (NP(opt)+NC(opt)+NRDCB(opt)+NIFCB(opt)+NT(opt)+NI(opt))
#define NB(opt)      (NF(opt)*MAXSAT)
#define NX(opt)      (NR(opt)+NB(opt))
#define IC(s,opt)    (NP(opt)+(s))
#define IRDCB(s,opt) (NP(opt)+NC(opt)+(s))
#define IIFCB(s,opt) (NP(opt)+NC(opt)+NRDCB(opt)+(s))
#define IT(opt)      (NP(opt)+NC(opt)+NRDCB(opt)+NIFCB(opt))
#define II(s,opt)    (NP(opt)+NC(opt)+NRDCB(opt)+NIFCB(opt)+NT(opt)+(s)-1)
#define IB(s,f,opt)  (NR(opt)+MAXSAT*(f)+(s)-1)

extern int iamb_ppp(const prcopt_t *opt,int sat,int f)
{
    return IB(sat,f,opt);
}

/* standard deviation of state -----------------------------------------------*/
static double STD(rtk_t *rtk, int i)
{
    if (rtk->sol.stat==SOLQ_FIX) return SQRT(rtk->Pa[i+i*rtk->na]);
    return SQRT(rtk->P[i+i*rtk->nx]);
}
/* write solution status for PPP ---------------------------------------------*/
extern int pppoutstat(rtk_t *rtk, char *buff)
{
    ssat_t *ssat;
    double tow,pos[3],vel[3],acc[3],*x;
    int i,j,week;
    char id[32],*p=buff;
    
    if (!rtk->sol.stat) return 0;
    
    trace(3,"pppoutstat:\n");
    
    tow=time2gpst(rtk->sol.time,&week);
    
    x=rtk->sol.stat==SOLQ_FIX?rtk->xa:rtk->x;
    
    /* receiver position */
    p+=sprintf(p,"$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
               rtk->sol.stat,x[0],x[1],x[2],STD(rtk,0),STD(rtk,1),STD(rtk,2));
    
    /* receiver velocity and acceleration */
    if (rtk->opt.dynamics) {
        ecef2pos(rtk->sol.rr,pos);
        ecef2enu(pos,rtk->x+3,vel);
        ecef2enu(pos,rtk->x+6,acc);
        p+=sprintf(p,"$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,"
                   "%.4f,%.5f,%.5f,%.5f\n",week,tow,rtk->sol.stat,vel[0],vel[1],
                   vel[2],acc[0],acc[1],acc[2],0.0,0.0,0.0,0.0,0.0,0.0);
    }
    /* receiver clocks */
    i=IC(0,&rtk->opt);
    p+=sprintf(p,"$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",
               week,tow,rtk->sol.stat,1,x[i]*1E9/CLIGHT,x[i+1]*1E9/CLIGHT,
               STD(rtk,i)*1E9/CLIGHT,STD(rtk,i+1)*1E9/CLIGHT);
    
    /* tropospheric parameters */
    if (rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG) {
        i=IT(&rtk->opt);
        p+=sprintf(p,"$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n",week,tow,rtk->sol.stat,
                   1,x[i],STD(rtk,i));
    }
    if (rtk->opt.tropopt==TROPOPT_ESTG) {
        i=IT(&rtk->opt);
        p+=sprintf(p,"$TRPG,%d,%.3f,%d,%d,%.5f,%.5f,%.5f,%.5f\n",week,tow,
                   rtk->sol.stat,1,x[i+1],x[i+2],STD(rtk,i+1),STD(rtk,i+2));
    }
    /* ionosphere parameters */
    if (rtk->opt.ionoopt==IONOOPT_UC||rtk->opt.ionoopt==IONOOPT_UC_CONS) {
        for (i=0;i<MAXSAT;i++) {
            ssat=rtk->ssat+i;
            if (!ssat->vs) continue;
            j=II(i+1,&rtk->opt);
            if (rtk->x[j]==0.0) continue;
            satno2id(i+1,id);
            p+=sprintf(p,"$ION,%d,%.3f,%d,%s,%.1f,%.1f,%.4f,%.4f\n",week,tow,
                       rtk->sol.stat,id,rtk->ssat[i].azel[0]*R2D,
                       rtk->ssat[i].azel[1]*R2D,x[j],STD(rtk,j));
        }
    }

    return (int)(p-buff);
}

/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/
static void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs,rtk_t *rtk)
{
    double rsun[3],esun[3],r,ang,erpv[5]={0},cosa;
    int i,j,week=0;
    const char *type;
    double sec,ex[3]={0},ep[6]={0};
    
    trace(4,"testeclipse:\n");
    
    /* unit vector of sun direction (ecef) */
    sunmoonpos(gpst2utc(obs[0].time),erpv,rsun,NULL,NULL);
    normv3(rsun,esun);
    
    for (i=0;i<n;i++) {
        sec=time2gpst(obs[0].time,&week);
        time2epoch(obs[0].time,ep);
        week=cal_eclips(obs[i].sat,rs+6*i,rs+i*6+3,rsun,sec,ex,nav,(int)ep[3]);

        if(week != 0){
            trace(3,"%s %s eclipse, type=%d\n",time_str(obs[0].time,3),sat_id(obs[i].sat),week);
            rtk->ssat[obs[i].sat-1].eclipse=4.0;
        }
        else{
            rtk->ssat[obs[i].sat-1].eclipse=1.0;
        }

        type=nav->pcvs[obs[i].sat-1].type;
        
        if ((r=norm(rs+i*6,3))<=0.0) continue;
        
        /* only block IIA */
        if (*type&&!strstr(type,"BLOCK IIA")) continue;
        
        /* sun-earth-satellite angle */
        cosa=dot(rs+i*6,esun,3)/r;
        cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
        ang=acos(cosa);
        
        /* test eclipse */
        if (ang<PI/2.0||r*sin(ang)>RE_WGS84) continue;
        
        trace(3,"eclipsing sat excluded %s sat=%2d\n",time_str(obs[0].time,0),
              obs[i].sat);
        
        for (j=0;j<3;j++) rs[j+i*6]=0.0;
    }
}
/* nominal yaw-angle ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
    if (fabs(beta)<1E-12&&fabs(mu)<1E-12) return PI;
    return atan2(-tan(beta),sin(mu))+PI;
}
/* yaw-angle of satellite ----------------------------------------------------*/
extern int yaw_angle(int sat, const char *type, int opt, double beta, double mu,
                     double *yaw)
{
    *yaw=yaw_nominal(beta,mu);
    return 1;
}
/* satellite attitude model --------------------------------------------------*/
static int sat_yaw(gtime_t time, int sat, const char *type, int opt,
                   const double *rs, double *exs, double *eys)
{
    double rsun[3],ri[6],es[3],esun[3],n[3],p[3],en[3],ep[3],ex[3],E,beta,mu;
    double yaw,cosy,siny,erpv[5]={0};
    int i;
    
    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,NULL);
    
    /* beta and orbit angle */
    matcpy(ri,rs,6,1);
    ri[3]-=OMGE*ri[1];
    ri[4]+=OMGE*ri[0];
    cross3(ri,ri+3,n);
    cross3(rsun,n,p);
    if (!normv3(rs,es)||!normv3(rsun,esun)||!normv3(n,en)||
        !normv3(p,ep)) return 0;
    beta=PI/2.0-acos(dot(esun,en,3));
    E=acos(dot(es,ep,3));
    mu=PI/2.0+(dot(es,esun,3)<=0?-E:E);
    if      (mu<-PI/2.0) mu+=2.0*PI;
    else if (mu>=PI/2.0) mu-=2.0*PI;
    
    /* yaw-angle of satellite */
    if (!yaw_angle(sat,type,opt,beta,mu,&yaw)) return 0;
    
    /* satellite fixed x,y-vector */
    cross3(en,es,ex);
    cosy=cos(yaw);
    siny=sin(yaw);
    for (i=0;i<3;i++) {
        exs[i]=-siny*en[i]+cosy*ex[i];
        eys[i]=-cosy*en[i]-siny*ex[i];
    }
    return 1;
}
/* phase windup model --------------------------------------------------------*/
static int model_phw(gtime_t time, int sat, const char *type, int opt,
                     const double *rs, const double *rr, double *phw)
{
    double exs[3],eys[3],ek[3],exr[3],eyr[3],eks[3],ekr[3],E[9];
    double dr[3],ds[3],drs[3],r[3],pos[3],cosp,ph;
    int i;
    
    if (opt<=0) return 1; /* no phase windup */
    
    /* satellite yaw attitude model */
    if (!sat_yaw(time,sat,type,opt,rs,exs,eys)) return 0;
    
    /* unit vector satellite to receiver */
    for (i=0;i<3;i++) r[i]=rr[i]-rs[i];
    if (!normv3(r,ek)) return 0;
    
    /* unit vectors of receiver antenna */
    ecef2pos(rr,pos);
    xyz2enu(pos,E);
    exr[0]= E[1]; exr[1]= E[4]; exr[2]= E[7]; /* x = north */
    eyr[0]=-E[0]; eyr[1]=-E[3]; eyr[2]=-E[6]; /* y = west  */
    
    /* phase windup effect */
    cross3(ek,eys,eks);
    cross3(ek,eyr,ekr);
    for (i=0;i<3;i++) {
        ds[i]=exs[i]-ek[i]*dot(ek,exs,3)-eks[i];
        dr[i]=exr[i]-ek[i]*dot(ek,exr,3)+ekr[i];
    }
    cosp=dot(ds,dr,3)/norm(ds,3)/norm(dr,3);
    if      (cosp<-1.0) cosp=-1.0;
    else if (cosp> 1.0) cosp= 1.0;
    ph=acos(cosp)/2.0/PI;
    cross3(ds,dr,drs);
    if (dot(ek,drs,3)<0.0) ph=-ph;
    
    *phw=ph+floor(*phw-ph+0.5); /* in cycle */
    return 1;
}
/* measurement error variance ------------------------------------------------*/
static double varerr(int sat, int sys, double el, double snr_rover, int freq, 
                     int type, const prcopt_t *opt)
{
    double a, b, snr_max;
    double fact = 1.0;
    double sinel = sin(el);
    int prn;
    satsys(sat,&prn);
    
    if (type == 1) fact *= opt->eratio[(freq == 0) ? 0 : 1];
    
    switch (sys) {
        case SYS_GPS: fact *= EFACT_GPS; break;
        case SYS_GLO: fact *= EFACT_GLO; break;
        case SYS_GAL: fact *= EFACT_GAL; break;
        case SYS_SBS: fact *= EFACT_SBS; break;
        default:      fact *= EFACT_GPS; break;
    }

#if 0
    if(sys==SYS_CMP){
        if(prn<=5) fact*=100.0;
//        else if((prn>=6&&prn<=10)||(prn==13||prn==16)) fact*=4.0;
    }
#endif

    a = fact * opt->err[1];
    b = fact * opt->err[2];
    snr_max = opt->err[5];

    if(opt->ionoopt==IONOOPT_IFLC||opt->ionoopt==IONOOPT_IF2){
        if(opt->nf==1) fact*=SQR(0.5);
        else if(opt->nf==2) fact*=SQR(3.0);
        else if(opt->nf==3){
            if(opt->ionoopt==IONOOPT_IFLC){

            }
            else if(opt->ionoopt==IONOOPT_IF2){

            }
        }
    }

    switch (opt->weightmode) {
        case WEIGHTOPT_ELEVATION: return fact * ( SQR(a) + SQR(b / sinel) );
        case WEIGHTOPT_SNR      : return fact * SQR(a) * pow(10, 0.1 * MAX(snr_max - snr_rover, 0)); 
                                                   ;
        default: return 0;
    }
}
/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t *rtk, double xi, double var, int i)
{
    int j;
    rtk->x[i]=xi;
    for (j=0;j<rtk->nx;j++) {
        rtk->P[i+j*rtk->nx]=rtk->P[j+i*rtk->nx]=i==j?var:0.0;
    }
}

/* geometry-free phase measurement (m)-------------------------------------------*/
static double gfmeas(const obsd_t *obs, const nav_t *nav)
{
    double freq1,freq2;

    freq1=sat2freq(obs->sat,obs->code[0],nav);
    freq2=sat2freq(obs->sat,obs->code[1],nav);
    if (freq1==0.0||freq2==0.0||obs->L[0]==0.0||obs->L[1]==0.0) return 0.0;
    return (obs->L[0]/freq1-obs->L[1]/freq2)*CLIGHT;
}
/* Melbourne-Wubbena linear combination (cycle)--------------------------------------*/
static double mwmeas(const prcopt_t *opt,const obsd_t *obs, const nav_t *nav,double *var,double el)
{
    int f1=0,f2=1,sys,prn,sys_idx=-1;
    double freq1,freq2,lam_wl,lam1,lam2,P1C1=0.0,P2C2=0.0,cbias[NFREQ]={0};
    double osb_L1=0.0,osb_L2=0.0,osb_P1=0.0,osb_P2=0.0;
    double mea_L1=0.0,mea_L2=0.0,mea_P1=0.0,mea_P2=0.0;

    sys=satsys(obs->sat,&prn);
    sys_idx=satsysidx(obs->sat);
    if(sys_idx==-1) return  0.0;
    f1=opt->gnss_frq_idx[sys_idx][0]-1;
    f2=opt->gnss_frq_idx[sys_idx][1]-1;

    freq1=sat2freq(obs->sat,obs->code[f1],nav);
    freq2=sat2freq(obs->sat,obs->code[f2],nav);
    if (freq1==0.0||freq2==0.0||obs->L[f1]==0.0||obs->L[f2]==0.0||
        obs->P[f1]==0.0||obs->P[f2]==0.0) return 0.0;

    mea_L1=obs->L[f1];mea_L2=obs->L[f2];
    mea_P1=obs->P[f1];mea_P2=obs->P[f2];

    if(opt->modear==ARMODE_PPPAR_ILS&&(opt->arprod>=AR_PROD_OSB_GRM&&opt->arprod<=AR_PROD_OSB_CNT)){
        if(opt->arprod==AR_PROD_OSB_WHU){
            mea_P1+=corrISC(opt,nav->cbias[obs->sat-1],obs->code[f1],obs->sat);
            mea_P2+=corrISC(opt,nav->cbias[obs->sat-1],obs->code[f2],obs->sat);
            if(sys==SYS_GPS){
                if(obs->code[f1]==CODE_L1C){
                    osb_P1=nav->osbs->sat_osb[0].code[obs->sat-1][CODE_L1C];
                    osb_L1=nav->osbs->sat_osb[0].phase[obs->sat-1][CODE_L1C];
                }
                if(obs->code[f2]==CODE_L2W){
                    osb_P2=nav->osbs->sat_osb[0].code[obs->sat-1][CODE_L2W];
                    osb_L2=nav->osbs->sat_osb[0].phase[obs->sat-1][CODE_L2W];
                }
            }
            else if(sys==SYS_GLO){
                osb_L1=osb_L2=0.0;
            }
            else if(sys==SYS_CMP){
                if(obs->code[f1]==CODE_L2I) osb_L1=nav->osbs->sat_osb[0].phase[obs->sat-1][CODE_L2I];
                if(obs->code[f2]==CODE_L6I) osb_L2=nav->osbs->sat_osb[0].phase[obs->sat-1][CODE_L6I];
            }
            else if(sys==SYS_GAL){
                if(obs->code[f1]==CODE_L1C) osb_L1=nav->osbs->sat_osb[0].phase[obs->sat-1][CODE_L1C];
                if(obs->code[f2]==CODE_L5Q) osb_L2=nav->osbs->sat_osb[0].phase[obs->sat-1][CODE_L5Q];
            }
        }
        else if(opt->arprod==AR_PROD_OSB_GRM||opt->arprod==AR_PROD_OSB_CNT||opt->arprod==AR_PROD_OSB_COM||opt->arprod==AR_PROD_OSB_SGG){
            matchcposb(opt->arprod,obs,nav,f1,&osb_P1,&osb_L1);
            matchcposb(opt->arprod,obs,nav,f2,&osb_P2,&osb_L2);
        }
    }
    else{
        mea_P1+=corrISC(opt,nav->cbias[obs->sat-1],obs->code[f1],obs->sat);
        mea_P2+=corrISC(opt,nav->cbias[obs->sat-1],obs->code[f2],obs->sat);
    }

    lam_wl=CLIGHT/(freq1-freq2);
    lam1=CLIGHT/freq1;
    lam2=CLIGHT/freq2;
#if 0
    return ((obs->L[0]-obs->L[1])*CLIGHT/(freq1-freq2)-
           (freq1*(obs->P[0]+P1C1)+freq2*(obs->P[1]+P2C2))/(freq1+freq2))/lam_wl;
#else
    double mw=0.0;
    mw=(freq1*(mea_L1*lam1-osb_L1)-freq2*(mea_L2*lam2-osb_L2))/(freq1-freq2)
       -(freq1*(mea_P1-osb_P1)+freq2*(mea_P2-osb_P2))/(freq1+freq2);
    mw/=lam_wl;

    if(var){
        *var=SQR(0.7*opt->eratio[0])*(SQR(opt->err[1])+SQR(opt->err[2]/sin(el)))/SQR(lam_wl);
    }
    return mw;
#endif
}
/* antenna corrected measurements --------------------------------------------*/
static void corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel,
                      const prcopt_t *opt, const double *dantr,
                      const double *dants, double phw, double *L, double *P,
                      double *Lc, double *Pc,double *dcb)
{
    double freq[NFREQ]={0},C1,C2;
    int i,sys=satsys(obs->sat,NULL);

    for (i=0;i<NFREQ;i++) {
        L[i]=P[i]=0.0;
        freq[i]=sat2freq(obs->sat,obs->code[i],nav);
        if (freq[i]==0.0||obs->L[i]==0.0||obs->P[i]==0.0) continue;
        if (testsnr(0,0,azel[1],obs->SNR[i]*SNR_UNIT,&opt->snrmask)) continue;

        /* antenna phase center and phase windup correction */
        L[i]=obs->L[i]*CLIGHT/freq[i]-dants[i]-dantr[i]-phw*CLIGHT/freq[i];
        P[i]=obs->P[i]-dants[i]-dantr[i];

        /* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) */
        if (sys==SYS_GPS||sys==SYS_GLO) {
            if (obs->code[i]==CODE_L1C) P[i]+=nav->cbias[obs->sat-1][1];
            if (obs->code[i]==CODE_L2C) P[i]+=nav->cbias[obs->sat-1][2];
        }
    }
    /* iono-free LC */
    *Lc=*Pc=0.0;
    if (freq[0]==0.0||freq[1]==0.0) return;
    C1= SQR(freq[0])/(SQR(freq[0])-SQR(freq[1]));
    C2=-SQR(freq[1])/(SQR(freq[0])-SQR(freq[1]));

    if (L[0]!=0.0&&L[1]!=0.0) *Lc=C1*L[0]+C2*L[1];
    if (P[0]!=0.0&&P[1]!=0.0) *Pc=C1*P[0]+C2*P[1];
}
/* detect cycle slip by LLI --------------------------------------------------*/
static void detslp_ll(rtk_t *rtk, const obsd_t *obs, int n)
{
#if 0
    int i,j,nf;
    trace(4,"detslp_ll: n=%d\n",n);
    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<rtk->opt.nf;j++) {
        if (obs[i].L[j]==0.0||!(obs[i].LLI[j]&3)) continue;
        trace(3,"detslp_ll: slip detected %s f=%d el=%4.2f\n",sat_id(obs[i].sat),j+1,rtk->ssat[obs[i].sat-1].azel[1]*R2D);
        nf=rtk->opt.nf;/* fixes gcc compiler warning */
        rtk->ssat[obs[i].sat-1].slip[j<nf?j:nf]=1;
    }
#else
    unsigned int slip,LLI;
    int i,f,sat;

    for(i=0;i<n;i++){
        sat=obs[i].sat;
        for (f=0;f<NFREQ;f++) {

            if ((obs[i].L[f]==0.0&&obs[i].LLI[f]==0)||
                fabs(timediff(obs[i].time,rtk->ssat[sat-1].pt[0][f]))<DTTOL) {
                continue;
            }
            /* restore previous LLI */
            LLI=getbitu(&rtk->ssat[sat-1].slip[f],0,2); /* rover */

            /* detect slip by cycle slip flag in LLI */
            if (rtk->tt>=0.0) { /* forward */
                if (obs[i].LLI[f]&1) {
                    trace(3,"%s(%d): slip detected forward (sat=%s F=%d LLI=%x)\n",
                          time_str(obs[i].time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(sat),f+1,obs[i].LLI[f]);
                }
                slip=obs[i].LLI[f];
            }
            else { /* backward */
                if (LLI&1) {
                    trace(3,"%s(%d): slip detected backward (sat=%s F=%d LLI=%x)\n",
                          time_str(obs[i].time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(sat),f+1,LLI);
                }
                slip=LLI;
            }
            /* detect slip by parity unknown flag transition in LLI */
            if (((LLI&2)&&!(obs[i].LLI[f]&2))||(!(LLI&2)&&(obs[i].LLI[f]&2))) {
                trace(3,"%s(%d): slip detected half-cyc (sat=%s F=%d LLI=%x->%x)\n",
                       time_str(obs[i].time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(sat),f+1,LLI,obs[i].LLI[f]);
                slip|=1;
            }
            /* save current LLI */
            setbitu(&rtk->ssat[sat-1].slip[f],0,2,obs[i].LLI[f]);

            /* save slip and half-cycle valid flag */
            rtk->ssat[sat-1].slip[f]|=(unsigned char)slip;
            rtk->ssat[sat-1].half[f]=(obs[i].LLI[f]&2)?0:1;
        }
    }
#endif
}

/* detect cycle slip by geometry free phase jump -----------------------------*/
static void detslp_gf(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    /*精度较高，可以探测小周跳，但对比值接近于lam2/lam1的(L1,L2)周跳不敏感*/
    double g0,g1,el,thres;
    int i,j;
    
    trace(4,"detslp_gf: n=%d\n",n);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        
        if ((g1=gfmeas(obs+i,nav))==0.0){
//            rtk->ssat[obs[i].sat-1].gf[0]=0.0;
            rtk->ssat[obs[i].sat-1].delta_gf[0]=rtk->ssat[obs[i].sat-1].delta_gf[1]=0.0;
            continue;
        }
        
        g0=rtk->ssat[obs[i].sat-1].gf[0];
        if(g0==0.0) continue;
        rtk->ssat[obs[i].sat-1].gf[0]=g1;

        trace(4,"detslip_gf: %s gf0=%8.3f gf1=%8.3f\n",sat_id(obs[i].sat),g0,g1);
        el=rtk->ssat[obs[i].sat-1].azel[1]*R2D;

        if(el<rtk->opt.elmin*R2D) continue;

        if(el>=15.0)  thres=rtk->opt.cs_gf;
        else thres=(-0.4*el+7.0)*rtk->opt.cs_gf;

        rtk->ssat[obs[i].sat-1].delta_gf[0]=fabs(g1-g0);
        rtk->ssat[obs[i].sat-1].delta_gf[1]=thres;
        if (g0!=0.0&&fabs(g1-g0)>thres) {
            trace(2,"%s(%d): gf slip detected %s gf=%8.3f-> %8.3f el=%4.2f thres=%4.2f\n",
                  time_str(obs[i].time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(obs[i].sat),g0,g1,rtk->ssat[obs[i].sat-1].azel[1]*R2D,thres);
            
            for (j=0;j<rtk->opt.nf;j++) rtk->ssat[obs[i].sat-1].slip[j]|=1;
        }
    }
}
/* detect slip by Melbourne-Wubbena linear combination jump ------------------*/
static void detslp_mw(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    /*只能探测不小于三周的周跳，且对L1\L2大小相等的周跳组合不敏感*/
    double w0,w1;
    int i,j,sat;
    double thres,el;
    
    trace(4,"detslp_mw: n=%d\n",n);


    for (i=0;i<n&&i<MAXOBS;i++) {
        sat=obs[i].sat;
        if (fabs(timediff(rtk->sol.time, rtk->ssat[sat - 1].ct)) > rtk->opt.sample+1.0) {
            rtk->ssat[sat-1].mw[1] = 0.0;
            rtk->ssat[sat-1].mw[2] = 0.0;
            rtk->ssat[sat-1].delta_mw[0]=rtk->ssat[sat-1].delta_mw[1]=0.0;
        }

        if ((w1=mwmeas(&rtk->opt,obs+i,nav,NULL,0.0))==0.0) continue;

        w0=rtk->ssat[obs[i].sat-1].mw[0];
        rtk->ssat[obs[i].sat-1].mw[0]=w1;
        
        trace(4,"detslip_mw: %s mw0=%8.3f mw1=%8.3f\n",sat_id(obs[i].sat),w0,w1);
        if(w0==0.0) continue;
        el=rtk->ssat[obs[i].sat-1].azel[1];

        if(el<rtk->opt.elmin) continue;

        if(el*R2D>15.0) thres=rtk->opt.cs_mw;
        else thres=(-0.2*el*R2D+4.0)*rtk->opt.cs_mw;

        rtk->ssat[obs[i].sat-1].delta_mw[0]=fabs(w1-w0);
        rtk->ssat[obs[i].sat-1].delta_mw[1]=thres;
        if (w0!=0.0&&fabs(w1-w0)>THRES_MW_JUMP) {
            trace(2,"%s(%d): mw slip detected %s mw=%8.3f->%8.3f el=%4.2f thres=%4.2f\n",
                  time_str(obs[i].time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(obs[i].sat),w0,w1,rtk->ssat[obs[i].sat-1].azel[1]*R2D,thres);
            
            for (j=0;j<rtk->opt.nf;j++) rtk->ssat[obs[i].sat-1].slip[j]|=1;
        }
    }
}

static void saveinfo(const obsd_t *obs,rtk_t *rtk,int n,const nav_t *nav)
{
    int i,j,sat;
    prcopt_t *popt=&rtk->opt;
    double w0,w1,var0,var1,gf;

    for(i=0;i<MAXSAT;i++){
        rtk->ssat[i].mw[0]=0.0;
        rtk->ssat[i].gf[0]=0.0;
    }

    for (i=0;i<n&&i<MAXOBS;i++) {
        sat=obs[i].sat;
        rtk->ssat[sat-1].ct=rtk->sol.time;

        if ((gf=gfmeas(obs+i,nav))!=0.0)
            rtk->ssat[sat-1].gf[0]=gf;

        if ((w1=mwmeas(&rtk->opt,obs+i,nav,&var1,rtk->ssat[sat-1].azel[1]))==0.0) {
            continue;
        }

        w0=rtk->ssat[sat-1].mw[1];
        rtk->ssat[sat-1].mw[0]=w1;

        if (rtk->ssat[sat-1].mw[2]>0){
            double K=0.0;
            K=rtk->ssat[sat-1].mw[3]/(rtk->ssat[sat-1].mw[3]+var1);
            rtk->ssat[sat-1].mw[1]+=K*(w1-w0);
            rtk->ssat[sat-1].mw[3]-=K*rtk->ssat[sat-1].mw[3];
            rtk->ssat[sat-1].mw[2]++;
//            j=rtk->ssat[sat-1].mw[2];
//            var0=rtk->ssat[sat-1].mw[3];
//            var1=(w1-w0)*(w1-w0)-var0;
//            var1=var0+var1/j;
//
//            rtk->ssat[sat-1].mw[1]=(w0*j+w1)/(j+1);
//            rtk->ssat[sat-1].mw[2]++;
//            rtk->ssat[sat-1].mw[3]=var1;
        }
        else{
            rtk->ssat[sat-1].mw[1]=w1;
            rtk->ssat[sat-1].mw[2]++;
            rtk->ssat[sat-1].mw[3]=SQR(5);
        }
    }
}

static void detecs_ppp(const obsd_t *obs, rtk_t *rtk, int n, const nav_t *nav)
{
    int i,j;

    if(rtk->opt.sample<=1){
        rtk->opt.cs_mw=1.0; /*cycle*/
        rtk->opt.cs_gf=0.05; /*m*/
    }
    else if(rtk->opt.sample<=15){
        rtk->opt.cs_mw=1.5;
        rtk->opt.cs_gf=0.10;
    }
    else{
        rtk->opt.cs_mw=2.0;
        rtk->opt.cs_gf=0.15;
    }

    for (i=0;i<MAXSAT;i++){
        for (j=0;j<NFREQ;j++){
            rtk->ssat[i].slip[j]=0;
        }
    }

    detslp_ll(rtk,obs,n);
    if(rtk->opt.nf>=2){
        /*MW GF方法消除了卫星至接收机的几何距离影响，即与接收机的运动状态无关，因此，
         * 这两种方法非常适合动态数据的周跳探测*/
        detslp_mw(rtk,obs,n,nav);
        detslp_gf(rtk,obs,n,nav);
        saveinfo(obs,rtk,n,nav);
    }
}

/* temporal update of position -----------------------------------------------*/
static void udpos_ppp(rtk_t *rtk)
{
    double *F,*P,*FP,*x,*xp,pos[3],Q[9]={0},Qv[9];
    int i,j,*ix,nx,restart=0;

    if(rtk->opt.restart>0&&timediff(rtk->sol.time,rtk->filter_start)>rtk->opt.restart*3600.0){
        restart=0;
    }
    
    trace(4,"udpos_ppp:\n");


    /* fixed mode */
    if (rtk->opt.mode==PMODE_PPP_FIXED) {
        for (i=0;i<3;i++) initx(rtk,rtk->opt.ru[i],1E-8,i);
        return;
    }
    /* initialize position for first epoch */
    if (norm(rtk->x,3)<=0.0||restart) {
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        if (rtk->opt.dynamics) {
            for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
            for (i=6;i<9;i++) initx(rtk,1E-6,VAR_ACC,i);
        }
    }
    /* static ppp mode */
    if (rtk->opt.mode==PMODE_PPP_STATIC) {
        for (i=0;i<3;i++) {
            rtk->P[i*(1+rtk->nx)]+=SQR(rtk->opt.prn[5])*fabs(rtk->tt);
        }
        return;
    }
    /* kinmatic mode without dynamics */
    if (!rtk->opt.dynamics) {
        for (i=0;i<3;i++) {
            initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        }
        return;
    }
    /* generate valid state index */
    ix=imat(rtk->nx,1);
    for (i=nx=0;i<rtk->nx;i++) {
        if (rtk->x[i]!=0.0&&rtk->P[i+i*rtk->nx]>0.0) ix[nx++]=i;
    }
    if (nx<9) {
        free(ix);
        return;
    }
    /* state transition of position/velocity/acceleration */
    F=eye(nx); P=mat(nx,nx); FP=mat(nx,nx); x=mat(nx,1); xp=mat(nx,1);


    for (i=0;i<6;i++) {
        F[i+(i+3)*nx]=rtk->tt;
    }
    for (i=0;i<3;i++) {
        F[i+(i+6)*nx]=SQR(rtk->tt)/2.0;
    }
    for (i=0;i<nx;i++) {
        x[i]=rtk->x[ix[i]];
        for (j=0;j<nx;j++) {
            P[i+j*nx]=rtk->P[ix[i]+ix[j]*rtk->nx];
        }
    }

    /* x=F*x, P=F*P*F+Q */
    matmul("NN",nx,1,nx,1.0,F,x,0.0,xp);
    matmul("NN",nx,nx,nx,1.0,F,P,0.0,FP);
    matmul("NT",nx,nx,nx,1.0,FP,F,0.0,P);
    
    for (i=0;i<nx;i++) {
        rtk->x[ix[i]]=xp[i];
        for (j=0;j<nx;j++) {
            rtk->P[ix[i]+ix[j]*rtk->nx]=P[i+j*nx];
        }
    }
    /* process noise added to only acceleration */
    Q[0]=Q[4]=SQR(rtk->opt.prn[3])*fabs(rtk->tt);
    Q[8]=SQR(rtk->opt.prn[4])*fabs(rtk->tt);
    ecef2pos(rtk->x,pos);
    covecef(pos,Q,Qv);
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
        rtk->P[i+6+(j+6)*rtk->nx]+=Qv[i+j*3];
    }
    free(ix); free(F); free(P); free(FP); free(x); free(xp);
}
/* temporal update of clock --------------------------------------------------*/
static void udclk_ppp(rtk_t *rtk,const obsd_t *obs,int n)
{
    double dtr;
    int i,sat,sat_idx,main_irc,irc,num_sys=NSYS;
    static int last_main_irc;

    trace(4,"udclk_ppp:\n");

    if((rtk->opt.navsys&SYS_BD3)&&rtk->opt.bd3opt>=BD3OPT_BD2_3) num_sys+=1;

    for(i=0;i<NSYS+1;i++){
        for(int j=0;j<n;j++){
            sat=obs[j].sat;
            sat_idx=satsysidx(sat);
            if(i==sat_idx){
                rtk->exist_sys[i]=1;
                break;
            }
        }
    }

    for(i=0;i<num_sys;i++) if(rtk->sol.dtr[i]!=0.0&&rtk->exist_sys[i]){
        main_irc=i;
            break;
    }

    /* initialize every epoch for clock (white noise) */
    for (i=0;i<num_sys;i++) {
        irc=IC(i,&rtk->opt);
        dtr=rtk->sol.dtr[i];
        if(rtk->x[irc]==0.0&&rtk->exist_sys[i]){
            /// init
            initx(rtk,CLIGHT*dtr,VAR_CLK,irc); /*初始化clk和isb*/
        }
        else{
            /// translate
            if((rtk->opt.navsys&SYS_GPS)&&rtk->exist_sys[0]){
                if(i==0) initx(rtk,CLIGHT*dtr,VAR_CLK,irc+i);  /*GPS为主系统*/
                else{
                    rtk->P[irc+irc*rtk->nx]+=SQR(0.0001)*fabs(rtk->tt);
                }
            }
            else if((rtk->opt.navsys&SYS_GLO)&&rtk->exist_sys[1]){
                if(i==0) initx(rtk,0,SQR(0.01),irc);
                else if(i==1) initx(rtk,CLIGHT*dtr,VAR_CLK,irc);
                else{
                    if(fabs(dtr-rtk->x[irc])>50) initx(rtk,CLIGHT*dtr,VAR_CLK,irc);
                    rtk->P[irc+irc*rtk->nx]+=SQR(0.0001)*fabs(rtk->tt);
                }
            }
            else if((rtk->opt.navsys&SYS_GAL)&&rtk->exist_sys[2]){
                if(i==0) initx(rtk,0.0,SQR(0.01),irc);
                if(i==1) initx(rtk,0.0,SQR(0.01),irc);
                else if(i==2) initx(rtk,CLIGHT*dtr,VAR_CLK,irc);
                else{
                    if(last_main_irc!=main_irc){
                        initx(rtk,CLIGHT*dtr,VAR_CLK,irc);
                    }
                    else{
                        rtk->P[irc+irc*rtk->nx]+=SQR(0.0001)*fabs(rtk->tt);
                    }
//                    if(fabs(dtr-rtk->x[irc])>50) initx(rtk,CLIGHT*dtr,VAR_CLK,irc);
                }
            }
            else if((rtk->opt.navsys&SYS_CMP)&&rtk->exist_sys[3]){
                if(i==0) initx(rtk,0.0,SQR(0.01),irc);
                else if(i==1) initx(rtk,0.0,SQR(0.01),irc);
                else if(i==2) initx(rtk,0.0,SQR(0.01),irc);
                else if(i==3) initx(rtk,CLIGHT*dtr,VAR_CLK,irc);
                else{
                    if(fabs(dtr-rtk->x[irc])>50) initx(rtk,CLIGHT*dtr,VAR_CLK,irc);
                    rtk->P[irc+irc*rtk->nx]+=SQR(0.0001)*fabs(rtk->tt);
                }
            }
        }
    }
    last_main_irc=main_irc;
}

/* temporal update of L5-receiver-dcb parameters -----------------------------*/
static void uddcb_ppp(rtk_t *rtk)
{
    /// TODO: tc
    int i,num_sys=NSYS;


    trace(4,"uddcb_ppp:\n");

    if((rtk->opt.navsys&SYS_BD3)&&rtk->opt.bd3opt>=BD3OPT_BD2_3) num_sys+=1;

    for(int j=0;j<num_sys;j++){
        i=IRDCB(j,&rtk->opt);
        if (rtk->x[i]==0.0) {
            if(j==0&&rtk->opt.navsys&SYS_GPS) initx(rtk,0.1,VAR_DCB,i);
            if(j==1&&rtk->opt.navsys&SYS_GLO) initx(rtk,0.1,VAR_DCB,i);
            if(j==2&&rtk->opt.navsys&SYS_GAL) initx(rtk,0.1,VAR_DCB,i);
            if(j==3&&rtk->opt.navsys&SYS_CMP) initx(rtk,0.1,VAR_DCB,i);
            if(j==4&&rtk->opt.navsys&SYS_QZS) initx(rtk,0.1,VAR_DCB,i);
            if(j==5&&(rtk->opt.navsys&SYS_GPS)&&rtk->opt.bd3opt>=BD3OPT_BD2_3) initx(rtk,0.1,VAR_DCB,i);
        }
        else{
#if 0
            rtk->P[i+i*rtk->nx]+=SQR(0.00001)*fabs(rtk->tt);
#else
//            initx(rtk,rtk->x[i],VAR_DCB,i);
#endif
        }
    }
}

static void udifcb_ppp(rtk_t *rtk)
{
    int i,num_sys=NSYS;


    trace(4,"uddcb_ppp:\n");

    if((rtk->opt.navsys&SYS_BD3)&&rtk->opt.bd3opt>=BD3OPT_BD2_3) num_sys+=1;

    for(int j=0;j<num_sys;j++){
        i=IIFCB(j,&rtk->opt);
        if (rtk->x[i]==0.0) {
            if(j==0&&rtk->opt.navsys&SYS_GPS) initx(rtk,0.1,VAR_IFCB,i);
            if(j==1&&rtk->opt.navsys&SYS_GLO) initx(rtk,0.1,VAR_IFCB,i);
            if(j==2&&rtk->opt.navsys&SYS_GAL) initx(rtk,0.1,VAR_IFCB,i);
            if(j==3&&rtk->opt.navsys&SYS_CMP) initx(rtk,0.1,VAR_IFCB,i);
            if(j==4&&rtk->opt.navsys&SYS_QZS) initx(rtk,0.1,VAR_IFCB,i);
            if(j==5&&(rtk->opt.navsys&SYS_GPS)&&rtk->opt.bd3opt>=BD3OPT_BD2_3) initx(rtk,0.1,VAR_IFCB,i);
        }
        else{
#if 0
            rtk->P[i+i*rtk->nx]+=SQR(0.00001)*fabs(rtk->tt);
#else
//            initx(rtk,rtk->x[i],VAR_DCB,i);
#endif
        }
    }
}

/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop_ppp(rtk_t *rtk)
{
    double pos[3],azel[]={0.0,PI/2.0},ztd,var;
    int i,j,restart=0;
    
    trace(4,"udtrop_ppp:\n");

    if(rtk->opt.restart>0&&timediff(rtk->sol.time,rtk->filter_start)>rtk->opt.restart*3600.0){
        restart=0;
    }

    i=IT(&rtk->opt);
    if (rtk->x[i]==0.0) {
        initx(rtk,0.15,rtk->opt.std[2]==0.0?VAR_ZWD:SQR(rtk->opt.std[2]),i);
        if (rtk->opt.tropopt>=TROPOPT_ESTG) {
            for (j=i+1;j<i+3;j++) initx(rtk,1E-6,VAR_GRA,j);
        }
    }
    else {
        rtk->P[i+i*rtk->nx]+=SQR(rtk->opt.prn[2])*fabs(rtk->tt);
        
        if (rtk->opt.tropopt>=TROPOPT_ESTG) {
            for (j=i+1;j<i+3;j++) {
                rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[2]*0.1)*fabs(rtk->tt);
            }
        }
    }
}
/* temporal update of ionospheric parameters ---------------------------------*/
static void udiono_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    double freq1,freq2,ion,sinel,pos[3],*azel,ion_var,beta=0.0,frq1=0.0,frq2=0.0,P1=0.0,P2=0.0;
    char *p;
    int i,j,sat,sys,sys_idx=-1,prn,gap_resion=GAP_RESION,restart=0;

    if(rtk->opt.restart>0&&timediff(rtk->sol.time,rtk->filter_start)>rtk->opt.restart*3600.0){
        restart=0;
    }
    
    trace(3,"udiono_ppp:\n");
    
    if ((p=strstr(rtk->opt.pppopt,"-GAP_RESION="))) {
        sscanf(p,"-GAP_RESION=%d",&gap_resion);
    }
    for (i=0;i<MAXSAT;i++) {
        j=II(i+1,&rtk->opt);
        if (rtk->x[j]!=0.0&&(int)rtk->ssat[i].outc[0]>gap_resion) {
            rtk->x[j]=0.0;
        }
    }

    ecef2pos(rtk->sol.rr,pos);
    for (i=0;i<n;i++) {
        sat=obs[i].sat;
        sys=satsys(sat,&prn);
        sys_idx=satsysidx(sat);
        if(sys_idx==-1) continue;
        if(rtk->ssat[sat-1].azel[1]<rtk->opt.elmin) continue;
        j=II(obs[i].sat,&rtk->opt);
        if (rtk->x[j]==0.0||restart) { /*init*/
            if(rtk->opt.ionoopt==IONOOPT_UC_CONS&&rtk->opt.nf==1){ /*single-freq or ionospheric-constraint*/
                iontec(obs[i].time,nav,pos,rtk->ssat[obs[i].sat-1].azel,1,&ion,&ion_var);
            }
            else{ /*dual-freq triple-freq case*/
                frq1=sat2freq(sat,obs[i].code[rtk->opt.gnss_frq_idx[sys_idx][0]-1],nav);
                frq2=sat2freq(sat,obs[i].code[rtk->opt.gnss_frq_idx[sys_idx][1]-1],nav);
                P1=obs[i].P[rtk->opt.gnss_frq_idx[sys_idx][0]-1];
                P2=obs[i].P[rtk->opt.gnss_frq_idx[sys_idx][1]-1];

                if(frq1==0.0||frq2==0.0||P1==0.0||P2==0.0) continue;

                beta=-SQR(frq2)/(SQR(frq1)-SQR(frq2));
                ion=(P1-P2)*beta;
            }
            initx(rtk,ion,VAR_IONO,j);
        }
        else {
#if 0
            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[1])*fabs(rtk->tt);
#else
            initx(rtk,rtk->x[j],VAR_IONO,j);
#endif
        }
    }
}

/* temporal update of phase biases -------------------------------------------*/
static void udbias_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    double L[NFREQ]={0},P[NFREQ]={0},Lc[NFREQ]={0},Pc[NFREQ]={0},freqs[NFREQ],bias[MAXOBS]={0},offset=0.0,pos[3]={0};
    int i,j,k,f,sat,slip[MAXOBS]={0},clk_jump=0,sys,sys_idx=-1,prn;
    int iamb=0,iion,restart=0;
    double dcb[3]={0};

    trace(4,"udbias  : n=%d\n",n);

    if(rtk->opt.restart>0&&timediff(rtk->sol.time,rtk->filter_start)>rtk->opt.restart*3600.0){
        restart=1;
    }

    /* handle day-boundary clock jump */
    if (rtk->opt.posopt[5]) {
        clk_jump=ROUND(time2gpst(obs[0].time,NULL)*10)%864000==0;
    }

    ecef2pos(rtk->sol.rr,pos);
    
    for (f=0;f<NF(&rtk->opt);f++) {
        
        /* reset phase-bias if expire obs outage counter */
        for (i=0;i<MAXSAT;i++) {
            if (++rtk->ssat[i].outc[f]>(unsigned int)rtk->opt.maxout||
                rtk->opt.modear==ARMODE_INST||clk_jump) {
                iamb=IB(i+1,f,&rtk->opt);
                initx(rtk,0.0,0.0,iamb);
            }
        }
        for (i=k=0;i<n&&i<MAXOBS;i++) {
            sat=obs[i].sat;
            sys=satsys(sat,&prn);
            sys_idx=satsysidx(sat);
            if(sys_idx==-1) continue;
            if(rtk->ssat[sat-1].azel[1]<rtk->opt.elmin) continue;
            j=IB(sat,f,&rtk->opt);
            iion=II(sat,&rtk->opt);

            getcorrobs(&rtk->opt,obs+i,nav,rtk->opt.gnss_frq_idx[sys_idx],NULL,NULL,0.0,L,P,Lc,Pc,freqs,NULL,NULL);
            bias[i]=0.0;
            if (rtk->opt.ionoopt==IONOOPT_IFLC||rtk->opt.ionoopt==IONOOPT_IF2) {  /*ionospheric-free*/
                if(rtk->opt.nf==1){ /* UofC */
                    bias[i]=(L[f]-P[f])*0.5;
                }
                else if(rtk->opt.nf==2){ /* dual-frequency ionospheric-free */
                    bias[i]=Lc[f]-Pc[f];
                    slip[i]=rtk->ssat[sat-1].slip[0]||rtk->ssat[sat-1].slip[1];
                }
                else if(rtk->opt.nf==3){
                    if(rtk->opt.ionoopt==IONOOPT_IFLC){

                    }
                    else if(rtk->opt.ionoopt==IONOOPT_IF2){

                    }
                }
            }
            else if (L[f]!=0.0&&P[f]!=0.0) {
                double frq1=0.0,frq2=0.0;
                if(sys==SYS_GPS){
                    frq1=FREQ1;
                    frq2=sat2freq(sat,obs[i].code[rtk->opt.gnss_frq_idx[sys_idx][f]-1],nav);
                }
                else if(sys==SYS_GLO){
                    frq1=sat2freq(sat,CODE_L1P,nav);
                    frq2=sat2freq(sat,obs[i].code[rtk->opt.gnss_frq_idx[sys_idx][f]-1],nav);
                }
                else if(sys==SYS_GAL){
                    frq1=FREQ1;
                    frq2=sat2freq(sat,obs[i].code[rtk->opt.gnss_frq_idx[sys_idx][f]-1],nav);
                }
                else if(sys==SYS_CMP){
                    if(prn<=18){
                        frq1=FREQ1_CMP;
                        frq2=sat2freq(sat,obs[i].code[rtk->opt.gnss_frq_idx[sys_idx][f]-1],nav);
                    }
                    else{
                        frq1=FREQ1_CMP;
                        frq2=sat2freq(sat,obs[i].code[rtk->opt.gnss_frq_idx[sys_idx][f]-1],nav);
                    }
                }
                else if(sys==SYS_QZS){
                    frq1=FREQ1;
                    frq2=sat2freq(sat,obs[i].code[rtk->opt.gnss_frq_idx[sys_idx][f]-1],nav);
                }
                slip[i]=rtk->ssat[sat-1].slip[0]||rtk->ssat[sat-1].slip[1];
                bias[i]=L[f]-P[f]+2.0*rtk->x[iion]*SQR(frq1/frq2);
            }
            if (rtk->x[j]==0.0||slip[i]||bias[i]==0.0) continue;
            
            offset+=bias[i]-rtk->x[j];
            k++;
        }
        /* correct phase-code jump to ensure phase-code coherency */
#if 0
        if (k>=2&&fabs(offset/k)>0.0005*CLIGHT) {
            for (i=0;i<MAXSAT;i++) {
                j=rtk->tc?xiAmb(&rtk->opt.insopt,i+1,f):IB(i+1,f,&rtk->opt);
                if (rtk->x[j]!=0.0) rtk->x[j]+=offset/k;
            }
            trace(2,"phase-code jump corrected: %s n=%2d dt=%12.9fs\n",
                  time_str(rtk->sol.time,0),k,offset/k/CLIGHT);
        }
#endif
        for (i=0;i<n&&i<MAXOBS;i++) {
            sat=obs[i].sat;
            j=IB(sat,f,&rtk->opt);

            if(rtk->opt.nf==1&&rtk->opt.ionoopt==IONOOPT_IFLC) rtk->opt.prn[0]=1E-4;

            if(rtk->opt.ionoopt==IONOOPT_UC&&rtk->opt.nf==3){
                if(sys==SYS_GPS&&f==2) rtk->P[j+j*rtk->nx]+=3E-5*fabs(rtk->tt);
                else if(sys==SYS_GAL&&f==1){
                    rtk->P[j+j*rtk->nx]+=SQR(3E-5)*fabs(rtk->tt);
                }
                else if(sys==SYS_CMP){
                    if(!strcmp(rtk->opt.ac_name,"com")){
                        if(f==2) rtk->P[j+j*rtk->nx]+=SQR(3E-5)*fabs(rtk->tt);
                        else{
                            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[0])*fabs(rtk->tt);
                        }
                    }
                    else{
                        if(f==1) rtk->P[j+j*rtk->nx]+=SQR(3E-5)*fabs(rtk->tt);
                        else{
                            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[0])*fabs(rtk->tt);
                        }
                    }
                }
                else{
                    rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[0])*fabs(rtk->tt);
                }
            }
            else{
                rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[0])*fabs(rtk->tt);
            }

            if(restart){
                iamb=IB(sat,f,&rtk->opt);
                initx(rtk,bias[i],SQR(rtk->opt.std[0]),iamb);
            }
            else{
                if (bias[i]==0.0||(rtk->x[j]!=0.0&&!slip[i])) continue;

                /* reinitialize phase-bias if detecting cycle slip */
                iamb=IB(sat,f,&rtk->opt);
                initx(rtk,bias[i],SQR(rtk->opt.std[0]),iamb);

                /* reset fix flags */
                for (k=0;k<MAXSAT;k++) rtk->ambc[sat-1].flags[k]=0;


                trace(3,"udbias_ppp: %s L%s_%d bias=%8.3f lock=%5d slip=%d el=%4.2f\n",
                      sat_id(sat),rtk->opt.ionoopt==IONOOPT_IFLC?"IF":"UC",f+1,bias[i],rtk->ssat[sat-1].lock[f],rtk->ssat[sat-1].slip[f],rtk->ssat[sat-1].azel[1]*R2D);
                rtk->ssat[sat-1].lock[f]=-rtk->opt.minlock;
            }
        }
    }
}
/* temporal update of states --------------------------------------------------*/
static void udstate_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    trace(4,"udstate_ppp: n=%d\n",n);
    
    /* temporal update of position */
    udpos_ppp(rtk);
    
    /* temporal update of clock */
    if(!rtk->opt.sdopt){
        udclk_ppp(rtk,obs,n);
    }

    /* temporal update of receiver-dcb parameters */
    if (rtk->opt.ionoopt==IONOOPT_UC_CONS&&rtk->opt.nf>=2) {
        uddcb_ppp(rtk);
    }

    /* temporal update of inter-frequency clock bias parameters */
    if(rtk->opt.nf>=3){
        udifcb_ppp(rtk);
    }
    
    /* temporal update of tropospheric parameters */
    if (rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG) {
        udtrop_ppp(rtk);
    }
    /* temporal update of ionospheric parameters */
    if (rtk->opt.ionoopt==IONOOPT_UC||rtk->opt.ionoopt==IONOOPT_UC_CONS) {
        udiono_ppp(rtk,obs,n,nav);
    }

    /* temporal update of phase-bias */
    udbias_ppp(rtk,obs,n,nav);

    if(rtk->opt.restart>0&&timediff(rtk->sol.time,rtk->filter_start)>rtk->opt.restart*3600.0){
        rtk->filter_start=rtk->sol.time;
    }
}

extern void clk_repair(obsd_t* obs,int n,rtk_t *rtk)
{
    int i,sat,vG=0,cjG=0;
    double delta0=0.0,delta1=0.0,d1,d2,d3,d4,ddd1,ddd2;
    double lam1,lam2;
    double CJ_F1,CJ_F2;
    int flag[MAXPRNGPS];
    sat_model_t smod={0};

    for(i=0;i<MAXPRNGPS;i++) flag[i]=0;

    vG=cjG=0;

    for(i=0;i<n;i++){
        sat=obs[i].sat;
        smod=rtk->smod[sat-1];
        if(sat>MAXPRNGPS) continue;

        if(obs[i].P[0]*obs[i].P[0]*obs[i].L[0]*obs[i].L[0]==0.0) continue;

        if(smod.raw_P[0]*smod.raw_P[1]*smod.raw_L[0]*smod.raw_L[1]==0.0) continue;

        lam1=CLIGHT/sat2freq(sat,obs[i].code[0],NULL);
        lam2=CLIGHT/sat2freq(sat,obs[i].code[1],NULL);

        vG++;
        d1=obs[i].P[0]-smod.raw_P[0];
        d2=obs[i].P[1]-smod.raw_P[1];
        d3=(obs[i].L[0]-smod.raw_L[0])*lam1;
        d4=(obs[i].L[1]-smod.raw_L[1])*lam2;

        if(fabs(d1-d3)>290000){
            delta0+=d1-d3;
            delta1+=d2-d4;
            cjG++;
        }
    }

    if(cjG!=0&&cjG==vG){
        d1=delta0/cjG;
        d2=delta1/cjG;

        CJ_F1=d1/CLIGHT*1000.0;
        CJ_F2=newround(CJ_F1);
        if(fabs(CJ_F1-CJ_F2)<2.5E-2){
            trace(2,"%s: receiver clock jump\n",CJ_F2);
        }
    }

    for(i=0;i<n;i++){
        sat=obs[i].sat;
        lam1=CLIGHT/sat2freq(sat,obs[i].code[0],NULL);
        lam2=CLIGHT/sat2freq(sat,obs[i].code[1],NULL);

        if(sat>MAXPRNGPS) continue;
        flag[sat-1]=1;

        rtk->smod[sat-1].raw_P[0]=obs[i].P[0];
        rtk->smod[sat-1].raw_P[1]=obs[i].P[1];
        rtk->smod[sat-1].raw_L[0]=obs[i].L[0];
        rtk->smod[sat-1].raw_L[1]=obs[i].L[1];

        ddd1=rtk->clk_jump*CLIGHT/1000.0;
        ddd2=rtk->clk_jump*CLIGHT/1000.0;

        if(obs[i].L[0]!=0.0) obs[i].L[0]+=ddd1/lam1;
        if(obs[i].L[1]!=0.0) obs[i].L[1]+=ddd2/lam2;
    }
}

extern void bd2_multipath(rtk_t *rtk,obsd_t *obs,int n)
{
#if 1
    int i,j,sat,prn,b;
    double dmp[3],elev,a;
    const static double IGSOCOEF[3][10]={		/* m */
            {-0.55,-0.40,-0.34,-0.23,-0.15,-0.04,0.09,0.19,0.27,0.35},	//B1
            {-0.71,-0.36,-0.33,-0.19,-0.14,-0.03,0.08,0.17,0.24,0.33},	//B2
            {-0.27,-0.23,-0.21,-0.15,-0.11,-0.04,0.05,0.14,0.19,0.32},	//B3
    };
    const static double MEOCOEF[3][10]={		/* m */
            {-0.47,-0.38,-0.32,-0.23,-0.11,0.06,0.34,0.69,0.97,1.05},	//B1
            {-0.40,-0.31,-0.26,-0.18,-0.06,0.09,0.28,0.48,0.64,0.69},	//B2
            {-0.22,-0.15,-0.13,-0.10,-0.04,0.05,0.14,0.27,0.36,0.47},	//B3
    };

    for(i=0;i<n;i++){
        sat=obs[i].sat;
        if (satsys(sat,&prn)!=SYS_CMP) continue;
        if (prn<=5||prn>=19) continue;
        elev=rtk->ssat[sat-1].azel[1]*R2D;
        if (elev<=rtk->opt.elmin*R2D) continue;

        for (j=0;j<3;j++) dmp[j]=0.0;

        a=elev*0.1;
        b=(int)a;
        if ((prn>=6&&prn<=10)||(prn==13||prn==16)) { // IGSO(C06, C07, C08, C09, C10,C13,C16)
            if (b<0) {
                for (j=0;j<3;j++) dmp[j]=IGSOCOEF[j][0];
            }
            else if (b>=9) {
                for (j=0;j<3;j++) dmp[j]=IGSOCOEF[j][9];
            }
            else {
                for (j=0;j<3;j++) dmp[j]=IGSOCOEF[j][b]*(1.0-a+b)+IGSOCOEF[j][b+1]*(a-b);
            }
        }
        else if (prn==11||prn==12||prn==14) {   // MEO(C11, C12, C13, C14)
            if (b<0) {
                for (j=0;j<3;j++) dmp[j]=MEOCOEF[j][0];
            }
            else if (b>=9) {
                for (j=0;j<3;j++) dmp[j]=MEOCOEF[j][9];
            }
            else {
                for (j=0;j<3;j++) dmp[j]=MEOCOEF[j][b]*(1.0-a+b)+MEOCOEF[j][b+1]*(a-b);
            }
        }
        if(obs[i].P[0]!=0.0) obs[i].P[0]+=dmp[0];
        if(obs[i].P[1]!=0.0) obs[i].P[1]+=dmp[1];
        if(obs[i].P[3]!=0.0) obs[i].P[3]+=dmp[2];
    }
#else
    int i,j,sat,sys,prn;
    double el,mp=0.0;
    const static double geo_coef[3][3]={
            {-0.436, 1.158, -0.333}, /*B1: a1 a2 a3*/
            {-0.275, 1.087, -0.452}, /*B2: */
            {-0.048, 0.566, -0.185}, /*B3: */
    };
    const static double igso_coef[3][3]={
            {-0.590, 1.624, -0.645},
            {-0.257, 0.995, -0.381},
            {-0.102, 0.748, -0.307},
    };
    const static double meo_coef[3][3]={
            {-0.946, 2.158, -0.642},
            { 0.598, 1.635, -0.556},
            {-0.177, 0.652, -0.178},
    };

    for(i=0;i<n;i++){
        sat=obs[i].sat;
        sys=satsys(sat,&prn);
        el=rtk->ssat[sat-1].azel[1];
        if(!(sys&SYS_CMP)) continue;
        if(prn>=18) continue;

        if(prn<=5){  /*GEO C01 C02 C03 C04 C05*/
            /*B1I*/
            if(obs[i].P[0]!=0.0){
                mp=geo_coef[0][0]*el+geo_coef[0][1]*SQR(el)+geo_coef[0][2]*SQR(el)*el;
                obs[i].P[0]+=mp;
            }
            /*B2I*/
            if(obs[i].P[1]!=0.0){
                mp=geo_coef[1][0]*el+geo_coef[1][1]*SQR(el)+geo_coef[1][2]*SQR(el)*el;
                obs[i].P[1]+=mp;
            }
            /*B3I*/
            if(obs[i].P[3]!=0.0){
                mp=geo_coef[2][0]*el+geo_coef[2][1]*SQR(el)+geo_coef[2][2]*SQR(el)*el;
                obs[i].P[3]+=mp;
            }
        }
        else if((prn>=6&&prn<=10)||(prn==13||prn==16)){ /*IGSO C06 C07 C08 C09 C10 C13 C14*/
            /*B1I*/
            if(obs[i].P[0]!=0.0){
                mp=igso_coef[0][0]*el+igso_coef[0][1]*SQR(el)+igso_coef[0][2]*SQR(el)*el;
                obs[i].P[0]+=mp;
            }
            /*B2I*/
            if(obs[i].P[1]!=0.0){
                mp=igso_coef[1][0]*el+igso_coef[1][1]*SQR(el)+igso_coef[1][2]*SQR(el)*el;
                obs[i].P[1]+=mp;
            }
            /*B3I*/
            if(obs[i].P[3]!=0.0){
                mp=igso_coef[2][0]*el+igso_coef[2][1]*SQR(el)+igso_coef[2][2]*SQR(el)*el;
                obs[i].P[3]+=mp;
            }
        }
        else if(prn==11||prn==12||prn==14){ /*MEO C11 C12 C14*/
            /*B1I*/
            if(obs[i].P[0]!=0.0){
                mp=meo_coef[0][0]*el+meo_coef[0][1]*SQR(el)+meo_coef[0][2]*SQR(el)*el;
                obs[i].P[0]+=mp;
            }
            /*B2I*/
            if(obs[i].P[1]!=0.0){
                mp=meo_coef[1][0]*el+meo_coef[1][1]*SQR(el)+meo_coef[1][2]*SQR(el)*el;
                obs[i].P[1]+=mp;
            }
            /*B3I*/
            if(obs[i].P[3]!=0.0){
                mp=meo_coef[2][0]*el+meo_coef[2][1]*SQR(el)+meo_coef[2][2]*SQR(el)*el;
                obs[i].P[3]+=mp;
            }
        }
    }
#endif


}

/* satellite antenna phase center variation ----------------------------------*/
static void satantpcv(const double *rs, const double *rr, const pcv_t *pcv,
                      double *dant)
{
    double ru[3],rz[3],eu[3],ez[3],nadir,cosa;
    int i;
    
    for (i=0;i<3;i++) {
        ru[i]=rr[i]-rs[i];
        rz[i]=-rs[i];
    }
    if (!normv3(ru,eu)||!normv3(rz,ez)) return;
    
    cosa=dot(eu,ez,3);
    cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
    nadir=acos(cosa);
    
    antmodel_s(pcv,nadir,dant);
}
/* ionospheric model ---------------------------------------------------------*/
extern int model_iono(gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, int sat, const double *x,
                      const nav_t *nav, double *dion, double *var)
{
    static double iono_p[MAXSAT]={0},std_p[MAXSAT]={0};
    static gtime_t time_p;
    
    if (opt->ionoopt==IONOOPT_SBAS) {
        return sbsioncorr(time,nav,pos,azel,dion,var);
    }
    if (opt->ionoopt==IONOOPT_TEC) {
        return iontec(time,nav,pos,azel,1,dion,var);
    }
    if (opt->ionoopt==IONOOPT_BRDC) {
        *dion=klobuchar_GPS(time,nav->ion_gps,pos,azel);
        *var=SQR(*dion*ERR_BRDCI);
        return 1;
    }
    if (opt->ionoopt==IONOOPT_UC||opt->ionoopt==IONOOPT_UC_CONS) {
        int i,tc;
        tc=opt->mode>=PMODE_TC_SPP&&opt->mode<=PMODE_TC_PPP;
        i=II(sat,opt);
        *dion=x[i];
        *var=0.0;
        return 1;
    }
    if (opt->ionoopt==IONOOPT_IFLC||opt->ionoopt==IONOOPT_IF2) {
        *dion=*var=0.0;
        return 1;
    }
    if (opt->ionoopt==IONOOPT_STEC) {
        if (timediff(time,time_p)!=0.0&&
            !pppcorr_stec(&nav->pppcorr,time,pos,iono_p,std_p)) return 0;
        if (iono_p[sat-1]==0.0||std_p[sat-1]>0.1) return 0;
        time_p=time;
        *dion=iono_p[sat-1];
        *var=SQR(std_p[sat-1]);
        return 1;
    }
    return 0;
}
/* precise tropospheric model ------------------------------------------------*/
static double trop_model_prec(const prcopt_t *opt,gtime_t time, const double *pos,
                              const double *azel, const double *x, double *dtdx,
                              double *var,double ztrp[2],double mtrp[2])
{
    const double zazel[]={0.0,PI/2.0};
    double zhd,m_h,m_w,cotz,grad_n,grad_e;

    /* zenith hydrostatic delay */
    zhd=saastamoinen(time,pos,zazel,0.0,1,NULL,NULL);

    /* mapping function */
    m_h=tropmapf(opt,time,pos,azel,&m_w);
    if(mtrp){
        mtrp[0]=m_h;
        mtrp[1]=m_w;
    }

    if (opt->tropopt==TROPOPT_ESTG&&azel[1]>0.0) {
        /* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
        cotz=1.0/tan(azel[1]);
        grad_n=m_w*cotz*cos(azel[0]);
        grad_e=m_w*cotz*sin(azel[0]);
        m_w+=grad_n*x[1]+grad_e*x[2];
        dtdx[1]=grad_n*(x[0]-zhd);
        dtdx[2]=grad_e*(x[0]-zhd);
    }
    dtdx[0]=m_w;
    *var=SQR(0.01);
    if(ztrp){
        ztrp[0]=zhd;
        ztrp[1]=x[0];
    }
    return m_h*zhd+m_w*(x[0]);
}
/* tropospheric model ---------------------------------------------------------*/
extern int model_trop(gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, const double *x, double *dtdx,
                      const nav_t *nav, double *dtrp,double *ztrp,double *mtrp, double *var,int it)
{
    double trp[3]={0},std[3];

    if (opt->tropopt==TROPOPT_SAAS) {
        *dtrp=saastamoinen(time,pos,azel,0.7,0,NULL,NULL);
        *var=SQR(0.3);
        return 1;
    }
    if (opt->tropopt==TROPOPT_SBAS) {
        *dtrp=sbstropcorr(time,pos,azel,var);
        return 1;
    }
    if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
        matcpy(trp,x+it,opt->tropopt==TROPOPT_EST?1:3,1);
        *dtrp=trop_model_prec(opt,time,pos,azel,trp,dtdx,var,ztrp,mtrp);
        return 1;
    }
    if (opt->tropopt==TROPOPT_ZTD) {
        if (pppcorr_trop(&nav->pppcorr,time,pos,trp,std)) {
            *dtrp=trop_model_prec(opt,time,pos,azel,trp,dtdx,var,ztrp,mtrp);
            *var=SQR(dtdx[0]*std[0]);
            return 1;
        }
    }
    return 0;
}

static double shapiro_corr(int sys,const double *rs,double *rr)
{
    double drr,drs,r,mu;

    drr=sqrt(rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2]);
    drs=sqrt(rs[0]*rs[0]+rs[1]*rs[1]+rs[2]*rs[2]);
    r=sqrt(SQR(rs[0]-rr[0])+SQR(rs[1]-rr[1])+SQR(rs[2]-rr[2]));

    switch(sys){
        case SYS_CMP:
            mu=3.986004418E14;break;
        case SYS_GLO:
            mu=3.9860044E14;break;
        case SYS_GAL:
            mu=3.986004418E14;break;
        default:
            mu=3.9860050E14;break;
    }

    return 2.0*mu/(CLIGHT*CLIGHT)*log((drs+drr+r)/(drs+drr-r));
}

/* constraint to local correction --------------------------------------------*/
static int const_corr(const obsd_t *obs, int n, const int *exc,
                      const nav_t *nav, const double *x, const double *pos,
                      const double *azel, rtk_t *rtk, double *v, double *H,
                      double *var)
{
    gtime_t time=obs[0].time;
    double trop[3],std_trop[3],iono[MAXSAT],std_iono[MAXSAT];
    int i,j,k,sat,nv=0;
    
    /* constraint to external troposphere correction */
    if ((rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG)&&
        pppcorr_trop(&nav->pppcorr,time,pos,trop,std_trop)) {
        
        for (i=0;i<(rtk->opt.tropopt==TROPOPT_EST?1:3);i++) {
            if (std_trop[i]==0.0) continue;
            j=IT(&rtk->opt)+i;
            v[nv]=trop[i]-x[j];
            for (k=0;k<rtk->nx;k++) H[k+nv*rtk->nx]=k==j?1.0:0.0;
            var[nv++]=SQR(std_trop[i]);
        }
    }
    /* constraint to external ionosphere correction */
    if (rtk->opt.ionoopt==IONOOPT_UC&&
        pppcorr_stec(&nav->pppcorr,time,pos,iono,std_iono)) {
        
        for (i=0;i<n;i++) {
            sat=obs[i].sat;
            if (exc[i]||iono[sat-1]==0.0||std_iono[sat-1]>0.5) continue;
            j=II(sat,&rtk->opt);
            v[nv]=iono[sat-1]-x[j];
            for (k=0;k<rtk->nx;k++) H[k+nv*rtk->nx]=k==j?1.0:0.0;
            var[nv++]=SQR(std_iono[sat-1]);
        }
    }
    return nv;
}

static void tracesatinfo(const prcopt_t *opt,int level,int f,const int *frq_idxs,int type,const rtk_t *rtk,const obsd_t *obs,double meas,double range,double dtr,double isb,double rdcb,double rifcb,
        double trp,double ion,double antr,double shapiro,const double *rs,double dts,double sdcb,double ants,double phw,double amb,double res,double sig)
{
    static int last_sat_no=0;
    char buff[1024]={'\0'};
    const ssat_t *sat_info=&rtk->ssat[obs->sat-1];
    int slip_flag;
    int sys=satsys(obs->sat,NULL);

//    for(int i=0;i<opt->nf;i++){
//        if(sys==SYS_G)
//    }


    if(last_sat_no==0){ /*first line*/
        last_sat_no=obs->sat;
        sprintf(buff,"%s lock=%4d outc=%4d slip=%d el=%3.1f dtr=%10.3f isb=%6.3f trp=%6.3f shapiro=%6.3f dts=%12.3f sx=%12.3f sy=%12.3f sz=%12.3f\n",
                obs->id,sat_info->lock[f],sat_info->outc[f],sat_info->slip[0]&3||sat_info->slip[1]&2||sat_info->slip[2]&3,sat_info->azel[1]*R2D,dtr,isb,trp,shapiro,dts*CLIGHT,rs[0],rs[1],rs[2]);
        trace(level,buff);
        /*obs info*/
        char *c1=NULL,*c2=NULL,*c3=NULL,type_str[20]={'\0'};

        sprintf(type_str,"%s",type==0?"L":"P");

        if(rtk->opt.ionoopt==IONOOPT_IFLC){
            if(rtk->opt.nf==1){ /*SF*/
                c1=code2obs(obs->code[frq_idxs[f]-1]);
                sprintf(type_str,"%s%s+%s%s","P",c1,"L",c1);
            }
            if(rtk->opt.nf>=2){ /*DF*/
                c1=code2obs(obs->code[frq_idxs[f]-1]);
                c2=code2obs(obs->code[frq_idxs[f+1]-1]);
                sprintf(type_str,"%s%s+%s",type_str,c1,c2);
            }
            if(rtk->opt.nf>=3){ /*TF*/
                c1=code2obs(obs->code[frq_idxs[f]-1]);
                c2=code2obs(obs->code[frq_idxs[f+1]-1]);
                c3=code2obs(obs->code[frq_idxs[f+2]-1]);
                sprintf(type_str,"%s%s+%s+%s",type_str,c1,c2,c3);
            }
        }
        else if(rtk->opt.ionoopt==IONOOPT_IF2){ /*TF with two dual-frequency IF*/
            if(f==0){ /*L1L2*/
                c1=code2obs(obs->code[frq_idxs[f]-1]);
                c2=code2obs(obs->code[frq_idxs[f+1]-1]);
                sprintf(type_str,"%s%s+%s",type_str,c1,c2);
            }
            else if(f==1){ /*L1L3*/
                c1=code2obs(obs->code[frq_idxs[f]-1]);
                c2=code2obs(obs->code[frq_idxs[f+2]-1]);
                sprintf(type_str,"%s%s+%s",type_str,c1,c2);
            }
        }
        else{ /*UC*/
            c1=code2obs(obs->code[frq_idxs[f]-1]);
            sprintf(type_str,"%s%s",type_str,c1);
        }
        buff[0]='\0';
        if(type==1){ /*P*/
            sprintf(buff,"%s res=%9.3f sig=%9.6f obs=%12.3f r=%12.3f cbias=%6.3f ion=%6.3f dantr=%6.3f dants=%6.3f rdcb=%6.3f rifcb=%6.3f\n",
                    type_str,res,sig,meas,range,sdcb,ion,antr,ants,rdcb,rifcb);
        }
        else{ /*L*/
            sprintf(buff,"%s res=%9.3f sig=%9.6f obs=%12.3f r=%12.3f ambig=%6.3f ion=%6.3f dantr=%6.3f dants=%6.3f phw=%6.3f\n",
                    type_str,res,sig,meas,range,amb,ion,antr,ants,phw);
        }
        trace(level,buff);
    }
    else if(last_sat_no!=0&&last_sat_no==obs->sat){
        /*obs info*/
        char *c1=NULL,*c2=NULL,*c3=NULL,type_str[20]={'\0'};

        sprintf(type_str,"%s",type==0?"L":"P");

        if(rtk->opt.ionoopt==IONOOPT_IFLC){
            if(rtk->opt.nf==1){ /*SF*/
                c1=code2obs(obs->code[frq_idxs[f]-1]);
                sprintf(type_str,"%s%s+%s%s","P",c1,"L",c1);
            }
            if(rtk->opt.nf>=2){ /*DF*/
                c1=code2obs(obs->code[frq_idxs[f]-1]);
                c2=code2obs(obs->code[frq_idxs[f+1]-1]);
                sprintf(type_str,"%s%s+%s",type_str,c1,c2);
            }
            if(rtk->opt.nf>=3){ /*TF*/
                c1=code2obs(obs->code[frq_idxs[f]-1]);
                c2=code2obs(obs->code[frq_idxs[f+1]-1]);
                c3=code2obs(obs->code[frq_idxs[f+2]-1]);
                sprintf(type_str,"%s%s+%s+%s",type_str,c1,c2,c3);
            }
        }
        else if(rtk->opt.ionoopt==IONOOPT_IF2){ /*TF with two dual-frequency IF*/
            if(f==0){ /*L1L2*/
                c1=code2obs(obs->code[frq_idxs[f]-1]);
                c2=code2obs(obs->code[frq_idxs[f+1]-1]);
                sprintf(type_str,"%s%s+%s",type_str,c1,c2);
            }
            else if(f==1){ /*L1L3*/
                c1=code2obs(obs->code[frq_idxs[f]-1]);
                c2=code2obs(obs->code[frq_idxs[f+2]-1]);
                sprintf(type_str,"%s%s+%s",type_str,c1,c2);
            }
        }
        else{ /*UC*/
            c1=code2obs(obs->code[frq_idxs[f]-1]);
            sprintf(type_str,"%s%s",type_str,c1);
        }

        buff[0]='\0';
        if(type==1){ /*P*/
            sprintf(buff,"%s res=%9.3f sig=%9.6f obs=%12.3f r=%12.3f cbias=%6.3f ion=%6.3f dantr=%6.3f dants=%6.3f rdcb=%6.3f rifcb=%6.3f\n",
                    type_str,res,sig,meas,range,sdcb,ion,antr,ants,rdcb,rifcb);
        }
        else{ /*L*/
            sprintf(buff,"%s res=%9.3f sig=%9.6f obs=%12.3f r=%12.3f ambig=%6.3f ion=%6.3f dantr=%6.3f dants=%6.3f phw=%6.3f\n",
                    type_str,res,sig,meas,range,amb,ion,antr,ants,phw);
        }
        trace(level,buff);
    }
    else if(last_sat_no!=0&&last_sat_no!=obs->sat){
        last_sat_no=obs->sat;
        sprintf(buff,"%s lock=%4d outc=%4d slip=%d el=%3.1f dtr=%10.3f isb=%6.3f rdcb=%6.3f rifcb=%6.3f trp=%6.3f dts=%12.3f sx=%12.3f sy=%12.3f sz=%12.3f\n",
                obs->id,sat_info->lock[f],sat_info->outc[f],sat_info->slip[0]&3||sat_info->slip[1]&2||sat_info->slip[2]&3,sat_info->azel[1]*R2D,dtr,isb,rdcb,rifcb,trp,dts*CLIGHT,rs[0],rs[1],rs[2]);
        trace(level,buff);
        /*obs info*/
        char *c1=NULL,*c2=NULL,*c3=NULL,type_str[20]={'\0'};

        sprintf(type_str,"%s",type==0?"L":"P");

        if(rtk->opt.ionoopt==IONOOPT_IFLC){
            if(rtk->opt.nf==1){ /*SF*/
                c1=code2obs(obs->code[frq_idxs[f]-1]);
                sprintf(type_str,"%s%s+%s%s","P",c1,"L",c1);
            }
            if(rtk->opt.nf>=2){ /*DF*/
                c1=code2obs(obs->code[frq_idxs[f]-1]);
                c2=code2obs(obs->code[frq_idxs[f+1]-1]);
                sprintf(type_str,"%s%s+%s",type_str,c1,c2);
            }
            if(rtk->opt.nf>=3){ /*TF*/
                c1=code2obs(obs->code[frq_idxs[f]-1]);
                c2=code2obs(obs->code[frq_idxs[f+1]-1]);
                c3=code2obs(obs->code[frq_idxs[f+2]-1]);
                sprintf(type_str,"%s%s+%s+%s",type_str,c1,c2,c3);
            }
        }
        else if(rtk->opt.ionoopt==IONOOPT_IF2){ /*TF with two dual-frequency IF*/
            if(f==0){ /*L1L2*/
                c1=code2obs(obs->code[frq_idxs[f]-1]);
                c2=code2obs(obs->code[frq_idxs[f+1]-1]);
                sprintf(type_str,"%s%s+%s",type_str,c1,c2);
            }
            else if(f==1){ /*L1L3*/
                c1=code2obs(obs->code[frq_idxs[f]-1]);
                c2=code2obs(obs->code[frq_idxs[f+2]-1]);
                sprintf(type_str,"%s%s+%s",type_str,c1,c2);
            }
        }
        else{ /*UC*/
            c1=code2obs(obs->code[frq_idxs[f]-1]);
            sprintf(type_str,"%s%s",type_str,c1);
        }
        buff[0]='\0';
        if(type==1){ /*P*/
            sprintf(buff,"%s res=%9.3f sig=%9.6f obs=%12.3f r=%12.3f cbias=%6.3f ion=%6.3f dantr=%6.3f dants=%6.3f rdcb=%6.3f rifcb=%6.3f\n",
                    type_str,res,sig,meas,range,sdcb,ion,antr,ants,rdcb,rifcb);
        }
        else{ /*L*/
            sprintf(buff,"%s res=%9.3f sig=%9.6f obs=%12.3f r=%12.3f ambig=%6.3f ion=%6.3f dantr=%6.3f dants=%6.3f phw=%6.3f\n",
                    type_str,res,sig,meas,range,amb,ion,antr,ants,phw);
        }
        trace(level,buff);
    }
}

static int zdres(int post, const obsd_t *obs, int n, const double *rs,
                 const double *dts, const double *var_rs, const int *svh,
                 const double *dr, int *exc, const nav_t *nav,
                 const double *x, rtk_t *rtk, double *v, double *los, double *mw,
                 double *gamma,double *azel,double *rpos,double *var_sat,int *vflg)
{
    prcopt_t *opt=&rtk->opt;
    double y,r,bias,C=1.0,rr[3],pos[3],dtdx[3],L[NFREQ]={0},P[NFREQ]={0},Lc[NFREQ]={0},Pc[NFREQ]={0},freqs[NFREQ]={0},freq_base=0.0;
    double var[MAXOBS*2],dtrp=0.0,dion=0.0,vart=0.0,vari=0.0;
    double dantr[NFREQ]={0},dants[NFREQ]={0};
    double ve[MAXOBS*2*NFREQ]={0},vare[MAXOBS*2*NFREQ]={0},cbias[NFREQ]={0},shapiro=0,isb=0,rdcb=0.0,rifcb=0.0;
    double tec_fact=1.0,freq;
    char str[32];
    int ne=0,obsi[MAXOBS*2*NFREQ]={0},frqi[MAXOBS*2*NFREQ];
    int i,j,k,sat,sys,sys_idx=-1,nv=0,it=0,prn;
    int iamb=0,iion=0;
    int *frq_idxs=NULL,level=3;
    int vaild_ns=0;

    time2str(obs[0].time,str,2);

    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].vsat[j]=0;

    for (i=0;i<3;i++) rr[i]=rpos[i]+dr[i];
    ecef2pos(rr,pos);

    if(rtk->epoch>=6113){
        level=20;
    }
    trace(level,"%s(%5d) %s residual: rr=%12.3f %12.3f %12.3f dr=%6.3f %6.3f %6.3f\n",time_str(obs->time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,post?"post":"prior",rr[0],rr[1],rr[2],dr[0],dr[1],dr[2]);

    for (i=0;i<n&&i<MAXOBS;i++) {
        sat=obs[i].sat;
        sys=satsys(sat,&prn);
        sys_idx=satsysidx(sat);
        if(sys_idx==-1) continue;

        frq_idxs=rtk->opt.gnss_frq_idx[sys_idx];

        if ((r=geodist(rs+i*6,rr,los+3*i))<=0.0||
            satazel(pos,los+3*i,azel+i*2)<opt->elmin) {
            exc[i]=1;
            continue;
        }
        if (!(sys=satsys(sat,&prn))||!rtk->ssat[sat-1].vs||
            satexclude(obs[i].sat,var_rs[i],svh[i],opt)||exc[sat-1]) {
            exc[i]=1;
            continue;
        }
        if(!snrmask(&obs[i],&rtk->ssat[sat-1].azel[1],&rtk->opt)){
            exc[i]=1;
            continue;
        }

        /* tropospheric and ionospheric model */
        it=IT(opt);
        if (!model_trop(obs[i].time,pos,azel+i*2,opt,x,dtdx,nav,&dtrp,rtk->sol.ztrp,rtk->ssat[sat-1].mtrp,&vart,it)||
            !model_iono(obs[i].time,pos,azel+i*2,opt,sat,x,nav,&dion,&vari)) {
            continue;
        }

        /* satellite and receiver antenna model */
        if (opt->posopt[0]) satantpcv(rs+i*6,rr,nav->pcvs+sat-1,dants);
        antmodel(obs[i].sat,opt->pcvr,opt->antdel[0],azel+i*2,opt->posopt[1],dantr);

        /* phase windup model */
        if (!model_phw(rtk->sol.time,sat,nav->pcvs[sat-1].type,
                       opt->posopt[2]?2:0,rs+i*6,rr,&rtk->ssat[sat-1].phw)) {
            continue;
        }

        shapiro=shapiro_corr(satsys(sat,NULL),rs+i*6,rr);

        getcorrobs(&rtk->opt,obs+i,nav,frq_idxs,dantr,dants,rtk->ssat[sat-1].phw,L,P,Lc,Pc,freqs,cbias,&rtk->ssat[sat-1]);

        tec_fact=40.30E16/freq_base/freq_base;

        /* stack phase and code residuals {L1,P1,L2,P2,...} */
        for (j=0;j<2*NF(opt);j++) {

            bias=0.0;
            if (opt->ionoopt==IONOOPT_IFLC||opt->ionoopt==IONOOPT_IF2) {
                if ((y=j%2==0?Lc[j/2]:Pc[j/2])==0.0) continue;
            }
            else {
                if ((y=j%2==0?L[j/2]:P[j/2])==0.0) continue;
            }

            /*troposphere delay*/
            if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                for (k=0;k<(opt->tropopt>=TROPOPT_ESTG?3:1);k++) {
                    mw[k+i*(opt->tropopt>=TROPOPT_ESTG?3:1)]=dtdx[k];
                }
            }

            /*ionospheric delay*/
            if (opt->ionoopt==IONOOPT_UC||opt->ionoopt==IONOOPT_UC_CONS) {
                if ((freq=freqs[frq_idxs[j/2]-1])==0.0) continue;
                C=SQR(freq_base/freq)*(j%2==0?-1.0:1.0);

                iion=II(sat,opt);
                if (rtk->x[iion]==0.0) continue;
                gamma[j+i*NF(opt)*2]=C;
                if(post) rtk->ssat[sat-1].tec=rtk->x[iion]/tec_fact;
            }

            /*phase bias*/
            if (j%2==0) {
                iamb=IB(sat,j/2,opt);
                if ((bias=x[iamb])==0.0) continue;

                if(post){
                    if(rtk->opt.ionoopt==IONOOPT_IFLC||rtk->opt.ionoopt==IONOOPT_IF2){
                        rtk->ssat[sat-1].amb[j/2]=x[iamb];
                    }
                    else{
                        rtk->ssat[sat-1].amb[j/2]=x[iamb]*(CLIGHT/freqs[j/2]);
                    }
                }
            }

            /* residual */
            double a=y-(r-CLIGHT*dts[i*2]+dtrp+C*dion+bias-shapiro);
            v[j+i*NF(opt)*2]=y-(r-CLIGHT*dts[i*2]+dtrp+C*dion+bias-shapiro);

            if (j%2==0) rtk->ssat[sat-1].resc[j/2]=v[nv];
            else        rtk->ssat[sat-1].resp[j/2]=v[nv];

            /* variance */
            var_sat[j+i*NF(opt)*2]=varerr(obs[i].sat,sys,azel[1+i*2],0.25*rtk->ssat[sat-1].snr_rover[j/2],j/2,(opt->nf==1&&opt->ionoopt==IONOOPT_IFLC)?1:j%2,opt)+
                                   vart+SQR(C)*vari+var_rs[i];
            var_sat[j+i*NF(opt)*2]*=rtk->ssat[sat-1].var_fact[j%2][j/2];

            /* down weight glonass ifb */
            if (sys==SYS_GLO&&j%2==1) var_sat[nv]+=VAR_GLO_IFB;
            /* down weight eclipse satellite */
            var_sat[j+i*NF(opt)*2]*=rtk->ssat[sat-1].eclipse;

            tracesatinfo(opt,level,j/2,frq_idxs,j%2,rtk,obs+i,y,r,0.0,isb,rdcb,rifcb,dtrp,dion,dantr[frq_idxs[j/2]-1],shapiro,
                         rs+6*i,dts[i*2],cbias[j/2],dants[frq_idxs[j/2]-1],rtk->ssat[sat-1].phw,bias,v[j+i*NF(opt)*2],sqrt( var_sat[j+i*NF(opt)*2]));

            /* reject satellite by pre-fit residuals */
            if (!post&&opt->maxinno>0.0&&fabs(v[j+i*NF(opt)*2])>opt->maxinno) {
                trace(2,"outlier (%d) rejected %s %s %s%d res=%9.4f el=%4.1f\n",
                      post,str,sat_id(sat),j%2?"P":"L",j/2+1,v[nv],azel[1+i*2]*R2D);
                exc[i]=1; rtk->ssat[sat-1].rejc[j%2]++;
                continue;
            }
            /* record large post-fit residuals */
            if (post&&fabs(v[j+i*NF(opt)*2])>sqrt(var_sat[j+i*NF(opt)*2])*THRES_REJECT) {
                obsi[ne]=i; frqi[ne]=j; ve[ne]=v[j+i*NF(opt)*2];vare[ne]=var_sat[j+i*NF(opt)*2]; ne++;
            }
            if (j%2==0){
                vaild_ns++;
                rtk->ssat[sat-1].vsat[j/2]=1;
            }
            nv++;
        }

//        if(opt->ionoopt==IONOOPT_UC_CONS){
//            double var_tec=1.0;
//            iontec(obs[i].time,nav,pos,azel+i*2,1,&tec_ion,&var_tec);
//            var_tec=iontecvar(rtk->epoch,obs[0].time,pos,rtk->ssat[sat-1].azel);
//            iion=II(sat,opt);
//            v[nv]=tec_ion-x[iion];
//            for(int kk=0;kk<rtk->nx;kk++) H[kk+rtk->nx*nv]=kk==iion?1.0:0.0;
//            var[nv]=var_tec;
//            trace(3,"ION res=%9.3f sig=%9.6f tec_ion=%12.3f ion=%12.3f\n",v[nv],sqrt(var[nv]),tec_ion,x[iion]);
//            nv++;
//        }
    }

    return vaild_ns;
}

/* phase and code residuals --------------------------------------------------*/
static int ppp_res(int post, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *var_rs, const int *svh,
                   const double *dr, int *exc, const nav_t *nav,
                   const double *x, rtk_t *rtk, double *v, double *H, double *R,
                   double *azel,double *rpos,int *v_flag,int *valid_ns)
{
    prcopt_t *opt=&rtk->opt;
    int nv=0;

    if(opt->sdopt){
        int i,j,k,m,f,sysi,nb[5]={0},b=0,nf=NF(&rtk->opt),frq,code,nx=rtk->nx,ntrp=opt->tropopt>=TROPOPT_ESTG?3:1;
        double *y,*var_sat,*e,*mw,*gamma,*lam,*Ri,*Rj;
        y=zeros(nf*2,n);var_sat=zeros(nf*2,n);e=zeros(3,n);mw=zeros(1*ntrp,n);lam=zeros(nf,n);gamma=zeros(nf*2,n);
        Ri=zeros(n,nf*2);Rj=zeros(n,nf*2);


        *valid_ns=zdres(post,obs,n,rs,dts,var_rs,svh,dr,exc,nav,x,rtk,y,e,mw,gamma,azel,rpos,var_sat,v_flag);

        for(m=0;m<5;m++){

            /*find reference satellite*/
            for(i=-1,j=0;j<n;j++){
                sysi=satsys(obs[j].sat,NULL);
                if(!test_sys(sysi,m)) continue;
                if(y[j*2]==0.0) continue;
                if(i<0||azel[1+j*2]>=azel[1+i*2]) i=j;
            }
            if(i<0) continue;

            for(f=0;f<nf*2;f++){
                code=f%2==0?0:1;frq=f/2;

                for(j=0;j<n;j++){
                    if(j==i) continue;
                    sysi=satsys(obs[j].sat,NULL);
                    if(!test_sys(sysi,m)) continue;
                    if(azel[j*2+1]==0.0) continue;
                    if(y[f+j*nf*2]==0.0||y[f+i*nf*2]==0.0) continue;

                    double a=y[f+i*nf*2];
                    double c=y[f+j*nf*2];
                    double d=a-c;
                    if(fabs(d)>10.0) continue;
                    v[nv]=y[f+i*nf*2]-y[f+j*nf*2];
                    Ri[nv]=var_sat[f+i*nf*2];
                    Rj[nv]=var_sat[f+j*nf*2];

                    if(rtk->tc){
                        /*position */
                        for(k=0;k<rtk->nx;k++) H[k+nv*nx]=k<3?(e[k+3*i]-e[k+3*j]):0.0;
                        /*attitude*/
                        v3_t dpda={0};
//                        jacobi_tc_meas2a_ECEF(&dpda,&rtk->ins_kf->insstate->dcm,&rtk->opt.insopt.imup.lever_arm_gps,&e[3*i],&e[3*j]);
                    }
                    else{
                        for(k=0;k<rtk->nx;k++) H[k+nv*nx]=k<3?(-e[k+3*i]+e[k+3*j]):0.0;
                    }

                    /*tropospheric delay*/
                    int itrp=IT(&rtk->opt);
                    for(k=itrp;k<itrp+ntrp;k++){
                        H[k+nv*nx]=mw[k-itrp+i*ntrp]-mw[k-itrp+j*ntrp];
                    }

                    /*ionospheric delay*/
                    if(opt->ionoopt==IONOOPT_UC||opt->ionoopt==IONOOPT_UC_CONS){
                        int iion=II(obs[i].sat,&rtk->opt);
                        int jion=II(obs[j].sat,&rtk->opt);
                        H[iion+nv*nx]=gamma[f+i*nf*2];
                        H[jion+nv*nx]=gamma[f+j*nf*2];
                    }

                    /*Ambiguity*/
                    if(!code){
                        int iamb=0,jamb=0;
                        iamb=IB(obs[i].sat,frq,&rtk->opt);
                        jamb=IB(obs[j].sat,frq,&rtk->opt);
                        H[iamb+nv*nx]=1.0;
                        H[jamb+nv*nx]=-1.0;
                        if(post){
                            rtk->sdamb[obs[j].sat-1].lc=x[jamb];
                            rtk->sdamb[obs[i].sat-1].lc=x[iamb];
                            rtk->sdamb[obs[j].sat-1].ref_sat_no=obs[i].sat;
                            rtk->sdamb[obs[i].sat-1].ref_sat_no=0;
                        }
                    }

                    v_flag[nv]=(obs[i].sat<<16)|(obs[j].sat<<8)|(code?1:0)<<4|frq;
                    nv++;
                    nb[b]++;
                }
                b++;
            }
        }
        ddcov(nb,b,Ri,Rj,nv,R);
        free(y);free(var_sat);free(e);free(mw);
        free(gamma);free(lam);free(Ri);free(Rj);
        return nv;
    }
    else{
        int level=3;
        double y,r,cdtr,bias,C=1.0,rr[3],pos[3],e[3],dtdx[3],L[NFREQ]={0},P[NFREQ]={0},Lc[NFREQ]={0},Pc[NFREQ]={0},freqs[NFREQ]={0},freq_base=0.0,freq_base2=0.0,freq=0.0;
        double var[MAXOBS*2],dtrp=0.0,dion=0.0,vart=0.0,vari=0.0,tec_ion=0.0;
        double dantr[NFREQ]={0},dants[NFREQ]={0},ztrp[2]={0};
        double ve[MAXOBS*2*NFREQ]={0},vare[MAXOBS*2*NFREQ]={0},vmax=0,varmax=0,cbias[NFREQ]={0},shapiro=0,isb=0,rdcb=0.0,rifcb=0.0;
        double alpha,beta,tec_fact=1.0;
        char str[32];
        int ne=0,obsi[MAXOBS*2*NFREQ]={0},frqi[MAXOBS*2*NFREQ],maxobs,maxfrq,rej;
        int i,j,k,sat,sys,sys_idx=-1,nx=rtk->nx,stat=1,it=0,prn;
        int isys=0,iamb=0,itrp=0,iion=0,main_irc=0,irc=0,mask[6]={0},bd3_flag=0;
        int *frq_idxs=NULL,vs=0;

        time2str(obs[0].time,str,2);

        for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].vsat[j]=0;

        for (i=0;i<3;i++) rr[i]=rpos[i]+dr[i];
        ecef2pos(rr,pos);

        if(rtk->stc&&rtk->ins_kf->couple_epoch>=10000){
            level=1;
        }
        trace(level,"%5d %s residual(%s): rr=%12.3f %12.3f %12.3f dr=%6.3f %6.3f %6.3f\n",
                rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,post?"post":"prior",time_str(obs[0].time,1),rr[0],rr[1],rr[2],dr[0],dr[1],dr[2]);
        for (i=0;i<n&&i<MAXOBS;i++) {
            sat=obs[i].sat;
            sys=satsys(sat,&prn);
            sys_idx=satsysidx(sat);
            if(sys_idx==-1) continue;

            if ((r=geodist(rs+i*6,rr,e))<=0.0||
                satazel(pos,e,azel+i*2)<opt->elmin) {
                exc[sat-1]=1;
                continue;
            }
            if (!(sys=satsys(sat,&prn))||!rtk->ssat[sat-1].vs||
                satexclude(obs[i].sat,var_rs[i],svh[i],opt)||exc[sat-1]) {
                exc[sat-1]=1;
                continue;
            }

            /* tropospheric and ionospheric model */
            it=IT(opt);
            if (!model_trop(obs[i].time,pos,azel+i*2,opt,x,dtdx,nav,&dtrp,rtk->sol.ztrp,rtk->ssat[sat-1].mtrp,&vart,it)||
                !model_iono(obs[i].time,pos,azel+i*2,opt,sat,x,nav,&dion,&vari)) {
                rtk->ssat[sat-1].vs=0;
                continue;
            }

            /* satellite and receiver antenna model */
            if (opt->posopt[0]) satantpcv(rs+i*6,rr,nav->pcvs+sat-1,dants);
            antmodel(obs[i].sat,opt->pcvr,opt->antdel[0],azel+i*2,opt->posopt[1],dantr);

            /* phase windup model */
            if (!model_phw(rtk->sol.time,sat,nav->pcvs[sat-1].type,
                           opt->posopt[2]?2:0,rs+i*6,rr,&rtk->ssat[sat-1].phw)) {
                continue;
            }

            shapiro=shapiro_corr(satsys(sat,NULL),rs+i*6,rr);

            getcorrobs(&rtk->opt,obs+i,nav,rtk->opt.gnss_frq_idx[sys_idx],dantr,dants,rtk->ssat[sat-1].phw,L,P,Lc,Pc,freqs,cbias,&rtk->ssat[sat-1]);
            frq_idxs=rtk->opt.gnss_frq_idx[sys_idx];

            tec_fact=40.30E16/freq_base/freq_base;
#if 0
            /*check new sat*/
            if(rtk->ssat[sat-1].lock[0]<0&&rtk->ssat[sat-1].outc[0]>1){
                rtk->ssat[sat-1].new_sat=1;
                if(!post){
                    trace(3,"%s(%d):%s is new satellite lock=%4d outc=%4d\n",time_str(obs->time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,
                          sat_id(sat),rtk->ssat[sat-1].lock[0],rtk->ssat[sat-1].outc[0]);
                }
            }
#endif
            /* stack phase and code residuals {L1,P1,L2,P2,...} */
            for (j=0;j<2*NF(opt);j++) {

                bias=0.0;
                if (opt->ionoopt==IONOOPT_IFLC||opt->ionoopt==IONOOPT_IF2) {
                    if ((y=j%2==0?Lc[j/2]:Pc[j/2])==0.0) continue;
                }
                else {
                    if ((y=j%2==0?L[j/2]:P[j/2])==0.0) continue;

                    /* receiver DCB correction for P2 */
//                if (j/2==1) dcb=-nav->rbias[0][sys==SYS_GLO?1:0][0];
                }

                for (k=0;k<nx;k++) H[k+nx*nv]=k<3?-e[k]:0.0;

                /* receiver clock*/
                int main_clk;
                for(main_clk=0;main_clk<NSYS+1;main_clk++) if(rtk->sol.dtr[main_clk]!=0.0&&rtk->exist_sys[main_clk]) break;
                main_irc=IC(main_clk,opt);
                cdtr=x[main_irc];

                /*Multi-GNSS receiver DCB*/
                if (opt->ionoopt==IONOOPT_UC_CONS&&opt->nf>=2){
                    freq=freq_base2;
                    alpha=SQR(freq_base)/(SQR(freq_base)-SQR(freq));
                    beta=-SQR(freq)/(SQR(freq_base)-SQR(freq));

                    int irdcb=IRDCB(isys,opt);
                    rdcb=x[irdcb]*(j/2==0?beta:-alpha)*(j%2==0?0.0:1.0);
                    H[irdcb+nx*nv]=(j/2==0?beta:-alpha)*(j%2==0?0.0:1.0);
                }

                /*Multi-GNSS inter-frequency clock bias*/
                if(opt->nf>=3){
                    int iifcb=IIFCB(isys,opt);
                    if(sys==SYS_GPS||sys==SYS_QZS){
                        if(j%2==1&&j/2==2){ /*L5*/
                            rifcb=x[iifcb];
                            H[iifcb+nx*nv]=1.0;
                        }
                        else rifcb=0.0;
                    }
                    else if(sys==SYS_GAL){ /*E5b*/
                        if(j%2==1&&j/2==1){
                            rifcb=x[iifcb];
                            H[iifcb+nx*nv]=1.0;
                        }
                        else rifcb=0.0;
                    }
                    else if(sys==SYS_CMP){
                        if(!strcmp(opt->ac_name,"com")){ /*B3I*/
                            if(j%2==1&&j/2==2){
                                rifcb=x[iifcb];
                                H[iifcb+nx*nv]=1.0;
                            }
                            else rifcb=0.0;
                        }
                        else{ /*B2I*/
                            if(j%2==1&&j/2==1){
                                rifcb=x[iifcb];
                                H[iifcb+nx*nv]=1.0;
                            }
                            else rifcb=0.0;
                        }
                    }
                    else rifcb=0.0;
                }

                /*troposphere delay*/
                if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                    itrp=IT(opt);
                    for (k=0;k<(opt->tropopt>=TROPOPT_ESTG?3:1);k++) {
                        H[itrp+k+nx*nv]=dtdx[k];
                    }
                }

                /*ionospheric delay*/
                if (opt->ionoopt==IONOOPT_UC||opt->ionoopt==IONOOPT_UC_CONS) {
                    if ((freq=freqs[j/2])==0.0) continue;
                    freq_base=freqs[frq_idxs[0]-1];
                    C=SQR(freq_base/freq)*(j%2==0?-1.0:1.0);

                    iion=II(sat,opt);
                    if (rtk->x[iion]==0.0) continue;
                    H[iion+nx*nv]=C;
                    if(post) rtk->ssat[sat-1].tec=rtk->x[iion]/tec_fact;
                }

                /*phase bias*/
                if (j%2==0) {
                    vs++;
                    iamb=IB(sat,j/2,opt);
                    if ((bias=x[iamb])==0.0) continue;
                    H[iamb+nx*nv]=1.0;
                    if(post){
                        if(rtk->opt.ionoopt==IONOOPT_IFLC||rtk->opt.ionoopt==IONOOPT_IF2){
                            rtk->ssat[sat-1].amb[j/2]=x[iamb];
                        }
                        else{
                            rtk->ssat[sat-1].amb[j/2]=x[iamb]*(CLIGHT/freqs[j/2]);
                        }
                    }
                }

                /* residual */
                v[nv]=y-(r+cdtr+rdcb+rifcb-CLIGHT*dts[i*2]+dtrp+C*dion+bias-shapiro);

                /* variance */
                var[nv]=varerr(obs[i].sat,sys,azel[1+i*2],0.25*rtk->ssat[sat-1].snr_rover[j/2],j/2,(opt->nf==1&&opt->ionoopt==IONOOPT_IFLC)?1:j%2,opt)+
                        vart+SQR(C)*vari+var_rs[i];
                var[nv]*=rtk->ssat[sat-1].var_fact[j%2][j/2];

                /* down weight glonass ifb */
                if (sys==SYS_GLO&&j%2==1) var[nv]+=VAR_GLO_IFB;
                /* down weight eclipse satellite */
                var[nv]*=rtk->ssat[sat-1].eclipse;

                /*inter-system bias*/
                if(sys==SYS_GLO)      {
                    irc=IC(1,opt);
                    H[irc+nx*nv]=1.0;
                    if(main_irc!=irc){
                        H[main_irc+nx*nv]=1.0;
                        isb=x[irc];v[nv]-=isb;
                    }
                    mask[1]=1;
                }
                else if(sys==SYS_GAL) {
                    irc=IC(2,opt);
                    H[irc+nx*nv]=1.0;
                    if(main_irc!=irc){
                        H[main_irc+nx*nv]=1.0;
                        isb=x[irc];v[nv]-=isb;
                    }
                    mask[2]=1;
                }
                else if(sys==SYS_CMP) {
                    if(opt->navsys&SYS_BD3&&opt->bd3opt>=BD3OPT_BD2_3) bd3_flag=1;
                    if(opt->bd3opt==BD3OPT_OFF||opt->bd3opt==BD3OPT_BD23||(opt->bd3opt==BD3OPT_BD2_3&&prn<=18)){
                        irc=IC(3,opt);
                        H[irc+nx*nv]=1.0;
                        if(main_irc!=irc){
                            H[main_irc+nx*nv]=1.0;
                            isb=x[irc];v[nv]-=isb;
                        }
                        mask[3]=1;
                    }
                    else if((opt->bd3opt==BD3OPT_BD2_3&&prn>18)||opt->bd3opt==BD3OPT_BD3){
                        irc=IC(5,opt);
                        H[irc+nx*nv]=1.0;
                        if(main_irc!=irc){
                            H[main_irc+nx*nv]=1.0;
                            isb=x[irc];v[nv]-=isb;
                        }
                        mask[5]=1;
                    }
                }
                else if(sys==SYS_QZS) {
                    irc=IC(4,opt);
                    H[irc+nx*nv]=1.0;
                    if(main_irc!=irc){
                        isb=x[irc];v[nv]-=isb;
                        H[main_irc+nx*nv]=1.0;
                    }
                    mask[4]=1;
                }
                else{
                    irc=IC(0,opt);
                    H[irc+nx*nv]=1.0;
                    mask[0]=1;
                }
                if (j%2==0) rtk->ssat[sat-1].resc[j/2]=v[nv];
                else        rtk->ssat[sat-1].resp[j/2]=v[nv];

                tracesatinfo(opt,level,j/2,frq_idxs,j%2,rtk,obs+i,y,r,cdtr,isb,rdcb,rifcb,dtrp,dion,dantr[frq_idxs[j/2]-1],shapiro,
                                 rs+6*i,dts[i*2],cbias[j/2],dants[frq_idxs[j/2]-1],rtk->ssat[sat-1].phw,bias,v[nv],sqrt(var[nv]));

                /* reject satellite by pre-fit residuals */
                if (!post&&opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) {
                    trace(2,"outlier (%d) rejected %s %s %s%d res=%9.4f el=%4.1f\n",
                          post,str,sat_id(sat),j%2?"P":"L",j/2+1,v[nv],azel[1+i*2]*R2D);
                    exc[sat-1]=1; rtk->ssat[sat-1].rejc[j%2]++;
                    continue;
                }
                /* record large post-fit residuals */
                if (post&&fabs(v[nv])>sqrt(var[nv])*THRES_REJECT) {
                    obsi[ne]=i; frqi[ne]=j; ve[ne]=v[nv];vare[ne]=var[nv]; ne++;
                }
                if (j%2==0){
                    rtk->ssat[sat-1].vsat[j/2]=1;
                }
                v_flag[nv]=(sat<<8)|(j%2==1?1:0)<<4|(j/2);
                nv++;
            }

            if(opt->ionoopt==IONOOPT_UC_CONS){
                double var_tec=1.0;
                iontec(obs[i].time,nav,pos,azel+i*2,1,&tec_ion,&var_tec);
                var_tec=iontecvar(rtk->epoch,obs[0].time,pos,rtk->ssat[sat-1].azel);
                iion=II(sat,opt);
                v[nv]=tec_ion-x[iion];
                for(int kk=0;kk<rtk->nx;kk++) H[kk+rtk->nx*nv]=kk==iion?1.0:0.0;
                var[nv]=var_tec;
                trace(3,"ION res=%9.3f sig=%9.6f tec_ion=%12.3f ion=%12.3f\n",v[nv],sqrt(var[nv]),tec_ion,x[iion]);
                nv++;
            }
        }

        /* constraint to avoid rank-deficient */
        for(i=0;i<(bd3_flag?5+1:5);i++){
            if(mask[i]) continue;
            v[nv]=0.0;
            irc=IC(i,opt);
            for(j=0;j<nx;j++) H[j+nv*nx]=j==irc?1.0:0.0;
        }

        /* reject satellite with large and max post-fit residual */
        if (post&&ne>0&&post!=9) {
            vmax=ve[0]; maxobs=obsi[0]; maxfrq=frqi[0]; rej=0;
            for (j=1;j<ne;j++) {
                if (fabs(vmax)>=fabs(ve[j])) continue;
                vmax=ve[j];varmax=var[nv]; maxobs=obsi[j]; maxfrq=frqi[j]; rej=j;
            }
            sat=obs[maxobs].sat;
            trace(2,"%s(%d): outlier (%d) rejected %s %s%d res=%9.4f el=%4.1f\n",
                  str,rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,post,sat_id(sat),maxfrq%2?"P":"L",maxfrq/2+1,vmax,azel[1+maxobs*2]*R2D);
            exc[sat-1]=1; rtk->ssat[sat-1].rejc[maxfrq/2]++; stat=0;
            ve[rej]=0;
        }
        /* constraint to local correction */
        nv+=const_corr(obs,n,exc,nav,x,pos,azel,rtk,v+nv,H+nv*rtk->nx,var+nv);

        for (i=0;i<nv;i++) for (j=0;j<nv;j++) {
                R[i+j*nv]=i==j?var[i]:0.0;
        }
        *valid_ns=vs;
        return post?stat:nv;
    }
}
/* number of estimated states ------------------------------------------------*/
extern int pppnx(const prcopt_t *opt)
{
    return NX(opt);
}

extern int pppna(const prcopt_t *opt)
{
    return NR(opt);
}

/* update solution status ----------------------------------------------------*/
static void update_stat(rtk_t *rtk, const obsd_t *obs, int n, int stat)
{
    const prcopt_t *opt=&rtk->opt;
    int i,j;
    
    /* test # of valid satellites */
    rtk->sol.ns=0;
    for (i=0;i<MAXSAT;i++) {
        for (j=0;j<opt->nf;j++) {
            if (!rtk->ssat[i].vsat[j]){
                rtk->ssat[i].resc[j]=rtk->ssat[i].resp[j]=0;
                continue;
            }
            rtk->ssat[i].outc[j]=0;

            if(rtk->ssat[i].lock[j]<0||(rtk->nfix>0&&rtk->ssat[i].fix[j]==2)){
                rtk->ssat[i].lock[j]++;
            }
            if (j==0) rtk->sol.ns++;
        }
    }
    rtk->sol.stat=rtk->sol.ns<MIN_NSAT_SOL?SOLQ_NONE:stat;
    
    if (rtk->sol.stat==SOLQ_FIX) {
        if(!rtk->tc){
            for (i=0;i<3;i++) {
                rtk->sol.rr[i]=rtk->xa[i];
                rtk->sol.qr[i]=(float)rtk->Pa[i+i*rtk->na];
            }
            rtk->sol.qr[3]=(float)rtk->Pa[1];
            rtk->sol.qr[4]=(float)rtk->Pa[1+2*rtk->na];
            rtk->sol.qr[5]=(float)rtk->Pa[2];
        }
    }
    else {
        if(!rtk->tc&&stat!=SOLQ_SINGLE){
            for (i=0;i<3;i++) {
                rtk->sol.rr[i]=rtk->x[i];
                rtk->sol.qr[i]=(float)rtk->P[i+i*rtk->nx];
            }
            rtk->sol.qr[3]=(float)rtk->P[1];
            rtk->sol.qr[4]=(float)rtk->P[2+rtk->nx];
            rtk->sol.qr[5]=(float)rtk->P[2];
        }
    }

    int irc;
    for(i=0;i<6;i++){
        irc=IC(i,opt);
        rtk->sol.dtr[i]=rtk->x[irc]/CLIGHT;
    }

    int it;
    it=IT(opt);
    rtk->sol.ztrp[1]=rtk->x[it];

    if(opt->ionoopt==IONOOPT_UC_CONS){
        int irdcb;
        for(i=0;i<6;i++){
            irdcb=IRDCB(i,opt);
            rtk->sol.rdcb[i]=rtk->x[irdcb]/CLIGHT;
        }
    }

    if(opt->nf>=3&&(opt->ionoopt==IONOOPT_UC||opt->ionoopt==IONOOPT_UC_CONS)){
        int iifcb;
        for (i=0;i<6;i++){
            iifcb=IIFCB(i,opt);
            rtk->sol.ifcb[i]=rtk->x[iifcb]/CLIGHT;
        }
    }
    
    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<opt->nf;j++) {
        rtk->ssat[obs[i].sat-1].snr_rover[j]=obs[i].SNR[j];
        rtk->ssat[obs[i].sat-1].snr_base[j] =0;
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) {
        if (rtk->ssat[i].slip[j]&3) rtk->ssat[i].slipc[j]++;
        else rtk->ssat[i].slipc[j]=0;

        if (rtk->ssat[i].fix[j]==2&&stat!=SOLQ_FIX) rtk->ssat[i].fix[j]=1;
    }
}
/* test hold ambiguity -------------------------------------------------------*/
static int test_hold_amb(rtk_t *rtk)
{
    int i,j,stat=0;
    
    /* no fix-and-hold mode */
    if (rtk->opt.modear!=ARMODE_FIXHOLD) return 0;
    
    /* reset # of continuous fixed if new ambiguity introduced */
    for (i=0;i<MAXSAT;i++) {
        if (rtk->ssat[i].fix[0]!=2&&rtk->ssat[i].fix[1]!=2) continue;
        for (j=0;j<MAXSAT;j++) {
            if (rtk->ssat[j].fix[0]!=2&&rtk->ssat[j].fix[1]!=2) continue;
            if (!rtk->ambc[j].flags[i]||!rtk->ambc[i].flags[j]) stat=1;
            rtk->ambc[j].flags[i]=rtk->ambc[i].flags[j]=1;
        }
    }
    if (stat) {
        rtk->nfix=0;
        return 0;
    }
    /* test # of continuous fixed */
    return ++rtk->nfix>=rtk->opt.minfix;
}

extern int scan_ppp(rtk_t *rtk,const prcopt_t *opt,obsd_t *obs,int n)
{
    int i,m,sat,f1=0,f2=1,f3=2,f4=3,f5=4,sys,sys_idx=-1,prn;
    double dC1P1=0.0,dP1P2=0.0;

    for(i=m=0;i<n&&i<MAXOBS;i++){
        sat=obs[i].sat;
        sys=satsys(sat,&prn);
        sys_idx=satsysidx(sat);
        if(sys_idx==-1) continue;

        if(!rtk->ssat[sat-1].vs) continue;
        if(sys==SYS_CMP&&prn<=5&&!opt->geo_opt){
            continue;
        }

        if(sys==SYS_GPS){
            f1=opt->gnss_frq_idx[sys_idx][0]-1;
            if(opt->nf>=2) f2=opt->gnss_frq_idx[sys_idx][1]-1;
            if(opt->nf>=3) f3=opt->gnss_frq_idx[sys_idx][2]-1;
        }
        else if(sys==SYS_GLO){
            f1=opt->gnss_frq_idx[sys_idx][0]-1;
            if(opt->nf>=2) f2=opt->gnss_frq_idx[sys_idx][1]-1;
        }
        else if(sys==SYS_GAL){
            f1=opt->gnss_frq_idx[sys_idx][0]-1;
            if(opt->nf>=2) f2=opt->gnss_frq_idx[sys_idx][1]-1;
            if(opt->nf>=3) f3=opt->gnss_frq_idx[sys_idx][2]-1;
            if(opt->nf>=4) f4=opt->gnss_frq_idx[sys_idx][3]-1;
            if(opt->nf>=5) f5=opt->gnss_frq_idx[sys_idx][4]-1;
        }
        else if(sys==SYS_CMP){
            if(prn<=18){
                f1=opt->gnss_frq_idx[sys_idx][0]-1;
                if(opt->nf>=2) f2=opt->gnss_frq_idx[sys_idx][1]-1;
                if(opt->nf>=3) f3=opt->gnss_frq_idx[sys_idx][2]-1;
            }
            else{
                f1=opt->gnss_frq_idx[sys_idx][0]-1;
                if(opt->nf>=2) f2=opt->gnss_frq_idx[sys_idx][1]-1;
                if(opt->nf>=3) f3=opt->gnss_frq_idx[sys_idx][2]-1;
                if(opt->nf>=4) f4=opt->gnss_frq_idx[sys_idx][3]-1;
                if(opt->nf>=5) f5=opt->gnss_frq_idx[sys_idx][4]-1;
            }
        }
        else if(sys==SYS_QZS){
            f1=opt->gnss_frq_idx[sys_idx][0]-1;
            if(opt->nf>=2) f2=opt->gnss_frq_idx[sys_idx][1]-1;
            if(opt->nf>=3) f3=opt->gnss_frq_idx[sys_idx][2]-1;
        }

        if((sys==SYS_GPS||sys==SYS_GLO)&&NEXOBS>1){ /* C1 P1 pseudorange check*/
            if(obs[i].P[f1]*obs[i].P[NFREQ]!=0.0){
                dC1P1=fabs(obs[i].P[f1]-obs[i].P[NFREQ]);
                if(dC1P1>30.0){
                    rtk->ssat[obs[i].sat-1].vs=0;
                    trace(2,"%s: exclude %s due to dP1C1=%7.3f>30.0\n",time_str(obs[0].time,1),sat_id(sat),dC1P1);
                    continue;
                }
            }
        }

        if(opt->nf>=2){
            if(obs[i].P[f1]*obs[i].P[f2]!=0.0){
                dP1P2=fabs(obs[i].P[f1]-obs[i].P[f2]);
                if(dP1P2>60.0){
                    rtk->ssat[obs[i].sat-1].vs=0;
                    trace(2,"%s: exclude %s due to no dP1P2=%7.3f>60.0\n",time_str(obs[0].time,1),sat_id(sat),dP1P2);
                    continue;
                }
                if(obs[i].L[f1]*obs[i].L[f2]==0.0){ /*check obs*/
                    rtk->ssat[obs[i].sat-1].vs=0;
                    trace(3,"%s: exclude %s due to no phase observation\n",time_str(obs[0].time,1),sat_id(sat));
                    continue;
                }
            }
            else{
                rtk->ssat[obs[i].sat-1].vs=0;
                trace(3,"%s: exclude %s due to no pseudorange observation\n",time_str(obs[0].time,1),sat_id(sat));
                continue;
            }
        }

        obs[m]=obs[i];
        m++;
    }

    return m;
}

/* validation of solution ----------------------------------------------------*/
static int valpos(rtk_t *rtk, const double *v, const double *R, const int *vflg,
                  int nv, double thres)
{
    prcopt_t *opt=&rtk->opt;
    double vv=0.0;

    int i,stat=1;

    trace(3,"valpos  : nv=%d thres=%.1f\n",nv,thres);

    if (stat&&nv>NP(opt)) {

        /* chi-square validation */
        for (i=0;i<nv;i++) vv+=v[i]*v[i]/R[i+i*nv];

        if (vv>chisqr[nv-NP(opt)-1]) {
            trace(2,"%s(%d):residuals validation failed (nv=%d np=%d vv=%.2f cs=%.2f)\n",
                   time_str(rtk->sol.time,1),!rtk->tc?rtk->epoch:rtk->ins_kf->couple_epoch,nv,NP(opt),vv,chisqr[nv-NP(opt)-1]);
            stat=0;
        }
        else {
            trace(3,"valpos : validation ok (%s nv=%d np=%d vv=%.2f cs=%.2f)\n",
                  time_str(rtk->sol.time,2),nv,NP(opt),vv,chisqr[nv-NP(opt)-1]);
        }
    }
    return stat;
}


/* precise point positioning -------------------------------------------------*/
extern int pppos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    const prcopt_t *opt=&rtk->opt;
    double *rs,*dts,*var,*v,*H,*R,*azel,*xp,*Pp,*xa,*Pa,*norm_v,*post_v,*bias,dr[3]={0},std[3],rr[3],var_pos=0.0;
    char str[32];
    int i,j,nv,info,svh[MAXOBS],exc[MAXSAT]={0},stat=SOLQ_SINGLE,vflg[MAXOBS*NFREQ*2+1],qc_flag=0,valid_ns=0;
    res_t res={0};

    int check_obs=0,check_pri=0;
    time2str(obs[0].time,str,2);
    trace(3,"pppos   : time=%s nx=%d n=%d\n",str,rtk->nx,n);

    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel=zeros(2,n);
    
    for (i=0;i<MAXSAT;i++){
        rtk->ssat[i].sys=satsys(i+1,NULL);                                        /* gps system */
        for (j=0;j<opt->nf;j++) {
            rtk->ssat[i].fix[j]=0;
            rtk->ssat[i].snr_rover[j]=obs[i].SNR[j];
            rtk->ssat[i].snr_base[j] =0;
            rtk->ssat[i].eclipse = 1.0;
            rtk->ssat[i].var_fact[0][j]=1.0;
            rtk->ssat[i].var_fact[1][j]=1.0;
        }
    }

    /* satellite positions and clocks */
    satposs(&rtk->opt,obs[0].time,obs,n,nav,rtk->opt.sateph,rs,dts,var,svh);

    /* detect cycle slip */
    detecs_ppp(obs,rtk,n,nav);

    /* temporal update of ekf states */
    udstate_ppp(rtk,obs,n,nav);

    /* exclude measurements of eclipsing satellite (block IIA) */
    if (rtk->opt.posopt[3]) {
        testeclipse(obs,n,nav,rs,rtk);
    }

    for(j=0;j<3;j++) rr[j]=rtk->x[j];

    /* earth tides correction */
    if (opt->tidecorr) {
        tidedisp(gpst2utc(obs[0].time),rr,opt->tidecorr==1?1:7,&nav->erp,
                 opt->odisp[0],dr);
    }
    nv=n*rtk->opt.nf*2+MAXSAT+3;
    xp=mat(rtk->nx,1); Pp=zeros(rtk->nx,rtk->nx);
    xa=mat(rtk->nx,1); Pa=zeros(rtk->nx,rtk->nx);
    v=mat(nv,1); H=mat(rtk->nx,nv); R=mat(nv,nv);
    post_v=mat(nv,1);
    norm_v=mat(nv,1);
    bias=mat(rtk->nx,1);
    solins_t sol_copy={0};

    for (i=0;i<MAX_ITER;i++) {
        qc_flag=0;
        for(j=0;j<3;j++) rr[j]=rtk->x[j];

        matcpy(xp,rtk->x,rtk->nx,1);
        matcpy(Pp,rtk->P,rtk->nx,rtk->nx);

        /* prefit residuals */
        if (!(nv=ppp_res(0,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel,rr,vflg,&valid_ns))) {
            trace(2,"%s(%d): ppp (%d) no valid obs data\n",str,rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,i+1);
            break;
        }

        if(pri_res_check(obs[0].time,rtk,v,vflg,nv,exc)){
            continue;
        }
        /* measurement update of ekf states */
        int tra=0;
        init_prires(v,vflg,nv,&res);
        if ((info=filter(xp,Pp,H,v,R,rtk->nx,nv,opt->robust,opt->kfopt,&res,tra))) {
            trace(2,"%s ppp (%d) filter error info=%d\n",str,i+1,info);
            break;
        }

        /* postfit residuals */
        for(j=0;j<3;j++) rr[j]=xp[j];

        if(rtk->opt.robust){
            ppp_res(i+1,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,post_v,H,R,azel,rr,vflg,&valid_ns);
            init_postres(rtk,res.post_v,&res,res.Qvv);
            qc_flag=resqc(obs[0].time,rtk,&res,exc,1);
            if(qc_flag==1&&i>=5) qc_flag=0;
        }

        freeres(&res);

        if(qc_flag&&i<5) continue;
        else{
            if(!rtk->opt.robust){
                if(!ppp_res(i+1,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,post_v,H,R,azel,rr,vflg,&valid_ns)){
                    continue;
                }
            }

            stat=SOLQ_FLOAT;
            break;
        }
    }

    if (i>=MAX_ITER) {
        trace(2,"%s ppp (%d) iteration overflows\n",str,i);
    }

    /*check solution*/
    if(valpos(rtk,post_v,R,vflg,nv,4.0)){
        matcpy(rtk->x,xp,rtk->nx,1);
        matcpy(rtk->P,Pp,rtk->nx,rtk->nx);
    }
    else{
        stat=SOLQ_NONE;
    }

    if (stat!=SOLQ_NONE&&(opt->modear==ARMODE_PPPAR||opt->modear==ARMODE_PPPAR_ILS)) {
        matcpy(xa,xp,rtk->nx,1);

        /* ambiguity resolution in ppp */
        if(manage_pppar(rtk,bias,xa,Pa,1,obs,n,nav,exc)){
            for(int k=0;k<3;k++) rr[k]=xa[k];

            if(ppp_res(9,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel,rr,vflg,&valid_ns)){
                stat=SOLQ_FIX;
                rtk->fix_epoch++;
                rtk->nfix++;
                /* hold fixed ambiguities */

                if(rtk->tc){ /*continues*/
                    *rtk->ins_kf->sol=sol_copy;
                    matcpy(rtk->x,xa,rtk->nx,1);
                    *rtk->ins_kf->insstate=sol_copy;
                }
                if(rtk->tc){
                    matcpy(rtk->x,xa,rtk->nx,1);
                    *rtk->ins_kf->insstate=sol_copy;
                }
            }
            else {
                rtk->nfix=0;
            }
        }
    }

    /* update solution status */
    update_stat(rtk,obs,n,stat);

    free(norm_v);free(post_v);free(bias);
    free(rs); free(dts); free(var); free(azel);
    free(xp); free(Pp); free(v); free(H); free(R);
    free(xa);free(Pa);

    return stat==SOLQ_NONE?0:1;
}
