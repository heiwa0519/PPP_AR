/*------------------------------------------------------------------------------
* pntpos.c : standard positioning
*
*          Copyright (C) 2007-2018 by T.TAKASU, All rights reserved.
*
* version : $Revision:$ $Date:$
* history : 2010/07/28 1.0  moved from rtkcmn.c
*                           changed api:
*                               pntpos()
*                           deleted api:
*                               pntvel()
*           2011/01/12 1.1  add option to include unhealthy satellite
*                           reject duplicated observation data
*                           changed api: ionocorr()
*           2011/11/08 1.2  enable snr mask for single-mode (rtklib_2.4.1_p3)
*           2012/12/25 1.3  add variable snr mask
*           2014/05/26 1.4  support galileo and beidou
*           2015/03/19 1.5  fix bug on ionosphere correction for GLO and BDS
*           2018/10/10 1.6  support api change of satexclude()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants -----------------------------------------------------------------*/

#define SQR(x)      ((x)*(x))
#define MAX(x,y)    ((x)>=(y)?(x):(y))

#define NUM_SYS      6
#define NX          (3+NUM_SYS)       /* # of estimated parameters x,y,z,GPS_ClK,GR_ISB,GE_ISB GB2_ISB GJ_ISB GB3_ISB*/

#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay std (m) */
#define ERR_TROP    3.0         /* tropspheric delay std (m) */
#define ERR_SAAS    0.3         /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5         /* broadcast iono model error factor */
#define ERR_CBIAS   0.3         /* code bias error std (m) */
#define REL_HUMI    0.7         /* relative humidity for saastamoinen model */

/* pseudorange measurement error variance ------------------------------------*/
static double varerr(int type,const prcopt_t *opt, double el, double snr_rover, int sys)
{
#if 1
    double a,b,snr_max;
    double fact=type?opt->err[0]:1.0;
    double sinel=sin(el);
    
    switch (sys) {
        case SYS_GPS: fact *= EFACT_GPS; break;
        case SYS_GLO: fact *= EFACT_GLO; break;
        case SYS_SBS: fact *= EFACT_SBS; break;
        default:      fact *= EFACT_GPS; break;
    }
        
    a=fact*opt->err[1];
    b=fact*opt->err[2];
    snr_max=opt->err[5];
    
    /* note: SQR(3.0) is approximated scale factor for error variance 
       in the case of iono-free combination */
    fact = (opt->ionoopt == IONOOPT_IFLC) ? SQR(3.0) : 1.0;
    switch (opt->weightmode) {
        case WEIGHTOPT_ELEVATION: return fact * ( SQR(a) + SQR(b / sinel) );
        case WEIGHTOPT_SNR      : return fact * SQR(a) * pow(10, 0.1 * MAX(snr_max - snr_rover, 0)); 
                                                   ;
        default: return 0;
    }
#else
    double fact,varr;
    fact=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
    varr=SQR(opt->err[0])*(SQR(opt->err[1])+SQR(opt->err[2])/SQR(sin(el)));
    if (opt->ionoopt==IONOOPT_IFLC) varr*=SQR(3.0); /* iono-free */
    return SQR(fact)*varr;
#endif
}
/* get tgd parameter (m) -----------------------------------------------------*/
static double gettgd(int sat, const nav_t *nav,int type)
{
    int i,sys=satsys(sat,NULL);

    if (sys==SYS_GLO) {
        for (i=0;i<nav->ng;i++) {
            if (nav->geph[i].sat==sat) break;
        }
        return (i>=nav->ng)?0.0:-nav->geph[i].dtaun*CLIGHT;
    }
    else {
        for (i=0;i<nav->n;i++) {
            if (nav->eph[i].sat==sat) break;
        }
        return (i>=nav->n)?0.0:nav->eph[i].tgd[type]*CLIGHT;
    }
}
/* test SNR mask -------------------------------------------------------------*/
extern int snrmask(const obsd_t *obs, const double *azel, const prcopt_t *opt)
{
    if (testsnr(0,0,azel[1],obs->SNR[0]*SNR_UNIT,&opt->snrmask)) {
        return 0;
    }
    if (opt->ionoopt==IONOOPT_IFLC) {
        if (testsnr(0,1,azel[1],obs->SNR[1]*SNR_UNIT,&opt->snrmask)) return 0;
    }
    return 1;
}
/* psendorange with code bias correction -------------------------------------*/
static double prange(const obsd_t *obs, const nav_t *nav, const double *azel,
                     int iter, const prcopt_t *opt, double *var)
{
#if 0
    double P1=0.0,P2=0.0,gamma,b1,b2,freq1,freq2,alpha=0.0,beta=0.0,beta13,bias=0.0;
    int sat,sys,frq1_idx=0,frq2_idx=1,prn=0;

    sat=obs->sat;
    sys=satsys(sat,&prn);

    if(sys==SYS_GPS){
        frq1_idx = opt->gps_frq_idx[0]-1;
        frq2_idx = opt->gps_frq_idx[1]-1;
    }
    else if(sys==SYS_GLO){
        frq1_idx = opt->glo_frq_idx[0]-1;
        frq2_idx = opt->glo_frq_idx[1]-1;
    }
    else if(sys==SYS_GAL){
        frq1_idx = opt->gal_frq_idx[0]-1;
        frq2_idx = opt->gal_frq_idx[1]-1;
    }
    else if(sys==SYS_CMP){
        if(prn>18){
            frq1_idx = opt->bd3_frq_idx[0]-1;
            frq2_idx = opt->bd3_frq_idx[1]-1;
        }
        else{
            frq1_idx = opt->bd2_frq_idx[0]-1;
            frq2_idx = opt->bd2_frq_idx[1]-1;
        }
    }
    else if(sys==SYS_QZS){
        frq1_idx = opt->gps_frq_idx[0]-1;
        frq2_idx = opt->gps_frq_idx[1]-1;
    }

    if(frq1_idx==-1||(opt->ionoopt==IONOOPT_IFLC&&frq2_idx==-1)) return 0.0;
    freq1=sat2freq(sat,obs->code[frq1_idx],nav);
    freq2=sat2freq(sat,obs->code[frq2_idx],nav);

    P1=obs->P[frq1_idx];
    if(frq2_idx!=-1) P2=obs->P[frq2_idx];
    *var=0.0;

    if (P1==0.0||(opt->ionoopt==IONOOPT_IFLC&&P2==0.0)) return 0.0;

    /* P1-C1,P2-C2 DCB correction*/
    if (sys==SYS_GPS||sys==SYS_GLO) {
        if (obs->code[0]==CODE_L1C) P1+=nav->cbias[sat-1][1]; /* C1->P1 */
        if (obs->code[1]==CODE_L2C) P2+=nav->cbias[sat-1][2]; /* C2->P2 */
    }

    if (opt->ionoopt==IONOOPT_IFLC) { /* dual-frequency, GREJ broadcast ephemeris base on ionospheric-free dual-frequency, while BDS base on B3 frequency */
        if(opt->cbiaopt==CBIAS_OPT_BRD_TGD){
            if (sys==SYS_GPS||sys==SYS_QZS) { /* L1-L2,G1-G2 */
                gamma=SQR(FREQ1/FREQ2);
                return (P2-gamma*P1)/(1.0-gamma);
            }
            else if (sys==SYS_GLO) { /* G1-G2 */
                gamma=SQR(FREQ1_GLO/FREQ2_GLO);
                return (P2-gamma*P1)/(1.0-gamma);
            }
            else if (sys==SYS_GAL) { /* E1-E5b */
                gamma=SQR(FREQ1/FREQ7);
                if (getseleph(SYS_GAL)) { /* F/NAV */
                    P2-=gettgd(sat,nav,0)-gettgd(sat,nav,1); /* BGD_E5aE5b */
                }
                return (P2-gamma*P1)/(1.0-gamma);
            }
            else if (sys==SYS_CMP) { /* B1-B2 */
                gamma=SQR(((obs->code[0]==CODE_L2I)?FREQ1_CMP:FREQ1)/FREQ2_CMP);
                if      (obs->code[0]==CODE_L2I) b1=gettgd(sat,nav,0); /* TGD_B1I */
                else if (obs->code[0]==CODE_L1P) b1=gettgd(sat,nav,2); /* TGD_B1Cp */
                else b1=gettgd(sat,nav,2)+gettgd(sat,nav,4); /* TGD_B1Cp+ISC_B1Cd */
                b2=gettgd(sat,nav,1); /* TGD_B2I/B2bI (m) */
                return ((P2-gamma*P1)-(b2-gamma*b1))/(1.0-gamma);
            }
            else if (sys==SYS_IRN) { /* L5-S */
                gamma=SQR(FREQ5/FREQ9);
                return (P2-gamma*P1)/(1.0-gamma);
            }
        }
        else{ /*DCB correction*/
            double cbias=0.0;
            cbias=corr_code_bias(opt,nav,obs,frq1_idx);
            P1+=cbias;
            cbias=corr_code_bias(opt,nav,obs,frq2_idx);
            P2+=cbias;
            if(sys==SYS_GPS){
                if(frq2_idx==1){ /*L1 L2 ionospheric-free*/
                    alpha= SQR(FREQ1)/(SQR(FREQ1)-SQR(FREQ2));
                    beta =-SQR(FREQ2)/(SQR(FREQ1)-SQR(FREQ2));
                }
                else if(frq2_idx==2){ /*L1 L5 ionospheric-free*/
                    alpha= SQR(FREQ1)/(SQR(FREQ1)-SQR(FREQ2));
                    beta =-SQR(FREQ5)/(SQR(FREQ1)-SQR(FREQ5));
                }
                return alpha*P1+beta*P2;
            }
            else if(sys==SYS_GLO){
                double G1=sat2freq(sat,obs->code[frq1_idx],nav);
                double G2=sat2freq(sat,obs->code[frq2_idx],nav);
                alpha= SQR(G1)/(SQR(G1)-SQR(G2));
                beta =-SQR(G2)/(SQR(G1)-SQR(G2));
                return alpha*P1+beta*P2;
            }
            else if(sys==SYS_GAL){
                if(frq2_idx==1){ /*E1 E5b ionospheric-free*/
                    alpha= SQR(FREQ1)/(SQR(FREQ1)-SQR(FREQ7));
                    beta =-SQR(FREQ7)/(SQR(FREQ1)-SQR(FREQ7));
                }
                else if(frq2_idx==2){ /*E1 E5a ionospheric-free*/
                    alpha= SQR(FREQ1)/(SQR(FREQ1)-SQR(FREQ5));
                    beta =-SQR(FREQ5)/(SQR(FREQ1)-SQR(FREQ5));
                }
                else if(frq2_idx==3){ /*E1 E6 ionospheric-free*/
                    alpha= SQR(FREQ1)/(SQR(FREQ1)-SQR(FREQ6));
                    beta =-SQR(FREQ6)/(SQR(FREQ1)-SQR(FREQ6));
                }
                else if(frq2_idx==4){ /*E1 E5ab ionospheric-free*/
                    alpha= SQR(FREQ1)/(SQR(FREQ1)-SQR(FREQ8));
                    beta =-SQR(FREQ8)/(SQR(FREQ1)-SQR(FREQ8));
                }
                return alpha*P1+beta*P2;
            }
            else if(sys==SYS_CMP){
                satsys(sat,&prn);
                if(prn<=18){ /*BD2*/
                    if(frq2_idx==1){ /*B1 B2 ionospheric-free*/
                        alpha= SQR(FREQ1_CMP)/(SQR(FREQ1_CMP)-SQR(FREQ2_CMP));
                        beta =-SQR(FREQ2_CMP)/(SQR(FREQ1_CMP)-SQR(FREQ2_CMP));
                    }
                    else if(frq2_idx==3){ /*B1 B3 ionospheric-free*/
                        alpha= SQR(FREQ1_CMP)/(SQR(FREQ1_CMP)-SQR(FREQ3_CMP));
                        beta =-SQR(FREQ3_CMP)/(SQR(FREQ1_CMP)-SQR(FREQ3_CMP));
                    }
                    return alpha*P1+beta*P2;
                }
                else{ /*BD3*/
                    if(frq2_idx==2){ /*B1 B2a ionospheric-free*/
                        alpha= SQR(FREQ1_CMP)/(SQR(FREQ1_CMP)-SQR(FREQ5));
                        beta =-SQR(FREQ5)/(SQR(FREQ1_CMP)-SQR(FREQ5));
                    }
                    else if(frq2_idx==3){ /*B1 B3 ionospheric-free*/
                        alpha= SQR(FREQ1_CMP)/(SQR(FREQ1_CMP)-SQR(FREQ3_CMP));
                        beta =-SQR(FREQ3_CMP)/(SQR(FREQ1_CMP)-SQR(FREQ3_CMP));
                    }
                    else if(frq2_idx==4){ /*B1 B2ab ionospheric-free*/
                        alpha= SQR(FREQ1_CMP)/(SQR(FREQ1_CMP)-SQR(FREQ8));
                        beta =-SQR(FREQ8)/(SQR(FREQ1_CMP)-SQR(FREQ8));
                    }
                    else if(frq2_idx==NFREQ){ /*B1 B1C ionospheric-free*/
                        alpha= SQR(FREQ1_CMP)/(SQR(FREQ1_CMP)-SQR(FREQ1));
                        beta =-SQR(FREQ1)/(SQR(FREQ1_CMP)-SQR(FREQ1));
                    }
                    else if(frq2_idx==1+NFREQ){ /*B1 B2b ionospheric-free*/
                        alpha= SQR(FREQ1_CMP)/(SQR(FREQ1_CMP)-SQR(FREQ2_CMP));
                        beta =-SQR(FREQ2_CMP)/(SQR(FREQ1_CMP)-SQR(FREQ2_CMP));
                    }
                    return alpha*P1+beta*P2;
                }
            }
            else if(sys==SYS_QZS){
                if(frq2_idx==1){ /*L1 L2 ionospheric-free*/
                    alpha= SQR(FREQ1)/(SQR(FREQ1)-SQR(FREQ2));
                    beta =-SQR(FREQ2)/(SQR(FREQ1)-SQR(FREQ2));
                }
                else if(frq2_idx==2){ /*L1 L5 ionospheric-free*/
                    alpha= SQR(FREQ1)/(SQR(FREQ1)-SQR(FREQ2));
                    beta =-SQR(FREQ5)/(SQR(FREQ1)-SQR(FREQ5));
                }
                else if(frq2_idx==3){ /*L1 L6 ionospheric-free*/
                    alpha= SQR(FREQ1)/(SQR(FREQ1)-SQR(FREQ2));
                    beta =-SQR(FREQ6)/(SQR(FREQ1)-SQR(FREQ6));
                }
                return alpha*P1+beta*P2;
            }
        }
    }
    else { /* single-freq (L1/E1/B1) */
        *var=SQR(ERR_CBIAS);
        if(opt->cbiaopt==CBIAS_OPT_BRD_TGD){   /* TGD correction */
            if (sys==SYS_GPS||sys==SYS_QZS) { /* L1 */
                b1=gettgd(sat,nav,0); /* TGD (m) */
                return P1-b1;
            }
            else if (sys==SYS_GLO) { /* G1 */
                gamma=SQR(FREQ1_GLO/FREQ2_GLO);
                b1=gettgd(sat,nav,0); /* -dtaun (m) */
                return P1-b1/(gamma-1.0);
            }
            else if (sys==SYS_GAL) { /* E1 */
                if (getseleph(SYS_GAL)) b1=gettgd(sat,nav,0); /* BGD_E1E5a */
                else                    b1=gettgd(sat,nav,1); /* BGD_E1E5b */
                return P1-b1;
            }
            else if (sys==SYS_CMP) { /* B1I/B1Cp/B1Cd */
                if      (obs->code[0]==CODE_L2I) b1=gettgd(sat,nav,0); /* TGD_B1I */
                else if (obs->code[0]==CODE_L1P) b1=gettgd(sat,nav,2); /* TGD_B1Cp */
                else b1=gettgd(sat,nav,2)+gettgd(sat,nav,4); /* TGD_B1Cp+ISC_B1Cd */
                return P1-b1;
            }
            else if (sys==SYS_IRN) { /* L5 */
                gamma=SQR(FREQ9/FREQ5);
                b1=gettgd(sat,nav,0); /* TGD (m) */
                return P1-gamma*b1;
            }
        }
        else{  /* DCB correction */
            double cbias=0.0;
            cbias=corr_code_bias(opt,nav,obs,frq1_idx);
            P1+=cbias;
        }
    }
    return P1;
#else
    double P[NFREQ]={0},L[NFREQ]={0},Pc[NFREQ]={0},Lc[NFREQ]={0};
    int sys,prn,sys_idx=-1;

    sys=satsys(obs->sat,&prn);
    sys_idx=satsysidx(obs->sat);
    if(sys_idx==-1) return -1;
    getcorrobs(opt,obs,nav,opt->gnss_frq_idx[sys_idx],NULL,NULL,0.0,L,P,Lc,Pc,NULL,NULL,NULL);

    return opt->ionoopt==IONOOPT_IFLC?Pc[0]:P[0];
#endif
}
/* ionospheric correction ------------------------------------------------------
* compute ionospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          int    sat       I   satellite number
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    ionoopt   I   ionospheric correction option (IONOOPT_???)
*          double *ion      O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric delay (L1) variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var)
{
    trace(4,"ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f\n",
          time_str(time,3),ionoopt,sat,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* broadcast model */
    if (ionoopt==IONOOPT_BRDC) {
        int sys;
        sys=satsys(sat,NULL);
        if(sys==SYS_CMP||sys==SYS_BD3){
            *ion=klobuchar_BDS(time,nav->ion_cmp,pos,azel);
        }
        else{
            *ion=klobuchar_GPS(time,nav->ion_gps,pos,azel);
        }
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* sbas ionosphere model */
    if (ionoopt==IONOOPT_SBAS) {
        return sbsioncorr(time,nav,pos,azel,ion,var);
    }
    /* ionex tec model */
    if (ionoopt==IONOOPT_TEC) {
        return iontec(time,nav,pos,azel,1,ion,var);
    }
    /* qzss broadcast model */
    if (ionoopt==IONOOPT_QZS&&norm(nav->ion_qzs,8)>0.0) {
        *ion=klobuchar_GPS(time,nav->ion_qzs,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* lex ionosphere model */
    if (ionoopt==IONOOPT_LEX) {
        return lexioncorr(time,nav,pos,azel,ion,var);
    }

    *ion=0.0;
    *var=ionoopt==IONOOPT_OFF?SQR(ERR_ION):0.0;
    return 1;
}
/* tropospheric correction -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    tropopt   I   tropospheric correction option (TROPOPT_???)
*          double *trp      O   tropospheric delay (m)
*          double *var      O   tropospheric delay variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, int tropopt, double *trp, double *var,double *ztrph,double *ztrpw)
{
    trace(4,"tropcorr: time=%s opt=%d pos=%.3f %.3f azel=%.3f %.3f\n",
          time_str(time,3),tropopt,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* saastamoinen model */
    if (tropopt==TROPOPT_SAAS||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp=saastamoinen(time,pos,azel,REL_HUMI,0,ztrph,ztrpw);
//        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        *var=SQR(ERR_SAAS);
        if(*trp>100.0) *trp=100.0;
        else if(*trp<0.05) *trp=0.05;
        return 1;
    }
    /* sbas troposphere model */
    if (tropopt==TROPOPT_SBAS) {
        *trp=sbstropcorr(time,pos,azel,var);
        return 1;
    }
    /* no correction */
    *trp=0.0;
    *var=tropopt==TROPOPT_OFF?SQR(ERR_TROP):0.0;
    return 1;
}

static void zdres(const prcopt_t *opt,const obsd_t *obs,int n,const double *rs,const double *dts,const double *vare,
                 const int *svh,const nav_t *nav,const ssat_t *ssat,double *rr,double *y,double *e,
                 double *azel,int iter,double *var_sat,int *vsat)
{
    int i,sys,prn,sat,sys_idx=-1,ns;
    double pos[3],r,freq,dion,vion,dtrp,vtrp,ztrp[2],P,vmeas;
    double snr_rover=(ssat)?0.25*ssat->snr_rover[0]:opt->err[5];
    gtime_t time;

    ecef2pos(rr,pos);

    trace(3,"%s SPP residual(iter=%d)\n",time_str(obs->time,1),iter);
    for(i=ns=0;i<n;i++){
        vsat[i]=0;azel[2*i]=azel[2*i+1]=0.0;
        time=obs[i].time;
        sat=obs[i].sat;
        sys_idx=satsysidx(sat);
        if(sys_idx==-1) continue;
        if (!(sys=satsys(obs[i].sat,&prn))) continue;

        if((opt->bd3opt==BD3OPT_BD2_3||opt->bd3opt==BD3OPT_BD3)&&(sys==SYS_CMP&&prn>18)){
            sys=SYS_BD3;
        }

        /* reject duplicated observation data */
        if (i<n-1&&i<MAXOBS-1&&obs[i].sat==obs[i+1].sat) {
            trace(2,"duplicated observation data %s sat=%2d\n",
                  time_str(obs[i].time,3),obs[i].sat);
            i++;
            continue;
        }

        /* excluded satellite? */
        if (satexclude(obs[i].sat,vare[i],svh[i],opt)) continue;

        /* geometric distance/azimuth/elevation angle */
        if ((r=geodist(rs+i*6,rr,e+i*3))<=0.0) continue;

        if(satazel(pos,e+i*3,azel+i*2)<opt->elmin) continue;

        /* test SNR mask */
        if (!snrmask(obs+i,azel+i*2,opt)) continue;

        /* ionospheric correction */
        if (!ionocorr(time,nav,sat,pos,azel+i*2,opt->ionoopt,&dion,&vion)) {
            continue;
        }

        if ((freq=sat2freq(sat,obs[i].code[opt->gnss_frq_idx[sys_idx][0]-1],nav))==0.0) continue;
        dion*=SQR(FREQ1/freq);
        vion*=SQR(FREQ1/freq);

        /* tropospheric correction */
        if (!tropcorr(time,nav,pos,azel+i*2,opt->tropopt,&dtrp,&vtrp,&ztrp[0],&ztrp[1])) {
            continue;
        }

        /* psudorange with code bias correction */
        if ((P=prange(obs+i,nav,azel+i*2,iter,opt,&vmeas))==0.0) continue;

        y[i]=P-(r-CLIGHT*dts[i*2]+dion+dtrp);
        var_sat[i]=varerr(1,opt,azel[1+i*2],snr_rover,sys)+vare[i]+vion+vtrp;

        trace(3,"%s y=%10.3f P=%14.3f r=%14.3f dts=%12.3f dion=%6.3f dtrp=%6.3f var=%6.3f el=%3.1f\n",
                sat_id(sat),y[i],P,r,CLIGHT*dts[i*2],dion,dtrp,var_sat[i],azel[i*2+1]*R2D);
    }

}

/* pseudorange residuals -----------------------------------------------------*/
static int  rescode(int iep,int iter, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const nav_t *nav, const double *x, const prcopt_t *opt,
                   const ssat_t *ssat, double *v, double *H, double *var, 
                   double *azel, int *vsat, double *resp, int *ns,double *ztrp,
                   const kf_t *ins_kf,int sdopt,double *Ri,double*Rj,int *nb,int *b,int *vflg,int *exc)
{
    gtime_t time;
    double r,dion=0.0,dtrp,vmeas,vion,vtrp=0.0,rr[3],pos[3],cdtr=0.0,isb=0.0,P;
    double snr_rover = (ssat) ? 0.25 * ssat->snr_rover[0] : opt->err[5],freq;
    int i,j,nv=0,sys,sys_idx=-1,mask[NUM_SYS]={0},sat=0,prn=0;
    insopt_t iopt=opt->insopt;
    int tc=0,nx,irc=3,main_dtr_idx=0,level=3;

    trace(4,"resprng : n=%d\n",n);

    tc=opt->mode==PMODE_TC_SPP;
    nx=sdopt?3:NX;
    for(i=0;i<3;i++) rr[i]=x[i];

    if(!sdopt){
        irc=3;
        if(opt->navsys&SYS_GPS)     {cdtr=x[3+0];main_dtr_idx=0;}
        else if(opt->navsys&SYS_GLO){cdtr=x[3+1];main_dtr_idx=1;}
        else if(opt->navsys&SYS_GAL){cdtr=x[3+2];main_dtr_idx=2;}
        else if(opt->navsys&SYS_CMP){cdtr=x[3+3];main_dtr_idx=3;}
        else if(opt->navsys&SYS_QZS){cdtr=x[3+4];main_dtr_idx=4;}
        else if(opt->navsys&SYS_BD3&&(opt->bd3opt==BD3OPT_BD2_3||opt->bd3opt==BD3OPT_BD3)){
            cdtr=x[3+5];main_dtr_idx=5;
        }
    }
    ecef2pos(rr,pos);

    if(sdopt){
        int m,sysi;
        double *y,*var_sat,*e;
        y=zeros(1,n);var_sat=zeros(1,n);
        e=zeros(3,n);

        *b=0;
        for(i=0;i<5;i++) nb[i]=0;
        /*zero-different residual*/
        zdres(opt,obs,n,rs,dts,vare,svh,nav,ssat,rr,y,e,azel,iter,var_sat,vsat);

        for(m=0;m<5;m++){
            /*find reference satellite*/
            for (i=-1,j=0;j<n;j++) {
                sysi=satsys(obs[j].sat,NULL);
                if (!test_sys(sysi,m)) continue;
                if (y[j]==0.0) continue;
                if (i<0||azel[1+j*2]>=azel[1+i*2]) i=j;
            }
            if (i<0) continue;

            for(j=0;j<n;j++){
                if(j==i) continue;
                if(y[j]==0.0) continue;
                sysi=satsys(obs[j].sat,NULL);
                if (!test_sys(sysi,m)) continue;

                v[nv]=y[i]-y[j];
                Ri[nv]=var_sat[i];
                Rj[nv]=var_sat[j];

                for(int k=0;k<3;k++) H[k+nv*nx]=-e[i*3+k]+e[j*3+k];

                vsat[i]=1,vsat[j]=1; resp[i]=v[nv]; (*ns)++;
                vflg[nv]=(sat<<8)|1<<4|1;
                nv++;
                nb[*b]++;
            }
            (*b)++;
        }

        free(y);free(var_sat);free(e);
    }
    else{
        double e[3];
        if(iep>=135100){
            level=1;
        }
        trace(level,"%s(%d): iter=%d preudorange residual\n",time_str(obs[0].time,1),iep,iter);

        for (i=*ns=0;i<n&&i<MAXOBS;i++) {
            vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0;
            time=obs[i].time;
            sat=obs[i].sat;
            sys_idx=satsysidx(sat);
            if (!(sys=satsys(obs[i].sat,&prn))) continue;
            if(exc&&exc[sat-1]){
                vsat[i]=0; resp[i]=0.0;
                continue;
            }

            if((opt->bd3opt==BD3OPT_BD2_3||opt->bd3opt==BD3OPT_BD3)&&(sys==SYS_CMP&&prn>18)){
                sys=SYS_BD3;
            }

            /* reject duplicated observation data */
            if (i<n-1&&i<MAXOBS-1&&obs[i].sat==obs[i+1].sat) {
                trace(2,"duplicated observation data %s sat=%2d\n",
                      time_str(obs[i].time,3),obs[i].sat);
                i++;
                continue;
            }

            /* excluded satellite? */
            if (satexclude(obs[i].sat,vare[i],svh[i],opt)) continue;

            /* geometric distance/azimuth/elevation angle */
            if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;

            if(satazel(pos,e,azel+i*2)<opt->elmin) continue;

            /* test SNR mask */
            if (!snrmask(obs+i,azel+i*2,opt)) continue;

            /* ionospheric correction */
            if (!ionocorr(time,nav,sat,pos,azel+i*2,opt->ionoopt,&dion,&vion)) {
                continue;
            }

            int sys_idx=satsysidx(sat);

            if(sys==SYS_GPS){
                if ((freq=sat2freq(sat,obs[i].code[opt->gnss_frq_idx[sys_idx][0]-1],nav))==0.0) continue;
                dion*=SQR(FREQ1/freq);
                vion*=SQR(FREQ1/freq);
            }
            else if(sys==SYS_GLO){
                if ((freq=sat2freq(sat,obs[i].code[opt->gnss_frq_idx[sys_idx][0]-1],nav))==0.0) continue;
                dion*=SQR(FREQ1/freq);
            }
            else if(sys==SYS_GAL){
                if ((freq=sat2freq(sat,obs[i].code[opt->gnss_frq_idx[sys_idx][0]-1],nav))==0.0) continue;
                dion*=SQR(FREQ1/freq);
                vion*=SQR(FREQ1/freq);
            }
            else if(sys==SYS_CMP){
                if(prn<=18){
                    if ((freq=sat2freq(sat,obs[i].code[opt->gnss_frq_idx[sys_idx][0]-1],nav))==0.0) continue;
                    dion*=SQR(FREQ1_CMP/freq);
                    vion*=SQR(FREQ1_CMP/freq);
                }
                else{
                    if ((freq=sat2freq(sat,obs[i].code[opt->gnss_frq_idx[sys_idx][0]-1],nav))==0.0) continue;
                    dion*=SQR(FREQ1_CMP/freq);
                    vion*=SQR(FREQ1_CMP/freq);
                }
            }
            else if(sys==SYS_BD3){
                if ((freq=sat2freq(sat,obs[i].code[opt->gnss_frq_idx[sys_idx][0]-1],nav))==0.0) continue;
                dion*=SQR(FREQ1_CMP/freq);
                vion*=SQR(FREQ1_CMP/freq);
            }
            else if(sys==SYS_QZS){
                if ((freq=sat2freq(sat,obs[i].code[opt->gnss_frq_idx[sys_idx][0]-1],nav))==0.0) continue;
                dion*=SQR(FREQ1/freq);
                vion*=SQR(FREQ1/freq);
            }

            /* tropospheric correction */
            if (!tropcorr(time,nav,pos,azel+i*2,opt->tropopt,&dtrp,&vtrp,&ztrp[0],&ztrp[1])) {
                continue;
            }

            /* psudorange with code bias correction */
            if ((P=prange(obs+i,nav,azel+i*2,iter,opt,&vmeas))==0.0) continue;

            /* pseudorange residual */
            v[nv]=P-(r+cdtr-CLIGHT*dts[i*2]+dion+dtrp);
//            if(fabs(v[nv])>opt->maxinno) continue;

            /* design matrix */
            for (j=0;j<nx;j++) H[j+nv*nx]=j<3?-e[j]:(j==3?1.0:0.0);

            /* time system and receiver bias offset correction */
            /* time system and receiver bias offset correction */
            if      (sys==SYS_GLO) {
                if(1==main_dtr_idx) isb=0.0;
                else isb=x[irc+1];
                v[nv]-=isb; H[irc+1+nv*nx]=1.0; mask[1]=1;
            }
            else if (sys==SYS_GAL) {
                if(2==main_dtr_idx) isb=0.0;
                else isb=x[irc+2];
                v[nv]-=isb; H[irc+2+nv*nx]=1.0; mask[2]=1;
            }
            else if (sys==SYS_CMP) {
                if(3==main_dtr_idx) isb=0.0;
                else isb=x[irc+3];
                v[nv]-=isb; H[irc+3+nv*nx]=1.0; mask[3]=1;
            }
            else if (sys==SYS_QZS) {
                if(4==main_dtr_idx) isb=0.0;
                else isb=x[irc+4];
                v[nv]-=isb;H[irc+4+nv*nx]=1.0; mask[4]=1;
            }
            else if (sys==SYS_BD3) {
                if(5==main_dtr_idx) isb=0.0;
                else isb=x[irc+5];
                v[nv]-=isb;H[irc+5+nv*nx]=1.0; mask[5]=1;
            }
            else mask[0]=1;

            vsat[i]=1; resp[i]=v[nv]; (*ns)++;

            /* error variance */
            vflg[nv]=(sat<<8)|1<<4|1;
            var[nv++]=varerr(1,opt,azel[1+i*2],snr_rover,sys)+vare[i]+vion+vtrp;

            trace(level,"%s azel=%5.1f %4.1f res=%8.4f obs=%13.3f r=%13.3f rec_clk=%8.3f isb=%5.3f sat_clk=%9.3f trp=%7.3f ion=%7.3f var=%5.3f\n",sat_id(obs[i].sat),
                  azel[i*2]*R2D,azel[1+i*2]*R2D,resp[i],P,r,cdtr,isb,CLIGHT*dts[i*2],dtrp,dion,var[nv-1]);
        }

        /* constraint to avoid rank-deficient */
        for (i=0;i<NUM_SYS;i++) {
            if (mask[i]) continue;
            v[nv]=0.0;
            for (j=0;j<nx;j++) H[j+nv*nx]=j==i+3?1.0:0.0;

            vflg[nv]=(sat<<8)|0<<4|1;
            var[nv++]=0.01;
        }
    }

    return nv;
}
/* validate solution ---------------------------------------------------------*/
static int valsol(int iep,gtime_t t,const double *azel, const int *vsat, int n,
                  const prcopt_t *opt, const double *v, int nv, int nx,
                  char *msg,double *dop)
{
    double azels[MAXOBS*2],vv;
    int i,ns;
    
    trace(4,"valsol  : n=%d nv=%d\n",n,nv);
    
    /* chi-square validation of residuals */
    vv=dot(v,v,nv);
    if (nv>nx&&vv>chisqr[nv-nx-1]) {
        trace(2,"%s(%d) large chi-square error nv=%d vv=%.1f thres=%.1f in SPP\n",time_str(t,1),iep,nv,vv,chisqr[nv-nx-1]);
        /* return 0; */ /* threshold too strict for all use cases, report error but continue on */
    }
    /* large gdop check */
    for (i=ns=0;i<n;i++) {
        if (!vsat[i]) continue;
        azels[  ns*2]=azel[  i*2];
        azels[1+ns*2]=azel[1+i*2];
        ns++;
    }
    dops(ns,azels,opt->elmin,dop);
    if (dop[0]<=0.0||dop[0]>opt->maxgdop) {
        sprintf(msg,"gdop error nv=%d gdop=%.1f",nv,dop[0]);
        return 0;
    }
    return 1;
}

static int cmpres(const void *p1, const void *p2)
{
    double *q1 = (double *)p1, *q2 = (double *)p2;
    double delta = *q1 - *q2;
    return delta<-0.0 ? -1 : (delta>0.0 ? 1 : 0);
}

static int prefit_res_check(int iep,gtime_t  t,double *v,int nv,int *vflag,int *exc)
{
    double v_SYS[NUM_SYS][MAXOBS]={0},v_copy[MAXOBS]={0},mean;
    int nv_SYS[NUM_SYS]={0},sat_SYS[NUM_SYS][MAXOBS]={0};
    int i,j,sat,sys_idx,qc_flag=0;

    for(j=0;j<nv;j++){
        sat=(vflag[j]>>8)&0xFF;
        sys_idx=satsysidx(sat);
        sat_SYS[sys_idx][nv_SYS[sys_idx]]=sat;
        v_SYS[sys_idx][nv_SYS[sys_idx]++]=v[j];
    }

    for(j=0;j<NUM_SYS;j++){
        if(nv_SYS[j]<=0) continue;
        matcpy(v_copy,v_SYS[j],nv_SYS[j],1);
        mean=median(v_copy,nv_SYS[j]);
        for(i=0;i<nv_SYS[j];i++){
            if(fabs(v_SYS[j][i]-mean)>15.0){
                sat=sat_SYS[j][i];
                exc[sat-1]=1;
                qc_flag=1;
                trace(3,"%s(%d): large prior pseudorange residual in %s, res=%10.3f, mean=%10.3f\n",time_str(t,1),iep,sat_id(sat),v_SYS[j][i],mean);
            }
        }
    }

    return qc_flag;
}

/* estimate receiver position ------------------------------------------------*/
static int estpos(int iep,const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav,
                  const prcopt_t *opt, const ssat_t *ssat, sol_t *sol, double *azel,
                  int *vsat, double *resp, char *msg,int sdopt)
{
    int i,j,k,info,stat,nv=0,ns,est_nx=sdopt?3:NX,nb[5]={0},b=0,vflg[MAXOBS]={0},exc[MAXSAT]={0},sat;
    double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*Ri,*Rj,*var,sig,*R,*norm_v;

    trace(5,"estpos  : n=%d\n",n);

    if(sdopt){
        v=mat(n,1); H=mat(est_nx,n); var=mat(n,1);norm_v=mat(n,1);
    }
    else{
        v=mat(n+NUM_SYS,1); H=mat(est_nx,n+NUM_SYS); var=mat(n+NUM_SYS,1);norm_v=mat(n+NUM_SYS,1);
    }
    Ri=mat(n,1);Rj=mat(n,1);
    
    for (i=0;i<3;i++) x[i]=sol->rr[i];
    
    for (i=0;i<MAXITR;i++) {
        /* pseudorange residuals */
        nv=rescode(iep,i,obs,n,rs,dts,vare,svh,nav,x,opt,ssat,v,H,var,azel,vsat,resp,
                   &ns,sol->ztrp,NULL,sdopt,sdopt?Ri:NULL,sdopt?Rj:NULL,nb,&b,vflg,exc);

        if (sdopt?(nv<4):(nv-NUM_SYS)<4) {
            sol->ns=ns;
            trace(3,"%s(%d): SPP lack of valid sats ns=%d nv=%d\n",time_str(obs[0].time,1),iep,n,opt->sdopt?nv:nv-NUM_SYS);
            return 0;
        }

#if 0
        if(norm(sol->rr,3)>0){
            if(prefit_res_check(iep,obs->time,v,sdopt?nv:nv-NUM_SYS,vflg,exc)){
                continue;
            }
        }
#endif

        if(sdopt){
            R=mat(nv,nv);
            ddcov(nb,b,Ri,Rj,nv,R);
            if((info=lsq_(H,R,v,3,nv,dx,Q))){
                sprintf(msg,"lsq error info=%d",info);
                break;
            }
            for(j=0;j<nv;j++) v[j]/=sqrt(Rj[j]);
            if(R) free(R); R=NULL;
        }else{
            /* weight by variance */
            for (j=0;j<nv;j++) {
                sig=sqrt(var[j]);
                v[j]/=sig;
                for (k=0;k<NX;k++) H[k+j*NX]/=sig;
            }
            /* least square estimation */
            if ((info=lsq(H,v,NX,nv,dx,Q))) {
                sprintf(msg,"lsq error info=%d",info);
                break;
            }
        }

        for (j=0;j<est_nx;j++) x[j]+=dx[j];

        if (norm(dx,est_nx)<1E-4) {
            /*check residual*/
            for(j=0;j<nv;j++) norm_v[j]=fabs(v[j]);

            int res_flag=0;
            for(j=0;j<nv;j++){
                if(fabs(norm_v[j])>5.0){
                    sat=(vflg[j]>>8)&0xFF;
                    exc[sat-1]=1;
                    res_flag=1;
                }
            }

            if(res_flag){
                for (i=0;i<3;i++) x[i]=sol->rr[i];
                continue;
            }

            sol->type=0;
            if(!sdopt){
                sol->time=timeadd(obs[0].time,-x[3]/CLIGHT);
                if(opt->navsys&SYS_GPS) sol->dtr[0]=x[3]/CLIGHT; /* receiver clock bias (s) */
                if(opt->navsys&SYS_GLO) sol->dtr[1]=x[4]/CLIGHT; /* glo-gps time offset (s) */
                if(opt->navsys&SYS_GAL) sol->dtr[2]=x[5]/CLIGHT; /* gal-gps time offset (s) */
                if(opt->navsys&SYS_CMP) sol->dtr[3]=x[6]/CLIGHT; /* bds-gps time offset (s) */
                if(opt->navsys&SYS_QZS) sol->dtr[4]=x[7]/CLIGHT; /* qzs-gps time offset (s) */
                if(opt->navsys&SYS_BD3&&(opt->bd3opt==BD3OPT_BD2_3||opt->bd3opt==BD3OPT_BD3)){
                    sol->dtr[5]=x[8]/CLIGHT;
                }
            }
            else{
                sol->time=obs[0].time;
            }
            for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX];
            sol->qr[3]=(float)Q[1];    /* cov xy */
            sol->qr[4]=(float)Q[2+NX]; /* cov yz */
            sol->qr[5]=(float)Q[2];    /* cov zx */
            sol->ns=(unsigned char)ns;
            sol->age=sol->ratio=0.0;
            
            /* validate solution */
            if ((stat=valsol(iep,obs[0].time,azel,vsat,n,opt,v,nv,est_nx,msg,sol->dop))) {
                sol->stat=opt->sateph==EPHOPT_SBAS?SOLQ_SBAS:SOLQ_SINGLE;
            }
            free(v); free(H); free(var);free(norm_v);
            if(Ri) free(Ri); Ri=NULL;
            if(Rj) free(Rj); Rj=NULL;
            return stat;
        }
    }
    if (i>=MAXITR) sprintf(msg,"iteration divergent i=%d",i);
    
    free(v); free(H); free(var);free(norm_v);
    if(Ri) free(Ri); Ri=NULL;
    if(Rj) free(Rj); Rj=NULL;
    return 0;
}
/* raim fde (failure detection and exclution) -------------------------------*/
static int raim_fde(int iep,const obsd_t *obs, int n, const double *rs,
                    const double *dts, const double *vare, const int *svh,
                    const nav_t *nav, const prcopt_t *opt, const ssat_t *ssat, 
                    sol_t *sol, double *azel, int *vsat, double *resp, char *msg)
{
    obsd_t *obs_e;
    sol_t sol_e={{0}};
    char tstr[32],name[16],msg_e[128];
    double *rs_e,*dts_e,*vare_e,*azel_e,*resp_e,rms_e,rms=100.0;
    int i,j,k,nvsat,stat=0,*svh_e,*vsat_e,sat=0;
    
    trace(4,"raim_fde: %s n=%2d\n",time_str(obs[0].time,0),n);
    
    if (!(obs_e=(obsd_t *)malloc(sizeof(obsd_t)*n))) return 0;
    rs_e = mat(6,n); dts_e = mat(2,n); vare_e=mat(1,n); azel_e=zeros(2,n);
    svh_e=imat(1,n); vsat_e=imat(1,n); resp_e=mat(1,n); 
    
    for (i=0;i<n;i++) {
        
        /* satellite exclution */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            obs_e[k]=obs[j];
            matcpy(rs_e +6*k,rs +6*j,6,1);
            matcpy(dts_e+2*k,dts+2*j,2,1);
            vare_e[k]=vare[j];
            svh_e[k++]=svh[j];
        }
        /* estimate receiver position without a satellite */
        if (!estpos(iep,obs_e,n-1,rs_e,dts_e,vare_e,svh_e,nav,opt,ssat,&sol_e,azel_e,
                    vsat_e,resp_e,msg_e,0)) {
            trace(3,"raim_fde: exsat=%2d (%s)\n",obs[i].sat,msg);
            continue;
        }
        for (j=nvsat=0,rms_e=0.0;j<n-1;j++) {
            if (!vsat_e[j]) continue;
            rms_e+=SQR(resp_e[j]);
            nvsat++;
        }
        if (nvsat<5) {
            trace(3,"raim_fde: exsat=%2d lack of satellites nvsat=%2d\n",
                  obs[i].sat,nvsat);
            continue;
        }
        rms_e=sqrt(rms_e/nvsat);
        
        trace(3,"raim_fde: exsat=%2d rms=%8.3f\n",obs[i].sat,rms_e);
        
        if (rms_e>rms) continue;
        
        /* save result */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            matcpy(azel+2*j,azel_e+2*k,2,1);
            vsat[j]=vsat_e[k];
            resp[j]=resp_e[k++];
        }
        stat=1;
        sol_e.eventime = sol->eventime;
        *sol=sol_e;
        sat=obs[i].sat;
        rms=rms_e;
        vsat[i]=0;
        strcpy(msg,msg_e);
    }
    if (stat) {
        time2str(obs[0].time,tstr,2); satno2id(sat,name);
        trace(2,"%s: %s excluded by raim\n",tstr,name);
    }
    free(obs_e);
    free(rs_e ); free(dts_e ); free(vare_e); free(azel_e);
    free(svh_e); free(vsat_e); free(resp_e);
    return stat;
}

static int zdres_vel(const prcopt_t *opt,const obsd_t *obs, int n, const double *rs, const double *dts,
                     const nav_t *nav, const double *rr, const double *x,
                     const double *azel, const int *vsat, double err, double *y,double *e,double *meas_var,ssat_t *ssat)
{
    double freq,rate,pos[3],E[9],a[3],vs[3],cosel;
    int i,j,nv=0,sys;

    ecef2pos(rr,pos); xyz2enu(pos,E);
    for (i=0;i<n&&i<MAXOBS;i++) {
        y[i]=0.0;
        sys=satsys(obs[i].sat,NULL);

        freq=sat2freq(obs[i].sat,obs[i].code[0],nav);

        if(opt->velopt==VELOPT_DOPPLER){
            if (obs[i].D[0]==0.0||freq==0.0||!vsat[i]||norm(rs+3+i*6,3)<=0.0) {
                continue;
            }
        }
        else if(opt->velopt==VELOPT_TDCP){
            if (obs[i].L[0]==0.0||ssat[obs[i].sat-1].ph[0][0]==0.0||freq==0.0||!vsat[i]||norm(rs+3+i*6,3)<=0.0) {
                continue;
            }
        }

        /* LOS (line-of-sight) vector in ECEF */
        cosel=cos(azel[1+i*2]);
        a[0]=sin(azel[i*2])*cosel;
        a[1]=cos(azel[i*2])*cosel;
        a[2]=sin(azel[1+i*2]);
        matmul("TN",3,1,3,1.0,E,a,0.0,e+3*i);

        /* satellite velocity relative to receiver in ECEF */
        for (j=0;j<3;j++) {
            vs[j]=rs[j+3+i*6]-x[j];
        }
        /* range rate with earth rotation correction */
        rate=dot(vs,e+3*i,3)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
                                      rs[3+i*6]*rr[1]-rs[i*6]*x[1]);

        /* range rate residual (m/s) */

        if(opt->velopt==VELOPT_DOPPLER){
            y[i]=-obs[i].D[0]*CLIGHT/freq-(rate-CLIGHT*dts[1+i*2]);
            meas_var[i]=varerr(1,opt,azel[1+i*2],0.0,sys)*CLIGHT/freq;
        }
        else if(opt->velopt==VELOPT_TDCP){
            y[i]=(obs[i].L[0]-ssat[obs[i].sat-1].ph[0][0])*CLIGHT/freq-(rate-CLIGHT*dts[1+i*2]);
            meas_var[i]=varerr(0,opt,azel[1+i*2],0.0,sys)*CLIGHT/freq;
        }

        nv++;
    }
    return nv;
}

/* doppler residuals ---------------------------------------------------------*/
static int resvel(const prcopt_t *opt,const obsd_t *obs, int n, const double *rs, const double *dts,
                  const nav_t *nav, const double *rr, const double *x,
                  const double *azel, const int *vsat, double err, double *v,
                  double *H,int sdopt,double *Ri,double *Rj,int *nb,int *b,ssat_t *ssat)
{
    double freq,rate,pos[3],E[9],a[3],vs[3],cosel,sig;
    int i,j,nv=0,m=0,sysi;

    trace(3,"resdop  : n=%d\n",n);

    if(sdopt){
        double *e,*y,*var_sat;
        y=mat(1,n);var_sat=mat(1,n);e=mat(3,n);

        *b=0;
        for(i=0;i<5;i++) nb[i]=0;

        zdres_vel(opt,obs,n,rs,dts,nav,rr,x,azel,vsat,err,y,e,var_sat,ssat);

        for(m=0;m<5;m++){
            /*find reference satellite*/
            for (i=-1,j=0;j<n;j++) {
                sysi=satsys(obs[j].sat,NULL);
                if (!test_sys(sysi,m)) continue;
                if (y[j]==0.0) continue;
                if (i<0||azel[1+j*2]>=azel[1+i*2]) i=j;
            }
            if (i<0) continue;

            for(j=0;j<n;j++){
                if(j==i) continue;
                if(y[j]==0.0) continue;
                sysi=satsys(obs[j].sat,NULL);
                if (!test_sys(sysi,m)) continue;

                v[nv]=y[i]-y[j];
                Ri[nv]=var_sat[i];
                Rj[nv]=var_sat[j];

                /* design matrix */
                for (int k=0;k<3;k++) {
                    H[k+nv*3]=-e[3*i+k]+e[3*j+k];
                }

                nv++;
                nb[*b]++;
            }
            (*b)++;
        }
        free(e);free(y);free(var_sat);
    }
    else{
        double e[3];
        ecef2pos(rr,pos); xyz2enu(pos,E);
        for (i=0;i<n&&i<MAXOBS;i++) {

            freq=sat2freq(obs[i].sat,obs[i].code[0],nav);

            if(opt->velopt==VELOPT_DOPPLER){
                if (obs[i].D[0]==0.0||freq==0.0||!vsat[i]||norm(rs+3+i*6,3)<=0.0) {
                    continue;
                }
            }
            else if(opt->velopt==VELOPT_TDCP){
                if (obs[i].L[0]==0.0||ssat[obs[i].sat-1].ph[0][0]==0.0||freq==0.0||!vsat[i]||norm(rs+3+i*6,3)<=0.0) {
                    continue;
                }
            }
            /* LOS (line-of-sight) vector in ECEF */
            cosel=cos(azel[1+i*2]);
            a[0]=sin(azel[i*2])*cosel;
            a[1]=cos(azel[i*2])*cosel;
            a[2]=sin(azel[1+i*2]);
            matmul("TN",3,1,3,1.0,E,a,0.0,e);

            /* satellite velocity relative to receiver in ECEF */
            for (j=0;j<3;j++) {
                vs[j]=rs[j+3+i*6]-x[j];
            }
            /* range rate with earth rotation correction */
            rate=dot(vs,e,3)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
                                          rs[3+i*6]*rr[1]-rs[  i*6]*x[1]);

            /* Std of range rate error (m/s) */
            int sys=satsys(obs[i].sat,NULL);

            if(opt->velopt==VELOPT_DOPPLER){
                sig=varerr(1,opt,azel[1+i*2],0.0,sys)*CLIGHT/freq;;
                /* range rate residual (m/s) */
                v[nv]=(-obs[i].D[0]*CLIGHT/freq-(rate+x[3]-CLIGHT*dts[1+i*2]))/sig;
            }
            else if(opt->velopt==VELOPT_TDCP){
                sig=varerr(0,opt,azel[1+i*2],0.0,sys);;
                /* range rate residual (m/s) */
                v[nv]=((obs[i].L[0]-ssat[obs[i].sat-1].ph[0][0])*CLIGHT/freq-(rate+x[3]-CLIGHT*dts[1+i*2]))/sig;
            }

            /* design matrix */
            for (j=0;j<4;j++) {
                H[j+nv*4]=((j<3)?-e[j]:1.0)/sig;
            }
            nv++;
        }
    }

    return nv;
}

/* estimate receiver velocity ------------------------------------------------*/
static void estvel(const obsd_t *obs, int n, const double *rs, const double *dts,
                   const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                   const double *azel, const int *vsat,ssat_t *ssat)
{
    double x[4]={0},dx[4],Q[16],*v,*H;
    double err=opt->err[4]; /* Doppler error (Hz) */
    int i,j,nv=0,m=0,nb[5]={0},b=0,est_nx=0,tc;

    trace(5,"estvel  : n=%d\n",n);

    tc=opt->mode>=PMODE_TC_SPP&&opt->mode<=PMODE_TC_PPP;
    est_nx=opt->sdopt?3:4;

    v=mat(n,1);
    H=mat(est_nx,n);

    for (i=0;i<MAXITR;i++) {
        if(opt->sdopt){
            double *Ri,*Rj,*R;
            Ri=mat(n,1);Rj=mat(n,1);
            /* range rate residuals (m/s) */
            if ((nv=resvel(opt,obs,n,rs,dts,nav,sol->rr,x,azel,vsat,err,v,H,1,Ri,Rj,nb,&b,ssat))<4) {
                break;
            }
            R=mat(nv,nv);
            ddcov(nb,b,Ri,Rj,nv,R);
            if(lsq_(H,R,v,est_nx,nv,dx,Q)){
                break;
            }
            if(R) free(R);
            if(Ri) free(Ri);
            if(Rj) free(Rj);
        }
        else{
            /* range rate residuals (m/s) */
            if ((nv=resvel(opt,obs,n,rs,dts,nav,sol->rr,x,azel,vsat,err,v,H,0,NULL,NULL,NULL,NULL,ssat))<4) {
                break;
            }
            /* least square estimation */
            if (lsq(H,v,est_nx,nv,dx,Q)) break;
        }
        for (j=0;j<est_nx;j++) x[j]+=dx[j];

        if (norm(dx,est_nx)<1E-6) {
            matcpy(sol->rr+3,x,3,1);
            sol->qv[0]=(float)Q[0];  /* xx */
            sol->qv[1]=(float)Q[5];  /* yy */
            sol->qv[2]=(float)Q[10]; /* zz */
            sol->qv[3]=(float)Q[1];  /* xy */
            sol->qv[4]=(float)Q[6];  /* yz */
            sol->qv[5]=(float)Q[2];  /* zx */
            break;
        }
    }
    free(v); free(H);

}

/* single-point positioning ----------------------------------------------------
* compute receiver position, velocity, clock bias by single-point positioning
* with pseudorange and doppler observables
* args   : obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          prcopt_t *opt    I   processing options
*          sol_t  *sol      IO  solution
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
*          ssat_t *ssat     IO  satellite status              (NULL: no output)
*          char   *msg      O   error message for error exit
* return : status(1:ok,0:error)
* notes  : assuming sbas-gps, galileo-gps, qzss-gps, compass-gps time offset and
*          receiver bias are negligible (only involving glonass-gps time offset
*          and receiver bias)
*-----------------------------------------------------------------------------*/
extern int pntpos(int iep,const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, ssat_t *ssat,
                  char *msg,const insopt_t *iopt,kf_t *ins_kf)
{
    prcopt_t opt_=*opt;
    double *rs,*dts,*var,*azel_,*resp;
    int i,stat,vsat[MAXOBS]={0},svh[MAXOBS];

    sol->stat=SOLQ_NONE;
    
    if (n<=0) {strcpy(msg,"no observation data"); return 0;}
    
    trace(4,"pntpos  : tobs=%s n=%d\n",time_str(obs[0].time,3),n);
    
    sol->time=obs[0].time; msg[0]='\0';sol->ns=0;
    sol->eventime = obs[0].eventime;
    
    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel_=zeros(2,n); resp=mat(1,n);
    
    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            ssat[i].snr_rover[0]=0;
            ssat[i].snr_base[0] =0;
        }
        for (i=0;i<n;i++)
            ssat[obs[i].sat-1].snr_rover[0]=obs[i].SNR[0];
    }
    
    if (opt_.mode!=PMODE_SINGLE) { /* for precise positioning */
        opt_.sateph =EPHOPT_BRDC;
        opt_.ionoopt=IONOOPT_BRDC;
        opt_.tropopt=TROPOPT_SAAS;
    }
    /* satellite positons, velocities and clocks */
    satposs(opt,sol->time,obs,n,nav,opt_.sateph,rs,dts,var,svh);

    stat=estpos(iep,obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,msg,0);


    /* raim fde */
    if (!stat&&n>=6&&opt->posopt[4]) {
        stat=raim_fde(iep,obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,msg);
    }

    /* estimate receiver velocity */
    if (stat){
        estvel(obs,n,rs,dts,nav,&opt_,sol,azel_,vsat,ssat);
    }
    
    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];
    }
    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            ssat[i].vs=0;
            ssat[i].azel[0]=ssat[i].azel[1]=0.0;
            ssat[i].resp[0]=ssat[i].resc[0]=0.0;
            for(int k=0;k<NFREQ;k++) ssat[i].ph[0][k]=0.0;
        }
        for (i=0;i<n;i++) {
            ssat[obs[i].sat-1].azel[0]=azel_[  i*2];
            ssat[obs[i].sat-1].azel[1]=azel_[1+i*2];
            if (!vsat[i]) continue;
            ssat[obs[i].sat-1].vs=1;
            ssat[obs[i].sat-1].resp[0]=resp[i];
            ssat[obs[i].sat-1].snr[0]=obs[i].SNR[0];
            if(opt_.mode==PMODE_TDCP||opt_.mode==PMODE_TC_TDCP||opt_.velopt==VELOPT_TDCP){
                for(int k=0;k<NFREQ;k++){
                    ssat[obs[i].sat-1].ph[0][k]=obs[i].L[k];
                }
            }
        }
    }
    free(rs); free(dts); free(var); free(azel_); free(resp);
    return stat;
}
