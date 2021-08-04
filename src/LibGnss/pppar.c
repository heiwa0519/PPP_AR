//
// Created by chenc on 2021/3/10.
//
#include "rtklib.h"

#define ROUND(x)    (int)floor((x)+0.5)

#define LOG_PI      1.14472988584940017 /* log(pi) */
#define SQRT2       1.41421356237309510 /* sqrt(2) */
#define SWAP_I(x,y) do {int    _tmp=x; x=y; y=_tmp;} while (0)
#define SWAP_D(x,y) do {double _tmp=x; x=y; y=_tmp;} while (0)
#define MIN_AMB_RES 4         /* min number of ambiguities for ILS-AR */
#define MIN_LOCK_AR 15

/* number and index of ekf states */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics?9:3)
#define NC(opt)     ((opt)->sdopt?0:NSYS)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))
#define NI(opt)     ((opt)->ionoopt==IONOOPT_UC?MAXSAT:0)
#define ND(opt)     ((opt)->nf>=3?1:0)
#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt)+ND(opt))
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)

/* complementaty error function (ref [1] p.227-229) --------------------------*/
static double q_gamma(double a, double x, double log_gamma_a);
static double p_gamma(double a, double x, double log_gamma_a)
{
    double y,w;
    int i;

    if (x==0.0) return 0.0;
    if (x>=a+1.0) return 1.0-q_gamma(a,x,log_gamma_a);
    y=w=exp(a*log(x)-x-log_gamma_a)/a;
    for (i=1;i<100;i++) {
        w*=x/(a+i);
        y+=w;
        if (fabs(w)<1E-15) break;
    }
    return y;
}
static double q_gamma(double a, double x, double log_gamma_a)
{
    double y,w,la=1.0,lb=x+1.0-a,lc;
    int i;

    if (x<a+1.0) return 1.0-p_gamma(a,x,log_gamma_a);
    w=exp(-x+a*log(x)-log_gamma_a);
    y=w/lb;
    for (i=2;i<100;i++) {
        lc=((i-1-a)*(lb-la)+(i+x)*lb)/i;
        la=lb; lb=lc;
        w*=(i-1-a)/i;
        y+=w/la/lb;
        if (fabs(w/la/lb)<1E-15) break;
    }
    return y;
}

static double f_erfc(double x)
{
    return x>=0.0?q_gamma(0.5,x*x,LOG_PI/2.0):1.0+p_gamma(0.5,x*x,LOG_PI/2.0);
}

/* confidence function of integer ambiguity ----------------------------------*/
static double conf_func(int N, double B, double var)
{
    double x,p=1.0,sig=sqrt(var);
    int i;

    x=fabs(B-N);
    for (i=1;i<10;i++) {
        p-=f_erfc((i-x)/(SQRT2*sig))-f_erfc((i+x)/(SQRT2*sig));
    }
    return p;
}

static int matchnlupd(const gtime_t obst,int sat1,int sat2,double *nl_upd1,double *nl_upd2,const nav_t *nav)
{
    double upd1=0.0,upd2=0.0;
    int i,stat=0;

    for(i=0;i<nav->upds->nls.n;i++){
        if(timediff(nav->upds->nls.data[i].ts,obst)<=0&&timediff(nav->upds->nls.data[i].te,obst)>=0){
            upd1=nav->upds->nls.data[i].nl[sat1-1];
            upd2=nav->upds->nls.data[i].nl[sat2-1];
            stat=1;
            break;
        }
    }

    if(nl_upd1) *nl_upd1=upd1;
    if(nl_upd2) *nl_upd2=upd2;

    return stat;
}

static int matchnlfcb(const gtime_t obst,int sat1,int sat2,double *nl_fcb1,double *nl_fcb2,const nav_t *nav)
{
    double fcb1=0.0,fcb2=0.0;
    int i,stat=0;

    for(i=0;i<nav->fcbs->n;i++){
        if(timediff(nav->fcbs->data[i].ts,obst)<=0&&timediff(nav->fcbs->data[i].te,obst)>=0){
            fcb1=nav->fcbs->data[i].bias[sat1-1];
            fcb2=nav->fcbs->data[i].bias[sat2-1];
            stat=1;
            break;
        }
    }

    if(nl_fcb1) *nl_fcb1=fcb1;
    if(nl_fcb2) *nl_fcb2=fcb2;

    return stat;
}

static int gen_sat_sd(rtk_t *rtk,const nav_t *nav,const obsd_t *obs,int n,const int *exc,int *sat1,int *sat2,int *iu,int *ir,int f,double *el)
{
    const int sat_sys[]={SYS_GPS,SYS_GAL,SYS_CMP,SYS_BD3,0};
    double elmask,el_temp[MAXOBS]={0};
    int i,j,k,m,ns=0,prn,idxs[MAXOBS]={0},sat_no[MAXOBS]={0},frq=f==-1?0:f;
    int ref_sat_idx=0,sys,nf=NF(&rtk->opt);

    elmask=MAX(rtk->opt.elmaskar,rtk->opt.elmin);

    for(i=0;sat_sys[i];i++){ /*loop each system*/
        if((sat_sys[i]&SYS_GPS)&&rtk->opt.gpsmodear==ARMODE_OFF) continue;
        if((sat_sys[i]&SYS_GAL)&&rtk->opt.galmodear==ARMODE_OFF) continue;
        if(((sat_sys[i]&SYS_CMP)||sat_sys[i]&SYS_BD3)&&rtk->opt.bdsmodear==ARMODE_OFF) continue;

        for(j=m=0;j<n;j++){
            sys=satsys(obs[j].sat,&prn);

            if(!rtk->ssat[obs[j].sat-1].vs) continue;

            if(sys==SYS_CMP){
                if(rtk->opt.arprod==AR_PROD_OSB_WHU){
                    /*no support BD2 GEO*/
                    if((prn>=1&&prn<=5)||prn>36) continue;
                }
                else if(rtk->opt.arprod==AR_PROD_FCB){
                    /*no support BD3*/
                    if(prn>18) continue;
                }
                else if(rtk->opt.arprod==AR_PROD_OSB_COM||((rtk->opt.arprod==AR_PROD_IRC)&&strcmp(rtk->opt.ac_name,"grm")==0)){
                    /*no support BDS*/
                    continue;
                }
                else if(rtk->opt.arprod==AR_PROD_IRC&&strcmp(rtk->opt.ac_name,"gbm")==0){
                    /*no support BD2 GEO*/
                    if((prn>=1&&prn<=5)||prn==15||prn==17||prn==18||prn==31||prn>60) continue;
                }
#if 1
                if(prn>18){
                    sys=SYS_BD3;
                }
#endif
            }

            if(sys!=sat_sys[i]) continue;

            if((rtk->ssat[obs[j].sat-1].slip[0]||rtk->ssat[obs[j].sat-1].slip[1])){
                rtk->ssat[obs[j].sat-1].lock[f]=-rtk->opt.minlock;
                continue;
            }
            if(exc[j]||!rtk->ssat[obs[j].sat-1].vsat[frq]||rtk->ssat[obs[j].sat-1].azel[1]<elmask
                ||rtk->ssat[obs[j].sat-1].lock[f]<0){
                continue;
            }

#if 0
            int iamb=0;
            iamb=rtk->tc?xiAmb(&rtk->opt.insopt,obs[j].sat,0):IB(obs[j].sat,0,&rtk->opt);
            if(SQRT(rtk->P[iamb+iamb*rtk->nx])>0.30){
                continue;
            }
#endif
            sat_no[m]=obs[j].sat;
            el_temp[m]=rtk->ssat[obs[j].sat-1].azel[1];
            idxs[m]=j;
            m++;
        }

        if(m<=0) continue;

        for(j=0;j<m-1;j++) for(k=j+1;k<m;k++){
            if(el_temp[j]>=el_temp[k]) continue;
            SWAP_I(sat_no[j],sat_no[k]);
            SWAP_I(idxs[j],idxs[k]);
            SWAP_D(el_temp[j],el_temp[k]);
        }
        /* generate SD referenced to max elevation angle (too simple) */
        ref_sat_idx=0;
        for(j=0;j<m;j++){
            if(j==0) continue;
            sat1[ns]=sat_no[j];
            iu[ns]=idxs[j];
            el[ns]=el_temp[j];
            sat2[ns]=sat_no[ref_sat_idx];
            ir[ns]=idxs[ref_sat_idx];
            ns++;
        }
    }

    return ns;
}

static void resetamb(rtk_t *rtk,double *xa,const double *Bc,const int *sat1,const int *sat2,int nb)
{
    int i,iamb,jamb;
    for(i=0;i<nb;i++){
        iamb=IB(sat1[i],0,&rtk->opt);
        jamb=IB(sat2[i],0,&rtk->opt);
        xa[iamb]=Bc[i]+rtk->x[jamb];
    }
}

static int SDmat(rtk_t *rtk,const obsd_t *obs,int ns,const nav_t *nav,double *H_nl,double *H_if,int *sat1,
        int *sat2,int *iu,int *ir,double *el,double *Nw,double *Bw,double *Nl,double *Nc,double *sd_nl_fcb)
{
    prcopt_t opt=rtk->opt;
    int i,j,sat,ref_sat,prn,sys_idx=-1,iamb,jamb,nb=0;
    double frq1,frq2,lam1,lam2,lam_nl,lam_wl,gamma;

    /* clear fix flag for all sats (1=float, 2=fix) */
    for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
            rtk->ssat[i].fix[j]=0;
    }

    for(i=0;i<ns;i++){
        sat=sat1[i];
        ref_sat=sat2[i];
        satsys(sat,&prn);
        sys_idx=satsysidx(sat);
        if(sys_idx==-1) continue;
        if(rtk->sdamb[sat-1].fix_nl_flag!=1) continue;

        frq1 = sat2freq(sat, obs[iu[i]].code[opt.gnss_frq_idx[sys_idx][0] - 1], nav);
        frq2 = sat2freq(sat, obs[iu[i]].code[opt.gnss_frq_idx[sys_idx][1] - 1], nav);

        lam1 = CLIGHT / frq1;
        lam2 = CLIGHT / frq2;
        lam_nl = lam1 * lam2 / (lam2 + lam1);

        iamb = IB(sat, 0, &opt);
        jamb =  IB(ref_sat, 0, &opt);
        double sd_if = rtk->x[iamb] - rtk->x[jamb];

        H_nl[iamb + nb * rtk->nx] = 1.0 / lam_nl;
        H_nl[jamb + nb * rtk->nx] = -1.0 / lam_nl;
        H_if[iamb + (rtk->na + nb) * rtk->nx] = 1.0;
        H_if[jamb + (rtk->na + nb) * rtk->nx] = -1.0;

        sat1[nb] = sat;
        sat2[nb] = ref_sat;
        iu[nb] = iu[i];
        ir[nb] = ir[i];
        el[nb] = el[i];
        Nw[nb] = rtk->sdamb[sat-1].wl;
        Bw[nb] = ROUND(rtk->sdamb[sat-1].wl);
        Nl[nb] = rtk->sdamb[sat-1].nl;
        Nc[nb] = sd_if;
        if (opt.arprod == AR_PROD_FCB||opt.arprod==AR_PROD_UPD){
            double nl_fcb1=0.0,nl_fcb2=0.0;
            if(opt.arprod==AR_PROD_FCB) matchnlfcb(obs[i].time,sat,ref_sat,&nl_fcb1,&nl_fcb2,nav);
            else if(opt.arprod==AR_PROD_UPD) matchnlupd(obs[i].time,sat,ref_sat,&nl_fcb1,&nl_fcb2,nav);

            sd_nl_fcb[nb] = nl_fcb1 - nl_fcb2;
        }
        rtk->ssat[sat1[nb]-1].fix[0]=2;
        rtk->ssat[sat2[nb]-1].fix[0]=2;
        nb++;
    }
    return nb;
}

static int resamb_nl(rtk_t *rtk,double *H_nl,double *nl_amb,int num_nl)
{
    int nb=num_nl,stat=0;
    double *Qnl,*HP,s[2]={0},*b,*pb;

    Qnl=mat(nb,nb);HP=mat(nb,rtk->nx);b=mat(nb,2);pb=zeros(nb,2);
    matmul("TN",nb,rtk->nx,rtk->nx,1.0,H_nl,rtk->P,0.0,HP);
    matmul("NN",nb,nb,rtk->nx,1.0,HP,H_nl,0.0,Qnl);

    if(!lambda(nb,2,nl_amb,Qnl,b,s)){

        rtk->sol.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
        if (rtk->sol.ratio>999.9) rtk->sol.ratio=999.9f;
        rtk->sol.thres=(float)rtk->opt.thresar[0];

        if(rtk->sol.ratio<rtk->sol.thres&&nb>MIN_AMB_RES){
            stat=0;
        }
        else if(rtk->sol.ratio>rtk->sol.thres){
            matcpy(nl_amb,b,nb,1);
            stat=1;
        }
        else stat=0;
    }
    else stat=0;

    free(Qnl);free(HP);free(b);free(pb);

    return stat?nb:0;
}

static int fix_sol(rtk_t *rtk,const obsd_t *obs,const nav_t *nav,const double *sd_nl_fcb,double *H_if,const double *Bl,
        const double *Bw,int nb,const int *sat1,const int *sat2,const int *iu,double *xa)
{
    prcopt_t opt=rtk->opt;
    int i,j,ny,na=rtk->na,sat,sys,sys_idx=-1,prn,stat=1;
    double *y,*db,*Qb_if,*Qab,*QQ,*Qy,*DP,*Bc,*dx;
    double frq1=0.0,frq2=0.0,lam1,lam2,lam_nl,gamma;

    ny = na + nb;y = mat(ny, 1);
    db = mat(nb, 1);Qb_if = mat(nb, nb);Qab = mat(na, nb);
    QQ = mat(na, nb);Qy = mat(ny, ny);DP = mat(ny, rtk->nx);
    Bc=mat(nb,1);
    dx=mat(na,1);
    matmul("TN", ny, 1, rtk->nx, 1.0, H_if, rtk->x, 0.0, y);      /*y=D'*x，  星间单差无电离层组合模糊度*/
    matmul("TN", ny, rtk->nx, rtk->nx, 1.0, H_if, rtk->P, 0.0, DP);  /*DP=D'*P*/
    matmul("NN", ny, ny, rtk->nx, 1.0, DP, H_if, 0.0, Qy);           /*Qy=DP'*D，星间单差无电离层组合模糊度方差*/

    for (i = 0; i < nb; i++) {
        for (j = 0; j < nb; j++) {
            Qb_if[i + j * nb] = Qy[na + i + (na + j) * ny];  /*SD ionosphere-free ambiguity covariance*/
        }
    }

    for (i = 0; i < na; i++) for (j = 0; j < nb; j++) Qab[i + j * na] = Qy[i + (na + j) * ny];

    for(i=0;i<na;i++){
        rtk->xa[i]=xa[i];
        for(j=0;j<na;j++) rtk->Pa[i+j*na]=rtk->P[i+j*rtk->nx];
    }

    /*fix ionospheric-free ambiguity*/
    for(i=0;i<nb;i++){
        sat=sat1[i];
        sys=satsys(sat,&prn);
        sys_idx=satsysidx(sat);
        if(sys_idx==-1) continue;
        frq1=sat2freq(sat,obs[iu[i]].code[opt.gnss_frq_idx[sys_idx][0]-1],nav);
        frq2=sat2freq(sat,obs[iu[i]].code[opt.gnss_frq_idx[sys_idx][1]-1],nav);
        lam1=CLIGHT/frq1;
        lam2=CLIGHT/frq2;
        lam_nl=lam1*lam2/(lam2+lam1);
        gamma=CLIGHT*frq2/(SQR(frq1)-SQR(frq2));

        if(opt.arprod==AR_PROD_IRC||opt.arprod==AR_PROD_OSB_WHU||opt.arprod==AR_PROD_OSB_GRM
           ||opt.arprod==AR_PROD_OSB_CNT||opt.arprod==AR_PROD_OSB_COM||opt.arprod==AR_PROD_OSB_SGG){
            if(Bl[i]!=0.0) Bc[i]=lam_nl*Bl[i]+gamma*Bw[i];
            else Bc[i]=y[na+i];

            rtk->sdamb[sat1[i]-1].lc_fix=Bc[i];
            rtk->sdamb[sat1[i]-1].lc_res=Bc[i]-rtk->sdamb[sat1[i]-1].lc;
        }
        else if(opt.arprod==AR_PROD_FCB||opt.arprod==AR_PROD_UPD){
            Bc[i]=lam_nl*(Bl[i]+sd_nl_fcb[i])+gamma*Bw[i];
        }
    }
    /*y=differences between float and fixed ionospheric-free bias*/
    for(i=0;i<nb;i++){
        y[na+i]-=Bc[i];
    }

    /*adjuset non phase-bias states and covariance using fixed slution*/
    if(!matinv(Qb_if,nb)){

        /* rtk->Pa=rtk->P-Qab*Qb^-1*Qab') */
        matmul("NN",na,nb,nb, 1.0,Qab,Qb_if ,0.0,QQ);  /* QQ = Qab*Qb^-1 */
        matmul("NT",na,na,nb,-1.0,QQ ,Qab,1.0,rtk->Pa); /* rtk->Pa = rtk->P-QQ*Qab' */

        /* rtk->xa = rtk->x-Qab*Qb^-1*(b0-b) */
        matmul("NN",nb,1,nb, 1.0,Qb_if ,y+na,0.0,db); /* db = Qb^-1*(b0-b) */
        matmul("NN",na,1,nb,-1.0,Qab,db  ,0.0,dx); /* rtk->xa = rtk->x-Qab*db */

//        matprint(1,y+na,nb,1,15,6);
//        matprint(0,Qab,na,nb,15,6);
//        matprint(1,dx,na,1,15,6);

        for(i=0;i<na;i++){
            rtk->xa[i]+=dx[i];
        }

        for(i=0;i<na;i++){
            xa[i]=rtk->xa[i];
        }

        /*reset ambiguity*/
        resetamb(rtk,xa,Bc,sat1,sat2,nb);

        if(nb>=rtk->opt.minholdsats){
            rtk->holdamb=1;
        }
    }
    else{
        trace(2,"%s(%d): Qb_if inv error, AR failed\n",time_str(obs[0].time,1),
              rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch);
        stat=0;
    }
    free(y);free(db);free(Qb_if);free(Qab);free(QQ);free(Qy);free(DP);free(Bc);free(dx);

    return stat;
}

static int pppar_IF_ILS(rtk_t *rtk,double *xa,double *bias,const obsd_t *obs,int n,const int *exc,const nav_t *nav)
{
    prcopt_t opt = rtk->opt;
    int ns = 0, nb = 0;
    int i, j, sat,prn, ref_sat, sys,sys_idx=-1, sat1[MAXOBS] = {0}, sat2[MAXOBS] = {0}, iu[MAXOBS] = {0}, ir[MAXOBS] = {0}, na = rtk->na, stat = 0;
    double frq1 = 0.0, frq2 = 0.0, lam1, lam2, lam_nl, gamma,el[MAXOBS]={0};
    double *H_nl,*H_if;
    double *Nw, *Nl, *Nc;   /*float ambiguity*/
    double *Bw, *Bl, *Bc;   /*inter ambiguity*/
    double *sd_nl_fcb, *Qnl;

    /*generate satellite SD*/
    if (!(ns = gen_sat_sd(rtk,nav, obs, n, exc, sat1, sat2, iu, ir, 0,el))) {
        return 0;
    }

    for (i = 0; i < MAXSAT; i++) {
        rtk->sdamb[i].nl =0.0;
        rtk->sdamb[i].wl = rtk->sdamb[i].nl_res = rtk->sdamb[i].wl_res = 0.0;
        rtk->sdamb[i].fix_nl_flag = rtk->sdamb[i].fix_wl_flag = 0;
        rtk->sdamb[i].nl_fix =0.0;
        rtk->sdamb[i].wl_fix = 0;
        rtk->sdamb[i].ref_sat_no=0;
    }

    H_nl = zeros(rtk->nx, ns);
    H_if = zeros(rtk->nx, rtk->nx);
    for (i = 0; i < na; i++) H_if[i + i * rtk->nx] = 1.0;

    Nw = zeros(ns, 1);
    Nl = zeros(ns, 1);
    Nc = zeros(ns, 1);
    Bw = zeros(ns, 1);
    Bl = zeros(ns, 1);
    Bc = zeros(ns, 1);
    sd_nl_fcb = zeros(ns, 1);
    Qnl=zeros(ns,ns);
    /*generate WL-IF-NL*/
    for (i = 0; i < ns; i++) {
        sat = sat1[i];
        ref_sat = sat2[i];
        satsys(sat,&prn);
        sys_idx=satsysidx(sat);
        if(sys_idx==-1) continue;

        rtk->sdamb[ref_sat - 1].ref_sat_no = 0;

        frq1 = sat2freq(sat, obs[iu[i]].code[opt.gnss_frq_idx[sys_idx][0] - 1], nav);
        frq2 = sat2freq(sat, obs[iu[i]].code[opt.gnss_frq_idx[sys_idx][1] - 1], nav);
        lam1 = CLIGHT / frq1;
        lam2 = CLIGHT / frq2;
        lam_nl = lam1 * lam2 / (lam2 + lam1);
        gamma = CLIGHT * frq2 / (SQR(frq1) - SQR(frq2));

        /*round wl ambiguity*/
        double wl_amb = 0.0;
        double sd_wl = rtk->ssat[sat - 1].mw[1] - rtk->ssat[ref_sat - 1].mw[1];
        double wl_fcb = 0.0;
        if(opt.arprod==AR_PROD_UPD){
            wl_fcb = nav->upds->wls.wl[sat-1]-nav->upds->wls.wl[ref_sat-1];
        }
        else{
            wl_fcb = nav->wlbias[sat - 1] - nav->wlbias[ref_sat - 1];
        }
        double wl_var = rtk->ssat[sat - 1].mw[3] + rtk->ssat[ref_sat - 1].mw[3];

        if (opt.arprod == AR_PROD_IRC) {
            if(strcmp(opt.ac_name,"gbm")==0){
                wl_amb = sd_wl - wl_fcb;
            }
            else if(strcmp(opt.ac_name,"grm")==0){
                wl_amb = sd_wl + wl_fcb;
            }
        } else if (opt.arprod == AR_PROD_FCB||opt.arprod==AR_PROD_UPD) {
            wl_amb = sd_wl - wl_fcb;
        } else if (opt.arprod == AR_PROD_OSB_WHU || opt.arprod == AR_PROD_OSB_GRM || opt.arprod == AR_PROD_OSB_CNT ||
                   opt.arprod == AR_PROD_OSB_COM||opt.arprod==AR_PROD_OSB_SGG) {
            wl_amb = sd_wl;
        }

        rtk->sdamb[sat - 1].wl = wl_amb;
        rtk->sdamb[sat - 1].wl_fix = newround(wl_amb);
        rtk->sdamb[sat - 1].wl_res = wl_amb - newround(wl_amb);

        /*check wl to round*/
        if (conf_func(ROUND(wl_amb), wl_amb, wl_var) < rtk->opt.thresar[1] || fabs(newround(wl_amb) - wl_amb) > 0.25) {
            rtk->sdamb[sat - 1].fix_wl_flag = 0;
            continue;
        }
        rtk->sdamb[sat - 1].fix_wl_flag = 1;
        rtk->sdamb[sat - 1].ref_sat_no = ref_sat;

        /*float if ambiguity*/
        int iamb = IB(sat, 0, &opt);
        int jamb = IB(ref_sat, 0, &opt);
        double sd_if = rtk->x[iamb] - rtk->x[jamb];

        /*nl ambiguity*/
        double nl_amb = 0.0;
        double nl_fcb1 = 0.0, nl_fcb2 = 0.0;
        if (opt.arprod == AR_PROD_IRC || opt.arprod == AR_PROD_OSB_WHU || opt.arprod == AR_PROD_OSB_GRM ||
            opt.arprod == AR_PROD_OSB_CNT || opt.arprod == AR_PROD_OSB_COM||opt.arprod==AR_PROD_OSB_SGG) {
            nl_amb = (sd_if - gamma * ROUND(wl_amb)) / lam_nl;
        } else if (opt.arprod == AR_PROD_FCB || opt.arprod==AR_PROD_UPD) {
            if (opt.arprod==AR_PROD_FCB&&!matchnlfcb(obs[i].time, sat, ref_sat, &nl_fcb1, &nl_fcb2, nav)) {
                continue;
            }
            else if (opt.arprod==AR_PROD_UPD&&!matchnlupd(obs[i].time, sat, ref_sat, &nl_fcb1, &nl_fcb2, nav)) {
                continue;
            }

            nl_amb = (sd_if - gamma * newround(wl_amb)) / lam_nl - (nl_fcb1 - nl_fcb2);
        }
        double var_nl=(rtk->P[iamb+iamb*rtk->nx]+rtk->P[jamb+jamb*rtk->nx])/SQR(lam_nl);

        rtk->sdamb[sat - 1].nl = nl_amb;
        rtk->sdamb[sat - 1].nl_fix = newround(nl_amb);
        rtk->sdamb[sat - 1].nl_res = nl_amb - newround(nl_amb);
        rtk->sdamb[sat - 1].lc = sd_if;
        rtk->sdamb[sat - 1].fix_nl_flag = 1;  /*若NL能固定，则说明这颗卫星可固定*/
    }

    nb=SDmat(rtk,obs,ns,nav,H_nl,H_if,sat1,sat2,iu,ir,el,Nw,Bw,Nl,Nc,sd_nl_fcb);
    rtk->nb_ar=nb;

    if (nb >= MIN_AMB_RES) {
        nb=resamb_nl(rtk,H_nl,Nl,nb);

        if(nb&&fix_sol(rtk,obs,nav,sd_nl_fcb,H_if,Nl,Bw,nb,sat1,sat2,iu,xa)){
            stat=1;
        }
        else{
            trace(2,"%s(%d): PPP AR failed ratio=%5.2f\n",time_str(obs[0].time,1),
                    rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,rtk->sol.ratio);
        }
    }
    else{
        trace(2,"%s(%d): no enough sd ambiguites to PPP AR nb=%d\n",time_str(obs[0].time,1),
              rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,nb);
    }

    free(H_nl);free(H_if);
    free(Nw);free(Nl);free(Nc);
    free(Bw);free(Bl);free(Bc);
    free(sd_nl_fcb);free(Qnl);

    return stat?nb:0;
}

static int pppar_UC_ILS(rtk_t *rtk,double *xa,double *bias,const obsd_t *obs,int n,const int *exc,const nav_t *nav)
{
    return resamb(rtk,bias,xa,rtk->Pa,rtk->opt.elmaskar,rtk->opt.gpsmodear,rtk->opt.glomodear,0);
}

static int arfilter(rtk_t *rtk,const obsd_t *obs,int ns,int nf)
{
    int rerun=0,dly=0,i,f;
    dly=2;

    for(i=0;i<ns;i++){
        for(f=0;f<nf;f++){
            if (rtk->ssat[obs[i].sat-1].fix[f]!=2) continue;
            if(rtk->ssat[obs[i].sat-1].lock[f]==0){
                trace(2,"%s(%d): %s is new satellite and removed in AR lock=%d\n",
                        time_str(obs[0].time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(obs[i].sat),rtk->ssat[obs[i].sat-1].lock[f]);
                rtk->ssat[obs[i].sat-1].lock[f]=-rtk->opt.minlock-dly;
                dly+=2;
                rerun=1;
            }
        }
    }
    return rerun;
}

static int pppar(rtk_t *rtk,double *bias,double *xa,double *Pa,int nf,const obsd_t *obs,int ns,const nav_t *nav,int *exc)
{
    int i,nb=0;
    double var=0.0;
    prcopt_t opt=rtk->opt;
    rtk->sol.ratio=0.0;
    float ratio1=0.0;

    if(opt.thresar[0]<1.0){
        rtk->nb_ar=0;
        return 0;
    }

    /*skip AR if position variance too high to avoid false fix*/
    int ipos=0,npos=3;
    for (i=ipos;i<ipos+npos;i++) var+=SQRT(rtk->P[i+i*rtk->nx]);
    var=var/3.0; /* maintain compatibility with previous code */
    if(var>opt.thresar[2]){
        trace(2,"%s(%d): position variance too large, var=%7.3f thres=%7.3f\n",
                time_str(obs[0].time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,var,opt.thresar[2]);
        return 0;
    }


    if(rtk->opt.ionoopt==IONOOPT_IFLC){
        nb=pppar_IF_ILS(rtk,xa,bias,obs,ns,exc,nav);
    }
    else if(rtk->opt.ionoopt==IONOOPT_UC&&rtk->opt.arprod>=AR_PROD_OSB_GRM){
        nb=pppar_UC_ILS(rtk,xa,bias,obs,ns,exc,nav);
    }

    if(rtk->opt.arfilter){
        ratio1=rtk->sol.ratio;
        if(nb>=0&&rtk->sol.prev_ratio2>=rtk->sol.thres&&((rtk->sol.ratio<rtk->sol.thres)||
                (rtk->sol.ratio<rtk->opt.thresar[0]*1.1&&rtk->sol.ratio<rtk->sol.prev_ratio1/2))){
            if(arfilter(rtk,obs,ns,nf)){
                nb=pppar_IF_ILS(rtk,xa,bias,obs,ns,exc,nav);

                if(nb==0){
                    trace(2,"%s(%d): ar filter falied\n",time_str(obs[0].time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch);
                }
                else{
                    trace(2,"%s(%d): ar filter ok\n",time_str(obs[0].time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch);
                }
            }
        }
    }
    rtk->sol.prev_ratio1=ratio1>0?ratio1:rtk->sol.ratio;
    rtk->sol.prev_ratio2=rtk->sol.ratio;

    return nb;
}

extern int manage_pppar(rtk_t *rtk,double *bias,double *xa,double *Pa,int nf,const obsd_t *obs,int ns,const nav_t *nav,int *exc)
{
    int i,f,lockc[NFREQ],ar=0,excflag=0,arsats[MAXOBS]={0},sat,nb=0;
    double ratio1=0.0;

    trace(3,"previous sample ratio=%5.2f %5.2f\n",rtk->sol.prev_ratio1,rtk->sol.prev_ratio2);
    trace(3,"num ambiguites used in last AR: %d\n",rtk->nb_ar);
    /*if no fix on previous sample and enough sats, exclude next sat in list*/
    if(rtk->sol.prev_ratio2<rtk->sol.thres&&rtk->nb_ar>=rtk->opt.mindropsats){
        for(f=0;f<nf;f++) for(i=0;i<ns;i++){
            sat=obs[i].sat;
            if(rtk->ssat[sat-1].vsat[f]&&rtk->ssat[sat-1].lock[f]>=0&&rtk->ssat[sat-1].azel[1]>=rtk->opt.elmaskar){
                arsats[ar++]=i;
            }
        }
        if(rtk->excsat<ar){
            sat=obs[arsats[rtk->excsat]].sat;
            for(f=0;f<nf;f++){
                lockc[f]=rtk->ssat[sat-1].lock[f];
                rtk->ssat[sat-1].lock[f]=-rtk->nb_ar;
            }
            excflag=1;
        }
        else rtk->excsat=0;/*exclude none and reset to beginning of list*/
    }

    if(rtk->opt.glomodear!=GLO_ARMODE_FIXHOLD||rtk->holdamb){
        nb=pppar(rtk,bias,xa,Pa,nf,obs,ns,nav,exc);
    }
    else nb=0;

    /*restore excluded sat if still no fix or significant increase in ar ratio*/
    if(excflag&&(rtk->sol.ratio<rtk->sol.thres)&&(rtk->sol.ratio<(1.5*rtk->sol.prev_ratio2))){
        sat=obs[arsats[rtk->excsat++]].sat;
        for(f=0;f<nf;f++) rtk->ssat[sat-1].lock[f]=lockc[f];
    }

    return nb;
}

extern void holdamb_ppp(rtk_t *rtk,const double *xa)
{
    double *v,*H,*R;
    int i,n,m,f,info,index[MAXSAT],nb=rtk->nx-rtk->na,nv=0,nf=NF(&rtk->opt);

    trace(3,"holdamb_ppp :\n");

    v=mat(nb,1); H=zeros(nb,rtk->nx);

    for (m=0;m<5;m++) for (f=0;f<nf;f++) {

            for (n=i=0;i<MAXSAT;i++) {
                if (!test_sys(rtk->ssat[i].sys,m)||rtk->sdamb[i].fix_nl_flag!=1||
                    rtk->ssat[i].azel[1]<rtk->opt.elmaskhold) {
                    continue;
                }
                index[n++]=IB(i+1,f,&rtk->opt);
                rtk->ssat[i].fix[f]=3; /* hold */
            }
            /* use ambiguity resolution results to generate a set of pseudo-innovations
                    to feed to kalman filter based on error between fixed and float solutions */
            for (i=1;i<n;i++) {
                /* phase-biases are single diff, so subtract errors to get
                     double diff: v(nv)=err(i)-err(0) */
                v[nv]=(xa[index[0]]-xa[index[i]])-(rtk->x[index[0]]-rtk->x[index[i]]);

                H[index[0]+nv*rtk->nx]= 1.0;
                H[index[i]+nv*rtk->nx]=-1.0;
                nv++;
            }
        }
    /* return if less than min sats for hold (skip if fix&hold for GLONASS only) */
    if (rtk->opt.modear==ARMODE_FIXHOLD&&nv<rtk->opt.minholdsats) {
        trace(3,"holdamb: not enough sats to hold ambiguity\n");
        free(v); free(H);
        return;
    }

    rtk->holdamb=1;  /* set flag to indicate hold has occurred */
    R=zeros(nv,nv);
    for (i=0;i<nv;i++) R[i+i*nv]=rtk->opt.varholdamb;

    /* update states with constraints */
    if ((info=filter(rtk->x,rtk->P,H,v,R,rtk->nx,nv,0,0,NULL,0))) {
        trace(2,"%s(%d): ppp hold filter error (info=%d)\n",time_str(rtk->sol.time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,info);
    }

    free(R);free(v); free(H);
}