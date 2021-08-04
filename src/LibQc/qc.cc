//
// Created by chenc on 2021/3/22.
//

#include "rtklib.h"

/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics==0?3:9)
#define NI(opt)     ((opt)->ionoopt!=IONOOPT_UC?0:MAXSAT)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL(opt)     ((opt)->glomodear!=GLO_ARMODE_AUTOCAL?0:NFREQGLO)
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
#define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt))
#define NX(opt)     (NR(opt)+NB(opt))

/* state variable index */
#define II(s,opt)   (NP(opt)+(s)-1)                 /* ionos (s:satellite no) */
#define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r)) /* tropos (r:0=rov,1:ref) */
#define IL(f,opt)   (NP(opt)+NI(opt)+NT(opt)+(f))   /* receiver h/w bias */
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1) /* phase bias (s:satno,f:freq) */

#define NUM_SYS 6

/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t *rtk, double xi,double *P,double var, int i)
{
    int j;

    rtk->x[i]=xi;
    for (j=0;j<rtk->nx;j++) {
        P[i+j*rtk->nx]=P[j+i*rtk->nx]=i==j?var:0.0;
    }
}

static void res_class(res_t *res,int pri)
{
    int i,j=0,k=0,type;
    if(pri){
        for(i=0;i<res->nv;i++){
            type=(res->vflag[i]>>4)&0xF;
            if(type==1){
                res->npr++;
            }
            else if(type==0){
                res->ncp++;
            }
        }
        res->pri_pr=mat(res->npr,1);
        res->pri_cp=mat(res->ncp,1);
        res->post_pr=mat(res->npr,1);
        res->post_cp=mat(res->ncp,1);
        res->pr_idx=imat(res->npr,1);
        res->cp_idx=imat(res->ncp,1);
        res->norm_pr=mat(res->npr,1);
        res->norm_cp=mat(res->ncp,1);
    }

    for(i=0;i<res->nv;i++){
        type=(res->vflag[i]>>4)&0xF;
        if(type==1){
            if(pri){
                res->pr_idx[j]=i;
                res->pri_pr[j++]=fabs(res->pri_v[i]);
            }
            else{
                res->norm_pr[j]=res->post_v[i]/(res->sigma0*SQRT(res->R[i+i*res->nv]));
                res->post_pr[j++]=res->post_v[i];
            }
        }
        else if(type==0){
            if(pri){
                res->cp_idx[k]=i;
                res->pri_cp[k++]=fabs(res->pri_v[i]);
            }
            else{
                res->norm_cp[k]=res->post_v[i]/(res->sigma0*SQRT(res->R[i+i*res->nv]));
                res->post_cp[k++]=res->post_v[i];
            }
        }
    }
}

extern void init_prires(const double *v,const int *vflag,int nv,res_t *res)
{
    res->pri_v=mat(nv,1);
    res->vflag=imat(nv,1);
    matcpy(res->pri_v,v,nv,1);
    for(int i=0;i<nv;i++){
        res->vflag[i]=vflag[i];
    }
    res->nv=nv;
    res_class(res,1);
}

extern void init_postres(rtk_t *rtk,const double *post_v,res_t *res,const double *R)
{
    res->post_v=mat(res->nv,1);
    res->R=mat(res->nv,res->nv);
    matcpy(res->post_v,post_v,res->nv,1);
    matcpy(res->R,R,res->nv,res->nv);

    res_class(res,0);

    int sat,frq;
    for(int i=0;i<res->ncp;i++){
        sat=(res->vflag[res->cp_idx[i]]>>8)&0xFF;
        frq=(res->vflag[res->cp_idx[i]]&0xF);
        rtk->ssat[sat-1].norm_v[0][frq]=res->norm_cp[i];
        rtk->ssat[sat-1].norm_v[1][frq]=res->norm_pr[i];
    }
}

extern void freeres(res_t *res)
{
    res->npr=res->ncp=res->nv=0;
    if(res->vflag){
        free(res->vflag);res->vflag=nullptr;
    }
    if(res->pri_v){
        free(res->pri_v);res->pri_v=nullptr;
    }
    if(res->post_v){
        free(res->post_v);res->post_v=nullptr;
    }
    if(res->pr_idx){
        free(res->pr_idx);res->pr_idx=nullptr;
    }
    if(res->cp_idx){
        free(res->cp_idx);res->cp_idx=nullptr;
    }
    if(res->pri_pr){
        free(res->pri_pr);res->pri_pr=nullptr;
    }
    if(res->pri_cp){
        free(res->pri_cp);res->pri_cp=nullptr;
    }
    if(res->post_pr){
        free(res->post_pr);res->post_pr=nullptr;
    }
    if(res->post_cp){
        free(res->post_cp);res->post_cp=nullptr;
    }
    if(res->norm_pr){
        free(res->norm_pr);res->norm_pr=nullptr;
    }
    if(res->norm_cp){
        free(res->norm_cp);res->norm_cp=nullptr;
    }
    if(res->R){
        free(res->R);res->R=nullptr;
    }
    if(res->Qvv){
        free(res->Qvv);res->Qvv=nullptr;
    }
}

static int resqc_igg_pr(rtk_t *rtk,res_t *res,int *exc,int ppp){
/*only check norm post pseudorange residual*/
    int frq,sat,qc_flag=0;
    double el,k0=2.80,k1=4.13,fact;

    if(rtk->opt.igg_k0!=0.0){
        k0=rtk->opt.igg_k0;
    }
    if(rtk->opt.igg_k1!=0.0){
        k1=rtk->opt.igg_k1;
    }

    double max_n_pr;
    int max_n_pr_idx= findmax(res->norm_pr, res->npr, &max_n_pr);
    sat=(res->vflag[res->pr_idx[max_n_pr_idx]]>>8)&0xFF;
    frq=(res->vflag[res->pr_idx[max_n_pr_idx]]&0xF);
    el=rtk->ssat[sat-1].azel[1]*R2D;
    if(max_n_pr>k1){
        rtk->ssat[sat-1].var_fact[1][frq]=100000.0;
//        exc[sat-1]=1;
        qc_flag=1;
        trace(2,"%s(%d): %s P%d norm residual in rejected segment el=%4.2f v=%7.3f norm_v=%7.3f var=%7.3f\n",
              time_str(rtk->sol.time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(sat),frq+1,el,
              res->post_v[res->pr_idx[max_n_pr_idx]],max_n_pr,SQRT(res->R[res->pr_idx[max_n_pr_idx]+res->pr_idx[max_n_pr_idx]*res->nv]));
    }
    else if(max_n_pr>=k0&&max_n_pr<=k1){
        fact=(max_n_pr/k0)*SQR((k1-k0)/(k1-max_n_pr));
        rtk->ssat[sat-1].var_fact[1][frq]=fact;
        qc_flag=1;
        trace(3,"%s(%d): %s P%d norm residual in reduced segment el=%4.2f v=%7.3f norm_v=%7.3f var=%7.3f fact=%7.3f\n",
              time_str(rtk->sol.time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(sat),frq+1,el,
              res->post_v[res->pr_idx[max_n_pr_idx]],max_n_pr,SQRT(res->R[res->pr_idx[max_n_pr_idx]+res->pr_idx[max_n_pr_idx]*res->nv]),fact);
    }
    else{
        rtk->ssat[sat-1].var_fact[1][frq]=1.0;
    }

    return qc_flag;
}

static int resqc_igg_cp(rtk_t *rtk,res_t *res,int *exc,int ppp){
/*only check norm post phase residual*/
    int frq,sat,qc_flag=0;
    double el,k0=2.80,k1=4.13,fact;

    if(rtk->opt.igg_k0!=0.0){
        k0=rtk->opt.igg_k0;
    }
    if(rtk->opt.igg_k1!=0.0){
        k1=rtk->opt.igg_k1;
    }

#if 1
    double max_n_cp;
    int max_n_cp_idx= findmax(res->norm_cp, res->ncp, &max_n_cp);
    sat=(res->vflag[res->cp_idx[max_n_cp_idx]]>>8)&0xFF;
    frq=(res->vflag[res->cp_idx[max_n_cp_idx]]&0xF);
    el=rtk->ssat[sat-1].azel[1]*R2D;
    if(fabs(max_n_cp)>k1){
        rtk->ssat[sat-1].var_fact[0][frq]=100000.0;
        exc[sat-1]=1;
        qc_flag=1;
        trace(2,"%s(%d): %s L%d norm residual in rejected segment el=%4.2f v=%7.3f norm_v=%7.3f var=%7.3f\n",
              time_str(rtk->sol.time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(sat),frq+1,el,
              res->post_v[res->cp_idx[max_n_cp_idx]],max_n_cp,SQRT(res->R[res->cp_idx[max_n_cp_idx]+res->cp_idx[max_n_cp_idx]*res->nv]));
    }
    else if(fabs(max_n_cp)>=k0&&fabs(max_n_cp)<=k1){
        fact=(max_n_cp/k0)*SQR((k1-k0)/(k1-max_n_cp));
        rtk->ssat[sat-1].var_fact[0][frq]=fact;
        qc_flag=1;
        trace(3,"%s(%d): %s L%d norm residual in reduced segment el=%4.2f v=%7.3f norm_v=%7.3f var=%7.3f fact=%7.3f\n",
              time_str(rtk->sol.time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(sat),frq+1,el,
              res->post_v[res->cp_idx[max_n_cp_idx]],max_n_cp,SQRT(res->R[res->cp_idx[max_n_cp_idx]+res->cp_idx[max_n_cp_idx]*res->nv]),fact);
    }
    else{
        rtk->ssat[sat-1].var_fact[0][frq]=1.0;
    }
#endif
    return qc_flag;
}

static int resqc_igg(rtk_t *rtk,res_t *res,int *exc,int ppp){
    int frq,sat,qc_flag=0;
    double el,k0=2.80,k1=4.13,fact;

    if(rtk->opt.igg_k0!=0.0){
        k0=rtk->opt.igg_k0;
    }
    if(rtk->opt.igg_k1!=0.0){
        k1=rtk->opt.igg_k1;
    }

    double max_n_pr;
    int max_n_pr_idx= findmax(res->norm_pr, res->npr, &max_n_pr);
    sat=(res->vflag[res->pr_idx[max_n_pr_idx]]>>8)&0xFF;
    frq=(res->vflag[res->pr_idx[max_n_pr_idx]]&0xF);
    el=rtk->ssat[sat-1].azel[1]*R2D;
    if(max_n_pr>k1){
        rtk->ssat[sat-1].var_fact[1][frq]=100000.0;
        qc_flag=1;
        trace(2,"%s(%d): %s P%d norm residual in rejected segment el=%4.2f v=%7.3f norm_v=%7.3f var=%7.3f\n",
              time_str(rtk->sol.time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(sat),frq+1,el,
              res->post_v[res->pr_idx[max_n_pr_idx]],max_n_pr,SQRT(res->R[res->pr_idx[max_n_pr_idx]+res->pr_idx[max_n_pr_idx]*res->nv]));
    }
    else{
        double max_n_cp;
        int max_n_cp_idx= findmax(res->norm_cp, res->ncp, &max_n_cp);
        sat=(res->vflag[res->cp_idx[max_n_cp_idx]]>>8)&0xFF;
        frq=(res->vflag[res->cp_idx[max_n_cp_idx]]&0xF);
        el=rtk->ssat[sat-1].azel[1]*R2D;
        if(fabs(max_n_cp)>k1){
            rtk->ssat[sat-1].var_fact[0][frq]=100000.0;
            qc_flag=1;
            trace(2,"%s(%d): %s L%d norm residual in rejected segment el=%4.2f v=%7.3f norm_v=%7.3f var=%7.3f\n",
                  time_str(rtk->sol.time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(sat),frq+1,el,
                  res->post_v[res->cp_idx[max_n_cp_idx]],max_n_cp,SQRT(res->R[res->cp_idx[max_n_cp_idx]+res->cp_idx[max_n_cp_idx]*res->nv]));
        }
        else if(fabs(max_n_cp)>=k0&&fabs(max_n_cp)<=k1){
            fact=(max_n_cp/k0)*SQR((k1-k0)/(k1-max_n_cp));
            rtk->ssat[sat-1].var_fact[0][frq]=fact;
            qc_flag=1;
            trace(3,"%s(%d): %s L%d norm residual in reduced segment el=%4.2f v=%7.3f norm_v=%7.3f var=%7.3f fact=%7.3f\n",
                  time_str(rtk->sol.time,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(sat),frq+1,el,
                  res->post_v[res->cp_idx[max_n_cp_idx]],max_n_cp,SQRT(res->R[res->cp_idx[max_n_cp_idx]+res->cp_idx[max_n_cp_idx]*res->nv]),fact);
        }
        else{
            rtk->ssat[sat-1].var_fact[0][frq]=1.0;
        }
    }
    return qc_flag;
}

static int resqc_shi(rtk_t *rtk,res_t *res,int *exc,int ppp){
    int frq,sat,qc_flag=0,iamb;
    double el;

    /*STEP 1: check post pseudorange residual*/
    double max_pr;
    int max_pr_idx=findmax(res->post_pr,res->npr,&max_pr);
    sat=(res->vflag[res->pr_idx[max_pr_idx]]>>8)&0xFF;
    frq=(res->vflag[res->pr_idx[max_pr_idx]]&0xF);
    el=rtk->ssat[sat-1].azel[1]*R2D;
    if(fabs(max_pr)>3.0/sin(el*D2R)){
        exc[sat-1]=1;
        qc_flag=1;
        return qc_flag;
    }

    /*STEP 2: check norm pseudorange residual*/
    double max_n_pr;
    int max_n_pr_idx=findmax(res->norm_pr,res->npr,&max_n_pr);
    sat=(res->vflag[res->pr_idx[max_n_pr_idx]]>>8)&0xFF;
    frq=(res->vflag[res->pr_idx[max_n_pr_idx]]&0xF);
    el=rtk->ssat[sat-1].azel[1]*R2D;
    if(fabs(max_n_pr)>2.0){
        exc[sat-1]=1;
        qc_flag=1;
        return qc_flag;
    }

    /*STEP 3: check post phase residual*/
    double max_cp;
    int max_cp_idx=findmax(res->post_cp,res->ncp,&max_cp);
    sat=(res->vflag[res->cp_idx[max_cp_idx]]>>8)&0xFF;
    frq=(res->vflag[res->cp_idx[max_cp_idx]]&0xF);
    el=rtk->ssat[sat-1].azel[1]*R2D;
    iamb=(ppp?iamb_ppp(&rtk->opt,sat,frq):iamb_ppk(&rtk->opt,sat,frq));
    if(fabs(max_cp)>0.03/sin(el*D2R)){
        initx(rtk,rtk->x[iamb],rtk->P,SQR(rtk->opt.std[0]),iamb);
        rtk->ssat[sat-1].init_amb[frq]=1;
        qc_flag=1;
        return qc_flag;
    }

    /*STEP 4: check norm post phase residual*/
    double max_n_cp;
    int max_n_cp_idx=findmax(res->norm_cp,res->ncp,&max_n_cp);
    sat=(res->vflag[res->cp_idx[max_n_cp_idx]]>>8)&0xFF;
    frq=(res->vflag[res->cp_idx[max_n_cp_idx]]&0xF);
    el=rtk->ssat[sat-1].azel[1]*R2D;
    iamb=(ppp?iamb_ppp(&rtk->opt,sat,frq):iamb_ppk(&rtk->opt,sat,frq));
    if(fabs(max_n_cp)>2.0){
        initx(rtk,rtk->x[iamb],rtk->P,SQR(rtk->opt.std[0]),iamb);
        rtk->ssat[sat-1].init_amb[frq]=1;
        qc_flag=1;
        return qc_flag;
    }
    return 0;
}

static int resqc_zhao(gtime_t t,rtk_t *rtk,const double *v,const double *norm_v,const int *v_flag,int nv,int *exc){
    int i,j,k,m=0,n=0,type,frq,nf=NF(&rtk->opt),satj;
    double vv_p[MAXOBS*2]={0},vv_c[MAXOBS*2]={0},vv_np[MAXOBS*2]={0},vv_nc[MAXOBS*2]={0};
    int vv_sat_p[MAXOBS]={0},vv_sat_c[MAXOBS]={0},vv_frq_p[MAXOBS]={0},vv_frq_c[MAXOBS]={0};
    double vv_el_p[MAXOBS]={0},vv_el_c[MAXOBS]={0};
    double max_vp,max_vc,max_np,max_nc,max;
    int max_idx_vp,max_idx_vc,max_idx_np,max_idx_nc;
    double K0=1.5,K1=2.5;
    int qc_flag=0;

    for(i=0,j=0,k=0;i<nv;i++){
        satj=(v_flag[i]>>8)&0xFF;
        type=(v_flag[i]>>4)&0xF;
        frq=(v_flag[i]&0xF);

        if(type==1){
            vv_p[j]=fabs(v[i]);
            vv_np[j]=fabs(norm_v[i]);
            vv_frq_p[j]=frq;
            vv_el_p[j]=rtk->ssat[satj-1].azel[1];
            vv_sat_p[j++]=satj;
        }
        else if(type==0){
            vv_c[k]=fabs(v[i]);
            vv_nc[k]=fabs(norm_v[i]);
            vv_frq_c[k]=frq;
            vv_el_c[k]=rtk->ssat[satj-1].azel[1];
            vv_sat_c[k++]=satj;
        }
    }
    int sat=0;
    double el=0.0;
    double fact=1.0;
    /* pseudorange t-test*/
    double avp_vnp=0.0,delta_avp_vnp=0.0;
    for(i=0;i<j;i++){
        avp_vnp+=vv_np[i];
    }
    avp_vnp/=j;
    for(i=0;i<j;i++){
        delta_avp_vnp+=SQR(vv_np[i]-avp_vnp);
    }

    double thres1=tdistb_0050[j],thres2=tdistb_0010[j];  ///0050 and 0010 for cpt
    double T_p[MAXOBS*2]={0};
    for(i=0;i<j;i++){
        sat=vv_sat_p[i];
        el=vv_el_p[i];
        T_p[i]=fabs(vv_np[i]-avp_vnp)/SQRT(delta_avp_vnp/j);
        if(T_p[i]>thres1&&T_p[i]<thres2){
            fact=thres1/T_p[i]*SQR((thres2-T_p[i])/(thres2-thres1));
            rtk->ssat[sat-1].var_fact[0][frq]=1.0/fact;
            qc_flag=1;
            trace(2,"%s(%d): %s down weight by larger pseudorange residual el=%4.1f Tp=%5.2f thres1=%5.2f thres2=%5.2f\n",
                  time_str(t,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(sat),el*R2D,T_p[i],thres1,thres2);
        }
        else if(T_p[i]>thres2){
//            rtk->ssat[sat-1].var_fact=10000.0;
            exc[sat-1]=1;
            qc_flag=1;
            trace(2,"%s(%d): %s exclude by larger pseudorange residual el=%4.1f Tp=%5.2f thres=%5.2f\n",
                  time_str(t,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(sat),el*R2D,T_p[i],thres2);
        }
        else{
            rtk->ssat[sat-1].var_fact[0][frq]=1.0;
            qc_flag=0;
        }
    }

    /* phase t-test */
    double avp_vnc=0.0,delta_avp_vnc=0.0;
    for(i=0;i<k;i++){
        avp_vnc+=vv_nc[i];
    }
    avp_vnc/=k;
    for(i=0;i<k;i++){
        delta_avp_vnc+=SQR(vv_nc[i]-avp_vnc);
    }
    thres1=tdistb_0050[k],thres2=tdistb_0010[k];
    double T_c[MAXOBS*2]={0};
    for(i=0;i<k;i++){
        sat=vv_sat_c[i];
        el=vv_el_c[i];
        T_c[i]=fabs(vv_nc[i]-avp_vnc)/SQRT(delta_avp_vnc/k);
        if(T_c[i]>thres1&&T_c[i]<thres2){
            fact=thres1/T_c[i]*SQR((thres2-T_c[i])/(thres2-thres1));
            rtk->ssat[sat-1].var_fact[0][frq]=1.0/fact;
            qc_flag=1;
            trace(2,"%s(%d): %s down weight by larger phase residual el=%4.1f Tc=%5.2f thres1=%5.2f thres2=%5.2f\n",
                  time_str(t,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(sat),el*R2D,T_c[i],thres1,thres2);
        }
        else if(T_c[i]>thres2){
            exc[sat-1]=1;
            qc_flag=1;
            trace(2,"%s(%d): %s exclude by larger phase residual el=%4.1f Tp=%5.2f thres=%5.2f\n",
                  time_str(t,1),rtk->tc?rtk->ins_kf->couple_epoch:rtk->epoch,sat_id(sat),el*R2D,T_c[i],thres2);
        }
    }
    return qc_flag;
}

extern int resqc(gtime_t t,rtk_t *rtk,res_t *res,int *exc,int ppp){
    if(rtk->opt.robust==ROBUST_QC_OFF){
        return 0;
    }
    else if(rtk->opt.robust==ROBUST_QC_IGG_PR){
        return resqc_igg_pr(rtk,res,exc,ppp);
    }
    else if(rtk->opt.robust==ROBUST_QC_IGG_CP){
        return resqc_igg_cp(rtk,res,exc,ppp);
    }
    else if(rtk->opt.robust==ROBUST_QC_IGG){
        return resqc_igg(rtk,res,exc,ppp);
    }
    else if(rtk->opt.robust==ROBUST_QC_SHI){
        return resqc_shi(rtk,res,exc,ppp);
    }
    else return 0;
}

extern int pri_res_check(gtime_t t,rtk_t *rtk,const double *pri_v,const int *vflag,int nv,int *exc)
{
    double v_SYS[NUM_SYS][MAXOBS]={0},v_copy[MAXOBS]={0},mean,thres=20.0;
    int nv_SYS[NUM_SYS]={0},sat_SYS[NUM_SYS][MAXOBS]={0};
    int i,j,sat,sys_idx,qc_flag=0;
    int ppk=0,ppp=0;
    prcopt_t *popt=&rtk->opt;

    ppk= ((popt->mode >= PMODE_DGPS && popt->mode <= PMODE_STATIC_START) ||
           (popt->mode == PMODE_TC_DGPS || popt->mode == PMODE_TC_PPK) ||
           (popt->mode == PMODE_LC_DGPS || popt->mode == PMODE_LC_PPK) ||
           (popt->insopt.imu_align == INS_ALIGN_GNSS_PPK || popt->insopt.imu_align == INS_ALIGN_GNSS_DGPS));
    ppp = ((popt->mode >= PMODE_PPP_KINEMA && popt->mode <= PMODE_PPP_FIXED)|| (popt->mode == PMODE_TC_PPP||popt->mode==PMODE_LC_PPP));

    if(ppk){
        thres=3.0;
    }
    if(ppp){
        thres=100.0;
    }

    for(j=0;j<nv;j++){
        sat=(vflag[j]>>8)&0xFF;
        sys_idx=satsysidx(sat);
        sat_SYS[sys_idx][nv_SYS[sys_idx]]=sat;
        v_SYS[sys_idx][nv_SYS[sys_idx]++]=pri_v[j];
    }

    for(j=0;j<NUM_SYS;j++){
        if(nv_SYS[j]<=0) continue;
        matcpy(v_copy,v_SYS[j],nv_SYS[j],1);
        mean=median(v_copy,nv_SYS[j]);
        for(i=0;i<nv_SYS[j];i++){
            if(fabs(v_SYS[j][i]-mean)>thres){
                sat=sat_SYS[j][i];
                exc[sat-1]=1;
                qc_flag=1;
            }
        }
    }

    return qc_flag;
}