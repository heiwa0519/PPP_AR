/*------------------------------------------------------------------------------
* preceph.c : precise ephemeris and clock functions
*
*          Copyright (C) 2007-2017 by T.TAKASU, All rights reserved.
*
* references :
*     [1] S.Hilla, The Extended Standard Product 3 Orbit Format (SP3-c),
*         12 February, 2007
*     [2] J.Ray, W.Gurtner, RINEX Extensions to Handle Clock Information,
*         27 August, 1998
*     [3] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*     [4] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*         Space Technology Library, 2004
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2009/01/18 1.0  new
*           2009/01/31 1.1  fix bug on numerical error to read sp3a ephemeris
*           2009/05/15 1.2  support glonass,galileo,qzs
*           2009/12/11 1.3  support wild-card expansion of file path
*           2010/07/21 1.4  added api:
*                               eci2ecef(),sunmoonpos(),peph2pos(),satantoff(),
*                               readdcb()
*                           changed api:
*                               readsp3()
*                           deleted api:
*                               eph2posp()
*           2010/09/09 1.5  fix problem when precise clock outage
*           2011/01/23 1.6  support qzss satellite code
*           2011/09/12 1.7  fix problem on precise clock outage
*                           move sunmmonpos() to rtkcmn.c
*           2011/12/01 1.8  modify api readsp3()
*                           precede later ephemeris if ephemeris is NULL
*                           move eci2ecef() to rtkcmn.c
*           2013/05/08 1.9  fix bug on computing std-dev of precise clocks
*           2013/11/20 1.10 modify option for api readsp3()
*           2014/04/03 1.11 accept extenstion including sp3,eph,SP3,EPH
*           2014/05/23 1.12 add function to read sp3 velocity records
*                           change api: satantoff()
*           2014/08/31 1.13 add member cov and vco in peph_t sturct
*           2014/10/13 1.14 fix bug on clock error variance in peph2pos()
*           2015/05/10 1.15 add api readfcb()
*                           modify api readdcb()
*           2017/04/11 1.16 fix bug on antenna offset correction in peph2pos()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define SQR(x)      ((x)*(x))

#define NMAX        10              /* order of polynomial interpolation */
#define MAXDTE      900.0           /* max time difference to ephem time (s) */
#define EXTERR_CLK  1E-3            /* extrapolation error for clock (m/s) */
#define EXTERR_EPH  5E-7            /* extrapolation error for ephem (m/s^2) */

typedef struct {
    int sat;
    gtime_t ts,te;
    int code;
    int type;
    double bia;
}bias_t;

typedef struct {
    int nb,nmax;
    bias_t *data;
}biases_t;

/* satellite code to satellite system ----------------------------------------*/
static int code2sys(char code)
{
    if (code=='G'||code==' ') return SYS_GPS;
    if (code=='R') return SYS_GLO;
    if (code=='E') return SYS_GAL; /* extension to sp3-c */
    if (code=='J') return SYS_QZS; /* extension to sp3-c */
    if (code=='C') return SYS_CMP; /* extension to sp3-c */
    if (code=='L') return SYS_LEO; /* extension to sp3-c */
    return SYS_NONE;
}
/* read sp3 header -----------------------------------------------------------*/
static int readsp3h(FILE *fp, gtime_t *time, char *type, int *sats,
                    double *bfact, char *tsys)
{
    int i=0,j,k=0,ns=0,nl=5,sys,prn;
    char buff[1024];
    
    trace(3,"readsp3h:\n");
    
    while (fgets(buff,sizeof(buff),fp)) {

        if (buff[0]=='#'&&(buff[1]=='c'||buff[1]=='d')) {
            *type=buff[2];
            if (str2time(buff,3,28,time)) return 0;
        }
        else if (buff[0]=='+'&&buff[1]==' ') {
            if (i==2) {
                ns=(int)str2num(buff,3,3);
                if (ns>85) nl=ns/17+(ns%17!=0);
            }
            for (j=0;j<17&&k<ns;j++) {
                sys=code2sys(buff[9+3*j]);
                prn=(int)str2num(buff,10+3*j,2);
                if (k<MAXSAT) sats[k++]=satno(sys,prn);
            }
        }
        else if (i==2*nl+2) {/* %c */
            strncpy(tsys,buff+9,3); tsys[3]='\0';
        }
        else if (i==2*nl+4) {/* %f */
            bfact[0]=str2num(buff, 3,10);
            bfact[1]=str2num(buff,14,12);
        }
        else if (i==2*nl+11){
            break; /* at end of header */
        }
        i=i+1; /* line counter */
    }
    return ns;
}
/* add precise ephemeris -----------------------------------------------------*/
static int addpeph(nav_t *nav, peph_t *peph)
{
    peph_t *nav_peph;
    
    if (nav->ne>=nav->nemax) {
        nav->nemax+=256;
        if (!(nav_peph=(peph_t *)realloc(nav->peph,sizeof(peph_t)*nav->nemax))) {
            trace(1,"readsp3b malloc error n=%d\n",nav->nemax);
            free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
            return 0;
        }
        nav->peph=nav_peph;
    }
    nav->peph[nav->ne++]=*peph;
    return 1;
}
/* read sp3 body -------------------------------------------------------------*/
static void readsp3b(FILE *fp, char type, int *sats, int ns, double *bfact,
                     char *tsys, int index, int opt, nav_t *nav)
{
    peph_t peph;
    gtime_t time;
    double val,std,base;
    int i,j,sat,sys,prn,n=ns*(type=='P'?1:2),pred_o,pred_c,v;
    char buff[1024];
    
    trace(3,"readsp3b: type=%c ns=%d index=%d opt=%d\n",type,ns,index,opt);
    
    while (fgets(buff,sizeof(buff),fp)) {
        if(buff[0]=='%'||buff[0]=='/') continue;
        if (!strncmp(buff,"EOF",3)) break;
        
        if (buff[0]!='*'||str2time(buff,3,28,&time)) {
            trace(2,"sp3 invalid epoch %31.31s\n",buff);
            continue;
        }
        if (!strcmp(tsys,"UTC")) time=utc2gpst(time); /* utc->gpst */
        peph.time =time;
        peph.index=index;
        
        for (i=0;i<MAXSAT;i++) {
            for (j=0;j<4;j++) {
                peph.pos[i][j]=0.0;
                peph.std[i][j]=0.0f;
                peph.vel[i][j]=0.0;
                peph.vst[i][j]=0.0f;
            }
            for (j=0;j<3;j++) {
                peph.cov[i][j]=0.0f;
                peph.vco[i][j]=0.0f;
            }
        }
        for (i=pred_o=pred_c=v=0;i<n&&fgets(buff,sizeof(buff),fp);i++) {
            
            if (strlen(buff)<4||(buff[0]!='P'&&buff[0]!='V')) continue;
            
            sys=buff[1]==' '?SYS_GPS:code2sys(buff[1]);
            prn=(int)str2num(buff,2,2);
            if      (sys==SYS_SBS) prn+=100;
            else if (sys==SYS_QZS) prn+=192; /* extension to sp3-c */
            
            if (!(sat=satno(sys,prn))) continue;
            
            if (buff[0]=='P') {
                pred_c=strlen(buff)>=76&&buff[75]=='P';
                pred_o=strlen(buff)>=80&&buff[79]=='P';
            }
            for (j=0;j<4;j++) {
                
                /* read option for predicted value */
                if (j< 3&&(opt&1)&& pred_o) continue;
                if (j< 3&&(opt&2)&&!pred_o) continue;
                if (j==3&&(opt&1)&& pred_c) continue;
                if (j==3&&(opt&2)&&!pred_c) continue;
                
                val=str2num(buff, 4+j*14,14);
                std=str2num(buff,61+j* 3,j<3?2:3);
                
                if (buff[0]=='P') { /* position */
                    if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
                        peph.pos[sat-1][j]=val*(j<3?1000.0:1E-6);
                        v=1; /* valid epoch */
                    }
                    if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
                        peph.std[sat-1][j]=(float)(pow(base,std)*(j<3?1E-3:1E-12));
                    }
                }
                else if (v) { /* velocity */
                    if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
                        peph.vel[sat-1][j]=val*(j<3?0.1:1E-10);
                    }
                    if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
                        peph.vst[sat-1][j]=(float)(pow(base,std)*(j<3?1E-7:1E-16));
                    }
                }
            }
        }
        if (v) {
            if (!addpeph(nav,&peph)) return;
        }
    }
}
/* compare precise ephemeris -------------------------------------------------*/
static int cmppeph(const void *p1, const void *p2)
{
    peph_t *q1=(peph_t *)p1,*q2=(peph_t *)p2;
    double tt=timediff(q1->time,q2->time);
    return tt<-1E-9?-1:(tt>1E-9?1:q1->index-q2->index);
}
/* combine precise ephemeris -------------------------------------------------*/
static void combpeph(nav_t *nav, int opt)
{
    int i,j,k,m;
    
    trace(3,"combpeph: ne=%d\n",nav->ne);
    
    qsort(nav->peph,nav->ne,sizeof(peph_t),cmppeph);
    
    if (opt&4) return;
    
    for (i=0,j=1;j<nav->ne;j++) {
        
        if (fabs(timediff(nav->peph[i].time,nav->peph[j].time))<1E-9) {
            
            for (k=0;k<MAXSAT;k++) {
                if (norm(nav->peph[j].pos[k],4)<=0.0) continue;
                for (m=0;m<4;m++) nav->peph[i].pos[k][m]=nav->peph[j].pos[k][m];
                for (m=0;m<4;m++) nav->peph[i].std[k][m]=nav->peph[j].std[k][m];
                for (m=0;m<4;m++) nav->peph[i].vel[k][m]=nav->peph[j].vel[k][m];
                for (m=0;m<4;m++) nav->peph[i].vst[k][m]=nav->peph[j].vst[k][m];
            }
        }
        else if (++i<j) nav->peph[i]=nav->peph[j];
    }
    nav->ne=i+1;
    
    trace(4,"combpeph: ne=%d\n",nav->ne);
}
/* read sp3 precise ephemeris file ---------------------------------------------
* read sp3 precise ephemeris/clock files and set them to navigation data
* args   : char   *file       I   sp3-c precise ephemeris file
*                                 (wind-card * is expanded)
*          nav_t  *nav        IO  navigation data
*          int    opt         I   options (1: only observed + 2: only predicted +
*                                 4: not combined)
* return : none
* notes  : see ref [1]
*          precise ephemeris is appended and combined
*          nav->peph and nav->ne must by properly initialized before calling the
*          function
*          only files with extensions of .sp3, .SP3, .eph* and .EPH* are read
*-----------------------------------------------------------------------------*/
extern void readsp3(const char *file, nav_t *nav, int opt)
{
    FILE *fp;
    gtime_t time={0};
    double bfact[2]={0};
    int i,j,n,ns,sats[MAXSAT]={0};
    char *efiles[MAXEXFILE],*ext,type=' ',tsys[4]="";
    
    trace(3,"readpephs: file=%s\n",file);
    
    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return;
        }
    }
    /* expand wild card in file path */
    n=expath(file,efiles,MAXEXFILE);
    
    for (i=j=0;i<n;i++) {
        if (!(ext=strrchr(efiles[i],'.'))) continue;
        
        if (!strstr(ext+1,"sp3")&&!strstr(ext+1,"SP3")&&
            !strstr(ext+1,"eph")&&!strstr(ext+1,"EPH")) continue;
        
        if (!(fp=fopen(efiles[i],"r"))) {
            trace(2,"sp3 file open error %s\n",efiles[i]);
            continue;
        }
        /* read sp3 header */
        ns=readsp3h(fp,&time,&type,sats,bfact,tsys);
        
        /* read sp3 body */
        readsp3b(fp,type,sats,ns,bfact,tsys,j++,opt,nav);
        
        fclose(fp);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);
    
    /* combine precise ephemeris */
    if (nav->ne>0) combpeph(nav,opt);
}
/* read satellite antenna parameters -------------------------------------------
* read satellite antenna parameters
* args   : char   *file       I   antenna parameter file
*          gtime_t time       I   time
*          nav_t  *nav        IO  navigation data
* return : status (1:ok,0:error)
* notes  : only support antex format for the antenna parameter file
*-----------------------------------------------------------------------------*/
extern int readsap(const char *file, gtime_t time, nav_t *nav)
{
    pcvs_t pcvs={0};
    pcv_t pcv0={0},*pcv;
    int i;
    
    trace(3,"readsap : file=%s time=%s\n",file,time_str(time,0));
    
    if (!readpcv(file,&pcvs)) return 0;
    
    for (i=0;i<MAXSAT;i++) {
        pcv=searchpcv(i+1,"",time,&pcvs);
        nav->pcvs[i]=pcv?*pcv:pcv0;
    }
    free(pcvs.pcv);
    return 1;
}

/* add satellite fcb ---------------------------------------------------------*/
static void addfcb(fcbs_t *fcbs,const fcbd_t *fcb)
{
    fcbd_t *temp_fcbs;

    if(fcbs->nmax<=fcbs->n){
        fcbs->nmax+=256;
        if(!(temp_fcbs=(fcbd_t *)realloc(fcbs->data, sizeof(fcbd_t)*fcbs->nmax))){
            fprintf(stdout,"addfcb: memory allocation error\n");
            free(fcbs->data);fcbs->data=NULL;
            fcbs->n=fcbs->nmax=0;
            return;
        }
        fcbs->data=temp_fcbs;
    }
    fcbs->data[fcbs->n++]=*fcb;
}

static void readfcbhead(FILE *fp,nav_t *nav)
{
    double bias;
    char buff[1024],satid[8];
    int satno;
    while(fgets(buff, sizeof(buff),fp)){
        if (strlen(buff)<=60) continue;
        if(!strncmp(buff,"WL",2)){
            strncpy(satid,buff+4,3);
            satno=satid2no(satid);
            if(satno>0){
                bias=str2num(buff,14,6);
                nav->wlbias[satno-1]=bias;
            }
        }
        else if(strstr(buff,"END OF HEADER")){
            return;
        }
    }
}


/* read satellite fcb file ---------------------------------------------------*/
static int readfcbf(const char *file, nav_t *nav)
{
    FILE *fp;
    gtime_t ts;
    char buff[1024];

    trace(3,"readfcbf: file=%s\n",file);
    
    if (!(fp=fopen(file,"r"))) {
        trace(2,"fcb parameters file open error: %s\n",file);
        return 0;
    }

    nav->fcbs=(void *)malloc(sizeof(fcbs_t));
    nav->fcbs->n=nav->fcbs->nmax=0;
    nav->fcbs->data=NULL;
    readfcbhead(fp,nav);

    char sat_id[8];
    int sat_no;
    static const  fcbd_t fcb0={0};
    while(fgets(buff, sizeof(buff),fp)){
        if(buff[0]=='*'&&!str2time(buff,1,28,&ts)){
            addfcb(nav->fcbs,&fcb0);
            while(fgets(buff, sizeof(buff),fp)&&buff[0]=='P'){
                strncpy(sat_id,buff+1,3);
                sat_no=satid2no(sat_id);
                if(sat_no<=0){
                    if(!fgets(buff, sizeof(buff),fp)){
                        break;
                    }
                    continue;
                }
                nav->fcbs->data[nav->fcbs->n-1].bias[sat_no-1]=str2num(buff,24,6);
            }
            nav->fcbs->data[nav->fcbs->n-1].ts=ts;
            nav->fcbs->data[nav->fcbs->n-1].te=timeadd(ts,60.0*15);
            str2time(buff,1,28,&ts);
        }
        else if(buff[0]=='P'){
            addfcb(nav->fcbs,&fcb0);
            while(buff[0]=='P'){
                strncpy(sat_id,buff+1,3);
                sat_no=satid2no(sat_id);
                if(sat_no<=0){
                    if(!fgets(buff, sizeof(buff),fp)){
                        break;
                    }
                    continue;
                }
                nav->fcbs->data[nav->fcbs->n-1].bias[sat_no-1]=str2num(buff,24,6);
                if(!fgets(buff, sizeof(buff),fp)){
                    break;
                }
            }
            nav->fcbs->data[nav->fcbs->n-1].ts=ts;
            nav->fcbs->data[nav->fcbs->n-1].te=timeadd(ts,60.0*15);
            str2time(buff,1,28,&ts);
        }
    }

    fclose(fp);
    return 1;
}

/* read satellite fcb data -----------------------------------------------------
* read satellite fractional cycle bias (dcb) parameters
* args   : char   *file       I   fcb parameters file (wild-card * expanded)
*          nav_t  *nav        IO  navigation data
* return : status (1:ok,0:error)
* notes  : fcb data appended to navigation data
*-----------------------------------------------------------------------------*/
extern int readfcb(const char *file, nav_t *nav)
{
    char *efiles[MAXEXFILE]={0};
    int i,n;
    
    trace(3,"readfcb : file=%s\n",file);
    
    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return 0;
        }
    }
    n=expath(file,efiles,MAXEXFILE);
    
    for (i=0;i<n;i++) {
        readfcbf(efiles[i],nav);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);
    
    return 1;
}

static int biasstr2time(const char *s, int i, int n, gtime_t *t) {
    double ep[6] = {};
    ep[1] = 1.0;
    ep[2] = 1.0;
    double day, sec;
    char str[256], *p = str;
    if (i < 0 || (int) strlen(s) < i || (int) sizeof(str) - 1 < i) return -1;
    for (s += i; *s && --n >= 0;) *p++ = *s++;
    *p = '\0';
    if (sscanf(str, "%lf:%lf:%lf", ep, &day, &sec) < 3)
        return -1;
    if (ep[0] < 100.0) ep[0] += ep[0] < 80.0 ? 2000.0 : 1900.0;
    *t = timeadd(epoch2time(ep), 86400.0 * (day - 1.0) + sec);
    return 0;
}

static int readosbf(const char *file,biases_t *sat_bias)
{
    FILE *fp;
    bias_t *bias_temp;
    char buff[200];
    int i,code,sat,sys;
    const char *obscodes[]={       /* observation code strings */
            ""  ,"1C","1P","1W","1Y", "1M","1N","1S","1L","1E", /*  0- 9 */
            "1A","1B","1X","1Z","2C", "2D","2S","2L","2X","2P", /* 10-19 */
            "2W","2Y","2M","2N","5I", "5Q","5X","7I","7Q","7X", /* 20-29 */
            "6A","6B","6C","6X","6Z", "6S","6L","8L","8Q","8X", /* 30-39 */
            "2I","2Q","6I","6Q","3I", "3Q","3X","1I","1Q","5A", /* 40-49 */
            "5B","5C","9A","9B","9C", "9X","1D","5D","5P","5Z", /* 50-59 */
            "6E","7D","7P","7Z","8D", "8P","4A","4B","4X",""    /* 60-69 */
    };

    if (!(fp=fopen(file,"r"))) {
        trace(2,"fcb parameters file open error: %s\n",file);
        return 0;
    }

    memset(sat_bias,0, sizeof(biases_t));
    while(fgets(buff, sizeof(buff),fp)){
        if ((!strncmp(buff + 1, "OSB", 3)) && (!strncmp(buff + 65, "ns", 2))) {
            sat = satid2no(buff + 11);
            if(sat<=0) continue;
            int range = (buff[25] == 'C' ? 1 : 0);
            for (i = 0; i < MAXCODE; i++){
                if (!strncmp(obscodes[i], buff + 26, 2)){
                    code=i;
                    break;
                }
            }
            gtime_t t1, t2;
            char st1[20], st2[20];
            if (biasstr2time(buff, 35, 14, &t1)) continue;
            if (biasstr2time(buff, 50, 14, &t2)) continue;
            time2str(t1, st1, 0);
            time2str(t2, st2, 0);
            double value = atof(buff + 70);
            if (sat_bias->nb >= sat_bias->nmax) {
                sat_bias->nmax += 1024;
                if (!(bias_temp=(bias_t *)realloc(sat_bias->data,sizeof(bias_t) * (sat_bias->nmax)))) {
                    free(sat_bias->data);
                    sat_bias->data = NULL;
                    sat_bias->nb = sat_bias->nmax = 0;
                    return -1;
                }
                sat_bias->data =bias_temp ;
            }
            sat_bias->nb++;
            sat_bias->data[sat_bias->nb - 1].ts = t1;
            sat_bias->data[sat_bias->nb - 1].te = t2;
            sat_bias->data[sat_bias->nb - 1].sat = sat;
            sat_bias->data[sat_bias->nb - 1].code = code;
            sat_bias->data[sat_bias->nb - 1].type = range;
            sat_bias->data[sat_bias->nb - 1].bia = value;
        }
    }

    fclose(fp);

    return 1;
}

extern int  readosb(const char *file, nav_t *nav)
{
    biases_t biases={0};

    readosbf(file,&biases);

    int nb,ii,i;
    gtime_t tmin={0},tmax={0};
    double dt=0.0;

    for(i=0;i<biases.nb;i++){
        if(i==0){
            tmin=biases.data[i].ts;
            tmax=biases.data[i].te;
            dt=timediff(tmax,tmin);
        }
        if (timediff(biases.data[i].ts, tmin) < 0.0) tmin = biases.data[i].ts;
        if (timediff(biases.data[i].te, tmax) > 0.0) tmax = biases.data[i].te;
        if (timediff(biases.data[i].te, biases.data[i].ts) < dt) dt = timediff(biases.data[i].te, biases.data[i].ts);
    }

    nav->osbs=(void *)malloc(sizeof(osbs_t));
    char st1[20], st2[20];
    time2str(tmin, st1, 0);
    time2str(tmax, st2, 0);
    if (!dt) {
        nav->osbs->dt = 0.0;
        nav->osbs->sat_osb = NULL;
        return 0;
    }

    nav->osbs->dt=dt;
    nav->osbs->tmin=tmin;
    nav->osbs->tmax=tmax;
    nb=(int)(timediff(tmax,tmin)/dt);
    nav->osbs->sat_osb=(osb_t *)calloc(nb, sizeof(osb_t));

    if(nav->osbs==NULL){
        return 0;
    }

    int sat=0,code=0;
    for(i=0;i<biases.nb;i++){
        sat=biases.data[i].sat-1;
        code=biases.data[i].code;
        int i1=(int)(timediff(biases.data[i].ts,tmin)/dt);
        int i2=(int)(timediff(biases.data[i].te,tmin)/dt);
        for(ii=i1;ii<i2;ii++){
            if(biases.data[i].type){
                nav->osbs->sat_osb[ii].code[sat][code]=biases.data[i].bia*1E-9*CLIGHT;
            }
            else{
                nav->osbs->sat_osb[ii].phase[sat][code]=biases.data[i].bia*1E-9*CLIGHT;
            }
        }
    }

    free(biases.data);biases.data=NULL;biases.nmax=biases.nb=0;
    return 1;
}

static int readupdf_ewl(const char *file, wl_upds_t *ewls)
{
    FILE *fp;
    char buff[200]={'\0'};
    int sat;

    if(!ewls) return 0;
    if (!(fp=fopen(file,"r"))) {
        trace(2,"ewl upd parameters file open error: %s\n",file);
        return 0;
    }

    while(fgets(buff, sizeof(buff),fp)){
        if(!strncmp(buff, "EOF",3)) break;
        if(!strncmp(buff, "%",1)) continue;
        sat=satid2no(buff+1);
        if(sat<=0) continue;
        ewls->ewl[sat-1]=str2num(buff,14,6);
    }

    fclose(fp);

    return 1;
}

static int readupdf_wl(const char *file, wl_upds_t *wls)
{
    FILE *fp;
    char buff[200]={'\0'};
    int sat;

    if(!wls) return 0;
    if (!(fp=fopen(file,"r"))) {
        trace(2,"wl upd parameters file open error: %s\n",file);
        return 0;
    }

    while(fgets(buff, sizeof(buff),fp)){
        if(!strncmp(buff, "EOF",3)) break;
        if(!strncmp(buff, "%",1)) continue;
        sat=satid2no(buff+1);
        if(sat<=0) continue;
        wls->wl[sat-1]=str2num(buff,14,6);
    }

    fclose(fp);

    return 1;
}

static int readupdf_nl(const char *file, nl_upds_t *nls)
{
    FILE *fp;
    char buff[200]={'\0'};
    int sat,ep=0;
    mjd_t mjd;
    gtime_t t={0};
    nl_upd_t *nl_data_temp;

    if(!nls) return 0;
    if (!(fp=fopen(file,"r"))) {
        trace(2,"nl upd parameters file open error: %s\n",file);
        return 0;
    }

    while(fgets(buff, sizeof(buff),fp)){
        if(!strncmp(buff, "EOF",3)) break;
        if(!strncmp(buff, "%",1)) continue;
        if(!strncmp(buff+1,"EPOCH-TIME",10)){
            mjd.day=(long) str2num(buff,12,8);
            mjd.ds.sn=(long) str2num(buff,20,12);
            mjd.ds.tos=0.0;
            mjd2time(&mjd,&t);
            ep++;
            nls->n=ep;
            continue;
        }
        sat=satid2no(buff+1);
        if(sat<=0) continue;
        if(nls->n>=nls->nmax){
            nls->nmax+=1024;
            if(!(nl_data_temp=(nl_upd_t *)realloc(nls->data, sizeof(nl_upd_t)*nls->nmax))){
                free(nls->data);
                nls->data=NULL;
                nls->n=nls->nmax=0;
                return -1;
            }
            nls->data=nl_data_temp;
        }
        nls->data[ep-1].ts=t;
        nls->data[ep-1].te=timeadd(t,30.0);
        nls->data[ep-1].nl[sat-1]=str2num(buff,15,8);
        nls->data[ep-1].std[sat-1]=str2num(buff,25,8);
    }
    fclose(fp);
    return 0;
}

extern int  readupd(const prcopt_t *opt,char *file_ewl,char *file_wl,char *file_nl, nav_t *nav)
{
    nav->upds=(void *)malloc(sizeof(upds_t));
    readupdf_ewl(file_ewl,&nav->upds->wls);
    readupdf_wl(file_wl,&nav->upds->wls);
    readupdf_nl(file_nl,&nav->upds->nls);

    return 0;
}

/* polynomial interpolation by Neville's algorithm ---------------------------*/
static double interppol(const double *x, double *y, int n)
{
    int i,j;
    
    for (j=1;j<n;j++) {
        for (i=0;i<n-j;i++) {
            y[i]=(x[i+j]*y[i]-x[i]*y[i+1])/(x[i+j]-x[i]);
        }
    }
    return y[0];
}
/* satellite position by precise ephemeris -----------------------------------*/
static int pephpos(gtime_t time, int sat, const nav_t *nav, double *rs,
                   double *dts, double *vare, double *varc)
{
    double t[NMAX+1],p[3][NMAX+1],c[2],*pos,std=0.0,s[3],sinl,cosl;
    int i,j,k,index;
    
    trace(4,"pephpos : time=%s sat=%2d\n",time_str(time,3),sat);
    
    rs[0]=rs[1]=rs[2]=dts[0]=0.0;
    
    if (nav->ne<NMAX+1||
        timediff(time,nav->peph[0].time)<-MAXDTE||
        timediff(time,nav->peph[nav->ne-1].time)>MAXDTE) {
        trace(3,"no prec ephem %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    /* binary search */
    for (i=0,j=nav->ne-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->peph[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;
    
    /* polynomial interpolation for orbit */
    i=index-(NMAX+1)/2;
    if (i<0) i=0; else if (i+NMAX>=nav->ne) i=nav->ne-NMAX-1;
    
    for (j=0;j<=NMAX;j++) {
        t[j]=timediff(nav->peph[i+j].time,time);
        if (norm(nav->peph[i+j].pos[sat-1],3)<=0.0) {
            trace(3,"prec ephem outage %s sat=%2d\n",time_str(time,0),sat);
            return 0;
        }
    }

    for (j=0;j<=NMAX;j++) {
        pos=nav->peph[i+j].pos[sat-1];


#if 0
        p[0][j]=pos[0];
        p[1][j]=pos[1];
#else
        /* correciton for earh rotation ver.2.4.0 */
        sinl=sin(OMGE*t[j]);
        cosl=cos(OMGE*t[j]);
        p[0][j]=cosl*pos[0]-sinl*pos[1];
        p[1][j]=sinl*pos[0]+cosl*pos[1];
#endif
        p[2][j]=pos[2];
    }
    for (i=0;i<3;i++) {
        rs[i]=interppol(t,p[i],NMAX+1);
    }
    if (vare) {
        for (i=0;i<3;i++) s[i]=nav->peph[index].std[sat-1][i];
        std=norm(s,3);
        
        /* extrapolation error for orbit */
        if      (t[0   ]>0.0) std+=EXTERR_EPH*SQR(t[0   ])/2.0;
        else if (t[NMAX]<0.0) std+=EXTERR_EPH*SQR(t[NMAX])/2.0;
        *vare=SQR(std);
    }
    /* linear interpolation for clock */
    t[0]=timediff(time,nav->peph[index  ].time);
    t[1]=timediff(time,nav->peph[index+1].time);
    c[0]=nav->peph[index  ].pos[sat-1][3];
    c[1]=nav->peph[index+1].pos[sat-1][3];
    
    if (t[0]<=0.0) {
        if ((dts[0]=c[0])!=0.0) {
            std=nav->peph[index].std[sat-1][3]*CLIGHT-EXTERR_CLK*t[0];
        }
    }
    else if (t[1]>=0.0) {
        if ((dts[0]=c[1])!=0.0) {
            std=nav->peph[index+1].std[sat-1][3]*CLIGHT+EXTERR_CLK*t[1];
        }
    }
    else if (c[0]!=0.0&&c[1]!=0.0) {
        dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
        i=t[0]<-t[1]?0:1;
        std=nav->peph[index+i].std[sat-1][3]+EXTERR_CLK*fabs(t[i]);
    }
    else {
        dts[0]=0.0;
    }
    if (varc) *varc=SQR(std);
    return 1;
}

static void earth_rotation(const int k,double *pos,double p[3][NMAX+1],double dt)
{
    double sinl,cosl;
    sinl=sin(OMGE*dt);
    cosl=cos(OMGE*dt);
    p[0][k]=cosl*pos[0]-sinl*pos[1];
    p[1][k]=sinl*pos[0]+cosl*pos[1];
    p[2][k]=pos[2];
}

/* satellite clock by precise clock ------------------------------------------*/
static int pephclk(gtime_t time, int sat, const nav_t *nav, double *dts,
                   double *varc)
{
    double t[2],c[2],std;
    int i,j,k,index;
    
    trace(4,"pephclk : time=%s sat=%2d\n",time_str(time,3),sat);
    
    if (nav->nc<2||
        timediff(time,nav->pclk[0].time)<-MAXDTE||
        timediff(time,nav->pclk[nav->nc-1].time)>MAXDTE) {
        trace(4,"no prec clock %s sat=%2d\n",time_str(time,0),sat);
        return 1;
    }
    /* binary search */
    for (i=0,j=nav->nc-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->pclk[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;

    /* linear interpolation for clock */
    t[0]=timediff(time,nav->pclk[index  ].time);
    t[1]=timediff(time,nav->pclk[index+1].time);
    c[0]=nav->pclk[index  ].clk[sat-1][0];
    c[1]=nav->pclk[index+1].clk[sat-1][0];

    for (i=index;i>=0;i--) {
        if (nav->pclk[i].clk[sat-1][0]!=0.0) {
            t[0]=timediff(time,nav->pclk[i].time);
            c[0]=nav->pclk[i].clk[sat-1][0];
            break;
        }
    }

    for (i=index+1;i<nav->nc;i++) {
        if (nav->pclk[i].clk[sat-1][0]!=0.0) {
            t[1]=timediff(time,nav->pclk[i].time);
            c[1]=nav->pclk[i].clk[sat-1][0];
            index=i-1;
            break;
        }
    }
    
    if (t[0]<=0.0) {
        if ((dts[0]=c[0])==0.0) return 0;
        std=nav->pclk[index].std[sat-1][0]*CLIGHT-EXTERR_CLK*t[0];
    }
    else if (t[1]>=0.0) {
        if ((dts[0]=c[1])==0.0) return 0;
        std=nav->pclk[index+1].std[sat-1][0]*CLIGHT+EXTERR_CLK*t[1];
    }
    else if (c[0]!=0.0&&c[1]!=0.0) {
        dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
        i=t[0]<-t[1]?0:1;
        std=nav->pclk[index+i].std[sat-1][0]*CLIGHT+EXTERR_CLK*fabs(t[i]);
    }
    else {
        trace(3,"prec clock outage %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    if (varc) *varc=SQR(std);
    return 1;
}
/* satellite antenna phase center offset ---------------------------------------
* compute satellite antenna phase center offset in ecef
* args   : gtime_t time       I   time (gpst)
*          double *rs         I   satellite position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          double *dant       I   satellite antenna phase center offset (ecef)
*                                 {dx,dy,dz} (m) (iono-free LC value)
* return : none
*-----------------------------------------------------------------------------*/
extern void satantoff(const prcopt_t *popt,gtime_t time, const double *rs, int sat, const nav_t *nav,
                      double *dant)
{
    const pcv_t *pcv=nav->pcvs+sat-1;
    double ex[3],ey[3],ez[3],es[3],r[3],rsun[3],gmst,erpv[5]={0},freq[2];
    double C1,C2,dant1,dant2;
    int i,j=0,k=1,sys,prn;

    trace(4,"satantoff: time=%s sat=%2d\n",time_str(time,3),sat);

    dant[0]=dant[1]=dant[2]=0.0;

    /* sun position in ecef */
    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,&gmst);

    /* unit vectors of satellite fixed coordinates */
    for (i=0;i<3;i++) r[i]=-rs[i];
    if (!normv3(r,ez)) return;
    for (i=0;i<3;i++) r[i]=rsun[i]-rs[i];
    if (!normv3(r,es)) return;
    cross3(ez,es,r);
    if (!normv3(r,ey)) return;
    cross3(ey,ez,ex);

    /* iono-free LC coefficients */
    sys=satsys(sat,&prn);
    if (sys==SYS_GPS||sys==SYS_QZS) { /* L1-L2 */
        freq[0]=FREQ1;
        freq[1]=FREQ2;
        if(sys==SYS_QZS){
            j=0+4*NFREQ;
            k=1+4*NFREQ;
        }
    }
    else if (sys==SYS_GLO) { /* G1-G2 */
        freq[0]=sat2freq(sat,CODE_L1C,nav);
        freq[1]=sat2freq(sat,CODE_L2C,nav);
        j=0+NFREQ;
        k=1+NFREQ;
    }
    else if (sys==SYS_GAL) { /* E1-E5a */
        freq[0]=FREQ1;
        freq[1]=FREQ5;
        j=0+2*NFREQ;
        k=2+2*NFREQ;
    }
    else if (sys==SYS_CMP) {
        if(!strcmp(popt->ac_name,"com")){/* B1I-B2I */
            freq[0]=FREQ1_CMP;
            freq[1]=FREQ2_CMP;
            j=0+3*NFREQ;
            k=1+3*NFREQ;
        }
        else{ /* B1I-B3I */
            freq[0]=FREQ1_CMP;
            freq[1]=FREQ3_CMP;
            j=0+3*NFREQ;
            k=2+3*NFREQ;
        }
    }
    else return;

    C1= SQR(freq[0])/(SQR(freq[0])-SQR(freq[1]));
    C2=-SQR(freq[1])/(SQR(freq[0])-SQR(freq[1]));

    /* iono-free LC */
    for (i=0;i<3;i++) {
        dant1=pcv->off[j][0]*ex[i]+pcv->off[j][1]*ey[i]+pcv->off[j][2]*ez[i];
        dant2=pcv->off[k][0]*ex[i]+pcv->off[k][1]*ey[i]+pcv->off[k][2]*ez[i];
        dant[i]=C1*dant1+C2*dant2;
    }
}
/* satellite position/clock by precise ephemeris/clock -------------------------
* compute satellite position/clock with precise ephemeris/clock
* args   : gtime_t time       I   time (gpst)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          int    opt         I   sat postion option
*                                 (0: center of mass, 1: antenna phase center)
*          double *rs         O   sat position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          double *dts        O   sat clock {bias,drift} (s|s/s)
*          double *var        IO  sat position and clock error variance (m)
*                                 (NULL: no output)
* return : status (1:ok,0:error or data outage)
* notes  : clock includes relativistic correction but does not contain code bias
*          before calling the function, nav->peph, nav->ne, nav->pclk and
*          nav->nc must be set by calling readsp3(), readrnx() or readrnxt()
*          if precise clocks are not set, clocks in sp3 are used instead
*-----------------------------------------------------------------------------*/
extern int peph2pos(const prcopt_t *popt,gtime_t time, int sat, const nav_t *nav, int opt,
                    double *rs, double *dts, double *var)
{
    gtime_t time_tt;
    double rss[3],rst[3],dtss[1],dtst[1],dant[3]={0},vare=0.0,varc=0.0,tt=1E-3;
    int i;
    
    trace(4,"peph2pos: time=%s sat=%2d opt=%d\n",time_str(time,3),sat,opt);
    
    if (sat<=0||MAXSAT<sat) return 0;
    
    /* satellite position and clock bias */
    if (!pephpos(time,sat,nav,rss,dtss,&vare,&varc)||
        !pephclk(time,sat,nav,dtss,&varc)) return 0;
    
    time_tt=timeadd(time,tt);
    if (!pephpos(time_tt,sat,nav,rst,dtst,NULL,NULL)||
        !pephclk(time_tt,sat,nav,dtst,NULL)) return 0;
    
    /* satellite antenna offset correction */
    if (opt) {
        satantoff(popt,time,rss,sat,nav,dant);
    }
    for (i=0;i<3;i++) {
        rs[i  ]=rss[i]+dant[i];
        rs[i+3]=(rst[i]-rss[i])/tt;
    }
    /* relativistic effect correction */
    if (dtss[0]!=0.0) {
        dts[0]=dtss[0]-2.0*dot(rs,rs+3,3)/CLIGHT/CLIGHT;
        dts[1]=(dtst[0]-dtss[0])/tt;
    }
    else { /* no precise clock */
        dts[0]=dts[1]=0.0;
    }
    if (var) *var=vare+varc;
    
    return 1;
}
