/*-----------------------------------------------------------------------------
* match-file.cc : match correspond stand files according to date
*
* version : reorganized form chen chao's code
* history : Created by lizhen on 2021/3/29.
*----------------------------------------------------------------------------*/

#include "rtklib.h"
#include <iostream>

using namespace std;

const string kpmodestr[] = {
        "SPP", "TDCP", "DGPS", "PPK-KINE", "PPK-STATIC", "PPK-S-START", "PPK-MOVEB", "PPK-FIXED", "PPP-KINE",
        "PPP-STATIC", "PPP-FIXED",
};

extern int loadprcfiles(const char *dir, const prcopt_t *popt, filopt_t *fopt, sta_t *sta, int *nsta) {
    gtime_t ts=popt->ts, te=popt->te;
    int no_gnss = 0;
    no_gnss=(popt->mode >= PMODE_INS_MECH && popt->mode <= PMODE_LC_POS)&&popt->insopt.imu_align<=INS_ALIGN_GNSS_PV;

    for (int i = 0; i < 3; i++) {
        if (!(fopt->navf[i] = (char *) malloc(MAXSTRPATH))) {
            for (i--; i >= 0; i--) free(fopt->navf[i]);
        } else fopt->navf[i][0] = '\0';
    }
    if (!no_gnss && !match_navfile(popt, ts, te, dir, fopt->navf)) return 0;

    if (!no_gnss && popt->sateph == EPHOPT_PREC) {
        for (int i = 0; i < 3; i++) {
            if (!(fopt->sp3f[i] = (char *) malloc(MAXSTRPATH)))
                for (i--; i >= 0; i--) free(fopt->sp3f[i]);
            else fopt->sp3f[i][0] = '\0';
            if (!(fopt->clkf[i] = (char *) malloc(MAXSTRPATH)))
                for (i--; i >= 0; i--) free(fopt->clkf[i]);
            else fopt->clkf[i][0] = '\0';
        }

        if (!match_sp3file(popt, ts, te, dir, fopt->sp3f)) return 0;
        if (!match_clkfile(popt, ts, te, dir, fopt->clkf)) return 0;
    }
    if (!no_gnss && popt->cbiaopt == CBIAS_OPT_BRD_TGD) {
        if (!match_dcbfile(popt, ts, te, dir, fopt->dcb)) return 0; ///should to match p1c1 dcb file
    } else if (!no_gnss && popt->cbiaopt == CBIAS_OPT_COD_DCB) {
        if (!match_dcbfile(popt, ts, te, dir, fopt->dcb)) return 0;
    } else if (!no_gnss && (popt->cbiaopt == CBIAS_OPT_IGG_DCB||popt->cbiaopt==CBIAS_OPT_GBM_DCB)) {
        if (!match_mgexdcbfile(ts, te, dir, fopt->mgexdcb,popt->cbiaopt)) return 0;
    } else if (!no_gnss && popt->cbiaopt == CBIAS_OPT_MIX_DCB) {
        if (!match_dcbfile(popt, ts, te, dir, fopt->dcb)) return 0;
        if (popt->navsys & SYS_CMP || popt->navsys & SYS_BD3 || popt->navsys & SYS_GAL || popt->navsys & SYS_QZS)
            if (!match_mgexdcbfile(ts, te, dir, fopt->mgexdcb,popt->cbiaopt)) return 0;
    }

    if (!no_gnss && (popt->ionoopt == IONOOPT_TEC || popt->ionoopt == IONOOPT_UC_CONS)) {
        if (!match_ionfile(ts, te, dir, fopt->iono)) return 0;
    }

    int ppp = 0;
    ppp = ((popt->mode >= PMODE_PPP_KINEMA && popt->mode <= PMODE_PPP_FIXED)
           || (popt->mode == PMODE_TC_PPP||popt->mode==PMODE_LC_PPP||popt->mode==PMODE_STC_PPP)) &&
          !no_gnss;
    if (!no_gnss && ppp) {
        if (!match_atxfile(ts, te, dir, fopt->atx, popt)) return 0;
        if ((popt->modear == ARMODE_PPPAR || popt->modear == ARMODE_PPPAR_ILS)) {
            if(popt->arprod==AR_PROD_FCB){
                if(!match_fcbfile(popt,ts,te,dir,fopt->fcb)) return 0;
            }
            else if(popt->arprod>=AR_PROD_OSB_GRM&&popt->arprod<=AR_PROD_OSB_CNT){
                if(!match_biafile(popt,ts,te,dir,fopt->bia)) return 0;
            }
            else if(popt->arprod==AR_PROD_UPD){
                for (int i = 0; i < 3; i++) {
                    if (!(fopt->updf[i] = (char *) malloc(MAXSTRPATH)))
                        for (i--; i >= 0; i--) free(fopt->updf[i]);
                    else
                        fopt->updf[i][0] = '\0';
                }
                if(!match_updfile(popt,ts,te,dir,fopt->updf)) return 0;
            }
        }
    }

    if (popt->tidecorr&2&&!no_gnss) {
        if (!match_blqfile(ts,te,dir,fopt->blq)) return 0;
    }
    if (popt->tidecorr&2&&!no_gnss) {
        if (!match_eopfile(ts, te, dir, fopt->eop, popt)) return 0;
    }

    int ins = popt->mode >= PMODE_INS_MECH;
    int ppk = ((popt->mode >= PMODE_DGPS && popt->mode <= PMODE_STATIC_START) ||
               (popt->mode == PMODE_TC_DGPS || popt->mode == PMODE_TC_PPK) ||
               (popt->mode == PMODE_LC_DGPS || popt->mode == PMODE_LC_PPK) ||popt->mode==PMODE_STC_PPK||
               (popt->insopt.imu_align == INS_ALIGN_GNSS_PPK || popt->insopt.imu_align == INS_ALIGN_GNSS_DGPS)) &&
              !no_gnss;
    if (!no_gnss && ppk) {
        if (!match_baseofile(popt,popt->ts,popt->prcdir, fopt->bobsf)){
            trace(2,"Miss base obs file\n");
            return 0;
        }
    }

    return 1;
}

extern int freeprcfiles(const prcopt_t *popt, filopt_t *fopt) {
    int no_gnss = popt->mode == PMODE_LC_POS || popt->mode==PMODE_INS_MECH;
    for (int i = 0; i < 3; i++) {
        if (fopt->navf[i] != nullptr) free(fopt->navf[i]);
    }

    if (!no_gnss && popt->sateph == EPHOPT_PREC) {
        for (int i = 0; i < 3; i++) {
            if (fopt->sp3f[i] != nullptr) free(fopt->sp3f[i]);
            if (fopt->clkf[i] != nullptr) free(fopt->clkf[i]);
        }
    }

    for(int i=0;i<3;i++){
        if(fopt->updf[i]!= nullptr) free(fopt->updf[i]);
    }

    return 0;
}

extern int loadconf(const char *conf_file, prcopt_t *popt, solopt_t *sopt, filopt_t *fopt) {
    resetsysopts();
    if (!loadopts(conf_file, sysopts) || !loadopts(conf_file, insopts)) {
        return 0;
    }
    getsysopts(popt, sopt, fopt);

    popt->insopt.is_imu_samenoise = 1;

    return 1;
}

extern int parsecmd(int arc, char *arv[], prcopt_t *popt, solopt_t *sopt, filopt_t *fopt, int *port) {
    string mode;
    const char *p;
    int mask = SYS_NONE;
    int level = 128,ar=ARMODE_OFF;
    string conf_file;

    for (int i = 0; i < arc; i++) {
        if (!strcmp(arv[i], "-C") && i + 1 < arc) {
            conf_file = arv[++i];
            if (conf_file.empty()) {
                fprintf(stderr, "OPEN CONFIGURATION FILE ERROR file=%s\n", conf_file.c_str());
                fflush(stderr);
                return 0;
            } else {
                if (!loadconf(conf_file.c_str(), popt, sopt, fopt)) {
                    fprintf(stderr, "LOAD CONFIGURATION INFORMATION ERROR\n");
                    fflush(stderr);
                    return 0;
                }
            }
        } else if (!strcmp(arv[i], "-M") && i + 1 < arc) {
            mode = arv[++i];
        } else if (!strcmp(arv[i], "-P") && i + 1 < arc) {
            *port = atoi(arv[++i]);
        } else if (!strcmp(arv[i], "-S") && i + 1 < arc) {
            p = arv[++i];
            for (; *p && *p != ' '; p++) {
                switch (*p) {
                    case 'G':
                        mask |= SYS_GPS;
                        break;
                    case 'R':
                        mask |= SYS_GLO;
                        break;
                    case 'E':
                        mask |= SYS_GAL;
                        break;
                    case 'C':
                        if (popt->bd3opt == BD3OPT_OFF) {
                            mask |= SYS_CMP;
                            break;
                        } else if (popt->bd3opt == BD3OPT_BD2_3 || popt->bd3opt == BD3OPT_BD23) {
                            mask |= SYS_CMP;
                            mask |= SYS_BD3;
                            break;
                        } else if (popt->bd3opt == BD3OPT_BD3) {
                            mask |= SYS_BD3;
                            break;
                        }
                    case 'J':
                        mask |= SYS_QZS;
                        break;
                }
            }
            if (mask == SYS_NONE) {
                fprintf(stderr, "SATELLITE SYSTEM SET ERROR: %s\n", p);
                fflush(stderr);
                return 0;
            }
        } else if (!strcmp(arv[i], "-L") && i + 1 < arc) {
            level = atoi(arv[++i]);
        }
        else if(!strcmp(arv[i],"-A")&&i+1<arc){
            ar=atoi(arv[++i]);
        }
    }

    /*single-frequency PPP should use GIM product to solve ionospheric delay*/
    if (popt->nf == 1 && popt->ionoopt == IONOOPT_UC) popt->ionoopt = IONOOPT_UC_CONS;

    int i;
    popt->navsys = mask;
    for (i = 0; kpmodestr; i++) {
        if (mode == kpmodestr[i]) {
            popt->mode = i;
            break;
        }
    }
    if (i >= 24) return 0;
    sopt->trace = level;
    popt->modear=ar;
    return 1;
}

extern int match_navfile(const prcopt_t *popt, gtime_t ts, gtime_t te, const char prc_dir[], char *navpaths[]) {
    char dir[1024];
    char sep = (char) FILEPATHSEP;


    /* compute day of year */
    gtime_t t2 = ts;
    gtime_t t1 = timeadd(t2, -86400.0);
    gtime_t t3 = timeadd(t2, 86400.0);
    int yyyy1, doy1, yyyy2, doy2, yyyy3, doy3;
    double ep[6] = {0};
    doy1 = (int) time2doy(t1);
    doy2 = (int) time2doy(t2);
    doy3 = (int) time2doy(t3);
    time2epoch(t1, ep);
    yyyy1 = (int) ep[0];
    time2epoch(t2, ep);
    yyyy2 = (int) ep[0];
    time2epoch(t3, ep);
    yyyy3 = (int) ep[0];
    int yy1 = yyyy1 >= 2000 ? yyyy1 - 2000 : yyyy1 - 1900;
    int yy2 = yyyy2 >= 2000 ? yyyy2 - 2000 : yyyy2 - 1900;
    int yy3 = yyyy3 >= 2000 ? yyyy3 - 2000 : yyyy3 - 1900;

    char tmp[MAXSTRPATH] = {'\0'};
    /* the day before current day */
    sprintf(dir, "%s%c%4d%c%03d%cproducts%cnav", prc_dir,sep,yyyy1,sep,doy1,sep,sep);
    sprintf(tmp, "%s%cbrdm%03d0.%02dp", dir, sep, doy1, yy1);
    if (access(tmp, 0) != -1) setstr(navpaths[0], tmp, MAXSTRPATH);
    else {
        sprintf(tmp, "%s%cbrdc%03d0.%02dn", dir, sep, doy1, yy1);
        if (access(tmp, 0) != -1) setstr(navpaths[0], tmp, MAXSTRPATH);
        else{
            sprintf(tmp, "%s%cBRDM00DLR_S_%04d%03d0000_01D_MN.rnx", dir, sep, yyyy1,doy1);
            if (access(tmp, 0) != -1) setstr(navpaths[0], tmp, MAXSTRPATH);
        }
    }

    sprintf(dir, "%s%c%4d%c%03d%cproducts%cnav", prc_dir,sep,yyyy2,sep,doy2,sep,sep);
    /* current day */
    sprintf(tmp, "%s%cbrdm%03d0.%02dp", dir, sep, doy2, yy2);
    if (access(tmp, 0) != -1) setstr(navpaths[1], tmp, MAXSTRPATH);
    else {
        sprintf(tmp, "%s%cbrdc%03d0.%02dn", dir, sep, doy2, yy2);
        if (access(tmp, 0) != -1) setstr(navpaths[1], tmp, MAXSTRPATH);
        else {
            sprintf(tmp, "%s%cBRDM00DLR_S_%04d%03d0000_01D_MN.rnx", dir, sep, yyyy2,doy2);
            if (access(tmp, 0) != -1) setstr(navpaths[1], tmp, MAXSTRPATH);
            else{
                sprintf(tmp, "%s%cBRDC00IGN_R_%04d%03d0000_01D_MN.rnx", dir, sep, yyyy2,doy2);
                if (access(tmp, 0) != -1) setstr(navpaths[1], tmp, MAXSTRPATH);
                else{
                    trace(2,"Miss navigation file: %s\n",tmp);
                    return 0;
                }
            }
        }
    }

    sprintf(dir, "%s%c%4d%c%03d%cproducts%cnav", prc_dir,sep,yyyy3,sep,doy3,sep,sep);
    /* the day after current day */
    sprintf(tmp, "%s%cbrdm%03d0.%02dp", dir, sep, doy3, yy3);
    if (access(tmp, 0) != -1) setstr(navpaths[2], tmp, MAXSTRPATH);
    else {
        sprintf(tmp, "%s%cbrdc%03d0.%02dn", dir, sep, doy3, yy3);
        if (access(tmp, 0) != -1) setstr(navpaths[2], tmp, MAXSTRPATH);
        else{
            sprintf(tmp, "%s%cBRDM00DLR_S_%04d%03d0000_01D_MN.rnx", dir, sep, yyyy3,doy3);
            if (access(tmp, 0) != -1) setstr(navpaths[2], tmp, MAXSTRPATH);
        }
    }

    return 1;
}

char* getupac(const char ac[]){
    static char s[5];
    if(!(strcmp(ac,"com"))){
        sprintf(s,"%s","COD");
    }
    else if(!strcmp(ac,"wum")){
        sprintf(s,"%s","WUM");
    }
    else if(!strcmp(ac,"grm")){
        sprintf(s,"%s","GRG");
    }
    else if(!strcmp(ac,"gfz")){
        sprintf(s,"%s","GFZ");
    }
    else if(!strcmp(ac,"gbm")){
        sprintf(s,"%s","GBM");
    }

    return s;
}

extern int match_clkfile(const prcopt_t *popt, gtime_t ts, gtime_t te, const char prc_dir[], char *clkpaths[]) {
    char dir[1024];
    char sep = (char) FILEPATHSEP;

    /* compute week, day of week */
    gtime_t t2 = ts;
    gtime_t t1 = timeadd(t2, -86400.0);
    gtime_t t3 = timeadd(t2, 86400.0);
    int wwww1, wwww2, wwww3;
    int dow1 = (int)(time2gpst(t1, &wwww1) / 86400.0);
    int dow2 = (int)(time2gpst(t2, &wwww2) / 86400.0);
    int dow3 = (int)(time2gpst(t3, &wwww3) / 86400.0);
    int yyyy1, doy1, yyyy2, doy2, yyyy3, doy3;
    double ep[6] = {0};
    doy1 = (int) time2doy(t1);
    doy2 = (int) time2doy(t2);
    doy3 = (int) time2doy(t3);
    time2epoch(t1, ep);
    yyyy1 = (int) ep[0];
    time2epoch(t2, ep);
    yyyy2 = (int) ep[0];
    time2epoch(t3, ep);
    yyyy3 = (int) ep[0];


    char tmp[MAXSTRPATH] = {'\0'},file_name[1024]={'\0'};
    /* the day before current day */
    sprintf(dir, "%s%c%4d%c%03d%cproducts%c%s", prc_dir,sep,yyyy1,sep,doy1,sep,sep,popt->ac_name);
    if(popt->aclong){
        if(popt->arprod==AR_PROD_OSB_WHU&&!strcmp(popt->ac_name,"wum")){
            sprintf(file_name,"%3s5MGXFIN_%04d%03d0000_01D_%3s_CLK.CLK","WHU",yyyy1,doy1,popt->clk_int);
        }
        else{
            sprintf(file_name,"%3s0MGX%3s_%04d%03d0000_01D_%3s_CLK.CLK",getupac(popt->ac_name),popt->prdtype,yyyy1,doy1,popt->clk_int);
        }
    }
    else{
        if(popt->arprod==AR_PROD_OSB_WHU&&!strcmp(popt->ac_name,"wum")){
            sprintf(file_name,"%s%04d%d.%s","whu",wwww1,dow1,"clk");
        }
        else{
            sprintf(file_name,"%s%04d%d.%s",popt->ac_name,wwww1,dow1,"clk");
        }
    }
    sprintf(tmp, "%s%c%s", dir,sep, file_name);
    if (access(tmp, 0) != -1) setstr(clkpaths[0], tmp, MAXSTRPATH);

    /* current day */
    sprintf(dir, "%s%c%4d%c%03d%cproducts%c%s", prc_dir,sep,yyyy2,sep,doy2,sep,sep,popt->ac_name);
    if(popt->aclong){
        if(popt->modear==ARMODE_PPPAR_ILS&&popt->arprod==AR_PROD_OSB_WHU&&!strcmp(popt->ac_name,"wum")){
            sprintf(file_name,"%3s5MGXFIN_%04d%03d0000_01D_30S_CLK.CLK","WHU",yyyy2,doy2);
        }
        else{
            sprintf(file_name,"%3s0MGX%3s_%04d%03d0000_01D_%3s_CLK.CLK",getupac(popt->ac_name),popt->prdtype,yyyy2,doy2,popt->clk_int);
        }
    }
    else{
        if(popt->modear==ARMODE_PPPAR_ILS&&popt->arprod==AR_PROD_OSB_WHU&&!strcmp(popt->ac_name,"wum")){
            sprintf(file_name,"%s%04d%d.%s","whu",wwww2,dow2,"clk");
        }
        else{
            sprintf(file_name,"%s%04d%d.%s",popt->ac_name,wwww2,dow2,"clk");
        }
    }
    sprintf(tmp, "%s%c%s", dir,sep, file_name);

    if (access(tmp, 0) != -1) setstr(clkpaths[1], tmp, MAXSTRPATH);
    else {
        trace(2,"Miss precise clock file: %s\n",tmp);
        return 0;
    }

    /* the day after current day */
    sprintf(dir, "%s%c%4d%c%03d%cproducts%c%s", prc_dir,sep,yyyy3,sep,doy3,sep,sep,popt->ac_name);
    if(popt->aclong){
        if(popt->modear==ARMODE_PPPAR_ILS&&popt->arprod==AR_PROD_OSB_WHU&&!strcmp(popt->ac_name,"wum")){
            sprintf(file_name,"%3s5MGXFIN_%04d%03d0000_01D_30S_CLK.CLK","WHU",yyyy3,doy3);
        }
        else{
            sprintf(file_name,"%3s0MGX%3s_%04d%03d0000_01D_%3s_CLK.CLK",getupac(popt->ac_name),popt->prdtype,yyyy3,doy3,popt->clk_int);
        }
    }
    else{
        if(popt->modear==ARMODE_PPPAR_ILS&&popt->arprod==AR_PROD_OSB_WHU&&!strcmp(popt->ac_name,"wum")){
            sprintf(file_name,"%s%04d%d.%s","whu",wwww3,dow3,"clk");
        }
        else{
            sprintf(file_name,"%s%04d%d.%s",popt->ac_name,wwww3,dow3,"clk");
        }
    }
    sprintf(tmp, "%s%c%s", dir,sep, file_name);

    if (access(tmp, 0) != -1) setstr(clkpaths[2], tmp, MAXSTRPATH);

    return 1;
}

extern int match_sp3file(const prcopt_t *popt, gtime_t ts, gtime_t te, const char prc_dir[], char *sp3paths[]) {
    char dir[1024];
    char sep = (char) FILEPATHSEP;

    /* compute week, day of week */
    gtime_t t2 = ts;
    gtime_t t1 = timeadd(t2, -86400.0);
    gtime_t t3 = timeadd(t2, 86400.0);
    int wwww1, wwww2, wwww3;
    int dow1 = (int)(time2gpst(t1, &wwww1) / 86400.0);
    int dow2 = (int)(time2gpst(t2, &wwww2) / 86400.0);
    int dow3 = (int)(time2gpst(t3, &wwww3) / 86400.0);
    int yyyy1, doy1, yyyy2, doy2, yyyy3, doy3;
    double ep[6] = {0};
    doy1 = (int) time2doy(t1);
    doy2 = (int) time2doy(t2);
    doy3 = (int) time2doy(t3);
    time2epoch(t1, ep);
    yyyy1 = (int) ep[0];
    time2epoch(t2, ep);
    yyyy2 = (int) ep[0];
    time2epoch(t3, ep);
    yyyy3 = (int) ep[0];

    char tmp[MAXSTRPATH] = {'\0'},file_name[1024]={'\0'};
    /* the day before current day */
    sprintf(dir, "%s%c%4d%c%03d%cproducts%c%s", prc_dir,sep,yyyy1,sep,doy1,sep,sep,popt->ac_name);
    if(popt->aclong){
        sprintf(file_name,"%3s0MGX%3s_%04d%03d0000_01D_%3s_ORB.SP3",getupac(popt->ac_name),popt->prdtype,yyyy1,doy1,popt->eph_int);
    }
    else{
        sprintf(file_name,"%s%04d%d.%s",popt->ac_name,wwww1,dow1,"sp3");
    }
    sprintf(tmp, "%s%c%s", dir,sep, file_name);
    if (access(tmp, 0) != -1) setstr(sp3paths[0], tmp, MAXSTRPATH);

    /* current day */
    sprintf(dir, "%s%c%4d%c%03d%cproducts%c%s", prc_dir,sep,yyyy2,sep,doy2,sep,sep,popt->ac_name);
    if(popt->aclong){
        sprintf(file_name,"%3s0MGX%3s_%04d%03d0000_01D_%3s_ORB.SP3",getupac(popt->ac_name),popt->prdtype,yyyy2,doy2,popt->eph_int);
    }
    else{
        sprintf(file_name,"%s%04d%d.%s",popt->ac_name,wwww2,dow2,"sp3");
    }
    sprintf(tmp, "%s%c%s", dir, sep, file_name);

    if (access(tmp, 0) != -1) setstr(sp3paths[1], tmp, MAXSTRPATH);
    else {
        trace(2,"Miss precise orbit: %s\n",tmp);
        return 0;
    }

    /* the day after current day */
    sprintf(dir, "%s%c%4d%c%03d%cproducts%c%s", prc_dir,sep,yyyy3,sep,doy3,sep,sep,popt->ac_name);
    if(popt->aclong){
        sprintf(file_name,"%3s0MGX%3s_%04d%03d0000_01D_%3s_ORB.SP3",getupac(popt->ac_name),popt->prdtype,yyyy3,doy3,popt->eph_int);
    }
    else{
        sprintf(file_name,"%s%04d%d.%s",popt->ac_name,wwww3,dow3,"sp3");
    }
    sprintf(tmp, "%s%c%s", dir, sep, file_name);
    if (access(tmp, 0) != -1) setstr(sp3paths[2], tmp, MAXSTRPATH);

    return 1;
}

extern int match_eopfile(gtime_t ts, gtime_t te, const char prc_dir[], char eoppath[], const prcopt_t *popt) {
    int wwww, yy,doy;
    double ep[6] = {0};
    time2gpst(ts, &wwww);
    time2epoch(ts, ep);
    doy=(int)(time2doy(ts));
    yy = (int) ep[0] >= 2000 ? ep[0] - 2000 : ep[0] - 1900;
    char sep = (char) FILEPATHSEP;
    char dir[1024]={'\0'};
    sprintf(dir,"%s%c%04d%c%03d%cproducts%c%s",prc_dir,sep,int(ep[0]),sep,doy,sep,sep,popt->ac_name);

    char tmp[MAXSTRPATH] = {'\0'},file_name[1024]={'\0'};
    if(popt->aclong){
        if(!strcmp(popt->ac_name,"com")){
            sprintf(file_name,"COD0MGXFIN_%4d%03d0000_03D_12H_ERP.ERP",int(ep[0]),doy);
        }
        else if(!strcmp(popt->ac_name,"gbm")){
            sprintf(file_name,"GBM0MGXRAP_%4d%03d0000_01D_01D_ERP.ERP",int(ep[0]),doy);
        }
        else{
            sprintf(file_name,"%3s0MGXFIN_%4d%03d0000_01D_01D_ERP.ERP",getupac(popt->ac_name),int(ep[0]),doy);
        }
    }
    else{
        sprintf(file_name,"%s%04d7.erp",popt->ac_name,wwww);
    }
    sprintf(tmp, "%s%c%s", dir, sep, file_name);
    if (access(tmp, 0) != -1) setstr(eoppath, tmp, MAXSTRPATH);
    else {
        sprintf(tmp, "%s%c%04d%cigs_erp%cigs%4d7.erp", prc_dir, sep,int(ep[0]),sep,sep,wwww);
        if (access(tmp, 0) != -1) setstr(eoppath, tmp, MAXSTRPATH);
        else {
            trace(2,"Miss erp file: %s\n",tmp);
            return 0;
        }
    }

    return 1;
}

extern int match_dcbfile(const prcopt_t *popt, gtime_t ts, gtime_t te, const char prc_dir[], char dcbpath[]) {
    char sep = (char) FILEPATHSEP;
    double ep[6] = {0};
    time2epoch(ts, ep);

    int yyyy=(int) ep[0];
    int mon = (int) ep[1];
    if (mon == 0) return 0;

    char dir[1024];
    sprintf(dir, "%s%c%04d%cdcb", prc_dir,sep,yyyy,sep);
    sprintf(dcbpath, "%s%c*%02d*.DCB", dir, sep, mon);

    return 1;
}

extern int match_ionfile(gtime_t ts, gtime_t te, const char *prc_dir, char *ionpath) {
    /* calculate day of year */
    int doy = 0;
    doy = (int)time2doy(ts);
    double ct[6] = {0};
    time2epoch(ts, ct);
    int yyyy = (int) ct[0];
    int yy = yyyy >= 2000 ? yyyy - 2000 : yyyy - 1900;

    char tmp[MAXSTRPATH] = {'\0'};
    char sep = (char) FILEPATHSEP;
    char dir[1024]={'\0'};
    sprintf(dir,"%s%c%04d%c%03dproducts%cigs_ion",prc_dir,sep,yyyy,sep,doy,sep);

    sprintf(tmp, "%s%cCODG%03d0.%02dI", dir, sep, doy, yy);
    if (access(tmp, 0) != -1) setstr(ionpath, tmp, MAXSTRPATH);
    return 1;
}

extern int match_fcbfile(const prcopt_t *opt, gtime_t ts, gtime_t te, const char prc_dir[], char fcbpath[]) {
    int week, wod;
    wod = (int) (time2gpst(ts, &week) / 86400.0);
    int doy=(int)(time2doy(ts));
    double ep[6];
    time2epoch(ts,ep);
    int yyyy=(int)ep[0];

    char tmp[MAXSTRPATH] = {'\0'};
    char sep = (char) FILEPATHSEP;
    char dir[1024]={'\0'};
    sprintf(dir,"%s%c%04d%c%03d%cproducts%cfcb",prc_dir,sep,yyyy,sep,doy,sep,sep);

    sprintf(tmp, "%s%c%s%4d%d_%3s0MGX%3s.fcb", dir, sep, "sgg", week, wod, getupac(opt->ac_name),opt->prdtype);

    if (access(tmp, 0) != -1) setstr(fcbpath, tmp, MAXSTRPATH);
    else {
        trace(2,"Miss fcb file: %s\n",tmp);
        return 0;
    }

    return 1;
}

extern int match_biafile(const prcopt_t *opt, gtime_t ts, gtime_t te, const char prc_dir[], char biapath[]) {
    int week, wod,doy,yyyy;
    wod = (int)(time2gpst(ts, &week)/86400.0);
    doy = (int)time2doy(ts);
    double ep[6];
    time2epoch(ts,ep);
    yyyy=(int)ep[0];

    char tmp[MAXSTRPATH] = {'\0'},file_name[1024]={'\0'};
    char sep = (char) FILEPATHSEP;
    char dir[1024]={'\0'};
    sprintf(dir,"%s%c%04d%c%03d%cproducts",prc_dir,sep,yyyy,sep,doy,sep);

    if(opt->arprod==AR_PROD_OSB_WHU&&!strcmp(opt->ac_name,"wum")){
        if(opt->aclong){
            sprintf(file_name,"WHU0MGXFIN_%04d%03d0000_01D_01D_ABS.BIA",yyyy,doy);
        }
        else{
            sprintf(file_name,"whu%04d%d.bia",yyyy,wod);
        }
        sprintf(tmp, "%s%c%s%c%s", dir, sep,opt->ac_name,sep, file_name);
    }
    else if(opt->arprod==AR_PROD_OSB_CNT&&!strcmp(opt->ac_name,"cnt")){
        sprintf(tmp, "%s%c%s%c%s%4d%d.bia", dir, sep,opt->ac_name,sep, "cnt", week, wod);
    }
    else if(opt->arprod==AR_PROD_OSB_GRM&&!strcmp(opt->ac_name,"gbm")){
        sprintf(tmp, "%s%c%s%c%s%4d%d.bia", dir, sep,opt->ac_name,sep, "gbm", week, wod);
    }
    else if(opt->arprod==AR_PROD_OSB_COM&&!strcmp(opt->ac_name,"com")){
        if(opt->aclong){
            sprintf(file_name,"COD0MGXFIN_%04d%03d0000_01D_01D_OSB.BIA",yyyy,doy);
        }
        else{
            sprintf(file_name,"com%04d%d.bia",yyyy,wod);
        }
        sprintf(tmp, "%s%c%s%c%s", dir, sep,opt->ac_name,sep, file_name);
    }
    else if(opt->arprod==AR_PROD_OSB_SGG&&!strcmp(opt->ac_name,"com")){
        sprintf(file_name,"SGG%04d%d.BIA",week,wod);
        sprintf(tmp, "%s%c%s%c%s", dir, sep,opt->ac_name,sep, file_name);
    }
    if (access(tmp, 0) != -1) setstr(biapath, tmp, MAXSTRPATH);
    else{
        trace(2,"Miss bias file: %s\n",tmp);
        return 0;
    }

    return 1;
}

extern int match_updfile(const prcopt_t *opt, gtime_t ts, gtime_t te, const char prc_dir[], char *updpath[]) {
    double ep[6]={0};
    int doy=0;
    char tmp[MAXSTRPATH] = {'\0'};
    char sep = (char) FILEPATHSEP;
    char dir[1024]={'\0'};
    sprintf(dir,"%s%cproducts",prc_dir,sep);

    time2epoch(ts,ep);
    doy=(int)time2doy(ts);
    sprintf(tmp, "%s%cupd%cupd_ewl_%4d%03d_GEC",dir,sep,sep,(int)(ep[0]),doy);
    if(access(tmp,0)!=-1) setstr(updpath[0],tmp,MAXSTRPATH);

    sprintf(tmp, "%s%cupd%cupd_wl_%4d%03d_GREC",dir,sep,sep,(int)(ep[0]),doy);
    if(access(tmp,0)!=-1) setstr(updpath[1],tmp,MAXSTRPATH);
    else return 0;

    sprintf(tmp, "%s%cupd%cupd_nl_%4d%03d_GREC",dir,sep,sep,(int)(ep[0]),doy);
    if(access(tmp,0)!=-1) setstr(updpath[2],tmp,MAXSTRPATH);
    else{
        trace(2,"Miss upd file: %s\n",tmp);
        return 0;
    }

    return 1;
}

extern int match_baseofile(const prcopt_t *opt,gtime_t ts,const char *prc_dir, char *basepath) {
    DIR *dir;
    struct dirent *file;
    int ok = 0,doy,yyyy;
    char prcdir[1024]={'\0'};
    double ep[6]={0};
    time2epoch(ts,ep);
    yyyy=int(ep[0]);
    doy=int(time2doy(ts));

    if(opt->obsdir[0]!='\0')
        sprintf(prcdir,"%s%c%04d%c%03d%c%s",opt->prcdir,FILEPATHSEP,yyyy,FILEPATHSEP,doy,FILEPATHSEP,opt->obsdir);
    else{
        sprintf(prcdir,"%s%c%04d%c%03d%c%s",opt->prcdir,FILEPATHSEP,yyyy,FILEPATHSEP,doy,FILEPATHSEP,"obs");
    }

    if (!(dir = opendir(prcdir))) {
        return 0;
    }

    while ((file = readdir(dir)) != nullptr) {
        if (strncmp(file->d_name, ".", 1) == 0) continue;
        else if (strstr(file->d_name, "base")) {
            sprintf(basepath, "%s%c%s", prcdir, FILEPATHSEP, file->d_name);
            ok = 1;
            break;
        } else continue;
    }
    closedir(dir);
    return ok;
}

extern int match_mgexdcbfile(gtime_t ts, gtime_t te, const char prc_dir[], char dcbpath[],int opt) {
    /* calculate day of year */
    int doy = 0;
    doy = (int)(time2doy(ts));
    double ct[6] = {0};
    time2epoch(ts, ct);
    int yyyy = (int) ct[0];

    char tmp[MAXSTRPATH] = {'\0'};
    char sep = (char) FILEPATHSEP;
    char dir[1024]={'\0'};

    if(opt==CBIAS_OPT_IGG_DCB||opt==CBIAS_OPT_MIX_DCB){
        sprintf(dir,"%s%c%04d%c%03d%cproducts%ccas",prc_dir,sep,yyyy,sep,doy,sep,sep);
        sprintf(tmp, "%s%cCAS0MGXRAP_%04d%03d0000_01D_01D_DCB.BSX", dir, sep, yyyy, doy);
    }
    else if(opt==CBIAS_OPT_GBM_DCB){
        sprintf(dir,"%s%c%04d%c%03d%cproducts%cgbm",prc_dir,sep,yyyy,sep,doy,sep,sep);
        sprintf(tmp, "%s%cGBM0MGXRAP_%04d%03d0000_01D_01D_REL.BIA", dir, sep, yyyy, doy);
    }

    if (access(tmp, 0) != -1) setstr(dcbpath, tmp, MAXSTRPATH);
    else {
        return 0;
    }

    return 1;
}

extern int match_reffile(const prcopt_t *opt,const char *obs_file,gtime_t ts, gtime_t te, const char prc_dir[], char refpath[],int dynamic) {
    int wwww,yyyy,doy,yy;
    double ep[6]={0};
    time2epoch(ts,ep);
    doy=int(time2doy(ts));
    time2gpst(ts, &wwww);
    yyyy=int(ep[0]);
    yy = yyyy >= 2000 ? yyyy - 2000 : yyyy - 1900;
    int gnss=0;

    gnss=opt->mode<PMODE_INS_MECH;
    char tmp[MAXSTRPATH] = {'\0'},dir[MAXSTRPATH];
    char sep = (char) FILEPATHSEP;

    if(dynamic){
        if(gnss){
            sprintf(tmp,"%s.gnss",obs_file);
        }
        else if(opt->mode==PMODE_INS_MECH){
            sprintf(tmp,"%s.mech",obs_file);
        }
        else{
            sprintf(tmp,"%s.tc",obs_file);
        }
    }
    else{
        sprintf(dir,"%s%c%04d%cigs_snx",prc_dir,sep,yyyy,sep);
        sprintf(tmp, "%s%cigs%2dP%4d.snx", dir, sep,yy, wwww);
    }

    if (access(tmp, 0) != -1) setstr(refpath, tmp, MAXSTRPATH);
    else return 0;

    return 1;
}

extern int match_blqfile(gtime_t ts, gtime_t te, const char prc_dir[], char blqpath[]) {
    char tmp[MAXSTRPATH] = {'\0'};
    char sep = (char) FILEPATHSEP;
    char dir[1024]={'\0'};
    sprintf(dir,"%s%cblq",prc_dir,sep);

    sprintf(tmp, "%s%c%s", dir, sep, "ocnload.blq\0");
    if (access(tmp, 0) != -1) setstr(blqpath, tmp, MAXSTRPATH);

    return 1;
}

extern int match_atxfile(gtime_t ts, gtime_t te, const char prc_dir[], char atxPath[], const prcopt_t *popt) {
    double ct_atx_1[6] = {2006, 11, 5, 0, 0, 0.00000000};
    double ct_atx_2[6] = {2011, 4, 17, 0, 0, 0.00000000};
    double ct_atx_3[6] = {2017, 1, 29, 0, 0, 0.00000000};
    gtime_t gt_atx_1 = epoch2time(ct_atx_1);
    gtime_t gt_atx_2 = epoch2time(ct_atx_2);
    gtime_t gt_atx_3 = epoch2time(ct_atx_3);
    double dt1 = timediff(ts, gt_atx_1);
    double dt2 = timediff(ts, gt_atx_2);
    double dt3 = timediff(ts, gt_atx_3);

    char tmp0[500] = {'\0'};
    if (strcmp(popt->ac_name, "gbm") == 0) {
        if (dt3 < 0) sprintf(tmp0, "igs08_gbm.atx");
        else sprintf(tmp0, "igs14_%4d.atx",popt->atx_week);
    } else if (strcmp(popt->ac_name, "wum") == 0) {
        if (dt3 < 0) sprintf(tmp0, "igs08_wum.atx");
        else sprintf(tmp0, "igs14_%4d.atx",popt->atx_week);
    } else {
        if (dt1 < 0) sprintf(tmp0, "igs01_igs.atx%c", '\0');
        else if (dt2 < 0) sprintf(tmp0, "igs05_igs.atx%c", '\0');
        else if (dt3 < 0) sprintf(tmp0, "igs08_igs.atx%c", '\0');
        else sprintf(tmp0, "igs14_%4d.atx",popt->atx_week);
    }

    char tmp[MAXSTRPATH] = {'\0'};
    char sep = (char) FILEPATHSEP;
    char dir[1024]={'\0'};
    sprintf(dir,"%s%catx",prc_dir,sep);
    sprintf(tmp, "%s%c%s", dir, sep, tmp0);
    if (access(tmp, 0) != -1) setstr(atxPath, tmp, MAXSTRPATH);
    else {
        return 0;
    }

    return 1;
}


extern void matchout(const prcopt_t *popt, const char *prcdir, filopt_t *fopt,const solopt_t *sopt) {
    char sep = (char) FILEPATHSEP, sys[8] = {'\0'};
    int no_gnss = 0;
    char outdir[1024] = {'\0'}, outname[50] = {'\0'};
    double ep[6];
    int yyyy,doy;
    time2epoch(popt->ts,ep);
    yyyy=(int)ep[0];
    doy=time2doy(popt->ts);
    no_gnss=(popt->mode >= PMODE_INS_MECH && popt->mode <= PMODE_LC_POS);

    sprintf(outdir, "%s%c%04d%c%03d%c%s%s%c", prcdir, sep,yyyy,sep,doy,sep, "result_",kpmodestr[popt->mode].c_str(),sep);

    if (!no_gnss) {
        if (popt->navsys & SYS_GPS) {
            sprintf(sys, "%s%s", sys, "G");
        }
        if (popt->navsys & SYS_GLO) {
            sprintf(sys, "%s%s", sys, "R");
        }
        if (popt->navsys & SYS_GAL) {
            sprintf(sys, "%s%s", sys, "E");
        }
        if ((popt->navsys & SYS_CMP) || (popt->navsys & SYS_BD3)) {
            if (popt->bd3opt == BD3OPT_OFF) sprintf(sys, "%s%s", sys, "B2");
            else if (popt->bd3opt == BD3OPT_BD23) sprintf(sys, "%s%s", sys, "C");
            else if (popt->bd3opt == BD3OPT_BD2_3) sprintf(sys, "%s%s", sys, "B2B3");
            else if (popt->bd3opt == BD3OPT_BD3) sprintf(sys, "%s%s", sys, "B3");
        }
        if (popt->navsys & SYS_QZS) {
            sprintf(sys, "%s%s", sys, "J");
        }
        sprintf(outname, "%s_%s", popt->site_name, sys);

        if (popt->nf == 1) {
            sprintf(outname, "%s_%s", outname, "SF");
        } else if (popt->nf == 2) {
            sprintf(outname, "%s_%s", outname, "DF");
        } else if (popt->nf == 3) {
            sprintf(outname, "%s_%s", outname, "TF");
        } else if (popt->nf == 4) {
            sprintf(outname, "%s_%s", outname, "QF");
        }

        if (popt->ionoopt == IONOOPT_IFLC) {
            sprintf(outname, "%s_%s", outname, "IF");
        } else if (popt->ionoopt == IONOOPT_IF2) {
            sprintf(outname, "%s_%s", outname, "IF2");

        } else if (popt->ionoopt == IONOOPT_UC) {
            sprintf(outname, "%s_%s", outname, "UC");
        } else if (popt->ionoopt == IONOOPT_UC_CONS) {
            sprintf(outname, "%s_%s", outname, "UC_CON");
        }

        int ar = popt->modear>ARMODE_OFF;
        if (ar) { /*PPP*/
            if(popt->modear==ARMODE_PPPAR_ILS){
                if(popt->arprod==AR_PROD_IRC){
                    sprintf(outname, "%s%s", outname, "_FIX_IRC");
                }
                else if(popt->arprod==AR_PROD_FCB){
                    sprintf(outname, "%s%s", outname, "_FIX_FCB");
                }
                else if(popt->arprod==AR_PROD_OSB_WHU||popt->arprod==AR_PROD_OSB_GRM||popt->arprod==AR_PROD_OSB_COM
                        ||popt->arprod==AR_PROD_OSB_COM||popt->arprod==AR_PROD_OSB_CNT||popt->arprod==AR_PROD_OSB_SGG){
                    sprintf(outname, "%s%s", outname, "_FIX_OSB");
                }
                else if(popt->arprod==AR_PROD_UPD){
                    sprintf(outname, "%s%s", outname, "_FIX_UPD");
                }
            }
            else{ /*PPK*/
                if(popt->modear==ARMODE_CONT){
                    sprintf(outname, "%s%s", outname, "_FIX");
                }
                else if(popt->modear==ARMODE_INST){
                    sprintf(outname, "%s%s", outname, "_INST");
                }
                else if(popt->modear==ARMODE_FIXHOLD){
                    sprintf(outname, "%s%s", outname, "_HOLD");
                }
            }
        } else {
            sprintf(outname, "%s%s", outname, "_FLOAT");
        }
        if(access(outdir,0)!=0){
            createdir(outdir);
        }
        sprintf(outdir,"%s%s%c",outdir,popt->obsdir,sep);
        if(access(outdir,0)!=0){
            createdir(outdir);
        }
        sprintf(fopt->solf, "%s%s%c%s.pos", outdir,popt->ac_name, sep, outname);

        if(popt->insopt.rts){
            sprintf(fopt->rts_ins_fw, "%s", "rts.bin");
            sprintf(fopt->rtsfile, "%s%s%c%s.rts", outdir,popt->ac_name, sep, outname);
        }

    } else {
        if(access(outdir,0)!=0) {
            createdir(outdir);
        }
        sprintf(outdir,"%s%s%c",outdir,popt->obsdir,sep);
        if(access(outdir,0)!=0){
            createdir(outdir);
        }
        if(popt->mode==PMODE_INS_MECH){
            sprintf(fopt->solf, "%s%s.pos", outdir, "INS_MECH");
        }
        else if (popt->mode == PMODE_LC_POS) {
            sprintf(fopt->solf, "%s%s.pos", outdir, "LC_POS");
        }
        if(popt->insopt.rts){
            sprintf(fopt->rts_ins_fw, "%s", "rts.bin");
            sprintf(fopt->rtsfile, "%s%c%s.rts", outdir, sep, "LC_SOL");
        }
    }

    if(sopt->ambres&&(popt->mode>=PMODE_PPP_KINEMA||popt->mode<=PMODE_PPP_FIXED)&&(popt->modear==ARMODE_PPPAR||popt->modear==ARMODE_PPPAR_ILS)){
        sprintf(fopt->wl_amb,"%s%c%s%c%s.wlamb",outdir,sep,popt->ac_name,sep,outname);
        sprintf(fopt->nl_amb,"%s%c%s%c%s.nlamb",outdir,sep,popt->ac_name,sep,outname);
        sprintf(fopt->lc_amb,"%s%c%s%c%s.lcamb",outdir,sep,popt->ac_name,sep,outname);
    }
}
