#include "rtklib.h"

extern int process(prcopt_t *popt, filopt_t *fopt, solopt_t *sopt) {
    int i;
    DIR *dir;
    struct dirent *file;
    char obs_dir[1024]={'\0'};
    char sep=(char)FILEPATHSEP;
    int yyyy,doy;
    double ep[6];
    time2epoch(popt->ts,ep);
    doy=int(time2doy(popt->ts));
    yyyy=int(ep[0]);
    if(popt->obsdir[0]!='\0')
        sprintf(obs_dir,"%s%c%04d%c%03d%c%s",popt->prcdir,sep,yyyy,sep,doy,sep,popt->obsdir);
    else{
        strcpy(popt->obsdir,"obs");
        sprintf(obs_dir,"%s%c%04d%c%03d%c%s",popt->prcdir,sep,yyyy,sep,doy,sep,"obs");
    }

    if (!(dir = opendir(obs_dir))) {
        return 0;
    }
    int ret = 0, n = 0;
    char *infile[MAXFILE];
    for (i = 0; i < MAXFILE; i++) {
        if (!(infile[i] = (char *) malloc(1024))) {
            for (; i >= 0; i--) free(infile[i]);
            return -1;
        }
        infile[i][0] = '\0';
    }

    char *ext;
    int ppk = (popt->mode >= PMODE_DGPS && popt->mode <= PMODE_STATIC_START) ||
              (popt->mode == PMODE_TC_DGPS || popt->mode == PMODE_TC_PPK||popt->mode==PMODE_STC_PPK) ||
              (popt->insopt.imu_align == INS_ALIGN_GNSS_PPK || popt->insopt.imu_align == INS_ALIGN_GNSS_DGPS);
    int ppp = ((popt->mode >= PMODE_PPP_KINEMA && popt->mode <= PMODE_PPP_FIXED) || (popt->mode == PMODE_TC_PPP||popt->mode==PMODE_LC_PPP||popt->mode==PMODE_STC_PPP));

    while ((file = readdir(dir)) != nullptr) {
        n=0;
        if (strncmp(file->d_name, ".", 1) == 0) continue;
        else if (strstr(file->d_name, "base")) continue;
        else if (strstr(file->d_name, "imu")) continue;
        else if (!(ext = strrchr(file->d_name, '.'))) continue;
        else if (!strstr(ext + 3, "o")) continue;
        char name[5]={'\0'};
        setstr(name,file->d_name,4);
        if(popt->site_list[0]!='\0'){
            if(!strstr(popt->site_list,name)){
                continue;
            }
        }
        setstr(popt->site_name,file->d_name,4);

        fprintf(stderr,"PROCESS %s %s\n",popt->site_name,time_str(popt->ts,0));
        fflush(stderr);
        sprintf(fopt->robsf, "%s%c%s", obs_dir, FILEPATHSEP, file->d_name);
        popt->site_idx++;
        for(i=0;i<4;i++){
            if(file->d_name[i]>='a'&&file->d_name[i]<='z'){
                popt->site_name[i]+='A'-'a';
            }
        }
        /*match output file*/
        matchout(popt, popt->prcdir, fopt,sopt);

        strcpy(infile[n], fopt->robsf);
        n++;

        /* base file*/
        if (ppk) {
            strcpy(infile[n], fopt->bobsf);
            n++;
        }

        /* brdc file*/
        for (i = 0; i < 3; i++) {
            if (strcmp(fopt->navf[i], "")) {
                strcpy(infile[n], fopt->navf[i]);
                n++;
            }
        }

        /* sp3 and clk file*/
        if(ppp){
            for (i=0; i<3; i++){
                if(strcmp(fopt->sp3f[i], "")){
                    strcpy(infile[n],fopt->sp3f[i]);
                    n++;
                }
            }
            for(i=0;i<3;i++){
                if(strcmp(fopt->clkf[i], "")) {
                    strcpy(infile[n], fopt->clkf[i]);
                    n++;
                }
            }
        }
        ret=!postpos(popt->ts,popt->te,0.0,0.0,popt,sopt,fopt,infile,n,fopt->solf,nullptr,nullptr);
        popt->site_name[0]='\0';
    }

    closedir(dir);
    for (i = 0; i < MAXFILE; i++) {
        free(infile[i]);
    }

    return ret;
}

int main(int argc,char **argv)
{
    long t1=clock();
    prcopt_t popt=prcopt_default,popt_;
    solopt_t sopt=solopt_default,sopt_;
    filopt_t fopt={0},fopt_;
    int nsta;
    int port = 0;

    if(!parsecmd(argc,argv,&popt,&sopt,&fopt,&port)) return 0;

    gtime_t ts=popt.ts,te=popt.te;
    int nday;

    if(ts.time!=0.0&&te.time!=0.0&&newround(timediff(popt.te,popt.ts)/86400.0)==0){
        nday=1;
    }
    else{
        nday=newround(timediff(popt.te,popt.ts)/86400.0);
        if(nday>1) popt.prctype=1;
    }

    for(int i=0;i<nday;i++){
        popt_=popt;
        sopt_=sopt;
        fopt_=fopt;
        if(popt.prctype){
            popt_.ts=timeadd(ts,86400.0*i);
            popt_.te=timeadd(ts,86400.0*(i+1));
        }

        if(!loadprcfiles(popt_.prcdir,&popt_,&fopt_, nullptr,&nsta)) return 0;

        if(process(&popt_,&fopt_,&sopt_)){
            long t2=clock();
            double t=(double)(t2-t1)/CLOCKS_PER_SEC;
            fprintf(stderr,"total sec: %5.2f\n",t);
            fflush(stderr);
        }
        else{
            fprintf(stderr,"%s PROCESS ERROR!\n",popt_.site_name);
            fflush(stderr);
        }

        freeprcfiles(&popt_,&fopt_);
    }
}