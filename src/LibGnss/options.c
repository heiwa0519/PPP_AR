/*------------------------------------------------------------------------------
* options.c : options functions
*
*          Copyright (C) 2010-2015 by T.TAKASU, All rights reserved.
*
* version : $Revision:$ $Date:$
* history : 2010/07/20  1.1  moved from postpos.c
*                            added api:
*                                searchopt(),str2opt(),opt2str(),opt2buf(),
*                                loadopts(),saveopts(),resetsysopts(),
*                                getsysopts(),setsysopts()
*           2010/09/11  1.2  add options
*                                pos2-elmaskhold,pos1->snrmaskena
*                                pos1-snrmask1,2,3
*           2013/03/11  1.3  add pos1-posopt1,2,3,4,5,pos2-syncsol
*                                misc-rnxopt1,2,pos1-snrmask_r,_b,_L1,_L2,_L5
*           2014/10/21  1.4  add pos2-bdsarmode
*           2015/02/20  1.4  add ppp-fixed as pos1-posmode option
*           2015/05/10  1.5  add pos2-arthres1,2,3,4
*           2015/05/31  1.6  add pos2-armaxiter, pos1-posopt6
*                            add selection precise for pos1-pospot3
*           2015/11/26  1.7  modify pos1-frequency 4:l1+l2+l5+l6 -> l1+l5
*           2015/12/05  1.8  add misc-pppopt
*           2016/06/10  1.9  add ant2-maxaveep,ant2-initrst
*           2016/07/31  1.10 add out-outsingle,out-maxsolstd
*           2017/06/14  1.11 add out-outvel
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* system options buffer -----------------------------------------------------*/
static prcopt_t prcopt_;
static solopt_t solopt_;
static filopt_t filopt_;
static int antpostype_[2];
static double elmask_,elmaskar_,elmaskhold_;
static double antpos_[2][3];
static char exsats_[1024];
static char snrmask_[NFREQ][1024];
static char prc_ts_[100];
static char prc_te_[100];
static char GPS_frq_[20];
static char GLO_frq_[20];
static char GAL_frq_[20];
static char BD2_frq_[20];
static char BD3_frq_[20];
static char QZS_frq_[20];
static char imucoord[10];

/* system options table ------------------------------------------------------*/
#define SWTOPT  "0:off,1:on"
#define MODOPT  "0:single,1:dgps,2:kinematic,3:static,4:static-start,5:movingbase,6:fixed,7:ppp-kine,8:ppp-static,9:ppp-fixed,10:ins_mech,11:lc,12:tc"
#define FRQOPT  "1:SF,2:DF,3:TF,4:QF"
#define TYPOPT  "0:forward,1:backward,2:combined"
#define IONOPT  "0:off,1:brdc,2:sbas,3:tec,4:IF,5:IF2,6:UC,7:UC-tec,8:qzs-brdc,9:qzs-lex,10:stec"
#define TRPOPT  "0:off,1:saas,2:sbas,3:est-ztd,4:est-ztdgrad,5:ztd"
#define EPHOPT  "0:brdc,1:precise,2:brdc+sbas,3:brdc+ssrapc,4:brdc+ssrcom"
#define NAVOPT  "1:gps+2:sbas+4:glo+8:gal+16:qzs+32:comp"
#define GAROPT  "0:off,1:on,2:autocal,3:fix-and-hold"
#define WEIGHTOPT "0:elevation,1:snr"
#define SOLOPT  "0:llh,1:xyz,2:enu,3:nmea,4:stat,5:gsif,6:ins_llh,7:ins_xyz,8:ins_ygm"
#define TSYOPT  "0:gpst,1:utc,2:jst"
#define TFTOPT  "0:tow,1:hms"
#define DFTOPT  "0:deg,1:dms"
#define HGTOPT  "0:ellipsoidal,1:geodetic"
#define GEOOPT  "0:internal,1:egm96,2:egm08_2.5,3:egm08_1,4:gsi2000"
#define STAOPT  "0:all,1:single"
#define STSOPT  "0:off,1:state,2:residual"
#define ARMOPT  "0:off,1:continuous,2:instantaneous,3:fix-and-hold,4:wl-nl,5:tcar,6:ppp-ar,7:ppp-ar-ils"
#define POSOPT  "0:llh,1:xyz,2:single,3:posfile,4:rinexhead,5:rtcm,6:raw"
#define TIDEOPT "0:off,1:on,2:otl,7:all"
#define PHWOPT  "0:off,1:on,2:precise"
/* receiver options table ----------------------------------------------------*/
#define TIMOPT  "0:gpst,1:utc,2:jst,3:tow"
#define CONOPT  "0:dms,1:deg,2:xyz,3:enu,4:pyl"
#define FLGOPT  "0:off,1:std+2:age/ratio/ns"
#define ISTOPT  "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,7:ntripcli,8:ftp,9:http"
#define OSTOPT  "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,6:ntripsvr,11:ntripc_c"
#define FMTOPT  "0:rtcm2,1:rtcm3,2:oem4,3:oem3,4:ubx,5:ss2,6:hemis,7:skytraq,8:gw10,9:javad,10:nvs,11:binex,12:rt17,13:sbf,14:cmr,15:tersus,17:rinex,18:sp3,19:rnxclk,20:sbas,21:nmea,22:imu_ygm_sim,24:imu_nvt_cpt,25:avp_ygm_sim,26:avp_pos,27:od_ygm_sim,28:imu_m39"
#define NMEOPT  "0:off,1:latlon,2:single"
#define MSGOPT  "0:all,1:rover,2:base,3:corr"
#define ARPROD  "0:off,1:irc,2:fcb,3:upd,4:osb_grm,5:osb_whu,6:osb_com,7:osb_sgg,8:osb_cnt"
#define CBIAOPT "0:off,1:brd_tgd,2:code_dcb,3:igg_dcb,4:gbm_dcb,5:mix_dcb,6:igg_bia,7:code_bia,8:grm_bia"
#define BD3OPT  "0:off,1:bd2-3,2:bd2/3,3:bd3_only"
#define VELOPT   "0:off,1:doppler,2:tdcp"
#define TRPMAP  "0:el,1:nmf,2:gmf,3:vmf"
#define MECHCOORD  "1:llh,2:ecef"
#define LOCALCOORD "1:enu-rfu,2:ned-frd"
#define IMUALIGN   "0:manual,1:corse,2:sol,3:spp,4:dgps,5:ppk,6:ppp"
#define IMUSTDET   "0:off,1:glrt,2:mv,3:mag,4:are,5:all"
#define RESQCOPT   "0:off,1:pr_igg,2:cp_igg,3:igg"
#define INSCYCLE   "0:off,1:li,2:han"
#define KFOPT      "0:kf,1:ada_inno,2:vbkf,3:sage"
#define REFOPT     "0:ie,1:a15,2:pos,3:tex,4:ygm_avp,5:ygm_avpde,6:sinex"
#define LCPOSFMT   "0:none,1:pos,2:ygm_pv"

const char kpmodestr[22][12]={
        "SPP","DGPS","PPK-KINE","PPK-STATIC","PPK-S-START","PPK-MOVEB","PPK-FIXED","PPP-KINE","PPP-STATIC","PPP-FIXED",
        "INS-MECH","IGLC-POS",
        "IGLC-SPP","IGLC-DGPS","IGLC-PPK","IGLC-PPP",
        "IGTC-SPP","IGTC-DGPS","IGTC-PPK","IGTC-PPP"
};

EXPORT opt_t sysopts[]={
    {"pos1-prcdir",    2,  (void *)&prcopt_.prcdir,       "" },
    {"pos1-obsdir",    2,  (void *)&prcopt_.obsdir,       "" },
    {"pos1-prcts",     2,  (void *)&prc_ts_,       "" },
    {"pos1-prcte",     2,  (void *)&prc_te_,       "" },
    {"pos1-site_list", 2,  (void *)&prcopt_.site_list,""},
    {"pos1-sample",    1,  (void *)&prcopt_.sample,       "" },
    {"pos1-acname",    2,  (void *)&prcopt_.ac_name,       "" },
    {"pos1-aclong",    0,  (void *)&prcopt_.aclong,       "" },
    {"pos1-prdtype",   2,  (void *)&prcopt_.prdtype,     "" },
    {"pos1-eph_int",   2,  (void *)&prcopt_.eph_int,     "" },
    {"pos1-clk_int",   2,  (void *)&prcopt_.clk_int,     "" },
    {"pos1-atx_week",  0,  (void *)&prcopt_.atx_week,     "" },
    {"pos1-nf",        3,  (void *)&prcopt_.nf,         FRQOPT },
    {"pos1-GPS_freq",  2,  (void *)&GPS_frq_,         ""},
    {"pos1-GLO_freq",  2,  (void *)&GLO_frq_,         "" },
    {"pos1-GAL_freq",  2,  (void *)&GAL_frq_,         "" },
    {"pos1-BD2_freq",  2,  (void *)&BD2_frq_,         "" },
    {"pos1-BD3_freq",  2,  (void *)&BD3_frq_,         "" },
    {"pos1-QZS_freq",  2,  (void *)&QZS_frq_,         "" },
    {"pos1-soltype",    3,  (void *)&prcopt_.soltype,    TYPOPT },
    {"pos1-elmask",     1,  (void *)&elmask_,            "deg"  },
    {"pos1-snrmask_r",  3,  (void *)&prcopt_.snrmask.ena[0],SWTOPT},
    {"pos1-snrmask_b",  3,  (void *)&prcopt_.snrmask.ena[1],SWTOPT},
    {"pos1-snrmask_L1", 2,  (void *)snrmask_[0],         ""     },
    {"pos1-snrmask_L2", 2,  (void *)snrmask_[1],         ""     },
    {"pos1-snrmask_L5", 2,  (void *)snrmask_[2],         ""     },
    {"pos1-dynamics",   3,  (void *)&prcopt_.dynamics,   SWTOPT },
    {"pos1-tidecorr",   3,  (void *)&prcopt_.tidecorr,   TIDEOPT},
    {"pos1-ionoopt",    3,  (void *)&prcopt_.ionoopt,    IONOPT },
    {"pos1-tropopt",    3,  (void *)&prcopt_.tropopt,    TRPOPT },
    {"pos1-sateph",     3,  (void *)&prcopt_.sateph,     EPHOPT },
    {"pos1-velopt",     3,  (void *)&prcopt_.velopt,     VELOPT },
    {"pos1-posopt1",    3,  (void *)&prcopt_.posopt[0],  SWTOPT },
    {"pos1-posopt2",    3,  (void *)&prcopt_.posopt[1],  SWTOPT },
    {"pos1-posopt3",    3,  (void *)&prcopt_.posopt[2],  PHWOPT },
    {"pos1-posopt4",    3,  (void *)&prcopt_.posopt[3],  SWTOPT },
    {"pos1-posopt5",    3,  (void *)&prcopt_.posopt[4],  SWTOPT },
    {"pos1-posopt6",    3,  (void *)&prcopt_.posopt[5],  SWTOPT },
    {"pos1-sdopt",      3,  (void *)&prcopt_.sdopt,  SWTOPT },
    {"pos1-exclsats",   2,  (void *)exsats_,             "prn ..."},
    {"pos1-navsys",     0,  (void *)&prcopt_.navsys,     NAVOPT },
    {"pos1-arprod",     3,  (void *)&prcopt_.arprod,     ARPROD },
    {"pos1-geo_opt",    3,  (void *)&prcopt_.geo_opt,    SWTOPT },
    {"pos1-bd3opt",     3,  (void *)&prcopt_.bd3opt,     BD3OPT},
    {"pos1-adjobs",     3,  (void *)&prcopt_.adjobs,     SWTOPT },
    {"pos1-cbiaopt",    3,  (void *)&prcopt_.cbiaopt,   CBIAOPT},
    {"pos1-trpmap",     3,  (void *)&prcopt_.tropmap,   TRPMAP},
    {"pos1-restart",     1,  (void *)&prcopt_.restart,   ""},

//    {"pos2-armode",     3,  (void *)&prcopt_.modear,     ARMOPT },
    {"pos2-paramb",     3,  (void *)&prcopt_.par,     SWTOPT },
    {"pos2-gloarmode",  3,  (void *)&prcopt_.glomodear,  GAROPT },
    {"pos2-bdsarmode",  3,  (void *)&prcopt_.bdsmodear,  SWTOPT },
    {"pos2-galarmode",  3,  (void *)&prcopt_.galmodear,  SWTOPT },
    {"pos2-arfilter",   3,  (void *)&prcopt_.arfilter,   SWTOPT },
    {"pos2-arthres",    1,  (void *)&prcopt_.thresar[0], ""     },
    {"pos2-arthres1",   1,  (void *)&prcopt_.thresar[1], ""     },
    {"pos2-arthres2",   1,  (void *)&prcopt_.thresar[2], ""     },
    {"pos2-arthres3",   1,  (void *)&prcopt_.thresar[3], ""     },
    {"pos2-arthres4",   1,  (void *)&prcopt_.thresar[4], ""     },
    {"pos2-varholdamb", 1,  (void *)&prcopt_.varholdamb, "cyc^2"},
    {"pos2-gainholdamb",1,  (void *)&prcopt_.gainholdamb,""     },
    {"pos2-arlockcnt",  0,  (void *)&prcopt_.minlock,    ""     },
    {"pos2-minfixsats", 0,  (void *)&prcopt_.minfixsats, ""     },
    {"pos2-minholdsats",0,  (void *)&prcopt_.minholdsats,""     },
    {"pos2-mindropsats",0,  (void *)&prcopt_.mindropsats,""     },
    {"pos2-rcvstds",    3,  (void *)&prcopt_.rcvstds,    SWTOPT },
    {"pos2-arelmask",   1,  (void *)&elmaskar_,          "deg"  },
    {"pos2-arminfix",   0,  (void *)&prcopt_.minfix,     ""     },
    {"pos2-armaxiter",  0,  (void *)&prcopt_.armaxiter,  ""     },
    {"pos2-elmaskhold", 1,  (void *)&elmaskhold_,        "deg"  },
    {"pos2-aroutcnt",   0,  (void *)&prcopt_.maxout,     ""     },
    {"pos2-maxage",     1,  (void *)&prcopt_.maxtdiff,   "s"    },
    {"pos2-syncsol",    3,  (void *)&prcopt_.syncsol,    SWTOPT },
    {"pos2-slipthres",  1,  (void *)&prcopt_.thresslip,  "m"    },
    {"pos2-robust",     3,  (void *)&prcopt_.robust,  RESQCOPT },
    {"pos2-igg_k0",      1,  (void *)&prcopt_.igg_k0,  "" },
    {"pos2-igg_k1",      1,  (void *)&prcopt_.igg_k1,  "" },
    {"pos2-kfopt",       3,  (void *)&prcopt_.kfopt,      KFOPT },
    {"pos2-rejionno",   1,  (void *)&prcopt_.maxinno,    "m"    },
    {"pos2-rejgdop",    1,  (void *)&prcopt_.maxgdop,    ""     },
    {"pos2-niter",      0,  (void *)&prcopt_.niter,      ""     },
    {"pos2-baselen",    1,  (void *)&prcopt_.baseline[0],"m"    },
    {"pos2-basesig",    1,  (void *)&prcopt_.baseline[1],"m"    },
    
    {"out-solformat",   3,  (void *)&solopt_.posf,       SOLOPT },
    {"out-outhead",     3,  (void *)&solopt_.outhead,    SWTOPT },
    {"out-outopt",      3,  (void *)&solopt_.outopt,     SWTOPT },
    {"out-outvel",      3,  (void *)&solopt_.outvel,     SWTOPT },
    {"out-outclk",      3,  (void *)&solopt_.outclk,     SWTOPT },
    {"out-outtrp",      3,  (void *)&solopt_.outtrp,     SWTOPT },
    {"out-outba",       3,  (void *)&solopt_.outba,     SWTOPT },
    {"out-outbg",       3,  (void *)&solopt_.outbg,     SWTOPT },
    {"out-outsa",       3,  (void *)&solopt_.outsa,     SWTOPT },
    {"out-outsg",       3,  (void *)&solopt_.outsg,     SWTOPT },
    {"out-outarm",      3,  (void *)&solopt_.outarm,     SWTOPT },
    {"out-timesys",     3,  (void *)&solopt_.times,      TSYOPT },
    {"out-timeform",    3,  (void *)&solopt_.timef,      TFTOPT },
    {"out-timendec",    0,  (void *)&solopt_.timeu,      ""     },
    {"out-degform",     3,  (void *)&solopt_.degf,       DFTOPT },
    {"out-fieldsep",    2,  (void *) solopt_.sep,        ""     },
    {"out-outsingle",   3,  (void *)&prcopt_.outsingle,  SWTOPT },
    {"out-maxsolstd",   1,  (void *)&solopt_.maxsolstd,  "m"    },
    {"out-height",      3,  (void *)&solopt_.height,     HGTOPT },
    {"out-geoid",       3,  (void *)&solopt_.geoid,      GEOOPT },
    {"out-solstatic",   3,  (void *)&solopt_.solstatic,  STAOPT },
    {"out-nmeaintv1",   1,  (void *)&solopt_.nmeaintv[0],"s"    },
    {"out-nmeaintv2",   1,  (void *)&solopt_.nmeaintv[1],"s"    },
    {"out-outstat",     3,  (void *)&solopt_.sstat,      STSOPT },
    {"out-ambres",      3,  (void *)&solopt_.ambres,      SWTOPT },
    {"out-solfrq",      0,  (void *)&solopt_.outfrq,      "Hz"},
    {"out-reffmt",      3,  (void *)&solopt_.ref_type,      REFOPT},

    {"stats-weightmode",3,  (void *)&prcopt_.weightmode, WEIGHTOPT},
    {"stats-eratio1",   1,  (void *)&prcopt_.eratio[0],  ""     },
    {"stats-eratio2",   1,  (void *)&prcopt_.eratio[1],  ""     },
    {"stats-eratio5",   1,  (void *)&prcopt_.eratio[2],  ""     },
    {"stats-errphase",  1,  (void *)&prcopt_.err[1],     "m"    },
    {"stats-errphaseel",1,  (void *)&prcopt_.err[2],     "m"    },
    {"stats-errphasebl",1,  (void *)&prcopt_.err[3],     "m/10km"},
    {"stats-errdoppler",1,  (void *)&prcopt_.err[4],     "Hz"   },
    {"stats-snrmax",    1,  (void *)&prcopt_.err[5],     "dB.Hz"},
    {"stats-stdbias",   1,  (void *)&prcopt_.std[0],     "m"    },
    {"stats-stdiono",   1,  (void *)&prcopt_.std[1],     "m"    },
    {"stats-stdtrop",   1,  (void *)&prcopt_.std[2],     "m"    },
    {"stats-prnaccelh", 1,  (void *)&prcopt_.prn[3],     "m/s^2"},
    {"stats-prnaccelv", 1,  (void *)&prcopt_.prn[4],     "m/s^2"},
    {"stats-prnbias",   1,  (void *)&prcopt_.prn[0],     "m"    },
    {"stats-prniono",   1,  (void *)&prcopt_.prn[1],     "m"    },
    {"stats-prntrop",   1,  (void *)&prcopt_.prn[2],     "m"    },
    {"stats-prnpos",    1,  (void *)&prcopt_.prn[5],     "m"    },
    {"stats-clkstab",   1,  (void *)&prcopt_.sclkstab,   "s/s"  },
    
    {"ant1-postype",    3,  (void *)&antpostype_[0],     POSOPT },
    {"ant1-pos1",       1,  (void *)&antpos_[0][0],      "deg|m"},
    {"ant1-pos2",       1,  (void *)&antpos_[0][1],      "deg|m"},
    {"ant1-pos3",       1,  (void *)&antpos_[0][2],      "m|m"  },
    {"ant1-anttype",    2,  (void *)prcopt_.anttype[0],  ""     },
    {"ant1-antdele",    1,  (void *)&prcopt_.antdel[0][0],"m"   },
    {"ant1-antdeln",    1,  (void *)&prcopt_.antdel[0][1],"m"   },
    {"ant1-antdelu",    1,  (void *)&prcopt_.antdel[0][2],"m"   },
    
    {"ant2-postype",    3,  (void *)&antpostype_[1],     POSOPT },
    {"ant2-pos1",       1,  (void *)&antpos_[1][0],      "deg|m"},
    {"ant2-pos2",       1,  (void *)&antpos_[1][1],      "deg|m"},
    {"ant2-pos3",       1,  (void *)&antpos_[1][2],      "m|m"  },
    {"ant2-anttype",    2,  (void *)prcopt_.anttype[1],  ""     },
    {"ant2-antdele",    1,  (void *)&prcopt_.antdel[1][0],"m"   },
    {"ant2-antdeln",    1,  (void *)&prcopt_.antdel[1][1],"m"   },
    {"ant2-antdelu",    1,  (void *)&prcopt_.antdel[1][2],"m"   },
    {"ant2-maxaveep",   0,  (void *)&prcopt_.maxaveep    ,""    },
    {"ant2-initrst",    3,  (void *)&prcopt_.initrst,    SWTOPT },
    
    {"misc-timeinterp", 3,  (void *)&prcopt_.intpref,    SWTOPT },
    {"misc-sbasatsel",  0,  (void *)&prcopt_.sbassatsel, "0:all"},
    {"misc-rnxopt1",    2,  (void *)prcopt_.rnxopt[0],   ""     },
    {"misc-rnxopt2",    2,  (void *)prcopt_.rnxopt[1],   ""     },
    {"misc-pppopt",     2,  (void *)prcopt_.pppopt,      ""     },
    
//    {"file-satantfile", 2,  (void *)&filopt_.satantp,    ""     },
//    {"file-atx",        2,  (void *)&filopt_.atx,    ""     },
//    {"file-staposfile", 2,  (void *)&filopt_.stapos,     ""     },
//    {"file-geoidfile",  2,  (void *)&filopt_.geoid,      ""     },
//    {"file-ionofile",   2,  (void *)&filopt_.iono,       ""     },
//    {"file-dcbfile",    2,  (void *)&filopt_.dcb,        ""     },
//    {"file-eopfile",    2,  (void *)&filopt_.eop,        ""     },
//    {"file-blqfile",    2,  (void *)&filopt_.blq,        ""     },
//    {"file-tempdir",    2,  (void *)&filopt_.tempdir,    ""     },
//    {"file-geexefile",  2,  (void *)&filopt_.geexe,      ""     },
//    {"file-solstatfile",2,  (void *)&filopt_.solstat,    ""     },
//    {"file-tracefile",  2,  (void *)&filopt_.trace,      ""     },
//    {"file-navfile",    2,  (void *)&filopt_.navfile,    ""     },
//    {"file-sp3file",    2,  (void *)&filopt_.sp3file,    ""     },
//    {"file-clkfile",    2,  (void *)&filopt_.clkfile,    ""     },
//    {"file-imupfile",   2,  (void *)&filopt_.imupf,      ""     },

    {"inpstr1-type",    3,  (void *)&filopt_.strtype[0],  ISTOPT},
    {"inpstr2-type",    3,  (void *)&filopt_.strtype[1],  ISTOPT },
    {"inpstr3-type",    3,  (void *)&filopt_.strtype[2],  ISTOPT },
    {"inpstr4-type",    3,  (void *)&filopt_.strtype[3],  ISTOPT },
    {"inpstr5-type",    3,  (void *)&filopt_.strtype[4],  ISTOPT },
    {"inpstr6-type",    3,  (void *)&filopt_.strtype[5],  ISTOPT },
    {"inpstr7-type",    3,  (void *)&filopt_.strtype[6],  ISTOPT },
    {"inpstr1-path",    2,  (void *)&filopt_.strpath[0],  ""     },
    {"inpstr2-path",    2,  (void *)&filopt_.strpath[1],  ""     },
    {"inpstr3-path",    2,  (void *)&filopt_.strpath[2],  ""     },
    {"inpstr4-path",    2,  (void *)&filopt_.strpath[3],  ""     },
    {"inpstr5-path",    2,  (void *)&filopt_.strpath[4],  ""     },
    {"inpstr6-path",    2,  (void *)&filopt_.strpath[5],  ""     },
    {"inpstr7-path",    2,  (void *)&filopt_.strpath[6],  ""     },
    {"inpstr1-format",  3,  (void *)&filopt_.strfmt[0],   FMTOPT },
    {"inpstr2-format",  3,  (void *)&filopt_.strfmt[1],   FMTOPT },
    {"inpstr3-format",  3,  (void *)&filopt_.strfmt[2],   FMTOPT },
    {"inpstr4-format",  3,  (void *)&filopt_.strfmt[3],   FMTOPT },
    {"inpstr5-format",  3,  (void *)&filopt_.strfmt[4],   FMTOPT },
    {"inpstr6-format",  3,  (void *)&filopt_.strfmt[5],   FMTOPT },
    {"inpstr7-format",  3,  (void *)&filopt_.strfmt[6],   FMTOPT },
    {"outstr1-type",    3,  (void *)&filopt_.strtype[7],  OSTOPT },
    {"outstr2-type",    3,  (void *)&filopt_.strtype[8],  OSTOPT },
    {"outstr1-path",    2,  (void *)&filopt_.strpath[7],  ""     },
    {"outstr2-path",    2,  (void *)&filopt_.strpath[8],  ""     },
    {"outstr1-format",  3,  (void *)&filopt_.strfmt[7],   SOLOPT },
    {"outstr2-format",  3,  (void *)&filopt_.strfmt[8],   SOLOPT },

    {"",0,NULL,""} /* terminator */
};

EXPORT opt_t insopts[]={
        {"ins-lcposfmt",  3,(void *)&prcopt_.insopt.posfmt,LCPOSFMT},
        {"ins-mech_coord",3,(void *)&prcopt_.insopt.mech_coord,MECHCOORD},
        {"ins-local_coord",2,  (void *)&imucoord,""},
        {"ins-align",      3,  (void *)&prcopt_.insopt.imu_align,   IMUALIGN},
        {"ins-multi_sample", 0,  (void *)&prcopt_.insopt.ms,   "" },
        {"ins-max_ny",    0,  (void *)&prcopt_.insopt.max_ny,   "" },
        {"ins-est_ba",     3,  (void *)&prcopt_.insopt.est_ba,     SWTOPT },
        {"ins-est_bg",     3,  (void *)&prcopt_.insopt.est_bg,     SWTOPT },
        {"ins-est_sa",     3,  (void *)&prcopt_.insopt.est_sa,     SWTOPT },
        {"ins-est_sg",     3,  (void *)&prcopt_.insopt.est_sg,     SWTOPT },
        {"ins-est_kod",    3,  (void *)&prcopt_.insopt.est_kod,    SWTOPT},
        {"ins-est_armgps", 3,  (void *)&prcopt_.insopt.est_armgps, SWTOPT},
        {"ins-est_igdt",   3,  (void *)&prcopt_.insopt.est_igdt, SWTOPT},

        {"ins-feedratio",  1,  (void *)&prcopt_.insopt.feedratio,    "" },
        {"ins-fb_pos",    3,  (void *)&prcopt_.insopt.fb_pos,    SWTOPT},
        {"ins-fb_vel",    3,  (void *)&prcopt_.insopt.fb_vel,    SWTOPT},
        {"ins-fb_att",    3,  (void *)&prcopt_.insopt.fb_att,    SWTOPT},
        {"ins-fb_ba",     3,  (void *)&prcopt_.insopt.fb_ba,    SWTOPT},
        {"ins-fb_bg",     3,  (void *)&prcopt_.insopt.fb_bg,    SWTOPT},
        {"ins-fb_sa",     3,  (void *)&prcopt_.insopt.fb_sa,    SWTOPT},
        {"ins-fb_sg",     3,  (void *)&prcopt_.insopt.fb_sg,    SWTOPT},
        {"ins-fb_armgps", 3,  (void *)&prcopt_.insopt.fb_armgps,SWTOPT},
        {"ins-fb_igdt",   3,  (void *)&prcopt_.insopt.fb_igdt,  SWTOPT},

        {"ins-zvu",     3,  (void *)&prcopt_.insopt.zupt,    SWTOPT},
        {"ins-zvu_var",     1,  (void *)&prcopt_.insopt.zupt_var,    ""},
        {"ins-zaru",     3,  (void *)&prcopt_.insopt.zaru,    SWTOPT},
        {"ins-zaru_var",     1,  (void *)&prcopt_.insopt.zaru_var,    ""},
        {"ins-nhc",     3,  (void *)&prcopt_.insopt.nhc,    SWTOPT},
        {"ins-nhc_var",     1,  (void *)&prcopt_.insopt.nhc_var,    ""},
        {"ins-aid_cs",     3,  (void *)&prcopt_.insopt.aid_cs,    INSCYCLE},
        {"ins-slip_scale",     1,  (void *)&prcopt_.insopt.scslp,    ""},
        {"ins-velopt",     3,  (void *)&prcopt_.insopt.velopt,     VELOPT },

        {"ins-rts",    0,  (void *)&prcopt_.insopt.rts,    "" },
        {"ins-zv_det",  3,  (void *)&prcopt_.insopt.zvopt.det,    IMUSTDET},
        {"ins-zv_mt",  0,  (void *)&prcopt_.insopt.zvopt.mt,    "" },
        {"ins-zv_ws",  0,  (void *)&prcopt_.insopt.zvopt.ws,    "" },
        {"ins-zv_gthres",  1,  (void *)&prcopt_.insopt.zvopt.gthres,    "" },
        {"ins-zv_athres1",  1,  (void *)&prcopt_.insopt.zvopt.athres[0],    "" },
        {"ins-zv_athres2",  1,  (void *)&prcopt_.insopt.zvopt.athres[1],    "" },
        {"ins-zv_athres3",  1,  (void *)&prcopt_.insopt.zvopt.athres[2],    "" },
        {"ins-zv_gyrothres1",  1,  (void *)&prcopt_.insopt.zvopt.gyrothres[0],    "" },
        {"ins-zv_gyrothres2",  1,  (void *)&prcopt_.insopt.zvopt.gyrothres[1],    "" },
        {"ins-zv_gyrothres3",  1,  (void *)&prcopt_.insopt.zvopt.gyrothres[2],    "" },
        {"ins-zv_siga",  1,  (void *)&prcopt_.insopt.zvopt.sig_a,    "" },
        {"ins-zv_sigg",  1,  (void *)&prcopt_.insopt.zvopt.sig_g,    "" },
        {"ins-zv_gamma1",1,  (void *)&prcopt_.insopt.zvopt.gamma[0],    "" },
        {"ins-zv_gamma2",1,  (void *)&prcopt_.insopt.zvopt.gamma[1],    "" },
        {"ins-zv_gamma3",1,  (void *)&prcopt_.insopt.zvopt.gamma[2],    "" },
        {"ins-zv_gamma4",1,  (void *)&prcopt_.insopt.zvopt.gamma[3],    "" },

        {"ins-gps_loss_s",1,  (void *)&prcopt_.insopt.gps_loss_s,    "" },
        {"ins-gps_loss_last",1,  (void *)&prcopt_.insopt.gps_loss_last,    "" },

        {"",0,NULL,""}, /* terminator */
};

/* discard space characters at tail ------------------------------------------*/
static void chop(char *str)
{
    char *p;
    if ((p=strchr(str,'#'))) *p='\0'; /* comment */
    for (p=str+strlen(str)-1;p>=str&&!isgraph((int)*p);p--) *p='\0';
}
/* enum to string ------------------------------------------------------------*/
static int enum2str(char *s, const char *comment, int val)
{
    char str[32],*p,*q;
    int n;
    
    n=sprintf(str,"%d:",val);
    if (!(p=strstr(comment,str))) {
        return sprintf(s,"%d",val);
    }
    if (!(q=strchr(p+n,','))&&!(q=strchr(p+n,')'))) {
        strcpy(s,p+n);
        return (int)strlen(p+n);
    }
    strncpy(s,p+n,q-p-n); s[q-p-n]='\0';
    return (int)(q-p-n);
}
/* string to enum ------------------------------------------------------------*/
static int str2enum(const char *str, const char *comment, int *val)
{
    const char *p;
    char s[32];
    
    for (p=comment;;p++) {
       if (!(p=strstr(p,str))) break;
       if (*(p-1)!=':') continue;
       for (p-=2;'0'<=*p&&*p<='9';p--) ;
       return sscanf(p+1,"%d",val)==1;
    }
    sprintf(s,"%30.30s:",str);
    if ((p=strstr(comment,s))) { /* number */
        return sscanf(p,"%d",val)==1;
    }
    return 0;
}
/* search option ---------------------------------------------------------------
* search option record
* args   : char   *name     I  option name
*          opt_t  *opts     I  options table
*                              (terminated with table[i].name="")
* return : option record (NULL: not found)
*-----------------------------------------------------------------------------*/
extern opt_t *searchopt(const char *name, const opt_t *opts)
{
    int i=0;
    
    trace(4,"searchopt: name=%s\n",name);
    
    for (i=0;*opts[i].name;i++) {
        if (strstr(opts[i].name,name)) return (opt_t *)(opts+i);
    }
    return NULL;
}
/* string to option value ------------------------------------------------------
* convert string to option value
* args   : opt_t  *opt      O  option
*          char   *str      I  option value string
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int str2opt(opt_t *opt, const char *str)
{
    switch (opt->format) {
        case 0: *(int    *)opt->var=atoi(str); break;
        case 1: *(double *)opt->var=atof(str); break;
        case 2: strcpy((char *)opt->var,str);  break;
        case 3: return str2enum(str,opt->comment,(int *)opt->var);
        default: return 0;
    }
    return 1;
}
/* option value to string ------------------------------------------------------
* convert option value to string
* args   : opt_t  *opt      I  option
*          char   *str      O  option value string
* return : length of output string
*-----------------------------------------------------------------------------*/
extern int opt2str(const opt_t *opt, char *str)
{
    char *p=str;
    
    trace(3,"opt2str : name=%s\n",opt->name);
    
    switch (opt->format) {
        case 0: p+=sprintf(p,"%d"   ,*(int   *)opt->var); break;
        case 1: p+=sprintf(p,"%.15g",*(double*)opt->var); break;
        case 2: p+=sprintf(p,"%s"   , (char  *)opt->var); break;
        case 3: p+=enum2str(p,opt->comment,*(int *)opt->var); break;
    }
    return (int)(p-str);
}
/* option to string -------------------------------------------------------------
* convert option to string (keyword=value # comment)
* args   : opt_t  *opt      I  option
*          char   *buff     O  option string
* return : length of output string
*-----------------------------------------------------------------------------*/
extern int opt2buf(const opt_t *opt, char *buff)
{
    char *p=buff;
    int n;
    
    trace(3,"opt2buf : name=%s\n",opt->name);
    
    p+=sprintf(p,"%-18s =",opt->name);
    p+=opt2str(opt,p);
    if (*opt->comment) {
        if ((n=(int)(buff+30-p))>0) p+=sprintf(p,"%*s",n,"");
        p+=sprintf(p," # (%s)",opt->comment);
    }
    return (int)(p-buff);
}
/* load options ----------------------------------------------------------------
* load options from file
* args   : char   *file     I  options file
*          opt_t  *opts     IO options table
*                              (terminated with table[i].name="")
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int loadopts(const char *file, opt_t *opts)
{
    FILE *fp;
    opt_t *opt;
    char buff[2048],*p;
    int n=0;
    
    trace(4,"loadopts: file=%s\n",file);
    
    if (!(fp=fopen(file,"r"))) {
        trace(1,"loadopts: options file open error (%s)\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        n++;
        chop(buff);
        
        if (buff[0]=='\0') continue;
        
        if (!(p=strstr(buff,"="))) {
            fprintf(stderr,"invalid option %s (%s:%d)\n",buff,file,n);
            continue;
        }
        *p++='\0';
        chop(buff);
        if (!(opt=searchopt(buff,opts))) continue;
        
        if (!str2opt(opt,p)) {
            fprintf(stderr,"invalid option value %s (%s:%d)\n",buff,file,n);
            continue;
        }
    }
    fclose(fp);
    
    return 1;
}

#ifdef MATLAB
extern int loadmatopts(MATFile *pmat)
{
    int char_len = 50,i=0;
    double *num_opt,*num_opt1;
    char   char_opt[50]={'\0'};

    mxArray *pMxArray = matGetVariable(pmat,"iRTKLIB_conf");
    mxArray *pStructOpt = NULL;

    char_opt[0]='\0';
    pStructOpt = mxGetField(pMxArray,0,"prc_dir");
    mxGetString(pStructOpt,char_opt,50);
    strcpy(prcopt_.prcdir,char_opt);

    /// process mode
    pStructOpt = mxGetField(pMxArray,0,"prc_mode");
    mxGetString(pStructOpt,char_opt,char_len);
    for(i=0;kpmodestr;i++){
        if(!strcmp(char_opt,kpmodestr[i])){
            prcopt_.mode=i;
            break;
        }
    }

    /// process frequency number
    pStructOpt = mxGetField(pMxArray,0,"prc_frq");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.nf = (int)*num_opt;

    /// process filter type
    pStructOpt = mxGetField(pMxArray,0,"prc_filtertype");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.soltype = (int)*num_opt-1;

    /// process trace level
    pStructOpt = mxGetField(pMxArray,0,"prc_trace");
    num_opt =  (double *)mxGetData(pStructOpt);
    solopt_.trace = (int)*num_opt;

    /// process time start
    char_opt[0]='\0';
    pStructOpt = mxGetField(pMxArray,0,"prc_ts");
    mxGetString(pStructOpt,char_opt,50);
    if(!strcmp(char_opt,"")){
        prcopt_.ts.time=0.0;prcopt_.ts.sec=0.0;
    }
    else{
        str2time1(char_opt,0,18,&prcopt_.ts);
    }

    /// process time end
    char_opt[0]='\0';
    pStructOpt = mxGetField(pMxArray,0,"prc_te");
    mxGetString(pStructOpt,char_opt,50);
    if(!strcmp(char_opt,"")){
        prcopt_.te.time=0.0;prcopt_.te.sec=0.0;
    }
    else{
        str2time1(char_opt,0,18,&prcopt_.te);
    }

    prcopt_.navsys=SYS_GPS;
    /// GPS frq idx
    pStructOpt = mxGetField(pMxArray,0,"gnss_GPS");
    num_opt =  (double *)mxGetData(pStructOpt);
    if(*num_opt==1){
        prcopt_.navsys|=SYS_GPS;
        char_opt[0]='\0';
        pStructOpt = mxGetField(pMxArray,0,"gnss_GPS_frq");
        mxGetString(pStructOpt,GPS_frq_,20);
    }

    /// GLO frq idx
    pStructOpt = mxGetField(pMxArray,0,"gnss_GLO");
    num_opt =  (double *)mxGetData(pStructOpt);
    if(*num_opt==1){
        prcopt_.navsys|=SYS_GLO;
        char_opt[0]='\0';
        pStructOpt = mxGetField(pMxArray,0,"gnss_GLO_frq");
        mxGetString(pStructOpt,GLO_frq_,20);
    }

    /// BD2 frq idx
    pStructOpt = mxGetField(pMxArray,0,"gnss_BD2");
    num_opt =  (double *)mxGetData(pStructOpt);
    if(*num_opt==1){
        prcopt_.navsys|=SYS_CMP;

        char_opt[0]='\0';
        pStructOpt = mxGetField(pMxArray,0,"gnss_BD2_frq");
        mxGetString(pStructOpt,BD2_frq_,20);
    }

    /// BD3 frq idx
    pStructOpt = mxGetField(pMxArray,0,"gnss_BD3");
    num_opt =  (double *)mxGetData(pStructOpt);
    if(*num_opt==1){
        prcopt_.navsys|=SYS_BD3;
        char_opt[0]='\0';
        pStructOpt = mxGetField(pMxArray,0,"gnss_BD3_frq");
        mxGetString(pStructOpt,BD3_frq_,20);
    }

    /// GAL frq idx
    pStructOpt = mxGetField(pMxArray,0,"gnss_GAL");
    num_opt =  (double *)mxGetData(pStructOpt);
    if(*num_opt==1){
        prcopt_.navsys|=SYS_GAL;

        char_opt[0]='\0';
        pStructOpt = mxGetField(pMxArray,0,"gnss_GAL_frq");
        mxGetString(pStructOpt,GAL_frq_,20);
    }

    /// QZS
    char_opt[0] = '\0';
    pStructOpt = mxGetField(pMxArray,0,"gnss_QZS");
    num_opt =  (double *)mxGetData(pStructOpt);
    if(*num_opt==1){
        prcopt_.navsys|=SYS_QZS;

        char_opt[0]='\0';
        pStructOpt = mxGetField(pMxArray,0,"gnss_QZS_frq");
        mxGetString(pStructOpt,QZS_frq_,20);
    }

    /// gnss sample rate
    pStructOpt = mxGetField(pMxArray,0,"gnss_sample_rate");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.sample = *num_opt;

    /// gnss eph opt
    pStructOpt = mxGetField(pMxArray,0,"gnss_pre_eph");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.sateph = (int)*num_opt;

    /// gnss rec dynamic
    pStructOpt = mxGetField(pMxArray,0,"gnss_rec_dym");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.dynamics = (int)*num_opt;

    /// gnss raim
    pStructOpt = mxGetField(pMxArray,0,"gnss_raim");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.posopt[4] = (int)*num_opt;

    /// gnss iono model
    pStructOpt = mxGetField(pMxArray,0,"gnss_ion_model");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.ionoopt = (int)*num_opt;

    /// gnss trop model
    pStructOpt = mxGetField(pMxArray,0,"gnss_trp_model");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.tropopt = (int)*num_opt;

    /// gnss trop map

    /// gnss tide model
    pStructOpt = mxGetField(pMxArray,0,"gnss_tid_model");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.tidecorr = (int)*num_opt;

    /// ant type
    prcopt_.anttype[0][0]='*';
    prcopt_.anttype[1][0]='*';
    /// gnss ant2
    pStructOpt = mxGetField(pMxArray,0,"gnss_ant2pos_type");
    num_opt =  (double *)mxGetData(pStructOpt);
    antpostype_[1]=(int)*num_opt;
    pStructOpt = mxGetField(pMxArray,0,"gnss_ant2pos1");
    antpos_[1][0] =  *(double *)mxGetData(pStructOpt);
    pStructOpt = mxGetField(pMxArray,0,"gnss_ant2pos2");
    antpos_[1][1] =  *(double *)mxGetData(pStructOpt);
    pStructOpt = mxGetField(pMxArray,0,"gnss_ant2pos3");
    antpos_[1][2] =  *(double *)mxGetData(pStructOpt);

    /// gnss sync sol
    pStructOpt = mxGetField(pMxArray,0,"gnss_ant2pos_type");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.syncsol=(int)*num_opt;

    /// gnss max age
    pStructOpt = mxGetField(pMxArray,0,"gnss_maxage");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.maxtdiff = *num_opt;

    /// gnss ac
    char_opt[0]='\0';
    pStructOpt = mxGetField(pMxArray,0,"gnss_ac_name");
    mxGetString(pStructOpt,char_opt,50);
    strcpy(prcopt_.ac_name,char_opt);

    /// gnss sat pcv
    pStructOpt = mxGetField(pMxArray,0,"gnss_sat_pcv");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.posopt[0] = (int)*num_opt;

    /// gnss rec pcv
    pStructOpt = mxGetField(pMxArray,0,"gnss_rec_pcv");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.posopt[1] = (int)*num_opt;

    /// gnss phw
    pStructOpt = mxGetField(pMxArray,0,"gnss_phw");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.posopt[2] = (int)*num_opt;

    /// gnss eclipsing
    prcopt_.posopt[3] = 1;

    /// gnss exclude sats
    char_opt[0]='\0';
    pStructOpt = mxGetField(pMxArray,0,"gnss_excsat");
    mxGetString(pStructOpt,char_opt,50);
    strcpy(exsats_,char_opt);

    /// gnss mw thres
    pStructOpt = mxGetField(pMxArray,0,"gnss_mwthres");
    num_opt =  (double *)mxGetData(pStructOpt);

    /// gnss gf thres
    pStructOpt = mxGetField(pMxArray,0,"gnss_gfthres");
    num_opt =  (double *)mxGetData(pStructOpt);

    /// gnss rejc ionno
    pStructOpt = mxGetField(pMxArray,0,"gnss_rejc_ionno");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.maxinno = *num_opt;

    /// gnss rejcdop
    pStructOpt = mxGetField(pMxArray,0,"gnss_rejc_dop");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.maxgdop = *num_opt;

    /// gnss error ratio1
    pStructOpt = mxGetField(pMxArray,0,"gnss_ratio1");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.eratio[0] = *num_opt;

    /// gnss error ratio1
    pStructOpt = mxGetField(pMxArray,0,"gnss_ratio2");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.eratio[1] = *num_opt;

    /// gnss error phase
    pStructOpt = mxGetField(pMxArray,0,"gnss_ratio3");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.err[1] = *num_opt;

    /// gnss error phaseel
    pStructOpt = mxGetField(pMxArray,0,"gnss_ratio4");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.err[2] = *num_opt;

    /// gnss error phasebl
    pStructOpt = mxGetField(pMxArray,0,"gnss_ratio5");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.err[3] = *num_opt;

    /// gnss error doppler
    pStructOpt = mxGetField(pMxArray,0,"gnss_ratio6");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.err[4] = *num_opt;

    /// gnss std_bias
    pStructOpt = mxGetField(pMxArray,0,"gnss_stdbias");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.std[0]=*num_opt;

    /// gnss std_iono
    pStructOpt = mxGetField(pMxArray,0,"gnss_stdIono");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.std[1]=*num_opt;

    /// gnss std_trop
    pStructOpt = mxGetField(pMxArray,0,"gnss_stdTrp");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.std[2]=*num_opt;

    /// gnss prn_bias
    pStructOpt = mxGetField(pMxArray,0,"gnss_prnbias");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.prn[0]=*num_opt;

    /// gnss prn_iono
    pStructOpt = mxGetField(pMxArray,0,"gnss_prnIono");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.prn[1]=*num_opt;

    /// gnss prn_trp
    pStructOpt = mxGetField(pMxArray,0,"gnss_prnTrp");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.prn[2]=*num_opt;

    /// gnss ar mode
    pStructOpt = mxGetField(pMxArray,0,"gnss_armode");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.modear=(int)*num_opt;

    /// gnss ar filter
    pStructOpt = mxGetField(pMxArray,0,"gnss_arfilter");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.arfilter=(int)*num_opt;

    /// gnss ar elmask
    pStructOpt = mxGetField(pMxArray,0,"gnss_ar_elmask");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.elmaskar=*num_opt;

    /// gnss glo ar
    pStructOpt = mxGetField(pMxArray,0,"gnss_ar_glo");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.glomodear = (int)*num_opt;

    /// gnss bds ar
    pStructOpt = mxGetField(pMxArray,0,"gnss_ar_bds");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.bdsmodear = (int)*num_opt;

    /// gnss ppp ar prod
    pStructOpt = mxGetField(pMxArray,0,"gnss_ar_pppar_prod");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.arprod = (int)*num_opt;

    /// gnss ar ratio
    pStructOpt = mxGetField(pMxArray,0,"gnss_ar_ratio");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.thresar[0] = *num_opt;

    /// gnss ar success rate
    pStructOpt = mxGetField(pMxArray,0,"gnss_ar_success_rate");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.thresar[1] = *num_opt;

    /// gnss ar wl res
    pStructOpt = mxGetField(pMxArray,0,"ar_wlres");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.thresar[2] = *num_opt;

    /// gnss ar nl res
    pStructOpt = mxGetField(pMxArray,0,"ar_nlres");
    num_opt =  (double *)mxGetData(pStructOpt);

    /// gnss min sat2fix
    pStructOpt = mxGetField(pMxArray,0,"ar_minsat");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.minfixsats = (int)*num_opt;

    /// gnss ar lock count
    pStructOpt = mxGetField(pMxArray,0,"ar_lockcnt");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.minlock = (int)*num_opt;

    /// gnss ar minsat2drop
    pStructOpt = mxGetField(pMxArray,0,"ar_minsat2drop");
    num_opt =  (double *)mxGetData(pStructOpt);
    prcopt_.mindropsats = (int)*num_opt;

    /// gnss solution opt
    pStructOpt = mxGetField(pMxArray,0,"sol_fmt");
    num_opt =  (double *)mxGetData(pStructOpt);
    solopt_.posf = (int)*num_opt;

    pStructOpt = mxGetField(pMxArray,0,"sol_time_fmt");
    num_opt =  (double *)mxGetData(pStructOpt);
    solopt_.timef = (int)*num_opt;

    pStructOpt = mxGetField(pMxArray,0,"sol_time_dec");
    num_opt =  (double *)mxGetData(pStructOpt);
    solopt_.timeu = (int)*num_opt;

    pStructOpt = mxGetField(pMxArray,0,"sol_height");
    num_opt =  (double *)mxGetData(pStructOpt);
    solopt_.height = (int)*num_opt;

    pStructOpt = mxGetField(pMxArray,0,"sol_deg_fmt");
    num_opt =  (double *)mxGetData(pStructOpt);
    solopt_.degf = (int)*num_opt;

    pStructOpt = mxGetField(pMxArray,0,"sol_time_sys");
    num_opt =  (double *)mxGetData(pStructOpt);
    solopt_.times = (int)*num_opt;

    char_opt[0]='\0';
    pStructOpt = mxGetField(pMxArray,0,"sol_sep");
    mxGetString(pStructOpt,char_opt,50);
    strcpy(solopt_.sep,char_opt);

    pStructOpt = mxGetField(pMxArray,0,"sol_out_head");
    num_opt =  (double *)mxGetData(pStructOpt);
    solopt_.outhead = (int)*num_opt;

    pStructOpt = mxGetField(pMxArray,0,"sol_out_opt");
    num_opt =  (double *)mxGetData(pStructOpt);
    solopt_.outopt = (int)*num_opt;

    pStructOpt = mxGetField(pMxArray,0,"sol_out_ba");
    num_opt =  (double *)mxGetData(pStructOpt);

    pStructOpt = mxGetField(pMxArray,0,"sol_out_bg");
    num_opt =  (double *)mxGetData(pStructOpt);

    /// ins options-imu_type
    if(prcopt_.mode>=PMODE_INS_MECH){
        char_opt[0]='\0';
        pStructOpt = mxGetField(pMxArray,0,"ins_imu_type");
        mxGetString(pStructOpt,char_opt,50);
        if(!strcmp(char_opt,"NVT_CPT")){
            prcopt_.insopt.imup.strfmt=STRFMT_IMU_NVT_CPT;
        }

        /// ins options-imu_coord_type
        pStructOpt = mxGetField(pMxArray,0,"ins_coord_type");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.imucoors = (int)*num_opt;

        /// ins options-ins rate
        pStructOpt = mxGetField(pMxArray,0,"ins_rate");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.freq_imu = (int)*num_opt;

        /// ins options-imu_dec_fmt
        pStructOpt = mxGetField(pMxArray,0,"ins_decfmt");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.imudecfmt = (int)*num_opt;

        /// ins options-imu_val_fmt
        pStructOpt = mxGetField(pMxArray,0,"ins_val_fmt");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.imuvalfmt = (int)*num_opt;

        /// ins options-imu_align method
        char_opt[0]='\0';
        pStructOpt = mxGetField(pMxArray,0,"ins_align");
        mxGetString(pStructOpt,char_opt,50);
        if(!strcmp(char_opt,"PPK")){
//        prcopt_.insopt.imup.strfmt=STRFMT_IMU_NVT_CPT;
        }

        /// ins options-imu_est_GpsArm
        pStructOpt = mxGetField(pMxArray,0,"ins_est_GpsArm");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.est_armgps = (int)*num_opt;

        /// ins options-arm
        if(prcopt_.insopt.est_armgps==0){
            pStructOpt = mxGetField(pMxArray,0,"ins_GpsArm1");
            num_opt =  (double *)mxGetData(pStructOpt);
            prcopt_.insopt.imup.lever_arm_gps.x = *num_opt;

            pStructOpt = mxGetField(pMxArray,0,"ins_GpsArm2");
            num_opt =  (double *)mxGetData(pStructOpt);
            prcopt_.insopt.imup.lever_arm_gps.y = *num_opt;

            pStructOpt = mxGetField(pMxArray,0,"ins_GpsArm3");
            num_opt =  (double *)mxGetData(pStructOpt);
            prcopt_.insopt.imup.lever_arm_gps.z = *num_opt;
        }

        pStructOpt = mxGetField(pMxArray,0,"ins_est_Ba");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.est_ba = (int)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_est_Bg");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.est_bg = (int)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_est_Ka");
        num_opt =  (double *)mxGetData(pStructOpt);
//        prcopt_.insopt.est_kax=prcopt_.insopt.est_kay=prcopt_.insopt.est_kaz = (int)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_est_Kg");
        num_opt =  (double *)mxGetData(pStructOpt);
//        prcopt_.insopt.est_kgx=prcopt_.insopt.est_kgy=prcopt_.insopt.est_kgz = (int)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_std_pos1");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.initr.x = (double)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_std_pos2");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.initr.y = (double)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_std_pos3");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.initr.z = (double)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_std_vel1");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.initv.x = (double)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_std_vel2");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.initv.y = (double)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_std_vel3");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.initv.z = (double)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_std_att1");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.inita.x = (double)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_std_att2");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.inita.y = (double)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_std_att3");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.inita.z = (double)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_std_Ba1");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.ba_err.x = (double)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_std_Ba2");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.ba_err.y = (double)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_std_Ba3");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.ba_err.z = (double)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_std_Bg1");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.bg_err.x = (double)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_std_Bg2");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.bg_err.y = (double)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_std_Bg3");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.bg_err.z = (double)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_arw");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.gyrnd.x=prcopt_.insopt.imup.gyrnd.y=prcopt_.insopt.imup.gyrnd.z = (double)*num_opt*DPH2RPS;

        pStructOpt = mxGetField(pMxArray,0,"ins_arrw");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.gbrw.x=prcopt_.insopt.imup.gbrw.y=prcopt_.insopt.imup.gbrw.z = (double)*num_opt;

        pStructOpt = mxGetField(pMxArray,0,"ins_vrw");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.accnd.x=prcopt_.insopt.imup.accnd.y=prcopt_.insopt.imup.accnd.z = (double)*num_opt*DPH2RPS;

        pStructOpt = mxGetField(pMxArray,0,"ins_vrrw");
        num_opt =  (double *)mxGetData(pStructOpt);
        prcopt_.insopt.imup.abrw.x=prcopt_.insopt.imup.abrw.y=prcopt_.insopt.imup.abrw.z = (double)*num_opt;
    }

    return 1;
}
#endif
/* save options to file --------------------------------------------------------
* save options to file
* args   : char   *file     I  options file
*          char   *mode     I  write mode ("w":overwrite,"a":append);
*          char   *comment  I  header comment (NULL: no comment)
*          opt_t  *opts     I  options table
*                              (terminated with table[i].name="")
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int saveopts(const char *file, const char *mode, const char *comment,
                    const opt_t *opts)
{
    FILE *fp;
    char buff[2048];
    int i;
    
    trace(3,"saveopts: file=%s mode=%s\n",file,mode);
    
    if (!(fp=fopen(file,mode))) {
        trace(1,"saveopts: options file open error (%s)\n",file);
        return 0;
    }
    if (comment) fprintf(fp,"# %s\n\n",comment);
    
    for (i=0;*opts[i].name;i++) {
        opt2buf(opts+i,buff);
        fprintf(fp,"%s\n",buff);
    }
    fclose(fp);
    return 1;
}
/* system options buffer to options ------------------------------------------*/
static void buff2sysopts(void)
{
    double pos[3],*rr;
    char buff[1024],*p,*id;
    int i,j,sat,*ps;
    
    prcopt_.elmin     =elmask_    *D2R;
    prcopt_.elmaskar  =elmaskar_  *D2R;
    prcopt_.elmaskhold=elmaskhold_*D2R;

    if(prc_ts_[0]){
        str2time(prc_ts_,0,25,&prcopt_.ts);
    }
    if(prc_te_[0]){
        str2time(prc_te_,0,25,&prcopt_.te);
    }

    if(GPS_frq_[0]) getobsfrqidx(GPS_frq_,SYS_GPS,prcopt_.nf,prcopt_.gnss_frq_idx[0]);
#ifdef ENAGLO
    if(GLO_frq_[0]) getobsfrqidx(GLO_frq_,SYS_GLO,prcopt_.nf,prcopt_.gnss_frq_idx[1]);
#endif
#ifdef ENAGAL
    if(GAL_frq_[0]) getobsfrqidx(GAL_frq_,SYS_GAL,prcopt_.nf,prcopt_.gnss_frq_idx[2]);
#endif
#ifdef ENACMP
    if(BD2_frq_[0]) getobsfrqidx(BD2_frq_,SYS_CMP,prcopt_.nf,prcopt_.gnss_frq_idx[3]);
    if(BD3_frq_[0]) getobsfrqidx(BD3_frq_,SYS_BD3,prcopt_.nf,prcopt_.gnss_frq_idx[NSYS]);
#endif
#ifdef ENAQZS
    if(QZS_frq_[0]) getobsfrqidx(QZS_frq_,SYS_QZS,prcopt_.nf,prcopt_.gnss_frq_idx[4]);
#endif

    if(imucoord[0]){
        if(!strcmp(imucoord,"ned-frd")){
            prcopt_.insopt.local_coord=INSLOCAL_NED;
            prcopt_.insopt.imu_coord=IMUCOOR_FRD;
        }
        if(!strcmp(imucoord,"enu-rfu")){
            prcopt_.insopt.local_coord=INSLOCAL_ENU;
            prcopt_.insopt.imu_coord=IMUCOOR_RFU;
        }
    }

    for (i=0;i<2;i++) {
        ps=i==0?&prcopt_.rovpos:&prcopt_.refpos;
        rr=i==0?prcopt_.ru:prcopt_.rb;
        
        if (antpostype_[i]==0) { /* lat/lon/hgt */
            *ps=0;
            pos[0]=antpos_[i][0]*D2R;
            pos[1]=antpos_[i][1]*D2R;
            pos[2]=antpos_[i][2];
            pos2ecef(pos,rr);
        }
        else if (antpostype_[i]==1) { /* xyz-ecef */
            *ps=0;
            rr[0]=antpos_[i][0];
            rr[1]=antpos_[i][1];
            rr[2]=antpos_[i][2];
        }
        else *ps=antpostype_[i]-1;
    }
    /* excluded satellites */
    for (i=0;i<MAXSAT;i++) prcopt_.exsats[i]=0;
    if (exsats_[0]!='\0') {
        strcpy(buff,exsats_);
        for (p=strtok(buff," ");p;p=strtok(NULL," ")) {
            if (*p=='+') id=p+1; else id=p;
            if (!(sat=satid2no(id))) continue;
            prcopt_.exsats[sat-1]=*p=='+'?2:1;
        }
    }
    /* snrmask */
    for (i=0;i<NFREQ;i++) {
        for (j=0;j<9;j++) prcopt_.snrmask.mask[i][j]=0.0;
        strcpy(buff,snrmask_[i]);
        for (p=strtok(buff,","),j=0;p&&j<9;p=strtok(NULL,",")) {
            prcopt_.snrmask.mask[i][j++]=atof(p);
        }
    }

}
/* options to system options buffer ------------------------------------------*/
static void sysopts2buff(void)
{
    double pos[3],*rr;
    char id[32],*p;
    int i,j,sat,*ps;
    
    elmask_    =prcopt_.elmin     *R2D;
    elmaskar_  =prcopt_.elmaskar  *R2D;
    elmaskhold_=prcopt_.elmaskhold*R2D;
    
    for (i=0;i<2;i++) {
        ps=i==0?&prcopt_.rovpos:&prcopt_.refpos;
        rr=i==0?prcopt_.ru:prcopt_.rb;
        
        if (*ps==0) {
            antpostype_[i]=0;
            ecef2pos(rr,pos);
            antpos_[i][0]=pos[0]*R2D;
            antpos_[i][1]=pos[1]*R2D;
            antpos_[i][2]=pos[2];
        }
        else antpostype_[i]=*ps+1;
    }
    /* excluded satellites */
    exsats_[0]='\0';
    for (sat=1,p=exsats_;sat<=MAXSAT&&p-exsats_<(int)sizeof(exsats_)-32;sat++) {
        if (prcopt_.exsats[sat-1]) {
            satno2id(sat,id);
            p+=sprintf(p,"%s%s%s",p==exsats_?"":" ",
                       prcopt_.exsats[sat-1]==2?"+":"",id);
        }
    }
    /* snrmask */
    for (i=0;i<NFREQ;i++) {
        snrmask_[i][0]='\0';
        p=snrmask_[i];
        for (j=0;j<9;j++) {
            p+=sprintf(p,"%s%.0f",j>0?",":"",prcopt_.snrmask.mask[i][j]);
        }
    }
}
/* reset system options to default ---------------------------------------------
* reset system options to default
* args   : none
* return : none
*-----------------------------------------------------------------------------*/
extern void resetsysopts(void)
{
    int i,j;
    
    trace(4,"resetsysopts:\n");
    
    prcopt_=prcopt_default;
    solopt_=solopt_default;
//    filopt_.satantp[0]='\0';
    filopt_.atx[0]='\0';
    filopt_.stapos [0]='\0';
    filopt_.geoid  [0]='\0';
    filopt_.dcb    [0]='\0';
    filopt_.blq    [0]='\0';
    filopt_.solstat[0]='\0';
    filopt_.trace  [0]='\0';
    for (i=0;i<2;i++) antpostype_[i]=0;
    elmask_=15.0;
    elmaskar_=0.0;
    elmaskhold_=0.0;
    for (i=0;i<2;i++) for (j=0;j<3;j++) {
        antpos_[i][j]=0.0;
    }
    exsats_[0] ='\0';
}
/* get system options ----------------------------------------------------------
* get system options
* args   : prcopt_t *popt   IO processing options (NULL: no output)
*          solopt_t *sopt   IO solution options   (NULL: no output)
*          folopt_t *fopt   IO file options       (NULL: no output)
* return : none
* notes  : to load system options, use loadopts() before calling the function
*-----------------------------------------------------------------------------*/
extern void getsysopts(prcopt_t *popt, solopt_t *sopt, filopt_t *fopt)
{
    trace(4,"getsysopts:\n");
    
    buff2sysopts();
    if (popt) *popt=prcopt_;
    if (sopt) *sopt=solopt_;
    if (fopt) *fopt=filopt_;
}
/* set system options ----------------------------------------------------------
* set system options
* args   : prcopt_t *prcopt I  processing options (NULL: default)
*          solopt_t *solopt I  solution options   (NULL: default)
*          filopt_t *filopt I  file options       (NULL: default)
* return : none
* notes  : to save system options, use saveopts() after calling the function
*-----------------------------------------------------------------------------*/
extern void setsysopts(const prcopt_t *prcopt, const solopt_t *solopt,
                       const filopt_t *filopt)
{
    trace(3,"setsysopts:\n");
    
    resetsysopts();
    if (prcopt) prcopt_=*prcopt;
    if (solopt) solopt_=*solopt;
    if (filopt) filopt_=*filopt;
    sysopts2buff();
}
