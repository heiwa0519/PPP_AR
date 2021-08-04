/*------------------------------------------------------------------------------
* rtklib.h : rtklib constants, types and function prototypes
*
*          Copyright (C) 2007-2019 by T.TAKASU, All rights reserved.
*
* options : -DENAGLO   enable GLONASS
*           -DENAGAL   enable Galileo
*           -DENAQZS   enable QZSS
*           -DENACMP   enable BeiDou
*           -DENAIRN   enable IRNSS
*           -DNFREQ=n  set number of obs codes/frequencies
*           -DNEXOBS=n set number of extended obs codes
*           -DMAXOBS=n set max number of obs data in an epoch
*           -DEXTLEX   enable QZSS LEX extension
*           -DWIN32    use WIN32 API
*           -DWIN_DLL  generate library as Windows DLL
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/01/13 1.0  rtklib ver.1.0.0
*           2007/03/20 1.1  rtklib ver.1.1.0
*           2008/07/15 1.2  rtklib ver.2.1.0
*           2008/10/19 1.3  rtklib ver.2.1.1
*           2009/01/31 1.4  rtklib ver.2.2.0
*           2009/04/30 1.5  rtklib ver.2.2.1
*           2009/07/30 1.6  rtklib ver.2.2.2
*           2009/12/25 1.7  rtklib ver.2.3.0
*           2010/07/29 1.8  rtklib ver.2.4.0
*           2011/05/27 1.9  rtklib ver.2.4.1
*           2013/03/28 1.10 rtklib ver.2.4.2
*           2016/01/26 1.11 rtklib ver.2.4.3
*-----------------------------------------------------------------------------*/
#ifndef RTKLIB_H
#define RTKLIB_H
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <dirent.h>
#include <stdint.h>
#ifdef WIN32
#include <winsock2.h>
#include <windows.h>
#include "unistd.h"
#ifdef MATLAB
#include "mat.h"
#endif

#define strcasecmp _stricmp
#define access     _access
#else
#include <pthread.h>
#include "unistd.h"
#include <dirent.h>
#include <limits.h>
#include  <sys/stat.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef WIN_DLL
#define EXPORT __declspec(dllexport) /* for Windows DLL */
#else
#define EXPORT
#endif
/* constants -----------------------------------------------------------------*/

#define VER_RTKLIB  "demo"             /* library version */

#define PATCH_LEVEL "b1"               /* patch level */


#define COPYRIGHT_RTKLIB \
            "Copyright (C) 2007-2019 T.Takasu\nAll rights reserved."
#define PI          3.1415926535897932  /* pi */
#define D2R         (PI/180.0)          /* deg to rad */
#define R2D         (180.0/PI)          /* rad to deg */
#define CLIGHT      299792458.0         /* speed of light (m/s) */
#define SC2RAD      3.1415926535898     /* semi-circle to radian (IS-GPS) */
#define AU          149597870691.0      /* 1 AU (m) */
#define AS2R        (D2R/3600.0)
#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))
#define MIN(x,y)    ((x)<=(y)?(x):(y))
#define MAX(x,y)    ((x)>=(y)?(x):(y))

#define OMGE        7.2921151467E-5     /* earth angular velocity (IS-GPS) (rad/s) */

#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */

#define HION        350000.0            /* ionosphere height (m) */

#define MAXFREQ     7                   /* max NFREQ */

#define FREQ0_GPS   1.023E9
#define FREQ1       1.57542E9           /* L1/E1/B1C  frequency (Hz) */
#define FREQ2       1.22760E9           /* L2         frequency (Hz) */
#define FREQ5       1.17645E9           /* L5/E5a/B2a frequency (Hz) */
#define FREQ6       1.27875E9           /* E6/L6  frequency (Hz) */
#define FREQ7       1.20714E9           /* E5b    frequency (Hz) */
#define FREQ8       1.191795E9          /* E5a+b  frequency (Hz) */
#define FREQ9       2.492028E9          /* S      frequency (Hz) */
#define FREQ1_GLO   1.60200E9           /* GLONASS G1 base frequency (Hz) */
#define DFRQ1_GLO   0.56250E6           /* GLONASS G1 bias frequency (Hz/n) */
#define FREQ2_GLO   1.24600E9           /* GLONASS G2 base frequency (Hz) */
#define DFRQ2_GLO   0.43750E6           /* GLONASS G2 bias frequency (Hz/n) */
#define FREQ3_GLO   1.202025E9          /* GLONASS G3 frequency (Hz) */
#define FREQ1a_GLO  1.600995E9          /* GLONASS G1a frequency (Hz) */
#define FREQ2a_GLO  1.248060E9          /* GLONASS G2a frequency (Hz) */
#define FREQ1_CMP   1.561098E9          /* BDS B1I     frequency (Hz) */
#define FREQ2_CMP   1.20714E9           /* BDS B2I/B2b frequency (Hz) */
#define FREQ3_CMP   1.26852E9           /* BDS B3      frequency (Hz) */

#define GPS_FREQ_L1   0x01
#define GPS_FREQ_L2   0x02
#define GPS_FREQ_L5   0x04
#define GLO_FREQ_G1   0x01
#define GLO_FREQ_G2   0x02
#define GLO_FREQ_G3   0x04
#define GLO_FREQ_G1a  0x08
#define GLO_FREQ_G2a  0x10
#define GAL_FREQ_E1   0x01
#define GAL_FREQ_E5b  0x02
#define GAL_FREQ_E5a  0x04
#define GAL_FREQ_E6   0x08
#define GAL_FREQ_E5ab 0x10
#define BDS_FREQ_B1C  0x01
#define BDS_FREQ_B1I  0x02
#define BDS_FREQ_B2I  0x04
#define BDS_FREQ_B2b  0x08
#define BDS_FREQ_B2a  0x10
#define BDS_FREQ_B3   0x20
#define BDS_FREQ_B2ab 0x40
#define QZS_FREQ_L1   0x01
#define QZS_FREQ_L2   0x02
#define QZS_FREQ_L5   0x04
#define QZS_FREQ_L6   0x08

#define EFACT_GPS   1.0                 /* error factor: GPS */
#define EFACT_GLO   2.0                 /* error factor: GLONASS */
#define EFACT_GAL   1.0                 /* error factor: Galileo */
#define EFACT_QZS   1.0                 /* error factor: QZSS */
#define EFACT_CMP   1.0                 /* error factor: BeiDou */
#define EFACT_IRN   1.5                 /* error factor: IRNSS */
#define EFACT_SBS   3.0                 /* error factor: SBAS */

#define SYS_NONE    0x00                /* navigation system: none */
#define SYS_GPS     0x01                /* navigation system: GPS */
#define SYS_SBS     0x02                /* navigation system: SBAS */
#define SYS_GLO     0x04                /* navigation system: GLONASS */
#define SYS_GAL     0x08                /* navigation system: Galileo */
#define SYS_QZS     0x10                /* navigation system: QZSS */
#define SYS_CMP     0x20                /* navigation system: BeiDou */
#define SYS_BD3     0x40
#define SYS_IRN     0x80                /* navigation system: IRNS */
#define SYS_LEO     0x100               /* navigation system: LEO */
#define SYS_ALL     0xFF                /* navigation system: all */

#define TSYS_GPS    0                   /* time system: GPS time */
#define TSYS_UTC    1                   /* time system: UTC */
#define TSYS_GLO    2                   /* time system: GLONASS time */
#define TSYS_GAL    3                   /* time system: Galileo time */
#define TSYS_QZS    4                   /* time system: QZSS time */
#define TSYS_CMP    5                   /* time system: BeiDou time */
#define TSYS_IRN    6                   /* time system: IRNSS time */

#ifndef NFREQ
#define NFREQ       3                   /* number of carrier frequencies */
#endif
#define NFREQGLO    2                   /* number of carrier frequencies of GLONASS */

#ifndef NEXOBS
#define NEXOBS      3                   /* number of extended obs codes */
#endif

#define SNR_UNIT    0.25               /* SNR unit (dBHz) */

#define MINPRNGPS   1                   /* min satellite PRN number of GPS */
#define MAXPRNGPS   32                  /* max satellite PRN number of GPS */
#define NSATGPS     (MAXPRNGPS-MINPRNGPS+1) /* number of GPS satellites */
#define NSYSGPS     1

#ifdef ENAGLO
#define MINPRNGLO   1                   /* min satellite slot number of GLONASS */
#define MAXPRNGLO   27                  /* max satellite slot number of GLONASS */
#define NSATGLO     (MAXPRNGLO-MINPRNGLO+1) /* number of GLONASS satellites */
#define NSYSGLO     1
#else
#define MINPRNGLO   0
#define MAXPRNGLO   0
#define NSATGLO     0
#define NSYSGLO     0
#endif
#ifdef ENAGAL
#define MINPRNGAL   1                   /* min satellite PRN number of Galileo */
#define MAXPRNGAL   36                  /* max satellite PRN number of Galileo */
#define NSATGAL    (MAXPRNGAL-MINPRNGAL+1) /* number of Galileo satellites */
#define NSYSGAL     1
#else
#define MINPRNGAL   0
#define MAXPRNGAL   0
#define NSATGAL     0
#define NSYSGAL     0
#endif
#ifdef ENAQZS
#define MINPRNQZS   193                 /* min satellite PRN number of QZSS */
#define MAXPRNQZS   202                 /* max satellite PRN number of QZSS */
#define MINPRNQZS_S 183                 /* min satellite PRN number of QZSS SAIF */
#define MAXPRNQZS_S 191                 /* max satellite PRN number of QZSS SAIF */
#define NSATQZS     (MAXPRNQZS-MINPRNQZS+1) /* number of QZSS satellites */
#define NSYSQZS     1
#else
#define MINPRNQZS   0
#define MAXPRNQZS   0
#define MINPRNQZS_S 0
#define MAXPRNQZS_S 0
#define NSATQZS     0
#define NSYSQZS     0
#endif
#ifdef ENACMP
#define MINPRNCMP   1                   /* min satellite sat number of BeiDou */
#define MAXPRNCMP   50                  /* max satellite sat number of BeiDou */
#define NSATCMP     (MAXPRNCMP-MINPRNCMP+1) /* number of BeiDou satellites */
#define NSYSCMP     1
#else
#define MINPRNCMP   0
#define MAXPRNCMP   0
#define NSATCMP     0
#define NSYSCMP     0
#endif
#ifdef ENAIRN
#define MINPRNIRN   1                   /* min satellite sat number of IRNSS */
#define MAXPRNIRN   7                   /* max satellite sat number of IRNSS */
#define NSATIRN     (MAXPRNIRN-MINPRNIRN+1) /* number of IRNSS satellites */
#define NSYSIRN     1
#else
#define MINPRNIRN   0
#define MAXPRNIRN   0
#define NSATIRN     0
#define NSYSIRN     0
#endif
#ifdef ENALEO
#define MINPRNLEO   1                   /* min satellite sat number of LEO */
#define MAXPRNLEO   10                  /* max satellite sat number of LEO */
#define NSATLEO     (MAXPRNLEO-MINPRNLEO+1) /* number of LEO satellites */
#define NSYSLEO     1
#else
#define MINPRNLEO   0
#define MAXPRNLEO   0
#define NSATLEO     0
#define NSYSLEO     0
#endif
#define NSYS        (NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS+NSYSCMP+NSYSIRN+NSYSLEO) /* number of systems */

#ifdef ENASBS
#define MINPRNSBS   120                 /* min satellite PRN number of SBAS */
#define MAXPRNSBS   142                 /* max satellite PRN number of SBAS */
#define NSATSBS     (MAXPRNSBS-MINPRNSBS+1) /* number of SBAS satellites */
#else
#define MINPRNSBS   0                 /* min satellite PRN number of SBAS */
#define MAXPRNSBS   0                 /* max satellite PRN number of SBAS */
#define NSATSBS     0                 /* number of SBAS satellites */
#endif

#define MAXSAT      (NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+NSATIRN+NSATSBS+NSATLEO)
                                        /* max satellite number (1 to MAXSAT) */
#define MAXSTA      255

#ifndef MAXOBS
#define MAXOBS      64                  /* max number of obs in an epoch */
#endif
#define MAXRCV      64                  /* max receiver number (1 to MAXRCV) */
#define MAXOBSTYPE  64                  /* max number of obs type in RINEX */
#ifdef OBS_100HZ
#define DTTOL       0.005               /* tolerance of time difference (s) */
#else
#define DTTOL       0.025               /* tolerance of time difference (s) */
#endif
#define INS_DTTOL   0.0025
#define MAXDTOE     7200.0              /* max time difference to GPS Toe (s) */
#define MAXDTOE_QZS 7200.0              /* max time difference to QZSS Toe (s) */
#define MAXDTOE_GAL 14400.0             /* max time difference to Galileo Toe (s) */
#define MAXDTOE_CMP 21600.0             /* max time difference to BeiDou Toe (s) */
#define MAXDTOE_GLO 1800.0              /* max time difference to GLONASS Toe (s) */
#define MAXDTOE_IRN 7200.0              /* max time difference to IRNSS Toe (s) */
#define MAXDTOE_SBS 360.0               /* max time difference to SBAS Toe (s) */
#define MAXDTOE_S   86400.0             /* max time difference to ephem toe (s) for other */
#define MAXGDOP     300.0               /* max GDOP */

#define INT_SWAP_TRAC 86400.0           /* swap interval of trace file (s) */
#define INT_SWAP_STAT 86400.0           /* swap interval of solution status file (s) */

#define MAXEXFILE   1024                /* max number of expanded files */
#define MAXSBSAGEF  30.0                /* max age of SBAS fast correction (s) */
#define MAXSBSAGEL  1800.0              /* max age of SBAS long term corr (s) */
#define MAXSBSURA   8                   /* max URA of SBAS satellite */
#define MAXBAND     10                  /* max SBAS band of IGP */
#define MAXNIGP     201                 /* max number of IGP in SBAS band */
#define MAXNGEO     4                   /* max number of GEO satellites */
#define MAXCOMMENT  10                  /* max number of RINEX comments */
#define MAXSTRPATH  300                /* max length of stream path */
#define MAXCHARS    300
#define MAXSTRMSG   1024                /* max length of stream message */
#define MAXSTRRTK   8                   /* max number of stream in RTK server */
#define MAXSBSMSG   32                  /* max number of SBAS msg in RTK server */
#define MAXSOLMSG   8191                /* max length of solution message */
#define MAXRAWLEN   4096                /* max length of receiver raw message */
#define MAXERRMSG   4096                /* max length of error/warning message */
#define MAXANT      64                  /* max length of station name/antenna type */
#define MAXSOLBUF   256                 /* max number of solution buffer */
#define MAXOBSBUF   128                 /* max number of observation data buffer */
#define MAXNRPOS    16                  /* max number of reference positions */
#define MAXLEAPS    64                  /* max number of leap seconds table */
#define MAXGISLAYER 32                  /* max number of GIS data layers */
#define MAXRCVCMD   4096                /* max length of receiver commands */
#define MAXFILE     12

#define RNX2VER     2.10                /* RINEX ver.2 default output version */
#define RNX3VER     3.00                /* RINEX ver.3 default output version */

#define OBSTYPE_PR  0x01                /* observation type: pseudorange */
#define OBSTYPE_CP  0x02                /* observation type: carrier-phase */
#define OBSTYPE_DOP 0x04                /* observation type: doppler-freq */
#define OBSTYPE_SNR 0x08                /* observation type: SNR */
#define OBSTYPE_ALL 0xFF                /* observation type: all */

#define FREQTYPE_L1 0x01                /* frequency type: L1/E1 */
#define FREQTYPE_L2 0x02                /* frequency type: L2/B1/E5b */
#define FREQTYPE_L5 0x04                /* frequency type: L5/E5a/L3 */
#define FREQTYPE_E6 0x08                /* frequency type: E6/LEX/B3 */
#define FREQTYPE_E5ab 0x10              /* frequency type: E5(a+b) */
#define FREQTYPE_S 0x20                 /* frequency type: S */
#define FREQTYPE_ALL 0xFF               /* frequency type: all */

#define CODE_NONE   0                   /* obs code: none or unknown */
#define CODE_L1C    1                   /* obs code: L1C/A,G1C/A,E1C (GPS,GLO,GAL,QZS,SBS) */
#define CODE_L1P    2                   /* obs code: L1P,G1P,B1P (GPS,GLO,BDS) */
#define CODE_L1W    3                   /* obs code: L1 Z-track (GPS) */
#define CODE_L1Y    4                   /* obs code: L1Y        (GPS) */
#define CODE_L1M    5                   /* obs code: L1M        (GPS) */
#define CODE_L1N    6                   /* obs code: L1codeless,B1codeless (GPS,BDS) */
#define CODE_L1S    7                   /* obs code: L1C(D)     (GPS,QZS) */
#define CODE_L1L    8                   /* obs code: L1C(P)     (GPS,QZS) */
#define CODE_L1E    9                   /* (not used) */
#define CODE_L1A    10                  /* obs code: E1A,B1A    (GAL,BDS) */
#define CODE_L1B    11                  /* obs code: E1B        (GAL) */
#define CODE_L1X    12                  /* obs code: E1B+C,L1C(D+P),B1D+P (GAL,QZS,BDS) */
#define CODE_L1Z    13                  /* obs code: E1A+B+C,L1S (GAL,QZS) */
#define CODE_L2C    14                  /* obs code: L2C/A,G1C/A (GPS,GLO) */
#define CODE_L2D    15                  /* obs code: L2 L1C/A-(P2-P1) (GPS) */
#define CODE_L2S    16                  /* obs code: L2C(M)     (GPS,QZS) */
#define CODE_L2L    17                  /* obs code: L2C(L)     (GPS,QZS) */
#define CODE_L2X    18                  /* obs code: L2C(M+L),B1_2I+Q (GPS,QZS,BDS) */
#define CODE_L2P    19                  /* obs code: L2P,G2P    (GPS,GLO) */
#define CODE_L2W    20                  /* obs code: L2 Z-track (GPS) */
#define CODE_L2Y    21                  /* obs code: L2Y        (GPS) */
#define CODE_L2M    22                  /* obs code: L2M        (GPS) */
#define CODE_L2N    23                  /* obs code: L2codeless (GPS) */
#define CODE_L5I    24                  /* obs code: L5I,E5aI   (GPS,GAL,QZS,SBS) */
#define CODE_L5Q    25                  /* obs code: L5Q,E5aQ   (GPS,GAL,QZS,SBS) */
#define CODE_L5X    26                  /* obs code: L5I+Q,E5aI+Q,L5B+C,B2aD+P (GPS,GAL,QZS,IRN,SBS,BDS) */
#define CODE_L7I    27                  /* obs code: E5bI,B2bI  (GAL,BDS) */
#define CODE_L7Q    28                  /* obs code: E5bQ,B2bQ  (GAL,BDS) */
#define CODE_L7X    29                  /* obs code: E5bI+Q,B2bI+Q (GAL,BDS) */
#define CODE_L6A    30                  /* obs code: E6A,B3A    (GAL,BDS) */
#define CODE_L6B    31                  /* obs code: E6B        (GAL) */
#define CODE_L6C    32                  /* obs code: E6C        (GAL) */
#define CODE_L6X    33                  /* obs code: E6B+C,LEXS+L,B3I+Q (GAL,QZS,BDS) */
#define CODE_L6Z    34                  /* obs code: E6A+B+C,L6D+E (GAL,QZS) */
#define CODE_L6S    35                  /* obs code: L6S        (QZS) */
#define CODE_L6L    36                  /* obs code: L6L        (QZS) */
#define CODE_L8I    37                  /* obs code: E5abI      (GAL) */
#define CODE_L8Q    38                  /* obs code: E5abQ      (GAL) */
#define CODE_L8X    39                  /* obs code: E5abI+Q,B2abD+P (GAL,BDS) */
#define CODE_L2I    40                  /* obs code: B1_2I      (BDS) */
#define CODE_L2Q    41                  /* obs code: B1_2Q      (BDS) */
#define CODE_L6I    42                  /* obs code: B3I        (BDS) */
#define CODE_L6Q    43                  /* obs code: B3Q        (BDS) */
#define CODE_L3I    44                  /* obs code: G3I        (GLO) */
#define CODE_L3Q    45                  /* obs code: G3Q        (GLO) */
#define CODE_L3X    46                  /* obs code: G3I+Q      (GLO) */
#define CODE_L1I    47                  /* obs code: B1I        (BDS) (obsolute) */
#define CODE_L1Q    48                  /* obs code: B1Q        (BDS) (obsolute) */
#define CODE_L5A    49                  /* obs code: L5A SPS    (IRN) */
#define CODE_L5B    50                  /* obs code: L5B RS(D)  (IRN) */
#define CODE_L5C    51                  /* obs code: L5C RS(P)  (IRN) */
#define CODE_L9A    52                  /* obs code: SA SPS     (IRN) */
#define CODE_L9B    53                  /* obs code: SB RS(D)   (IRN) */
#define CODE_L9C    54                  /* obs code: SC RS(P)   (IRN) */
#define CODE_L9X    55                  /* obs code: SB+C       (IRN) */
#define CODE_L1D    56                  /* obs code: B1D        (BDS) */
#define CODE_L5D    57                  /* obs code: L5D(L5S),B2aD (QZS,BDS) */
#define CODE_L5P    58                  /* obs code: L5P(L5S),B2aP (QZS,BDS) */
#define CODE_L5Z    59                  /* obs code: L5D+P(L5S) (QZS) */
#define CODE_L6E    60                  /* obs code: L6E        (QZS) */
#define CODE_L7D    61                  /* obs code: B2bD       (BDS) */
#define CODE_L7P    62                  /* obs code: B2bP       (BDS) */
#define CODE_L7Z    63                  /* obs code: B2bD+P     (BDS) */
#define CODE_L8D    64                  /* obs code: B2abD      (BDS) */
#define CODE_L8P    65                  /* obs code: B2abP      (BDS) */
#define CODE_L4A    66                  /* obs code: G1aL1OCd   (GLO) */
#define CODE_L4B    67                  /* obs code: G1aL1OCd   (GLO) */
#define CODE_L4X    68                  /* obs code: G1al1OCd+p (GLO) */
#define MAXCODE     68                  /* max number of obs code */

#define G1W2W       0
#define G1C1W       1
#define G2C2W       2
#define G1C5Q       3
#define G1C2W       4
#define G1C5X       5
#define G2W2S       6
#define G2W2L       7
#define G2W2X       8
#define R1P2P       0
#define R1C1P       1
#define R2C2P       2
#define R1C2P       3
#define R1C2C       4
#define E1C5Q       0
#define E1C6C       1
#define E1C7Q       2
#define E1C8Q       3
#define E1X5X       4
#define E1X7X       5
#define E1X8X       6
#define B2I7I       0
#define B2I6I       1
#define C1X5X       2
#define C1P5P       3
#define C1D5D       4
#define C1X6I       5
#define C1P6I       6
#define C1D6I       7
#define C2I6I       8
#define C1X7Z       9
#define C1X8X       10
#define J1C2L       0
#define J1C5X       1
#define J1C5Q       2
#define J1X2X       3
#define J1X5X       4
#define J1C1X       5

#define PMODE_SINGLE       0              /* positioning mode: single */
#define PMODE_TDCP         1              /* positioning mode: tdcp */
#define PMODE_DGPS         2              /* positioning mode: DGPS/DGNSS */
#define PMODE_KINEMA       3              /* positioning mode: kinematic */
#define PMODE_STATIC       4              /* positioning mode: static */
#define PMODE_STATIC_START 5              /* positioning mode: static */
#define PMODE_MOVEB        6              /* positioning mode: moving-base */
#define PMODE_FIXED        7              /* positioning mode: fixed */
#define PMODE_PPP_KINEMA   8              /* positioning mode: PPP-kinemaric */
#define PMODE_PPP_STATIC   9              /* positioning mode: PPP-static */
#define PMODE_PPP_FIXED    10             /* positioning mode: PPP-fixed */
#define PMODE_INS_MECH     11             /* positioning mode: INS mech*/
#define PMODE_LC_POS       12             /* positioning mode: gnss/ins loosely coupled with gnss solution(.pos) file*/
#define PMODE_LC_SPP       13             /* positioning mode: gnss/ins loosely coupled using spp */
#define PMODE_LC_DGPS      14             /* positioning mode: gnss/ins loosely coupled using dgps*/
#define PMODE_LC_PPK       15             /* positioning mode: gnss/ins loosely coupled using ppk */
#define PMODE_LC_PPP       16             /* positioning mode: gnss/ins loosely coupled using ppp */
#define PMODE_TC_SPP       17             /* positioning mode: gnss/ins tightly coupled using spp */
#define PMODE_TC_TDCP      18             /* positioning mode: gnss/ins tightly coupled using tdcp */
#define PMODE_TC_DGPS      19             /* positioning mode: gnss/ins tightly coupled using dgps*/
#define PMODE_TC_PPK       20             /* positioning mode: gnss/ins tightly coupled using ppk */
#define PMODE_TC_PPP       21             /* positioning mode: gnss/ins tightly coupled using ppp */
#define PMODE_STC_PPK      22
#define PMODE_STC_PPP      23

#define AR_PROD_IRC       1
#define AR_PROD_FCB       2
#define AR_PROD_UPD       3
#define AR_PROD_OSB_GRM   4
#define AR_PROD_OSB_WHU   5
#define AR_PROD_OSB_COM   6
#define AR_PROD_OSB_SGG   7
#define AR_PROD_OSB_CNT   8

#define SOLT_XYZ    0                   /* solution type: x/y/z-ecef */
#define SOLT_ENU    1                   /* solution type: e/n/u-baseline */

#define SOLF_LLH     0                   /* solution format: lat/lon/height */
#define SOLF_XYZ     1                   /* solution format: x/y/z-ecef */
#define SOLF_ENU     2                   /* solution format: e/n/u-baseline */
#define SOLF_NMEA    3                   /* solution format: NMEA-183 */
#define SOLF_STAT    4                   /* solution format: solution status */
#define SOLF_GSIF    5                   /* solution format: GSI F1/F2 */
#define SOLF_INS_LLH 6
#define SOLF_INS_XYZ 7
#define SOLF_INS_YGM 8

#define SOLQ_NONE       0                   /* solution status: no solution */
#define SOLQ_FIX        1                   /* solution status: fix */
#define SOLQ_FLOAT      2                   /* solution status: float */
#define SOLQ_SBAS       3                   /* solution status: SBAS */
#define SOLQ_DGPS       4                   /* solution status: DGPS/DGNSS */
#define SOLQ_SINGLE     5                   /* solution status: single */
#define SOLQ_PPP        6                   /* solution status: PPP */
#define SOLQ_PPP_AR     7
#define SOLQ_INS_ALIGN  8
#define SOLQ_INS        9
#define SOLQ_IGLC      10
#define SOLQ_IGTC      11
#define SOLQ_STC_PR    12
#define SOLQ_STC_CP    13
#define SOLQ_NHC_AID   14
#define SOLQ_ZV_AID    15
#define SOLQ_DOP_AID   16
#define SOLQ_VEL_AID   17
#define MAXSOLQ        17                   /* max number of solution status */

#define GNSS_STATUS_NONE    0
#define GNSS_STATUS_POS     1
#define GNSS_STATUS_SPP     2
#define GNSS_STATUS_DGPS    3
#define GNSS_STATUS_PPK     4
#define GNSS_STATUS_PPK_AR  5
#define GNSS_STATUS_PPP     6
#define GNSS_STATUS_PPP_AR  7


#define TIMES_GPST  0                   /* time system: gps time */
#define TIMES_UTC   1                   /* time system: utc */
#define TIMES_JST   2                   /* time system: jst */

#define BD3OPT_OFF   0
#define BD3OPT_BD23  1
#define BD3OPT_BD2_3 2
#define BD3OPT_BD3   3

#define CBIAS_OPT_OFF      0
#define CBIAS_OPT_BRD_TGD  1
#define CBIAS_OPT_COD_DCB  2
#define CBIAS_OPT_IGG_DCB  3
#define CBIAS_OPT_GBM_DCB  4
#define CBIAS_OPT_MIX_DCB  5
#define CBIAS_OPT_IGG_BIA  6
#define CBIAS_OPT_COD_BIA  7
#define CBIAS_OPT_GRG_BIA  8

#define IONOOPT_OFF     0                   /* ionosphere option: correction off */
#define IONOOPT_BRDC    1                  /* ionosphere option: broadcast model */
#define IONOOPT_SBAS    2                  /* ionosphere option: SBAS model */
#define IONOOPT_TEC     3                  /* ionosphere option: IONEX TEC model */
#define IONOOPT_IFLC    4                  /* ionosphere option: L1/L2 or L1/L5 iono-free LC */
#define IONOOPT_IF2     5                  /* ionosphere option: two dual-frequency IF */
#define IONOOPT_UC      6                  /* ionosphere option: estimation */
#define IONOOPT_UC_CONS 7
#define IONOOPT_QZS     8                  /* ionosphere option: QZSS broadcast model */
#define IONOOPT_LEX     9                  /* ionosphere option: QZSS LEX ionospehre */
#define IONOOPT_STEC    10                  /* ionosphere option: SLANT TEC model */

#define TROPOPT_OFF 0                   /* troposphere option: correction off */
#define TROPOPT_SAAS 1                  /* troposphere option: Saastamoinen model */
#define TROPOPT_SBAS 2                  /* troposphere option: SBAS model */
#define TROPOPT_EST 3                   /* troposphere option: ZTD estimation */
#define TROPOPT_ESTG 4                  /* troposphere option: ZTD+grad estimation */
#define TROPOPT_ZTD 5                   /* troposphere option: ZTD correction */

#define TROPMAP_SAAS 0
#define TROPMAP_NMF  1
#define TROPMAP_GMF  2
#define TROPMAP_VMF  3

#define EPHOPT_BRDC 0                   /* ephemeris option: broadcast ephemeris */
#define EPHOPT_PREC 1                   /* ephemeris option: precise ephemeris */
#define EPHOPT_SBAS 2                   /* ephemeris option: broadcast + SBAS */
#define EPHOPT_SSRAPC 3                 /* ephemeris option: broadcast + SSR_APC */
#define EPHOPT_SSRCOM 4                 /* ephemeris option: broadcast + SSR_COM */
#define EPHOPT_LEX  5                   /* ephemeris option: QZSS LEX ephemeris */

#define VELOPT_OFF     0
#define VELOPT_DOPPLER 1
#define VELOPT_TDCP    2

#define WEIGHTOPT_ELEVATION 0           /* weighting option: elevation */
#define WEIGHTOPT_SNR       1           /* weighting option: snr */  

#define ARMODE_OFF  0                   /* AR mode: off */
#define ARMODE_CONT 1                   /* AR mode: continuous */
#define ARMODE_INST 2                   /* AR mode: instantaneous */
#define ARMODE_FIXHOLD 3                /* AR mode: fix and hold */
#define ARMODE_WLNL 4                   /* AR mode: wide lane/narrow lane */
#define ARMODE_TCAR 5                   /* AR mode: triple carrier ar */
#define ARMODE_PPPAR 6
#define ARMODE_PPPAR_ILS 7

#define ROBUST_QC_OFF     0
#define ROBUST_QC_IGG_PR  1
#define ROBUST_QC_IGG_CP  2
#define ROBUST_QC_IGG     3
#define ROBUST_QC_LIU     4
#define ROBUST_QC_SHI     5
#define ROBUST_QC_ZHAO    6
#define ROBUST_QC_ZHAO_X  7
#define ROBUST_QC_ZHAO_XX 8

#define INSAID_OFF    0
#define INSAID_LT     1
#define INSAID_HHZ    2

#define KFOPT_OFF         0
#define KFOPT_ADA_INNO    1
#define KFOPT_VBKF        2
#define KFOPT_SAGE_HUSA   3

#define GLO_ARMODE_OFF     0            /* GLO AR mode: off */
#define GLO_ARMODE_ON      1            /* GLO AR mode: on */
#define GLO_ARMODE_AUTOCAL 2            /* GLO AR mode: autocal */
#define GLO_ARMODE_FIXHOLD 3            /* GLO AR mode: fix and hold */

#define SBSOPT_LCORR 1                  /* SBAS option: long term correction */
#define SBSOPT_FCORR 2                  /* SBAS option: fast correction */
#define SBSOPT_ICORR 4                  /* SBAS option: ionosphere correction */
#define SBSOPT_RANGE 8                  /* SBAS option: ranging */

#define POSOPT_POS   0                  /* pos option: LLH/XYZ */
#define POSOPT_SINGLE 1                 /* pos option: average of single pos */
#define POSOPT_FILE  2                  /* pos option: read from pos file */
#define POSOPT_RINEX 3                  /* pos option: rinex header pos */
#define POSOPT_RTCM  4                  /* pos option: rtcm station pos */
#define POSOPT_RAW   5                  /* pos option: raw station pos */

#define STR_NONE     0                  /* stream type: none */
#define STR_SERIAL   1                  /* stream type: serial */
#define STR_FILE     2                  /* stream type: file */
#define STR_TCPSVR   3                  /* stream type: TCP server */
#define STR_TCPCLI   4                  /* stream type: TCP client */
#define STR_NTRIPSVR 6                  /* stream type: NTRIP server */
#define STR_NTRIPCLI 7                  /* stream type: NTRIP client */
#define STR_FTP      8                  /* stream type: ftp */
#define STR_HTTP     9                  /* stream type: http */
#define STR_NTRIPC_S 10                 /* stream type: NTRIP caster server */
#define STR_NTRIPC_C 11                 /* stream type: NTRIP caster client */
#define STR_UDPSVR   12                 /* stream type: UDP server */
#define STR_UDPCLI   13                 /* stream type: UDP server */
#define STR_MEMBUF   14                 /* stream type: memory buffer */

#define STRFMT_RTCM2 0                  /* stream format: RTCM 2 */
#define STRFMT_RTCM3 1                  /* stream format: RTCM 3 */
#define STRFMT_OEM4  2                  /* stream format: NovAtel OEMV/4 */
#define STRFMT_CNAV  3                  /* stream format: ComNav */
#define STRFMT_UBX   4                  /* stream format: u-blox LEA-*T */
#define STRFMT_SBP   5                  /* stream format: Swift Navigation SBP */
#define STRFMT_CRES  6                  /* stream format: Hemisphere */
#define STRFMT_STQ   7                  /* stream format: SkyTraq S1315F */
#define STRFMT_GW10  8                  /* stream format: Furuno GW10 */
#define STRFMT_JAVAD 9                  /* stream format: JAVAD GRIL/GREIS */
#define STRFMT_NVS   10                 /* stream format: NVS NVC08C */
#define STRFMT_BINEX 11                 /* stream format: BINEX */
#define STRFMT_RT17  12                 /* stream format: Trimble RT17 */
#define STRFMT_SEPT  13                 /* stream format: Septentrio */
#define STRFMT_CMR   14                 /* stream format: CMR/CMR+ */
#define STRFMT_TERSUS 15                /* stream format: TERSUS */
#define STRFMT_LEXR  16                 /* stream format: Furuno LPY-10000 */
#define STRFMT_RINEX 17                 /* stream format: RINEX */
#define STRFMT_SP3   18                 /* stream format: SP3 */
#define STRFMT_RNXCLK 19                /* stream format: RINEX CLK */
#define STRFMT_SBAS    20               /* stream format: SBAS messages */
#define STRFMT_NMEA    21               /* stream format: NMEA 0183 */
#define STRFMT_IMU_YGM_SIM   22
#define STRFMT_IMU_NVT_CPT   23
#define STRFMT_IMU_M39       24
#define STRFMT_IMU_NVT_KVH   25
#define STRFMT_IMU_NVT_IGM   26
#define STRFMT_IMU_CPT_TOKYO 27
#define STRFMT_IMU_A15       28
#define STRFMT_IMU_BY        29
#define STRFMT_IMU_ATLAN     30
#define STRFMT_IMU_BOSCH     31
#define STRFMT_IMU_LORD      32
#define STRFMT_IMU_PROPERTY  33
#define STRFMT_YGM_PV        34
#define STRFMT_POS           35
#define STRFMT_OD_YGM_SIM    36
#ifndef EXTLEX
#define MAXRCVFMT    15                 /* max number of receiver format */
#else
#define MAXRCVFMT    16
#endif

#define STR_MODE_R  0x1                 /* stream mode: read */
#define STR_MODE_W  0x2                 /* stream mode: write */
#define STR_MODE_RW 0x3                 /* stream mode: read/write */

#define REF_SOL_FMT_IE        0
#define REF_SOL_FMT_A15       1
#define REF_SOL_FMT_POS       2
#define REF_SOL_FMT_TEX       3
#define REF_SOL_FMT_YGM_AVP   4
#define REF_SOL_FMT_YGM_AVPED 5
#define REF_SOL_FMT_SINEX     6


#define FT_IDX_ROVER   0
#define FT_IDX_BASE    1
#define FT_IDX_IMU     2
#define FT_IDX_AVP     3
#define FT_IDX_OD      4
#define FT_IDX_CAR     5
#define FT_IDX_IN1     6
#define FT_IDX_OUT1    7
#define FT_IDX_OUT2    8

#define GEOID_EMBEDDED    0             /* geoid model: embedded geoid */
#define GEOID_EGM96_M150  1             /* geoid model: EGM96 15x15" */
#define GEOID_EGM2008_M25 2             /* geoid model: EGM2008 2.5x2.5" */
#define GEOID_EGM2008_M10 3             /* geoid model: EGM2008 1.0x1.0" */
#define GEOID_GSI2000_M15 4             /* geoid model: GSI geoid 2000 1.0x1.5" */
#define GEOID_RAF09       5             /* geoid model: IGN RAF09 for France 1.5"x2" */

#define COMMENTH    "%"                 /* comment line indicator for solution */
#define MSG_DISCONN "$_DISCONNECT\r\n"  /* disconnect message */

#define DLOPT_FORCE   0x01              /* download option: force download existing */
#define DLOPT_KEEPCMP 0x02              /* download option: keep compressed file */
#define DLOPT_HOLDERR 0x04              /* download option: hold on error file */
#define DLOPT_HOLDLST 0x08              /* download option: hold on listing file */

#define LLI_SLIP    0x01                /* LLI: cycle-slip */
#define LLI_HALFC   0x02                /* LLI: half-cycle not resovled */
#define LLI_BOCTRK  0x04                /* LLI: boc tracking of mboc signal */
#define LLI_HALFA   0x40                /* LLI: half-cycle added */
#define LLI_HALFS   0x80                /* LLI: half-cycle subtracted */

#define IMUFMT_KVH  1                   /* imu data format KVH */
#define IMUFMT_GI310  2                 /* imu data format GI310 */
#define IMUFMT_UBX    3                 /* imu data format ublox: EVK-M8U-0-00 */

#define IMUCOOR_FRD        1            /* imu body coordinate frame: frd */
#define IMUCOOR_RFU        2            /* imu body coordinate frame: rfu */

#define IMUVALFMT_DEG      1            /* imu gyro measurement data value format: degree */
#define IMUVALFMT_RAD      2            /* imu gyro measurement data value format: rad */

#define IMUDECFMT_RATE     1            /* imu measurements data format: angular rate/acceleration */
#define IMUDECFMT_INCR     2            /* imu measurements data format: angular/velocity increment */

#define INSMECH_LLH        1
#define INSMECH_ECEF       2

#define INSLOCAL_ENU       1
#define INSLOCAL_NED       2

#define LCPOSFMT_NONE      0
#define LCPOSFMT_POS       1
#define LCPOSFMT_YGM_PV    2

#define IMUDETST_GLRT      1            /* runs the generalized likelihood test for detect static imu measurement */
#define IMUDETST_MV        2            /* runs the acceleration moving variance detector */
#define IMUDETST_MAG       3            /* runs the acceleration magnitude detector */
#define IMUDETST_ARE       4            /* runs the angular rate energy detector */
#define IMUDETST_ALL       5            /* runs all static detector */

#define P2_5        0.03125             /* 2^-5 */
#define P2_6        0.015625            /* 2^-6 */
#define P2_11       4.882812500000000E-04 /* 2^-11 */
#define P2_15       3.051757812500000E-05 /* 2^-15 */
#define P2_17       7.629394531250000E-06 /* 2^-17 */
#define P2_19       1.907348632812500E-06 /* 2^-19 */
#define P2_20       9.536743164062500E-07 /* 2^-20 */
#define P2_21       4.768371582031250E-07 /* 2^-21 */
#define P2_23       1.192092895507810E-07 /* 2^-23 */
#define P2_24       5.960464477539063E-08 /* 2^-24 */
#define P2_27       7.450580596923828E-09 /* 2^-27 */
#define P2_29       1.862645149230957E-09 /* 2^-29 */
#define P2_30       9.313225746154785E-10 /* 2^-30 */
#define P2_31       4.656612873077393E-10 /* 2^-31 */
#define P2_32       2.328306436538696E-10 /* 2^-32 */
#define P2_33       1.164153218269348E-10 /* 2^-33 */
#define P2_35       2.910383045673370E-11 /* 2^-35 */
#define P2_38       3.637978807091710E-12 /* 2^-38 */
#define P2_39       1.818989403545856E-12 /* 2^-39 */
#define P2_40       9.094947017729280E-13 /* 2^-40 */
#define P2_43       1.136868377216160E-13 /* 2^-43 */
#define P2_48       3.552713678800501E-15 /* 2^-48 */
#define P2_50       8.881784197001252E-16 /* 2^-50 */
#define P2_55       2.775557561562891E-17 /* 2^-55 */

#ifdef WIN32
#define thread_t    HANDLE
#define lock_t      CRITICAL_SECTION
#define initlock(f) InitializeCriticalSection(f)
#define lock(f)     EnterCriticalSection(f)
#define unlock(f)   LeaveCriticalSection(f)
#define FILEPATHSEP '\\'
#else
#define thread_t    pthread_t
#define lock_t      pthread_mutex_t
#define initlock(f) pthread_mutex_init(f,NULL)
#define lock(f)     pthread_mutex_lock(f)
#define unlock(f)   pthread_mutex_unlock(f)
#define FILEPATHSEP '/'
#endif

/* type definitions ----------------------------------------------------------*/

#ifdef WIN32
typedef struct {        /* time struct */
    time_t time;        /* time (s) expressed by standard time_t */
    double sec;         /* fraction of second under 1 s */
} gtime_t;
#else
typedef struct {        /* time struct */
    time_t time;        /* time (s) expressed by standard time_t */
    __attribute__ ((aligned (8)))double sec; /* fraction of second under 1 s */
} gtime_t;
#endif /*WIN32*/

typedef struct {
    long sn;
    double tos;
}sod_t;

typedef struct{
    long day;
    sod_t ds;
}mjd_t;

typedef union {
    struct {
        double x,y,z;
    };
    double v[3];
}v3_t;

typedef struct {
    double q0,q1,q2,q3;
}quat_t;

typedef union {
    struct{
        double m11,m12,m13,m21,m22,m23,m31,m32,m33;
    };
    double v[9];
}m3_t;

/* type definition -----------------------------------------------------------*/
typedef struct {                        /* signal index type */
    int n;                              /* number of index */
    int idx[MAXOBSTYPE];                /* signal freq-index */
    int frq[MAXOBSTYPE];                /* signal frequency (1:L1,2:L2,...) */
    int pos[MAXOBSTYPE];                /* signal index in obs data (-1:no) */
    uint8_t pri [MAXOBSTYPE];     /* signal priority (15-0) */
    uint8_t type[MAXOBSTYPE];     /* type (0:C,1:L,2:D,3:S) */
    uint8_t code[MAXOBSTYPE];     /* obs code (CODE_L??) */
    double shift[MAXOBSTYPE];           /* phase shift (cycle) */
} sigind_t;

typedef struct {        /* observation data record */
    gtime_t time;       /* receiver sampling time (GPST) */
    gtime_t eventime;   /* time of event (GPST) */
    int timevalid;      /* time is valid (Valid GNSS fix) for time mark */
    uint8_t sat,rcv; /* satellite/receiver number */
    char id[5];
    uint16_t SNR [NFREQ+NEXOBS]; /* signal strength (0.25 dBHz) */
    uint8_t LLI [NFREQ+NEXOBS]; /* loss of lock indicator */
    uint8_t code[NFREQ+NEXOBS]; /* code indicator (CODE_???) */
    uint8_t qualL[NFREQ+NEXOBS]; /* quality of carrier phase measurement */
    uint8_t qualP[NFREQ+NEXOBS]; /* quality of pseudorange measurement */
    uint8_t freq; /* GLONASS frequency channel (0-13) */
    double L[NFREQ+NEXOBS]; /* observation data carrier-phase (cycle) */
    double P[NFREQ+NEXOBS]; /* observation data pseudorange (m) */
    float  D[NFREQ+NEXOBS]; /* observation data doppler frequency (Hz) */
} obsd_t;

typedef struct {        /* observation data */
    int n,nmax;         /* number of obervation data/allocated */
    int flag;           /* epoch flag (0:ok,1:power failure,>1:event flag) */
    int rcvcount;       /* count of rcv event */
    int tmcount;        /* time mark count */
    obsd_t *data;       /* observation data records */

    sigind_t sind[2][8];
} obs_t;

typedef struct {
    gtime_t time;   /**< current time */
    v3_t gyro;      /**< gyro output, angular increment */
    v3_t accel;     /**< accelermeter output, velocity increment */

    unsigned int pps;
    unsigned int imuc;

    short odoc;
} imud_t;

typedef struct {
    int strfmt;
    unsigned int imudecfmt;
    unsigned int imucoors;
    unsigned int imuvalfmt;
    unsigned int freq_imu;  /**< IMU sample rate[Hz] */
    unsigned int freq_od;   /**< Odometer sample rate[Hz] */
    int week;
    double sow;
    gtime_t tstart;     /**< first epoch */
    v3_t gyro_noise;    /**< Gyro output noise [rad/s] */
    v3_t accel_noise;   /**< Accelermeter output notput noise [m/s^2] */
    v3_t gyrnd;         /**< gyroscope noise density [rad/sqrt(s)==rad/s/sqrt(Hz)] */
    v3_t gbrw;          /**< gyroscope (bias) random walk [rad/s/sqrt(s)==rad/s^2/sqrt(Hz)] */
    v3_t accnd;         /**< accelerometer noise density [m/s/sqrt(s)==m/s^2/sqrt(Hz)] */
    v3_t abrw;          /**< accelerometer (bias) random walk[m/s^2/sqrt(s)==m/s^3/sqrt(Hz)] */
    v3_t Ta;            /**< Accel bias correlation time(1st order Markov) [s] */
    v3_t Tg;            /**< Gryo bias correlation time(1st order Markov) [s] */
    gtime_t init_tag;
    v3_t inita;         /**< initial attitude, Enb [rad] */
    v3_t inita_err;
    v3_t initv;         /**< initial velocity, veb_e [m/s] */
    v3_t initv_err;
    v3_t initr;         /**< initial position, reb_e [m] */
    v3_t initr_err;
    v3_t ba;            /**< (initial) accel bias [m/s^2] */
    v3_t ba_err;        /**< (initial) accel bias stanadard error[m/s^2] */
    v3_t bg;            /**< (initial) gryo bias [rad/s] */
    v3_t bg_err;        /**< (initial) gryo bias stanadard error[m/s^2] */
    v3_t sa;            /**< Accelermeter scalar factor(or initial value) */
    v3_t sa_err;        /**< Standard error of accelermeter scalar factor(initial value) */
    v3_t sg;            /**< Gyro scalar factor(or initial value) */
    v3_t sg_err;        /**< Standard error of Gyro scalar factor(initial value) */
    double kod;         /**< odometer scalar factor, true/output */
    double kod_err;         /**< initial odometer scalar factor uncertainty [^2] */
    v3_t lever_arm_gps;     /**< gnss phase center position under imu frame[m] */
    v3_t lever_arm_gps_std; /**< gnss lever arm uncertainty [m] */
    v3_t lever_arm_od;      /**< odometer reference center under imu frame[m]*/
    v3_t lever_arm_od_std;  /**< odometer lever arm uncertainty [m] */
    v3_t lever_arm_car;     /**< car tailing wheel center position under imu frame[m] */
    v3_t lever_arm_car_std; /**< car tailing whell center uncertainty [m] */
    v3_t err_angle_imu;     /**< IMU install error angle(car-imu, roll, pitch, yaw)[rad] */
    v3_t err_angle_imu_std; /**< IMU install error angle uncertainty[rad] */
    v3_t err_angle_imu_rw;  /**< IMU install error angle randon walk [rad/sqrt(s)] */
    v3_t Terr_angle_imu;    /**< IMU install error angle  correlation time(1st order Markov)[s] */
    v3_t err_angle_gps;     /**< GPS install error angle(gps-imu, roll, pitch, yaw)[rad] */
    v3_t err_angle_gps_std; /**< GPS install error angle uncertainty[rad] */
    v3_t ref_point;         /**< reference point under b-frame, use for solution output [m]*/
    unsigned char gyro_axis[3];
    unsigned char accel_axis[3];
} imup_t;

typedef struct {
    unsigned int n, nmax;   /**< number of data/allocated */
    imud_t* data;           /**< IMU observation data record */
    imup_t *property;       /**< IMU property */
} imu_t;

typedef struct {        /* earth rotation parameter data type */
    double mjd;         /* mjd (days) */
    double xp,yp;       /* pole offset (rad) */
    double xpr,ypr;     /* pole offset rate (rad/day) */
    double ut1_utc;     /* ut1-utc (s) */
    double lod;         /* length of day (s/day) */
} erpd_t;

typedef struct {        /* earth rotation parameter type */
    int n,nmax;         /* number and max number of data */
    erpd_t *data;       /* earth rotation parameter data */
} erp_t;

typedef struct {        /* antenna parameter type */
    int sat;            /* satellite number (0:receiver) */
    char type[MAXANT];  /* antenna type */
    char code[MAXANT];  /* serial number or satellite code */
    gtime_t ts,te;      /* valid time start and end */
    double off[NSYS*NFREQ][ 3]; /* phase center offset e/n/u or x/y/z (m) */
    double var[NSYS*NFREQ][80*50]; /* phase center variation (m) */
                        /* el=90,85,...,0 or nadir=0,1,2,3,... (deg) */
    double dazi;
    double zen1,zen2,dzen;
} pcv_t;

typedef struct {        /* antenna parameters type */
    int n,nmax;         /* number of data/allocated */
    pcv_t *pcv;         /* antenna parameters data */
} pcvs_t;

typedef struct {        /* almanac type */
    int sat;            /* satellite number */
    int svh;            /* sv health (0:ok) */
    int svconf;         /* as and sv config */
    int week;           /* GPS/QZS: gps week, GAL: galileo week */
    gtime_t toa;        /* Toa */
                        /* SV orbit parameters */
    double A,e,i0,OMG0,omg,M0,OMGd;
    double toas;        /* Toa (s) in week */
    double f0,f1;       /* SV clock parameters (af0,af1) */
} alm_t;

typedef struct {        /* GPS/QZS/GAL broadcast ephemeris type */
    int sat;            /* satellite number */
    int iode,iodc;      /* IODE,IODC */
    int sva;            /* SV accuracy (URA index) */
    int svh;            /* SV health (0:ok) */
    int week;           /* GPS/QZS: gps week, GAL: galileo week */
    int code;           /* GPS/QZS: code on L2 */
    /* GAL: data source defined as rinex 3.03 */
    /* BDS: data source (0:unknown,1:B1I,2:B1Q,3:B2I,4:B2Q,5:B3I,6:B3Q) */
    int flag;           /* GPS/QZS: L2 P data flag */
    /* BDS: nav type (0:unknown,1:IGSO/MEO,2:GEO) */
    gtime_t toe,toc,ttr; /* Toe,Toc,T_trans */
    /* SV orbit parameters */
    double A,e,i0,OMG0,omg,M0,deln,OMGd,idot;
    double crc,crs,cuc,cus,cic,cis;
    double toes;        /* Toe (s) in week */
    double fit;         /* fit interval (h) */
    double f0,f1,f2;    /* SV clock parameters (af0,af1,af2) */
    double tgd[6];      /* group delay parameters */
                        /* GPS/QZS:tgd[0]=TGD */
                        /* GAL:tgd[0]=BGD_E1E5a,tgd[1]=BGD_E1E5b */
                        /* CMP:tgd[0]=TGD_B1I ,tgd[1]=TGD_B2I/B2b,tgd[2]=TGD_B1Cp */
                        /*     tgd[3]=TGD_B2ap,tgd[4]=ISC_B1Cd   ,tgd[5]=ISC_B2ad */
    double Adot,ndot;   /* Adot,ndot for CNAV */
} eph_t;

typedef struct {        /* GLONASS broadcast ephemeris type */
    int sat;            /* satellite number */
    int iode;           /* IODE (0-6 bit of tb field) */
    int frq;            /* satellite frequency number */
    int svh,sva,age;    /* satellite health, accuracy, age of operation */
    gtime_t toe;        /* epoch of epherides (gpst) */
    gtime_t tof;        /* message frame time (gpst) */
    double pos[3];      /* satellite position (ecef) (m) */
    double vel[3];      /* satellite velocity (ecef) (m/s) */
    double acc[3];      /* satellite acceleration (ecef) (m/s^2) */
    double taun,gamn;   /* SV clock bias (s)/relative freq bias */
    double dtaun;       /* delay between L1 and L2 (s) */
} geph_t;

typedef struct {        /* precise ephemeris type */
    gtime_t time;       /* time (GPST) */
    int index;          /* ephemeris index for multiple files */
    double pos[MAXSAT][4]; /* satellite position/clock (ecef) (m|s) */
    float  std[MAXSAT][4]; /* satellite position/clock std (m|s) */
    double vel[MAXSAT][4]; /* satellite velocity/clk-rate (m/s|s/s) */
    float  vst[MAXSAT][4]; /* satellite velocity/clk-rate std (m/s|s/s) */
    float  cov[MAXSAT][3]; /* satellite position covariance (m^2) */
    float  vco[MAXSAT][3]; /* satellite velocity covariance (m^2) */
} peph_t;

typedef struct {        /* precise clock type */
    gtime_t time;       /* time (GPST) */
    int index;          /* clock index for multiple files */
    double clk[MAXSAT][1]; /* satellite clock (s) */
    float  std[MAXSAT][1]; /* satellite clock std (s) */
} pclk_t;

typedef struct {        /* SBAS ephemeris type */
    int sat;            /* satellite number */
    gtime_t t0;         /* reference epoch time (GPST) */
    gtime_t tof;        /* time of message frame (GPST) */
    int sva;            /* SV accuracy (URA index) */
    int svh;            /* SV health (0:ok) */
    double pos[3];      /* satellite position (m) (ecef) */
    double vel[3];      /* satellite velocity (m/s) (ecef) */
    double acc[3];      /* satellite acceleration (m/s^2) (ecef) */
    double af0,af1;     /* satellite clock-offset/drift (s,s/s) */
} seph_t;

typedef struct {        /* norad two line element data type */
    char name [32];     /* common name */
    char alias[32];     /* alias name */
    char satno[16];     /* satellilte catalog number */
    char satclass;      /* classification */
    char desig[16];     /* international designator */
    gtime_t epoch;      /* element set epoch (UTC) */
    double ndot;        /* 1st derivative of mean motion */
    double nddot;       /* 2st derivative of mean motion */
    double bstar;       /* B* drag term */
    int etype;          /* element set type */
    int eleno;          /* element number */
    double inc;         /* orbit inclination (deg) */
    double OMG;         /* right ascension of ascending node (deg) */
    double ecc;         /* eccentricity */
    double omg;         /* argument of perigee (deg) */
    double M;           /* mean anomaly (deg) */
    double n;           /* mean motion (rev/day) */
    int rev;            /* revolution number at epoch */
} tled_t;

typedef struct {        /* norad two line element type */
    int n,nmax;         /* number/max number of two line element data */
    tled_t *data;       /* norad two line element data */
} tle_t;

typedef struct {        /* TEC grid type */
    gtime_t time;       /* epoch time (GPST) */
    int ndata[3];       /* TEC grid data size {nlat,nlon,nhgt} */
    double rb;          /* earth radius (km) */
    double lats[3];     /* latitude start/interval (deg) */
    double lons[3];     /* longitude start/interval (deg) */
    double hgts[3];     /* heights start/interval (km) */
    double *data;       /* TEC grid data (tecu) */
    float *rms;         /* RMS values (tecu) */
} tec_t;

typedef struct {        /* satellite fcb data type */
    gtime_t ts,te;      /* start/end time (GPST) */
    double bias[MAXSAT]; /* fcb value   (cyc) */
    double std [MAXSAT]; /* fcb std-dev (cyc) */
} fcbd_t;

typedef struct {
    int n,nmax;
    fcbd_t *data;
}fcbs_t;

typedef struct {
    gtime_t ts,te;
    double nl[MAXSAT];
    double std[MAXSAT];
}nl_upd_t;

typedef struct {
    int n,nmax;
    nl_upd_t *data;
}nl_upds_t;

typedef struct {
    double ewl[MAXSAT];
    double wl[MAXSAT];
}wl_upds_t;

typedef struct {
    nl_upds_t nls;
    wl_upds_t wls;
}upds_t;

typedef struct {        /* SBAS message type */
    int week,tow;       /* receiption time */
    uint8_t prn,rcv;            /* SBAS satellite PRN number */
    uint8_t msg[29]; /* SBAS message (226bit) padded by 0 */
} sbsmsg_t;

typedef struct {        /* SBAS messages type */
    int n,nmax;         /* number of SBAS messages/allocated */
    sbsmsg_t *msgs;     /* SBAS messages */
} sbs_t;

typedef struct {        /* SBAS fast correction type */
    gtime_t t0;         /* time of applicability (TOF) */
    double prc;         /* pseudorange correction (PRC) (m) */
    double rrc;         /* range-rate correction (RRC) (m/s) */
    double dt;          /* range-rate correction delta-time (s) */
    int iodf;           /* IODF (issue of date fast corr) */
    short udre;         /* UDRE+1 */
    short ai;           /* degradation factor indicator */
} sbsfcorr_t;

typedef struct {        /* SBAS long term satellite error correction type */
    gtime_t t0;         /* correction time */
    int iode;           /* IODE (issue of date ephemeris) */
    double dpos[3];     /* delta position (m) (ecef) */
    double dvel[3];     /* delta velocity (m/s) (ecef) */
    double daf0,daf1;   /* delta clock-offset/drift (s,s/s) */
} sbslcorr_t;

typedef struct {        /* SBAS satellite correction type */
    int sat;            /* satellite number */
    sbsfcorr_t fcorr;   /* fast correction */
    sbslcorr_t lcorr;   /* long term correction */
} sbssatp_t;

typedef struct {        /* SBAS satellite corrections type */
    int iodp;           /* IODP (issue of date mask) */
    int nsat;           /* number of satellites */
    int tlat;           /* system latency (s) */
    sbssatp_t sat[MAXSAT]; /* satellite correction */
} sbssat_t;

typedef struct {        /* SBAS ionospheric correction type */
    gtime_t t0;         /* correction time */
    short lat,lon;      /* latitude/longitude (deg) */
    short give;         /* GIVI+1 */
    float delay;        /* vertical delay estimate (m) */
} sbsigp_t;

typedef struct {        /* IGP band type */
    short x;            /* longitude/latitude (deg) */
    const short *y;     /* latitudes/longitudes (deg) */
    unsigned char bits; /* IGP mask start bit */
    unsigned char bite; /* IGP mask end bit */
} sbsigpband_t;

typedef struct {        /* SBAS ionospheric corrections type */
    int iodi;           /* IODI (issue of date ionos corr) */
    int nigp;           /* number of igps */
    sbsigp_t igp[MAXNIGP]; /* ionospheric correction */
} sbsion_t;

typedef struct {        /* DGPS/GNSS correction type */
    gtime_t t0;         /* correction time */
    double prc;         /* pseudorange correction (PRC) (m) */
    double rrc;         /* range rate correction (RRC) (m/s) */
    int iod;            /* issue of data (IOD) */
    double udre;        /* UDRE */
} dgps_t;

typedef struct {        /* SSR correction type */
    gtime_t t0[6];      /* epoch time (GPST) {eph,clk,hrclk,ura,bias,pbias} */
    double udi[6];      /* SSR update interval (s) */
    int iod[6];         /* iod ssr {eph,clk,hrclk,ura,bias,pbias} */
    int iode;           /* issue of data */
    int iodcrc;         /* issue of data crc for beidou/sbas */
    int ura;            /* URA indicator */
    int refd;           /* sat ref datum (0:ITRF,1:regional) */
    double deph [3];    /* delta orbit {radial,along,cross} (m) */
    double ddeph[3];    /* dot delta orbit {radial,along,cross} (m/s) */
    double dclk [3];    /* delta clock {c0,c1,c2} (m,m/s,m/s^2) */
    double hrclk;       /* high-rate clock corection (m) */
    float  cbias[MAXCODE]; /* code biases (m) */
    double pbias[MAXCODE]; /* phase biases (m) */
    float  stdpb[MAXCODE]; /* std-dev of phase biases (m) */
    double yaw_ang,yaw_rate; /* yaw angle and yaw rate (deg,deg/s) */
    unsigned char update; /* update flag (0:no update,1:update) */
} ssr_t;

typedef struct {        /* QZSS LEX message type */
    int prn;            /* satellite PRN number */
    int type;           /* message type */
    int alert;          /* alert flag */
    unsigned char stat; /* signal tracking status */
    unsigned char snr;  /* signal C/N0 (0.25 dBHz) */
    unsigned int ttt;   /* tracking time (ms) */
    unsigned char msg[212]; /* LEX message data part 1695 bits */
} lexmsg_t;

typedef struct {        /* QZSS LEX messages type */
    int n,nmax;         /* number of LEX messages and allocated */
    lexmsg_t *msgs;     /* LEX messages */
} lex_t;

typedef struct {        /* QZSS LEX ephemeris type */
    gtime_t toe;        /* epoch time (GPST) */
    gtime_t tof;        /* message frame time (GPST) */
    int sat;            /* satellite number */
    unsigned char health; /* signal health (L1,L2,L1C,L5,LEX) */
    unsigned char ura;  /* URA index */
    double pos[3];      /* satellite position (m) */
    double vel[3];      /* satellite velocity (m/s) */
    double acc[3];      /* satellite acceleration (m/s2) */
    double jerk[3];     /* satellite jerk (m/s3) */
    double af0,af1;     /* satellite clock bias and drift (s,s/s) */
    double tgd;         /* TGD */
    double isc[8];      /* ISC */
} lexeph_t;

typedef struct {        /* QZSS LEX ionosphere correction type */
    gtime_t t0;         /* epoch time (GPST) */
    double tspan;       /* valid time span (s) */
    double pos0[2];     /* reference position {lat,lon} (rad) */
    double coef[3][2];  /* coefficients lat x lon (3 x 2) */
} lexion_t;

typedef struct {        /* stec data type */
    gtime_t time;       /* time (GPST) */
    unsigned char sat;  /* satellite number */
    double ion;         /* slant ionos delay (m) */
    float std;          /* std-dev (m) */
    float azel[2];      /* azimuth/elevation (rad) */
    unsigned char flag; /* fix flag */
} stec_t;

typedef struct {        /* trop data type */
    gtime_t time;       /* time (GPST) */
    double trp[3];      /* zenith tropos delay/gradient (m) */
    float std[3];       /* std-dev (m) */
} trop_t;

typedef struct {        /* ppp corrections type */
    int nsta;           /* number of stations */
    char stas[MAXSTA][8]; /* station names */
    double rr[MAXSTA][3]; /* station ecef positions (m) */
    int ns[MAXSTA],nsmax[MAXSTA]; /* number of stec data */
    int nt[MAXSTA],ntmax[MAXSTA]; /* number of trop data */
    stec_t *stec[MAXSTA]; /* stec data */
    trop_t *trop[MAXSTA]; /* trop data */
} pppcorr_t;

typedef struct {
    double code[MAXSAT][MAXCODE];
    double phase[MAXSAT][MAXCODE];
}osb_t;

typedef struct {
    gtime_t tmin,tmax;
    double dt;
    osb_t *sat_osb;
}osbs_t;

typedef struct {        /* navigation data type */
    int n,nmax;         /* number of broadcast ephemeris */
    int ng,ngmax;       /* number of glonass ephemeris */
    int ns,nsmax;       /* number of sbas ephemeris */
    int ne,nemax;       /* number of precise ephemeris */
    int nc,ncmax;       /* number of precise clock */
    int na,namax;       /* number of almanac data */
    int nt,ntmax;       /* number of tec grid data */
    int nf,nfmax;       /* number of satellite fcb data */
    eph_t *eph;         /* GPS/QZS/GAL ephemeris */
    geph_t *geph;       /* GLONASS ephemeris */
    seph_t *seph;       /* SBAS ephemeris */
    peph_t *peph;       /* precise ephemeris */
    pclk_t *pclk;       /* precise clock */
    alm_t *alm;         /* almanac data */
    tec_t *tec;         /* tec grid data */
    fcbs_t *fcbs;
    osbs_t *osbs;
    upds_t *upds;
    erp_t  erp;         /* earth rotation parameters */
    double utc_gps[4];  /* GPS delta-UTC parameters {A0,A1,T,W} */
    double utc_glo[4];  /* GLONASS UTC GPS time parameters */
    double utc_gal[4];  /* Galileo UTC GPS time parameters */
    double utc_qzs[4];  /* QZS UTC GPS time parameters */
    double utc_cmp[4];  /* BeiDou UTC parameters */
    double utc_irn[4];  /* IRNSS UTC parameters */
    double utc_sbs[4];  /* SBAS UTC parameters */
    double ion_gps[8];  /* GPS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_gal[4];  /* Galileo iono model parameters {ai0,ai1,ai2,0} */
    double ion_qzs[8];  /* QZSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_cmp[8];  /* BeiDou iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_irn[8];  /* IRNSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    int leaps;          /* leap seconds (s) */
    double cbias[MAXSAT][12]; /* satellite dcb (0:p1-p2,1:p1-c1,2:p2-c2) (m) */
    double rbias[MAXRCV][2][3]; /* receiver dcb (0:p1-p2,1:p1-c1,2:p2-c2) (m) */
    double wlbias[MAXSAT];   /* wide-lane bias (cycle) */
    double glo_cpbias[4];    /* glonass code-phase bias {1C,1P,2C,2P} (m) */
    char glo_fcn[MAXPRNGLO+1]; /* glonass frequency channel number + 8 */
    pcv_t pcvs[MAXSAT]; /* satellite antenna pcv */
    sbssat_t sbssat;    /* SBAS satellite corrections */
    sbsion_t sbsion[MAXBAND+1]; /* SBAS ionosphere corrections */
    dgps_t dgps[MAXSAT]; /* DGPS corrections */
    ssr_t ssr[MAXSAT];  /* SSR corrections */
    lexeph_t lexeph[MAXSAT]; /* LEX ephemeris */
    lexion_t lexion;    /* LEX ionosphere correction */
    pppcorr_t pppcorr;  /* ppp corrections */
} nav_t;

typedef struct {        /* station parameter type */
    char name   [MAXANT]; /* marker name */
    char marker [MAXANT]; /* marker number */
    char antdes [MAXANT]; /* antenna descriptor */
    char antsno [MAXANT]; /* antenna serial number */
    char rectype[MAXANT]; /* receiver type descriptor */
    char recver [MAXANT]; /* receiver firmware version */
    char recsno [MAXANT]; /* receiver serial number */
    int antsetup;       /* antenna setup id */
    int itrf;           /* ITRF realization year */
    int deltype;        /* antenna delta type (0:enu,1:xyz) */
    double pos[3];      /* station position (ecef) (m) */
    double del[3];      /* antenna position delta (e/n/u or x/y/z) (m) */
    double hgt;         /* antenna height (m) */
    int glo_cp_align;   /* GLONASS code-phase alignment (0:no,1:yes) */
    double glo_cp_bias[4]; /* GLONASS code-phase biases {1C,1P,2C,2P} (m) */
} sta_t;

typedef struct {        /* solution type */
    gtime_t time;       /* time (GPST) */
    gtime_t eventime;   /* time of event (GPST) */
    int err_fmt;
    double rr[6];       /* position/velocity (m|m/s) */
                        /* {x,y,z,vx,vy,vz} or {e,n,u,ve,vn,vu} */
    double rpy[3];
    quat_t quat;
    double ba[3];
    double bg[3];
    double sa[3];
    double sg[3];
    double arm_gps[3];
    double  qr[6];       /* position variance/covariance (m^2) */
                        /* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or */
                        /* {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} */
    double  qv[6];       /* velocity variance/covariance (m^2/s^2) */
    double  qa[3];
    double  qba[3];
    double  qbg[3];
    double  qsa[3];
    double  qsg[3];
    double  qarm_gps[3];
    double dtr[6];      /* receiver clock bias to time systems (s) G R E B2 J B3*/
    double rdcb[6];
    double ifcb[6];
    double ztrp[2];     /* tropospheric delay (m) */
    double dop[4];
    unsigned char type; /* type (0:xyz-ecef,1:enu-baseline) */
    unsigned char stat; /* solution status (SOLQ_???) */
    unsigned char ns;   /* number of valid satellites */
    unsigned char ar_wl_ns;
    unsigned char ar_nl_ns;
    float age;          /* age of differential (s) */
    float ratio;        /* AR ratio factor for valiation */
    float prev_ratio1;   /* previous initial AR ratio factor for validation */
    float prev_ratio2;   /* previous final AR ratio factor for validation */
    float thres;        /* AR ratio threshold for valiation */
} sol_t;

typedef struct {        /* solution buffer type */
    int n,nmax;         /* number of solution/max number of buffer */
    int cyclic;         /* cyclic buffer flag */
    int start,end;      /* start/end index */
    gtime_t time;       /* current solution time */
    sol_t *data;        /* solution data */
    double rb[3];       /* reference position {x,y,z} (ecef) (m) */
    unsigned char buff[MAXSOLMSG+1]; /* message buffer */
    int nb;             /* number of byte in message buffer */
} solbuf_t;

typedef struct {        /* solution status type */
    gtime_t time;       /* time (GPST) */
    unsigned char sat;  /* satellite number */
    unsigned char frq;  /* frequency (1:L1,2:L2,...) */
    float az,el;        /* azimuth/elevation angle (rad) */
    float resp;         /* pseudorange residual (m) */
    float resc;         /* carrier-phase residual (m) */
    unsigned char flag; /* flags: (vsat<<5)+(slip<<3)+fix */
    unsigned char snr;  /* signal strength (0.25 dBHz) */
    unsigned short lock;  /* lock counter */
    unsigned short outc;  /* outage counter */
    unsigned short slipc; /* slip counter */
    unsigned short rejc;  /* reject counter */
    float amb;
    float var_amb;
    float nv_cp;
    float nv_pr;
} solstat_t;

typedef struct {        /* solution status buffer type */
    int n,nmax;         /* number of solution/max number of buffer */
    solstat_t *data;    /* solution status data */
} solstatbuf_t;

typedef struct {        /* RTCM control struct type */
    int staid;          /* station id */
    int stah;           /* station health */
    int seqno;          /* sequence number for rtcm 2 or iods msm */
    int outtype;        /* output message type */
    gtime_t time;       /* message time */
    gtime_t time_s;     /* message start time */
    obs_t obs;          /* observation data (uncorrected) */
    nav_t nav;          /* satellite ephemerides */
    sta_t sta;          /* station parameters */
    dgps_t *dgps;       /* output of dgps corrections */
    ssr_t ssr[MAXSAT];  /* output of ssr corrections */
    char msg[128];      /* special message */
    char msgtype[256];  /* last message type */
    char msmtype[6][128]; /* msm signal types */
    int obsflag;        /* obs data complete flag (1:ok,0:not complete) */
    int ephsat;         /* update satellite of ephemeris */
    int ephset;
    double cp[MAXSAT][NFREQ+NEXOBS]; /* carrier-phase measurement */
    unsigned short lock[MAXSAT][NFREQ+NEXOBS]; /* lock time */
    unsigned short loss[MAXSAT][NFREQ+NEXOBS]; /* loss of lock count */
    gtime_t lltime[MAXSAT][NFREQ+NEXOBS]; /* last lock time */
    int nbyte;          /* number of bytes in message buffer */ 
    int nbit;           /* number of bits in word buffer */ 
    int len;            /* message length (bytes) */
    unsigned char buff[1200]; /* message buffer */
    unsigned int word;  /* word buffer for rtcm 2 */
    unsigned int nmsg2[100]; /* message count of RTCM 2 (1-99:1-99,0:other) */
    unsigned int nmsg3[400]; /* message count of RTCM 3 (1-299:1001-1299,300-399:2000-2099,0:ohter) */
    char opt[256];      /* RTCM dependent options */
} rtcm_t;

typedef struct {        /* rinex control struct type */
    gtime_t time;       /* message time */
    double ver;         /* rinex version */
    char   type;        /* rinex file type ('O','N',...) */
    int    sys;         /* navigation system */
    int    tsys;        /* time system */
    char   tobs[7][MAXOBSTYPE][4]; /* rinex obs types */
    obs_t  obs;         /* observation data */
    nav_t  nav;         /* navigation data */
    sta_t  sta;         /* station info */
    int    ephsat;      /* ephemeris satellite number */
    int    ephset;
    char   opt[256];    /* rinex dependent options */
} rnxctr_t;

typedef struct {        /* download url type */
    char type[32];      /* data type */
    char path[1024];    /* url path */
    char dir [1024];    /* local directory */
    double tint;        /* time interval (s) */
} url_t;

typedef struct {        /* option type */
    const char *name;   /* option name */
    int format;         /* option format (0:int,1:double,2:string,3:enum) */
    void *var;          /* pointer to option variable */
    const char *comment; /* option comment/enum labels/unit */
} opt_t;

typedef struct {        /* extended receiver error model */
    int ena[4];         /* model enabled */
    double cerr[4][NFREQ*2]; /* code errors (m) */
    double perr[4][NFREQ*2]; /* carrier-phase errors (m) */
    double gpsglob[NFREQ]; /* gps-glonass h/w bias (m) */
    double gloicb [NFREQ]; /* glonass interchannel bias (m/fn) */
} exterr_t;

typedef struct {        /* SNR mask type */
    int ena[2];         /* enable flag {rover,base} */
    double mask[NFREQ][9]; /* mask (dBHz) at 5,10,...85 deg */
} snrmask_t;

typedef struct {
    int det;
    int ws;
    double mt;
    double sp;
    double gthres;
    double athres[3];
    double gyrothres[3];
    double sig_a,sig_g;
    double gamma[4];
}zvopt_t;

typedef struct {
    void *gnss_opt;
    int posfmt;
    int mech_coord;
    int local_coord;
    int imu_coord;
    int imu_align;
    int ms;
    int max_ny;
    int est_ba;
    int est_bg;
    int est_sa;
    int est_sg;
    int est_kod;
    int est_armgps;
    int est_igdt;
    int fb_pos;
    int fb_vel;
    int fb_att;
    int fb_ba;
    int fb_bg;
    int fb_sa;
    int fb_sg;
    int fb_armgps;
    int fb_igdt;
    int aid_cs;   /* ins aid cycle slip detection */
    int is_imu_samenoise;
    int odo;
    int rts;
    int zupt;
    int zaru;
    int nhc;

    double scslp;
    double feedratio;
    double zupt_var;
    double zaru_var;
    double nhc_var;

    imup_t imup;
    zvopt_t zvopt;

    double gps_loss_s;
    double gps_loss_last;
    int velopt;
    int back;
}insopt_t;

typedef struct {        /* processing options type */
    int prctype;              /* process type (0:batch, 1:single)*/
    char prcdir[MAXSTRPATH];  /* process dir */
    char site_list[1024];
    int site_idx;
    gtime_t ts;               /* process start time */
    gtime_t te;               /* process end time */
    char ac_name[6];
    double ti;                /* gnss obs sample rate*/
    int mode;           /* positioning mode (PMODE_???) */
    int soltype;        /* solution type (0:forward,1:backward,2:combined) */
    int nf;             /* number of frequencies (1:L1,2:L1+L2,3:L1+L2+L3,4:L1+L2+L3+L4) */
    int navsys;         /* navigation system */
    double elmin;       /* elevation mask angle (rad) */
    snrmask_t snrmask;  /* SNR mask */
    int sateph;         /* satellite ephemeris/clock (EPHOPT_???) */
    int modear;         /* AR mode (0:off,1:continuous,2:instantaneous,3:fix and hold,4:ppp-ar) */
    int glomodear;      /* GLONASS AR mode (0:off,1:on,2:auto cal,3:ext cal) */
    int gpsmodear;      /* GPS AR mode, debug/learning only (0:off,1:on) */
    int bdsmodear;      /* BeiDou AR mode (0:off,1:on) */
    int galmodear;      /* Galileo AR mode (0:off,1:on)*/
    int arfilter;       /* AR filtering to reject bad sats (0:off,1:on) */
    int maxout;         /* obs outage count to reset bias */
    int minlock;        /* min lock count to fix ambiguity */
    int minfixsats;     /* min sats to fix integer ambiguities */
    int minholdsats;    /* min sats to hold integer ambiguities */
    int mindropsats;    /* min sats to drop sats in AR */
    int minfix;         /* min fix count to hold ambiguity */
    int rcvstds;        /* use stdev estimates from receiver to adjust measurement variances */
    int armaxiter;      /* max iteration to resolve ambiguity */
    int ionoopt;        /* ionosphere option (IONOOPT_???) */
    int tropopt;        /* troposphere option (TROPOPT_???) */
    int dynamics;       /* dynamics model (0:none,1:velociy,2:accel) */
    int tidecorr;       /* earth tide correction (0:off,1:solid,2:solid+otl+pole) */
    int niter;          /* number of filter iteration */
    int codesmooth;     /* code smoothing window size (0:none) */
    int intpref;        /* interpolate reference obs (for post mission) */
    int sbascorr;       /* SBAS correction options */
    int sbassatsel;     /* SBAS satellite selection (0:all) */
    int rovpos;         /* rover position for fixed mode */
    int refpos;         /* base position for relative mode */
                        /* (0:pos in prcopt,  1:average of single pos, */
                        /*  2:read from file, 3:rinex header, 4:rtcm pos) */
    int weightmode;     /* weighting option (WEIGHTOPT_??) */
    double eratio[NFREQ]; /* code/phase error ratio */
    double err[6];      /* measurement error factor */
                        /* [0]:reserved */
                        /* [1-3]:error factor a/b/c of phase (m) */
                        /* [4]:doppler frequency (hz) */
                        /* [5]: snr max value (dB.Hz) */
    double std[3];      /* initial-state std [0]bias,[1]iono [2]trop */
    double prn[6];      /* process-noise std [0]bias,[1]iono [2]trop [3]acch [4]accv [5] pos */
    double sclkstab;    /* satellite clock stability (sec/sec) */
    double thresar[8];  /* AR validation threshold */
    double elmaskar;    /* elevation mask of AR for rising satellite (deg) */
    double elmaskhold;  /* elevation mask to hold ambiguity (deg) */
    double thresslip;   /* slip threshold of geometry-free phase (m) */
    double varholdamb;  /* variance for fix-and-hold psuedo measurements (cycle^2) */
    double gainholdamb; /* gain used for GLO and SBAS sats to adjust ambiguity */
    double maxtdiff;    /* max difference of time (sec) */
    double maxinno;     /* reject threshold of innovation (m) */
    double maxgdop;     /* reject threshold of gdop */
    double baseline[2]; /* baseline length constraint {const,sigma} (m) */
    double ru[3];       /* rover position for fixed mode {x,y,z} (ecef) (m) */
    double rb[3];       /* base position for relative mode {x,y,z} (ecef) (m) */
    char anttype[2][MAXANT]; /* antenna types {rover,base} */
    double antdel[2][3]; /* antenna delta {{rov_e,rov_n,rov_u},{ref_e,ref_n,ref_u}} */
    pcv_t pcvr[2];      /* receiver antenna parameters {rov,base} */
    unsigned char exsats[MAXSAT]; /* excluded satellites (1:excluded,2:included) */
    int  maxaveep;      /* max averaging epoches */
    int  initrst;       /* initialize by restart */
    int  outsingle;     /* output single by dgps/float/fix/ppp outage */
    char rnxopt[2][256]; /* rinex options {rover,base} */
    int  posopt[6];     /* positioning options */
    int  syncsol;       /* solution sync mode (0:off,1:on) */
    double odisp[2][6*11]; /* ocean tide loading parameters {rov,base} */
    exterr_t exterr;    /* extended receiver error model */
    int freqopt;        /* disable L2-AR */
    char pppopt[256];   /* ppp option */
    sigind_t sind[2][7];

    double sample;
    int adjobs;
    int bd3opt;
    int cbiaopt;
    int arprod;
    double cs_mw;
    double cs_gf;
    int gnss_frq_idx[NSYS+1][NFREQ];
    int tropmap;
    int robust;
    double igg_k0;
    double igg_k1;
    int kfopt;
    int sdopt;
    int velopt;
    int par;
    double restart;
    int aclong;
    char prdtype[6];
    char eph_int[6];
    char clk_int[6];
    char obsdir[100];
    int atx_week;
    char site_name[5];
    int geo_opt;
    insopt_t insopt;
} prcopt_t;

enum REFPOS{
    REFPOS_IMU,
    REFPOS_MANUAL,
    REFPOS_GPS,
    REFPOS_CAR
};

typedef struct {        /* solution options type */
    int posf;           /* solution format (SOLF_???) */
    int times;          /* time system (TIMES_???) */
    int timef;          /* time format (0:sssss.s,1:yyyy/mm/dd hh:mm:ss.s) */
    int timeu;          /* time digits under decimal point */
    int degf;           /* latitude/longitude format (0:ddd.ddd,1:ddd mm ss) */
    int outhead;        /* output header (0:no,1:yes) */
    int outopt;         /* output processing options (0:no,1:yes) */
    int outvel;         /* output velocity options (0:no,1:yes) */
    int outclk;
    int outtrp;
    int outba;
    int outbg;
    int outsa;
    int outsg;
    int outarm;
    int outigdt;
    int datum;          /* datum (0:WGS84,1:Tokyo) */
    int height;         /* height (0:ellipsoidal,1:geodetic) */
    int geoid;          /* geoid model (0:EGM96,1:JGD2000) */
    int solstatic;      /* solution of static mode (0:all,1:single) */
    int sstat;          /* solution statistics level (0:off,1:states,2:residuals) */
    int trace;          /* debug trace level (0:off,1-5:debug) */
    double nmeaintv[2]; /* nmea output interval (s) (<0:no,0:all) */
                        /* nmeaintv[0]:gprmc,gpgga,nmeaintv[1]:gpgsv */
    char sep[64];       /* field separator */
    char prog[64];      /* program name */
    double maxsolstd;   /* max std-dev for solution output (m) (0:all) */

    enum REFPOS sol_refpos;
    int outfrq;
    int ambres;
    int ref_type;
} solopt_t;

typedef struct {        /* file options type */
//    char satantp[MAXSTRPATH]; /* satellite antenna parameters file */
//    char rcvantp[MAXSTRPATH]; /* receiver antenna parameters file */
    char stapos [MAXSTRPATH]; /* station positions file */
    char geoid  [MAXSTRPATH]; /* external geoid data file */
    char tempdir[MAXSTRPATH]; /* ftp/http temporaly directory */
    char geexe  [MAXSTRPATH]; /* google earth exec file */
    char solstat[MAXSTRPATH]; /* solution statistics file */
    char trace  [MAXSTRPATH]; /* debug trace file */
    char navfile[MAXSTRPATH]; /* brdc file */
    char sp3file[MAXSTRPATH]; /* sp3 file  */
    char clkfile[MAXSTRPATH]; /* clk file  */
    char dcb    [MAXSTRPATH]; /* dcb data file */
    char rtsfile[MAXSTRPATH]; /* rts smoothing result file*/

    char robsf[MAXSTRPATH];
    char bobsf[MAXSTRPATH];
    char imuf[MAXSTRPATH];
    char *navf[3];
    char *sp3f[3];
    char *clkf[3];
    char *updf[3];
    char fcb[MAXSTRPATH];
    char bia[MAXSTRPATH];
    char atx[MAXSTRPATH];     /* antex file */
    char mgexdcb[MAXSTRPATH];
    char iono   [MAXSTRPATH]; /* ionosphere data file */
    char eop    [MAXSTRPATH]; /* eop data file */
    char blq    [MAXSTRPATH]; /* ocean tide loading blq file */
    char imupf[MAXSTRPATH];
    char posf[MAXSTRPATH];    /* avp pos file */

    char satf[MAXSTRPATH];
    char solf[MAXSTRPATH];    /* solution file */

    char rts_ins_fw[MAXSTRPATH]; /* file to save ins state for rts */
    char nl_amb[MAXSTRPATH];
    char wl_amb[MAXSTRPATH];
    char lc_amb[MAXSTRPATH];

    char strpath[MAXFILE][MAXSTRPATH];
    int  strtype[MAXFILE];
    int  strfmt[MAXFILE];
} filopt_t;

typedef struct {        /* RINEX options type */
    gtime_t ts,te;      /* time start/end */
    double tint;        /* time interval (s) */
    double ttol;        /* time tolerance (s) */
    double tunit;       /* time unit for multiple-session (s) */
    double rnxver;      /* RINEX version */
    int navsys;         /* navigation system */
    int obstype;        /* observation type */
    int freqtype;       /* frequency type */
    char mask[7][64];   /* code mask {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
    char staid [32];    /* station id for rinex file name */
    char prog  [32];    /* program */
    char runby [32];    /* run-by */
    char marker[64];    /* marker name */
    char markerno[32];  /* marker number */
    char markertype[32]; /* marker type (ver.3) */
    char name[2][32];   /* observer/agency */
    char rec [3][32];   /* receiver #/type/vers */
    char ant [3][32];   /* antenna #/type */
    double apppos[3];   /* approx position x/y/z */
    double antdel[3];   /* antenna delta h/e/n */
    double glo_cp_bias[4]; /* GLONASS code-phase biases (m) */
    char comment[MAXCOMMENT][64]; /* comments */
    char rcvopt[256];   /* receiver dependent options */
    uint8_t exsats[MAXSAT]; /* excluded satellites */
    int glofcn[32];     /* glonass fcn+8 */
    int outiono;        /* output iono correction */
    int outtime;        /* output time system correction */
    int outleaps;       /* output leap seconds */
    int autopos;        /* auto approx position */
    int phshift;        /* phase shift correction */
    int halfcyc;        /* half cycle correction */
    int sep_nav;        /* separated nav files */
    gtime_t tstart;     /* first obs time */
    gtime_t tend;       /* last obs time */
    gtime_t trtcm;      /* approx log start time for rtcm */
    char tobs[7][MAXOBSTYPE][4]; /* obs types {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
    double shift[7][MAXOBSTYPE]; /* phase shift (cyc) {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
    int nobs[7];        /* number of obs types {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
} rnxopt_t;

typedef struct {        /* satellite status type */
    uint8_t sys;  /* navigation system */
    uint8_t vs;   /* valid satellite flag single */
    double azel[2];     /* azimuth/elevation angles {az,el} (rad) */
    double resp[NFREQ]; /* residuals of pseudorange (m) */
    double resc[NFREQ]; /* residuals of carrier-phase (m) */
    double icbias[NFREQ];  /* glonass IC bias (cycles) */
    uint8_t vsat[NFREQ]; /* valid satellite flag */
    uint16_t snr[NFREQ]; /* signal strength (*SNR_UNIT dBHz) */
    uint16_t snr_rover [NFREQ]; /* rover signal strength (0.25 dBHz) */
    uint16_t snr_base  [NFREQ]; /* base signal strength (0.25 dBHz) */
    uint8_t fix [NFREQ]; /* ambiguity fix flag (1:float,2:fix,3:hold) */
    uint8_t slip[NFREQ]; /* cycle-slip flag */
    uint8_t half[NFREQ]; /* half-cycle valid flag */
    int lock [NFREQ];   /* lock counter of phase */
    uint32_t outc [NFREQ]; /* obs outage counter of phase */
    uint32_t slipc[NFREQ]; /* cycle-slip counter */
    uint32_t rejc [NFREQ]; /* reject counter */
    double eclipse;
    double  gf[NFREQ-1];   /* geometry-free phase L1-L2 (m) */
    double  mw[4];         /* MW-LC (m) mw,smw,mw_idx,mw_var*/
    double  phw;           /* phase windup (cycle) */
    double delta_mw[2];
    double delta_gf[2];
    double mtrp[2];         /* mtrph mtrpw*/
    double amb[NFREQ];
    double fix_amb[NFREQ];
    gtime_t ct;
    gtime_t pt[2][NFREQ]; /* previous carrier-phase time */
    double  ph[2][NFREQ]; /* previous carrier-phase observable (cycle) */
    double tec;
    double L[NFREQ];
    double P[NFREQ];
    double cor_L[NFREQ];
    double cor_P[NFREQ];
    double lam[NFREQ];
    double norm_v[2][NFREQ];
    double detect[NFREQ][2];      /* previous and now omc, using in ins aid cycle slip detection */
    double var_fact[2][NFREQ];
    int init_amb[NFREQ];
    int new_sat;
} ssat_t;

typedef struct {        /* ambiguity control type */
    gtime_t epoch[4];   /* last epoch */
    int n[4];           /* number of epochs */
    double LC [4];      /* linear combination average */
    double LCv[4];      /* linear combination variance */
    int fixcnt;         /* fix count */
    char flags[MAXSAT]; /* fix flags */
} ambc_t;

typedef struct {
    gtime_t t;
    int fix_wl_flag;
    int wl_fail_c;
    int fix_nl_flag;
    int ref_sat_no;
    int wl_fix;
    int nl_fix;
    double wl;
    double nl;
    double lc;
    double lc_fix;
    double lc_res;
    double wl_res;
    double nl_res;

    double count_fix;
}sdamb_t;

typedef struct {
    double wie; /**< rotation rate(rad s^-1) */
    double R0;  /**< Equatorial radius(m) */
    double RP;  /**< Polar radius(m) */
    double mu;  /**< gravitational constant, GM(m^3 s^-2) */
    double J2;  /**< 2nd-order gravitational Spherical Harmonics Function coefficient */
    double e;   /**< Eccentricity */
    double f;   /**< Flattening */
    v3_t pos;
    v3_t vel;
    v3_t wnie;
    v3_t wnen;
    v3_t wnin;
    v3_t wnien;
    double g;
    v3_t gn;
    v3_t gcc;
    double sl;
    double cl;
    double tl;
    double sl2;
    double RN;
    double RNh;
    double clRNh;
    double RM;
    double RMh;
    m3_t Mpv;
} earth_t;

/**
 * @brief ins solution struct
 */
typedef struct{
    gtime_t time;       /**< current solution time */
    unsigned int status;/**< solution status, see macro SOL_* */
    m3_t dcm;           /**< attitude in DCM */
    quat_t quat;        /**< attitude in quaternion */
    v3_t rpy;
    m3_t Qatt;          /**< var-covariance matrix of attitude */
    v3_t vel;           /**< velocity */
    m3_t Qvel;          /**< var-covariance matrix of velocity */
    v3_t acc;
    v3_t pos;           /**< position */
    m3_t Qpos;          /**< var-covariance matrix of postion  */
    v3_t ba;            /**< accelermeter bias */
    v3_t ba_std;        /**< standard error of accelermeter bais */
    v3_t bg;            /**< gryo bias */
    v3_t bg_std;        /**< standard error of gyro bias */
    v3_t sa;            /**< accelermeter scalar factor */
    v3_t sa_std;        /**< standard error of accelermeter scalar factor */
    v3_t sg;            /**< gyro scalar scalar factor */
    v3_t sg_std;        /**< standard error of gyro scalar factor */
    v3_t arm_gps;
    v3_t arm_gps_std;
    double t_delay;
    v3_t delay_pos;
    double kod;         /**< scalar factor of odometer, true/output */
    double std_kod;     /**< standard error of odometer scalar factor */
    m3_t Cbc;           /**< install error angle */
    v3_t std_Cbc;       /**< standard error of install error angle */
    double dtr[NSYS];
    int ns;
    int g_status;
    v3_t wib,wnb,web;
    v3_t fb,fn,an;
    v3_t Mpvvn;
    m3_t CW;
    m3_t MpvCnb;
    earth_t eth;
    int zero_flag;
} solins_t;             /**< ins solution struct */

/**
 * @brief kalman filter struct
 */
typedef struct{
    gtime_t time;       /**< current imu time */
    gtime_t last_couple_time;
    double idt;         /**< time interval of imu */
    double odt;         /**< time interval of od */
    int nx;             /**< length of full state of x */
    int ny;             /**< length of full state of y */
    int nix;
    int ngx;
    int na;
    double *x;          /**< state vector */
    double *P;          /**< var-covariance matrix */
    double *xa;
    double *Pa;
    double *F;          /**< transition matrix */
    double *Q;          /**< System noise covariance matrix */
    double *H;          /**< measurement matrix(transpose) */
    solins_t  *insstate;     /**< solution of kalman fileter */
    solins_t  *sol;
    gtime_t itg_start;  /**< intergral start time */
    double *itg;        /**< integral variables */
    double *R;          /**< necessary measurement noise(not for normal) */
    imud_t *imud;               /**< imu data list */
    unsigned short nimud;       /**< number of imud_t struct in kf_t.imud */
    unsigned short imudend;     /**< last imu data  */
    unsigned int ZST_count;     /**< zero speed test count */

    imud_t *imu_obs;
    imup_t *imup;
    int nsample;
    int couple_epoch;
    int ins_epoch;
    v3_t dthetap;
    v3_t dvp;
    v3_t omgb, fb;

} kf_t;     /**< kalman filter status struct */

typedef struct {
    gtime_t t;
    double raw_P[NFREQ];
    double raw_L[NFREQ];
    double cor_P[NFREQ];
    double cor_L[NFREQ];
    double LC_P[NFREQ];
    double LC_L[NFREQ];
    double code[NFREQ];

    double range;
    double dtr;
    double shaprio;
    double tide;
    double antr;
    double ants;

    double dts;
    double rs[3];
    double dcb[NFREQ];
    double trp[4];
    double ion[2];

    double phw;
    double amb[NFREQ];
    double LC_amb[NFREQ];
}sat_model_t;

typedef struct {        /* RTK control/result type */
    sol_t  sol;         /* RTK solution */
    double rb[6];       /* base position/velocity (ecef) (m|m/s) */
    int nx,na;          /* number of float states/fixed states */
    double tt;          /* time difference between current and previous (s) */
    int del_ep;
    double *x, *P;      /* float states and their covariance */
    double *xa,*Pa;     /* fixed states and their covariance */
    int nfix;           /* number of continuous fixes of ambiguity */
    int excsat;         /* index of next satellite to be excluded for partial ambiguity resolution */
    int nb_ar;          /* number of ambiguities used for AR last epoch */
	double com_bias;    /* phase bias common between all sats (used to be distributed to all sats */
    char holdamb;       /* set if fix-and-hold has occurred at least once */
    ambc_t ambc[MAXSAT]; /* ambiguity control */
    ssat_t ssat[MAXSAT]; /* satellite status */
    sat_model_t smod[MAXSAT];
    sdamb_t sdamb[MAXSAT];
    int neb;            /* bytes in error message buffer */
    char errbuf[MAXERRMSG]; /* error message buffer */
    prcopt_t opt;       /* processing options */
    int initial_mode;   /* initial positioning mode */
    int epoch;
    int tc;
    int stc;
    kf_t *ins_kf;
    int clk_jump;
    int fix_epoch;
    int pfix;         /* whether last epoch fixed or not 1:fix, 0:float*/
    int ins_aid;
    int ins_repair;
    int ins_rerun;
    double e[MAXSAT][3];
    int refsat[NSYS][NFREQ][2]; /* reference satellite of double-difference residuals* (0:gps/qzs/sbs,1:glo,2:gal,3:bds) */
    gtime_t filter_start;
    int exist_sys[NSYS+1];
    double dtrr;
} rtk_t;

typedef struct {
    int nv;            /* number of observation             */
    int npr;           /* number of pseudorange residual    */
    int ncp;           /* number of carrier phase residual  */

    double *pri_v;     /* priori residual include pseudorange and phase */
    double *post_v;    /* post residual include pseudorange and phase   */
    int   *vflag;      /* observation vaild flag              */
    int *pr_idx;       /* priori pseudorange residual index   */
    int *cp_idx;       /* priori carrier phase residual index */

    double sigma0;
    double *R;         /* covariance using for residual normalize */
    double *Qvv;       /* post covariance get from filter fun     */

    double *pri_pr;    /* priori pseudorange residual index   */
    double *pri_cp;    /* priori carrier phase residual       */
    double *post_pr;   /* post pseudorange residual           */
    double *post_cp;   /* post carrier phase residual         */

    double *norm_pr;   /* normalized post pseudorange residual    */
    double *norm_cp;   /* normalized post carrier phase residual  */
}res_t;

typedef struct half_cyc_tag {  /* half-cycle correction list type */
    unsigned char sat;  /* satellite number */
    unsigned char freq; /* frequency number (0:L1,1:L2,2:L5) */
    unsigned char valid; /* half-cycle valid flag */
    char corr;          /* half-cycle corrected (x 0.5 cyc) */
    gtime_t ts,te;      /* time start, time end */
    struct half_cyc_tag *next; /* pointer to next correction */
} half_cyc_t;

typedef struct {        /* receiver raw data control type */
    gtime_t time;       /* message time */
    gtime_t tobs[MAXSAT][NFREQ+NEXOBS]; /* observation data time */
    obs_t obs;          /* observation data */
    obs_t obuf;         /* observation data buffer */
    nav_t nav;          /* satellite ephemerides */
    sta_t sta;          /* station parameters */
    int ephsat;         /* sat number of update ephemeris (0:no satellite) */
    int ephset;         /* update set of ephemeris (0-1) */
    sbsmsg_t sbsmsg;    /* SBAS message */
    char msgtype[256];  /* last message type */
    unsigned char subfrm[MAXSAT][380];  /* subframe buffer */
    lexmsg_t lexmsg;    /* LEX message */
    double lockt[MAXSAT][NFREQ+NEXOBS]; /* lock time (s) */
    unsigned char lockflag[MAXSAT][NFREQ+NEXOBS]; /* used for carrying forward cycle slip */
    double icpp[MAXSAT],off[MAXSAT],icpc; /* carrier params for ss2 */
    double prCA[MAXSAT],dpCA[MAXSAT]; /* L1/CA pseudrange/doppler for javad */
    unsigned char halfc[MAXSAT][NFREQ+NEXOBS]; /* half-cycle add flag */
    char freqn[MAXOBS]; /* frequency number for javad */
    int nbyte;          /* number of bytes in message buffer */ 
    int len;            /* message length (bytes) */
    int iod;            /* issue of data */
    int tod;            /* time of day (ms) */
    int tbase;          /* time base (0:gpst,1:utc(usno),2:glonass,3:utc(su) */
    int flag;           /* general purpose flag */
    int outtype;        /* output message type */
    unsigned char buff[MAXRAWLEN]; /* message buffer */
    char opt[256];      /* receiver dependent options */
    half_cyc_t *half_cyc; /* half-cycle correction list */
    
    int format;         /* receiver stream format */
    void *rcv_data;     /* receiver dependent data */

    int dire;
    imud_t imu;
    imu_t imub;
    imu_t imut;
    int curb;
    void *optp;
    void *strp;
} raw_t;

typedef struct {        /* stream type */
    int type;           /* type (STR_???) */
    int mode;           /* mode (STR_MODE_?) */
    int state;          /* state (-1:error,0:close,1:open) */
    unsigned int inb,inr;   /* input bytes/rate */
    unsigned int outb,outr; /* output bytes/rate */
    unsigned int tick_i; /* input tick tick */
    unsigned int tick_o; /* output tick */
    unsigned int tact;  /* active tick */
    unsigned int inbt,outbt; /* input/output bytes at tick */
    lock_t lock;        /* lock flag */
    void *port;         /* type dependent port control struct */
    char path[MAXSTRPATH]; /* stream path */
    char msg [MAXSTRMSG];  /* stream message */
} stream_t;

typedef struct {        /* stream converter type */
    int itype,otype;    /* input and output stream type */
    int nmsg;           /* number of output messages */
    int msgs[32];       /* output message types */
    double tint[32];    /* output message intervals (s) */
    unsigned int tick[32]; /* cycle tick of output message */
    int ephsat[32];     /* satellites of output ephemeris */
    int stasel;         /* station info selection (0:remote,1:local) */
    rtcm_t rtcm;        /* rtcm input data buffer */
    raw_t raw;          /* raw  input data buffer */
    rtcm_t out;         /* rtcm output data buffer */
} strconv_t;

typedef struct {        /* stream server type */
    int state;          /* server state (0:stop,1:running) */
    int cycle;          /* server cycle (ms) */
    int buffsize;       /* input/monitor buffer size (bytes) */
    int nmeacycle;      /* NMEA request cycle (ms) (0:no) */
    int relayback;      /* relay back of output streams (0:no) */
    int nstr;           /* number of streams (1 input + (nstr-1) outputs */
    int npb;            /* data length in peek buffer (bytes) */
    char cmds_periodic[16][MAXRCVCMD]; /* periodic commands */
    double nmeapos[3];  /* NMEA request position (ecef) (m) */
    uint8_t *buff;      /* input buffers */
    uint8_t *pbuf;      /* peek buffer */
    uint32_t tick;      /* start tick */
    stream_t stream[16]; /* input/output streams */
    stream_t strlog[16]; /* return log streams */
    strconv_t *conv[16]; /* stream converter */
    thread_t thread;    /* server thread */
    lock_t lock;        /* lock flag */
} strsvr_t;

typedef struct {        /* RTK server type */
    int state;          /* server state (0:stop,1:running) */
    int cycle;          /* processing cycle (ms) */
    int nmeacycle;      /* NMEA request cycle (ms) (0:no req) */
    int nmeareq;        /* NMEA request (0:no,1:nmeapos,2:single sol) */
    double nmeapos[3];  /* NMEA request position (ecef) (m) */
    int buffsize;       /* input buffer size (bytes) */
    int format[3];      /* input format {rov,base,corr} */
    solopt_t solopt[2]; /* output solution options {sol1,sol2} */
    int navsel;         /* ephemeris select (0:all,1:rover,2:base,3:corr) */
    int nsbs;           /* number of sbas message */
    int nsol;           /* number of solution buffer */
    rtk_t rtk;          /* RTK control/result struct */
    int nb [3];         /* bytes in input buffers {rov,base} */
    int nsb[2];         /* bytes in soulution buffers */
    int npb[3];         /* bytes in input peek buffers */
    unsigned char *buff[3]; /* input buffers {rov,base,corr} */
    unsigned char *sbuf[2]; /* output buffers {sol1,sol2} */
    unsigned char *pbuf[3]; /* peek buffers {rov,base,corr} */
    sol_t solbuf[MAXSOLBUF]; /* solution buffer */
    unsigned int nmsg[3][10]; /* input message counts */
    raw_t  raw [3];     /* receiver raw control {rov,base,corr} */
    rtcm_t rtcm[3];     /* RTCM control {rov,base,corr} */
    gtime_t ftime[3];   /* download time {rov,base,corr} */
    char files[3][MAXSTRPATH]; /* download paths {rov,base,corr} */
    obs_t obs[3][MAXOBSBUF]; /* observation data {rov,base,corr} */
    nav_t nav;          /* navigation data */
    sbsmsg_t sbsmsg[MAXSBSMSG]; /* SBAS message buffer */
    stream_t stream[8]; /* streams {rov,base,corr,sol1,sol2,logr,logb,logc} */
    stream_t *moni;     /* monitor stream */
    unsigned int tick;  /* start tick */
    thread_t thread;    /* server thread */
    int cputime;        /* CPU time (ms) for a processing cycle */
    int prcout;         /* missing observation data count */
    int nave;           /* number of averaging base pos */
    double rb_ave[3];   /* averaging base pos */
    char cmds_periodic[3][MAXRCVCMD]; /* periodic commands */
    char cmd_reset[MAXRCVCMD]; /* reset command */
    double bl_reset;    /* baseline length to reset (km) */
    lock_t lock;        /* lock flag */
} rtksvr_t;

typedef struct {        /* gis data point type */
    double pos[3];      /* point data {lat,lon,height} (rad,m) */
} gis_pnt_t;

typedef struct {        /* gis data polyline type */
    int npnt;           /* number of points */
    double bound[4];    /* boundary {lat0,lat1,lon0,lon1} */
    double *pos;        /* position data (3 x npnt) */
} gis_poly_t;

typedef struct {        /* gis data polygon type */
    int npnt;           /* number of points */
    double bound[4];    /* boundary {lat0,lat1,lon0,lon1} */
    double *pos;        /* position data (3 x npnt) */
} gis_polygon_t;

typedef struct gisd_tag { /* gis data list type */
    int type;           /* data type (1:point,2:polyline,3:polygon) */
    void *data;         /* data body */
    struct gisd_tag *next; /* pointer to next */
} gisd_t;

typedef struct {        /* gis type */
    char name[MAXGISLAYER][256]; /* name */
    int flag[MAXGISLAYER];     /* flag */
    gisd_t *data[MAXGISLAYER]; /* gis data list */
    double bound[4];    /* boundary {lat0,lat1,lon0,lon1} */
} gis_t;

typedef void fatalfunc_t(const char *); /* fatal callback function type */

/* global variables ----------------------------------------------------------*/
extern const double tdistb_0250[];
extern const double tdistb_0100[];
extern const double tdistb_0050[];
extern const double tdistb_0025[];
extern const double tdistb_0010[];
extern const double tdistb_0005[];
extern const double tdistb_0001[];
extern const double chisqr[];        /* chi-sqr(n) table (alpha=0.001) */
extern const prcopt_t prcopt_default; /* default positioning options */
extern const solopt_t solopt_default; /* default solution output options */
extern const sbsigpband_t igpband1[9][8]; /* SBAS IGP band 0-8 */
extern const sbsigpband_t igpband2[2][5]; /* SBAS IGP band 9-10 */
extern const char *formatstrs[];     /* stream format strings */
extern opt_t sysopts[];              /* system options table */
extern opt_t insopts[];
extern const double Crf[9];               /* transform matrix of rfu-frame convert to frd-frame */

/* satellites, systems, codes functions --------------------------------------*/
EXPORT int  satno   (int sys, int prn);
EXPORT int  satsys  (int sat, int *prn);
EXPORT int  satsysidx(int sat);
EXPORT int  satid2no(const char *id);
EXPORT void satno2id(int sat, char *id);
EXPORT char *sat_id(int sat);
EXPORT uint8_t obs2code(const char *obs);
EXPORT char *code2obs(uint8_t code);
EXPORT int code2idx(int sys, uint8_t code);
EXPORT double code2freq(int sys, uint8_t code, int fcn);
EXPORT double sat2freq(int sat, uint8_t code, const nav_t *nav);
EXPORT void getobsfrqidx(char* frq_str,int sys,int nf,int *idxs);
EXPORT int  satexclude(int sat, double var, int svh, const prcopt_t *opt);
EXPORT int  testsnr(int base, int freq, double el, double snr, const snrmask_t *mask);
EXPORT int snrmask(const obsd_t *obs, const double *azel, const prcopt_t *opt);
EXPORT void setcodepri(int sys, int freq, const char *pri);
EXPORT int getcodepri(int sys, uint8_t code, const char *opt);

/* matrix and vector functions -----------------------------------------------*/
EXPORT int newround(double d);
EXPORT int findmax(const double *d,int n,double *max);
EXPORT int median(double *a, int n);
EXPORT double std_vec(double *a,int n);
EXPORT double rmse_vec(double *a,int n);
EXPORT double *mat  (int n, int m);
EXPORT double *mat_scale(int n,double a);
EXPORT int    *imat (int n, int m);
EXPORT double *zeros(int n, int m);
EXPORT double *eye  (int n);
EXPORT double dot (const double *a, const double *b, int n);
EXPORT double norm(const double *a, int n);
EXPORT void cross3(const double *a, const double *b, double *c);
EXPORT int  normv3(const double *a, double *b);
EXPORT void matcpy(double *A, const double *B, int n, int m);
EXPORT void asi_blk_mat(double *A,int m,int n,const double *B,int p ,int q,
                        int isr,int isc);
EXPORT void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C);
EXPORT int  matinv(double *A, int n);
EXPORT void matblock(const double *A,int r,int c,double *B,int p,int q,int isr,int isc);
EXPORT void asignmat(double *A,int r,int c,const double *B,int p,int q,int isr,int isc);
EXPORT void matmul33(const char *tr,const double *A,const double *B,const double *C,
                     int n,int p,int q,int m,double *D);
EXPORT void matmul3v(const char *tr, const double *A, const double *b, double *c);
EXPORT int  solve (const char *tr, const double *A, const double *Y, int n,
                   int m, double *X);
EXPORT int  lsq(const double *A, const double *y, int n, int m, double *x,
                   double *Q);
EXPORT int  lsq_(const double *H,const double *R, const double *y, int n, int m, double *x,
                   double *Q);
EXPORT int  filter(double *x, double *P, const double *H, const double *v,
                   const double *R, int n, int m,int qc, int kf_type,res_t *res,int tc);
EXPORT int  smoother(const double *xf, const double *Qf, const double *xb,
                     const double *Qb, int n, double *xs, double *Qs);
EXPORT void matprint (int trans,const double *A, int n, int m, int p, int q);
EXPORT void matfprint(const double *A, int n, int m, int p, int q, FILE *fp);

EXPORT void add_fatal(fatalfunc_t *func);

/* time and string functions -------------------------------------------------*/
EXPORT void trimspace(char *dstsrc);
EXPORT void strmid(char *dst, const char *src, int nPos, int nCount);
EXPORT void cutfilepathsep(char *strPath,int n,int m,int k,char *other_str);
EXPORT void upstr(const char *s,char *s1,int n);
EXPORT void setstr(char *dst, const char *src, int n);
EXPORT void strmcpy(const char* src,const char* str_flag,int start_idx,int m,char *s);
EXPORT double  str2num(const char *s, int i, int n);
EXPORT int     str2time(const char *s, int i, int n, gtime_t *t);
EXPORT int     str2time1(const char *s, int i, int n, gtime_t *t);
EXPORT void    time2str(gtime_t t, char *str, int n);
EXPORT gtime_t epoch2time(const double *ep);
EXPORT void    time2epoch(gtime_t t, double *ep);
EXPORT gtime_t gpst2time(int week, double sec);
EXPORT double  time2gpst(gtime_t t, int *week);
EXPORT gtime_t gst2time(int week, double sec);
EXPORT double  time2gst(gtime_t t, int *week);
EXPORT gtime_t bdt2time(int week, double sec);
EXPORT double  time2bdt(gtime_t t, int *week);
EXPORT double  jd2mjd(const double jd);
EXPORT double  time2mjday(const gtime_t t);
EXPORT double  epoch2mjday(const double *ep);
EXPORT void time2mjd(gtime_t t,mjd_t *mjd);
EXPORT void mjd2time(const mjd_t *mjd,gtime_t *t);
EXPORT char    *time_str(gtime_t t, int n);

EXPORT gtime_t timeadd  (gtime_t t, double sec);
EXPORT double  timediff (gtime_t t1, gtime_t t2);
EXPORT gtime_t gpst2utc (gtime_t t);
EXPORT gtime_t utc2gpst (gtime_t t);
EXPORT gtime_t gpst2bdt (gtime_t t);
EXPORT gtime_t bdt2gpst (gtime_t t);
EXPORT gtime_t timeget  (void);
EXPORT void    timeset  (gtime_t t);
EXPORT void    timereset(void);
EXPORT double  time2doy (gtime_t t);
EXPORT double  utc2gmst (gtime_t t, double ut1_utc);
EXPORT int read_leaps(const char *file);

EXPORT int adjgpsweek(int week);
EXPORT unsigned int tickget(void);
EXPORT void sleepms(int ms);

EXPORT int reppath(const char *path, char *rpath, gtime_t time, const char *rov,
                   const char *base);
EXPORT int reppaths(const char *path, char *rpaths[], int nmax, gtime_t ts,
                    gtime_t te, const char *rov, const char *base);

/* coordinates transformation ------------------------------------------------*/
EXPORT void ecef2pos(const double *r, double *pos);
EXPORT void pos2ecef(const double *pos, double *r);
EXPORT void dpos2decef(const double *pos, double *C);
EXPORT void ecef2enu(const double *pos, const double *r, double *e);
EXPORT void enu2ecef(const double *pos, const double *e, double *r);
EXPORT void covenu  (const double *pos, const double *P, double *Q);
EXPORT void covecef (const double *pos, const double *Q, double *P);
EXPORT void xyz2enu (const double *pos, double *E);
EXPORT void ned2xyz(const double *pos,double *Cne);
EXPORT int llh2ecef(v3_t* pos, v3_t* vel, m3_t* Cbn,int opt);
EXPORT int ecef2llh(v3_t* xyz, v3_t* vel, m3_t* Cbe,int opt);
EXPORT int llh2ecefQ(const v3_t *pos, m3_t *Qpos, m3_t *Qvel, m3_t *Qatt,int opt);
EXPORT int ecef2llhQ(const v3_t *xyz, m3_t *Qxyz, m3_t *Qvel, m3_t *Qatt,int opt);
EXPORT void eci2ecef(gtime_t tutc, const double *erpv, double *U, double *gmst);
EXPORT void deg2dms (double deg, double *dms, int ndec);
EXPORT double dms2deg(const double *dms);

/* input and output functions ------------------------------------------------*/
EXPORT int readobsnav(gtime_t ts, gtime_t te, double ti, char **infile,
               const int *index, int n, prcopt_t *prcopt,
               obs_t *obs, nav_t *nav, sta_t *sta);
EXPORT void freeobsnav(obs_t *obs, nav_t *nav);
EXPORT void adjustobs(const prcopt_t *popt,const obsd_t *obss,obsd_t *adj_obss,int n);
EXPORT void readpreceph(char **infile, int n, const prcopt_t *prcopt,const filopt_t *fopt,
                        nav_t *nav, sbs_t *sbs, lex_t *lex);
EXPORT void freepreceph(nav_t *nav, sbs_t *sbs, lex_t *lex);
EXPORT void readpos(const char *file, const char *rcv, double *pos);
EXPORT int  sortobs(obs_t *obs);
EXPORT void uniqnav(nav_t *nav);
EXPORT int  screent(gtime_t time, gtime_t ts, gtime_t te, double tint);
EXPORT int  readnav(const char *file, nav_t *nav);
EXPORT int  savenav(const char *file, const nav_t *nav);
EXPORT void freeobs(obs_t *obs);
EXPORT void freenav(nav_t *nav, int opt);
EXPORT int  readblq(const char *file, const char *sta, double *odisp);
EXPORT int  readerp(const char *file, erp_t *erp);
EXPORT int  geterp (const erp_t *erp, gtime_t time, double *val);
EXPORT int nextobsf(const obs_t *obs, int *i, int rcv);
EXPORT int nextobsb(const obs_t *obs, int *i, int rcv);
EXPORT int inputgnss(obsd_t *obs, int solq, const prcopt_t *popt,const obs_t *obss,int *rover_idx,int *base_idx,int revs);
EXPORT int read_gnss_file(prcopt_t *popt, const filopt_t *fopt, int mode,int need_base, gtime_t ts, gtime_t te, obs_t *obss, nav_t *navs);

/* debug trace functions -----------------------------------------------------*/
EXPORT void traceopen(const char *file);
EXPORT void traceclose(void);
EXPORT void tracelevel(int level);
EXPORT void trace    (int level, const char *format, ...);
EXPORT void tracet   (int level, const char *format, ...);
EXPORT void tracemat (int level, const double *A, int n, int m, int p, int q);
EXPORT void traceobs (int level, const obsd_t *obs, int n);
EXPORT void tracenav (int level, const nav_t *nav);
EXPORT void tracegnav(int level, const nav_t *nav);
EXPORT void tracehnav(int level, const nav_t *nav);
EXPORT void tracepeph(int level, const nav_t *nav);
EXPORT void tracepclk(int level, const nav_t *nav);
EXPORT void traceb   (int level, const unsigned char *p, int n);
EXPORT int gettracelevel(void);
EXPORT void tracestdout(void);
EXPORT void traceins(const solins_t *sol_ins,int post,const insopt_t *opt);

/* platform dependent functions ----------------------------------------------*/
EXPORT int execcmd(const char *cmd);
EXPORT int expath (const char *path, char *paths[], int nmax);
EXPORT void createdir(const char *path);

/* positioning models --------------------------------------------------------*/
EXPORT double satazel(const double *pos, const double *e, double *azel);
EXPORT double geodist(const double *rs, const double *rr, double *e);
EXPORT void dops(int ns, const double *azel, double elmin, double *dop);
EXPORT void csmooth(obs_t *obs, int ns);

/* code bias model*/
EXPORT double corrISC(const prcopt_t *popt,const double *cbias,uint8_t code,int sat);
EXPORT double corrDCB(const prcopt_t *popt,const nav_t *nav, const double *cbias,uint8_t code,int frq,int sat);
EXPORT double corr_code_bias(const prcopt_t *popt,const nav_t *nav,const obsd_t *obs,int frq);
EXPORT int cal_eclips(int prn, double *satp, const double *satv, double *sunp,
                      double TTAG, double SANTXYZ[3], const nav_t *nav,int hour);

/* atmosphere models ---------------------------------------------------------*/
EXPORT double ionmapf(const double *pos, const double *azel);
EXPORT double ionppp(const double *pos, const double *azel, double re,
                     double hion, double *pppos);
EXPORT int iontec(gtime_t time, const nav_t *nav, const double *pos,
                  const double *azel, int opt, double *delay, double *var);
EXPORT double iontecvar(int ep,gtime_t t,const double *pos,const double *azel);
EXPORT void readtec(const char *file, nav_t *nav, int opt);
EXPORT double klobuchar_GPS(gtime_t t, const double *ion, const double *pos,
                       const double *azel);
extern double klobuchar_BDS(gtime_t t, const double *ion, const double *pos,
                            const double *azel);
extern int model_iono(gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, int sat, const double *x,
                      const nav_t *nav, double *dion, double *var);

EXPORT double saastamoinen(gtime_t time, const double *pos, const double *azel,
                        double humi,int gpt,double *ztrph,double *ztrpw);
EXPORT double trop_UNB3(gtime_t t,double *pos,double el,double *ztrpw);
EXPORT double tropmapf(const prcopt_t *popt,gtime_t time, const double *pos, const double *azel,
                       double *mapfw);
EXPORT int model_trop(gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, const double *x, double *dtdx,
                      const nav_t *nav, double *dtrp,double *ztrp,double *mtrp, double *var,int it);

EXPORT int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var);
EXPORT int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, int tropopt, double *trp, double *var,double *ztrph,double *ztrpw);

/* antenna models ------------------------------------------------------------*/
EXPORT int  readpcv(const char *file, pcvs_t *pcvs);
EXPORT pcv_t *searchpcv(int sat, const char *type, gtime_t time,
                        const pcvs_t *pcvs);
EXPORT void antmodel(int sat,const pcv_t *pcv, const double *del, const double *azel,
                     int opt, double *dant);
EXPORT void antmodel_s(const pcv_t *pcv, double nadir, double *dant);

/* earth tide models ---------------------------------------------------------*/
EXPORT void sunmoonpos(gtime_t tutc, const double *erpv, double *rsun,
                       double *rmoon, double *gmst);
EXPORT void tidedisp(gtime_t tutc, const double *rr, int opt, const erp_t *erp,
                     const double *odisp, double *dr);

/* geiod models --------------------------------------------------------------*/
EXPORT int opengeoid(int model, const char *file);
EXPORT void closegeoid(void);
EXPORT double geoidh(const double *pos);

/* datum transformation ------------------------------------------------------*/
EXPORT int loaddatump(const char *file);
EXPORT int tokyo2jgd(double *pos);
EXPORT int jgd2tokyo(double *pos);

/* rinex functions -----------------------------------------------------------*/
EXPORT int readrnx (const prcopt_t *popt,const char *file, int rcv, const char *opt, obs_t *obs,
                    nav_t *nav, sta_t *sta);
EXPORT int readrnxt(const prcopt_t *popt,const char *file, int rcv, gtime_t ts, gtime_t te,
                    double tint, const char *opt, obs_t *obs, nav_t *nav,
                    sta_t *sta);
EXPORT int readrnxc(const prcopt_t *popt,const char *file, nav_t *nav);
EXPORT int outrnxobsh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxobsb(FILE *fp, const rnxopt_t *opt, const obsd_t *obsd, int n, int flag);
EXPORT int outrnxnavh (FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxgnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxhnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxlnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxqnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxcnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxinavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxnavb (FILE *fp, const rnxopt_t *opt, const eph_t *eph);
EXPORT int outrnxgnavb(FILE *fp, const rnxopt_t *opt, const geph_t *geph);
EXPORT int outrnxhnavb(FILE *fp, const rnxopt_t *opt, const seph_t *seph);
EXPORT int rtk_uncompress(const char *file, char *uncfile);
EXPORT int convrnx(const prcopt_t *popt,int format, rnxopt_t *opt, const char *file, char **ofile);
EXPORT int  init_rnxctr (rnxctr_t *rnx);
EXPORT void free_rnxctr (rnxctr_t *rnx);
EXPORT int  open_rnxctr (rnxctr_t *rnx, FILE *fp);
EXPORT int  input_rnxctr(const prcopt_t *popt,rnxctr_t *rnx, FILE *fp);

/* ephemeris and clock functions ---------------------------------------------*/
EXPORT double eph2clk (gtime_t time, const eph_t  *eph);
EXPORT double geph2clk(gtime_t time, const geph_t *geph);
EXPORT double seph2clk(gtime_t time, const seph_t *seph);
EXPORT void eph2pos (gtime_t time, const eph_t  *eph,  double *rs, double *dts,
                     double *var);
EXPORT void geph2pos(gtime_t time, const geph_t *geph, double *rs, double *dts,
                     double *var);
EXPORT void seph2pos(gtime_t time, const seph_t *seph, double *rs, double *dts,
                     double *var);
EXPORT int  peph2pos(const prcopt_t *popt,gtime_t time, int sat, const nav_t *nav, int opt,
                     double *rs, double *dts, double *var);
EXPORT void satantoff(const prcopt_t *popt,gtime_t time, const double *rs, int sat, const nav_t *nav,
                      double *dant);
EXPORT int  satpos(const prcopt_t *popt,gtime_t time, gtime_t teph, int sat, int ephopt,
                   const nav_t *nav, double *rs, double *dts, double *var,
                   int *svh);
EXPORT void satposs(const prcopt_t *popt,gtime_t time, const obsd_t *obs, int n, const nav_t *nav,
                    int sateph, double *rs, double *dts, double *var, int *svh);
EXPORT void satseleph(int sys, int sel);
EXPORT int  getseleph(int sys);
EXPORT void readsp3(const char *file, nav_t *nav, int opt);
EXPORT int  readsap(const char *file, gtime_t time, nav_t *nav);
EXPORT int  readdcb(const char *file, nav_t *nav, const sta_t *sta);
EXPORT int  readdcb_mgex(const char *file,const prcopt_t *popt, nav_t *nav);
EXPORT int  readfcb(const char *file, nav_t *nav);
EXPORT int  readosb(const char *file, nav_t *nav);
EXPORT int  readupd(const prcopt_t *opt,char *file_ewl,char *file_wl,char *file_nl, nav_t *nav);
EXPORT void readotl(prcopt_t *popt, const char *file, const sta_t *sta);
EXPORT void setpcv(gtime_t time, prcopt_t *popt, nav_t *nav, const pcvs_t *pcvs,
                   const pcvs_t *pcvr, const sta_t *sta);
EXPORT void alm2pos(gtime_t time, const alm_t *alm, double *rs, double *dts);

EXPORT int tle_read(const char *file, tle_t *tle);
EXPORT int tle_name_read(const char *file, tle_t *tle);
EXPORT int tle_pos(gtime_t time, const char *name, const char *satno,
                   const char *desig, const tle_t *tle, const erp_t *erp,
                   double *rs);

/* receiver raw data functions -----------------------------------------------*/
EXPORT uint32_t getbitu(const uint8_t *buff, int pos, int len);
EXPORT int32_t          getbits(const uint8_t *buff, int pos, int len);
EXPORT void setbitu(uint8_t *buff, int pos, int len, unsigned int data);
EXPORT void setbits(uint8_t *buff, int pos, int len, int data);
EXPORT uint32_t rtk_crc32  (const uint8_t *buff, int len);
EXPORT uint32_t rtk_crc24q (const uint8_t *buff, int len);
EXPORT uint16_t rtk_crc16(const uint8_t *buff, int len);
EXPORT int decode_word (unsigned int word, unsigned char *data);
EXPORT int decode_frame(const uint8_t *buff, eph_t *eph, alm_t *alm,
                        double *ion, double *utc);
EXPORT int test_glostr(const uint8_t *buff);
EXPORT int decode_glostr(const uint8_t *buff, geph_t *geph, double *utc);
EXPORT int decode_bds_d1(const uint8_t *buff, eph_t *eph, double *ion,
                         double *utc);
EXPORT int decode_bds_d2(const uint8_t *buff, eph_t *eph, double *utc);
EXPORT int decode_gal_inav(const uint8_t *buff, eph_t *eph, double *ion,
                           double *utc);
EXPORT int decode_gal_fnav(const uint8_t *buff, eph_t *eph, double *ion,
                           double *utc);
EXPORT int decode_irn_nav(const uint8_t *buff, eph_t *eph, double *ion,
                          double *utc);

EXPORT int init_raw   (raw_t *raw, int format);
EXPORT void free_raw  (raw_t *raw);
EXPORT int input_raw  (raw_t *raw, int format, unsigned char data);
EXPORT int input_rawf (raw_t *raw, int format, FILE *fp);

EXPORT int init_rt17  (raw_t *raw);
EXPORT int init_cmr   (raw_t *raw);
EXPORT void free_rt17 (raw_t *raw);
EXPORT void free_cmr  (raw_t *raw);
EXPORT int update_cmr (raw_t *raw, rtksvr_t *svr, obs_t *obs);

EXPORT int input_oem4  (raw_t *raw, unsigned char data);
EXPORT int input_cnav  (raw_t *raw, unsigned char data);
EXPORT int input_ubx   (raw_t *raw, unsigned char data);
EXPORT int input_sbp   (raw_t *raw, unsigned char data);
EXPORT int input_cres  (raw_t *raw, unsigned char data);
EXPORT int input_stq   (raw_t *raw, unsigned char data);
EXPORT int input_gw10  (raw_t *raw, unsigned char data);
EXPORT int input_javad (raw_t *raw, unsigned char data);
EXPORT int input_nvs   (raw_t *raw, unsigned char data);
EXPORT int input_bnx   (raw_t *raw, unsigned char data);
EXPORT int input_rt17  (raw_t *raw, unsigned char data);
EXPORT int input_sbf   (raw_t *raw, unsigned char data);
EXPORT int input_cmr   (raw_t *raw, unsigned char data);
EXPORT int input_tersus(raw_t *raw, unsigned char data);
EXPORT int input_lexr  (raw_t *raw, unsigned char data);
EXPORT int input_oem4f (raw_t *raw, FILE *fp);
EXPORT int input_cnavf (raw_t *raw, FILE *fp);
EXPORT int input_ubxf  (raw_t *raw, FILE *fp);
EXPORT int input_sbpf  (raw_t *raw, FILE *fp);
EXPORT int input_cresf (raw_t *raw, FILE *fp);
EXPORT int input_stqf  (raw_t *raw, FILE *fp);
EXPORT int input_gw10f (raw_t *raw, FILE *fp);
EXPORT int input_javadf(raw_t *raw, FILE *fp);
EXPORT int input_nvsf  (raw_t *raw, FILE *fp);
EXPORT int input_bnxf  (raw_t *raw, FILE *fp);
EXPORT int input_rt17f (raw_t *raw, FILE *fp);
EXPORT int input_sbff  (raw_t *raw, FILE *fp);
EXPORT int input_cmrf  (raw_t *raw, FILE *fp);
EXPORT int input_tersusf(raw_t *raw, FILE *fp);
EXPORT int input_lexrf (raw_t *raw, FILE *fp);

EXPORT int gen_ubx (const char *msg, unsigned char *buff);
EXPORT int gen_stq (const char *msg, unsigned char *buff);
EXPORT int gen_nvs (const char *msg, unsigned char *buff);
EXPORT int gen_lexr(const char *msg, unsigned char *buff);

/* rtcm functions ------------------------------------------------------------*/
EXPORT int init_rtcm   (rtcm_t *rtcm);
EXPORT void free_rtcm  (rtcm_t *rtcm);
EXPORT int input_rtcm2 (rtcm_t *rtcm, unsigned char data);
EXPORT int input_rtcm3 (rtcm_t *rtcm, unsigned char data);
EXPORT int input_rtcm2f(rtcm_t *rtcm, FILE *fp);
EXPORT int input_rtcm3f(rtcm_t *rtcm, FILE *fp);
EXPORT int gen_rtcm2   (rtcm_t *rtcm, int type, int sync);
EXPORT int gen_rtcm3   (rtcm_t *rtcm, int type, int subtype, int sync);

/* solution functions --------------------------------------------------------*/
EXPORT void initsolbuf(solbuf_t *solbuf, int cyclic, int nmax);
EXPORT void freesolbuf(solbuf_t *solbuf);
EXPORT void freesolstatbuf(solstatbuf_t *solstatbuf);
EXPORT sol_t *getsol(solbuf_t *solbuf, int index);
EXPORT int addsol(solbuf_t *solbuf, const sol_t *sol);
EXPORT int readsol (char *files[], int nfile, solbuf_t *sol);
EXPORT int readsolt(char *files[], int nfile, gtime_t ts, gtime_t te,
                    double tint, int qflag, solbuf_t *sol,int coord,const solopt_t *sopt);
EXPORT int readref(const char *file,solbuf_t *sol,const solopt_t *sopt,int week,double sow);
EXPORT int readsolstat(char *files[], int nfile, solstatbuf_t *statbuf);
EXPORT int readsolstatt(char *files[], int nfile, gtime_t ts, gtime_t te,
                        double tint, solstatbuf_t *statbuf);
EXPORT int inputsol(unsigned char data, gtime_t ts, gtime_t te, double tint,
                    int qflag, const solopt_t *opt, solbuf_t *solbuf,int coord);
EXPORT void readsolopt(FILE *fp, solopt_t *opt);
EXPORT int readsoldata(FILE *fp, gtime_t ts, gtime_t te, double tint, int qflag,
                       const solopt_t *opt, solbuf_t *solbuf,int coord);
EXPORT int outprcopts(unsigned char *buff, const prcopt_t *opt);
EXPORT int outhead(const char *outfile, char **infile, int n,
                   const prcopt_t *popt, const solopt_t *sopt);
EXPORT int outsolheads(unsigned char *buff, const solopt_t *opt,const prcopt_t *popt);
EXPORT int outsols  (unsigned char *buff, const sol_t *sol, const double *rb,
                     const solopt_t *opt,const prcopt_t *popt,const solins_t *ins_sol);
EXPORT int outsolexs(unsigned char *buff, const sol_t *sol, const ssat_t *ssat, const solopt_t *opt);
EXPORT void outprcopt(FILE *fp, const prcopt_t *opt);
EXPORT void outsolhead(FILE *fp, const solopt_t *opt,const prcopt_t *popt);
EXPORT void outsol  (FILE *fp, const sol_t *sol, const double *rb,
                     const solopt_t *opt,const prcopt_t *popt, const solins_t *ins_sol);
EXPORT void outsolex(FILE *fp, const sol_t *sol, const ssat_t *ssat,
                     const solopt_t *opt);
EXPORT int outnmea_rmc(unsigned char *buff, const sol_t *sol);
EXPORT int outnmea_gga(unsigned char *buff, const sol_t *sol);
EXPORT int outnmea_gsa(unsigned char *buff, const sol_t *sol, const ssat_t *ssat);
EXPORT int outnmea_gsv(unsigned char *buff, const sol_t *sol, const ssat_t *ssat);
EXPORT void outsatinfo(FILE *fp, const prcopt_t *popt,
        const rtk_t *rtk,const obsd_t *obs,int nobs);
/* google earth kml converter ------------------------------------------------*/
EXPORT int convkml(const char *infile, const char *outfile, gtime_t ts,
                   gtime_t te, double tint, int qflg, double *offset,
                   int tcolor, int pcolor, int outalt, int outtime);

/* gpx converter -------------------------------------------------------------*/
EXPORT int convgpx(const char *infile, const char *outfile, gtime_t ts,
                   gtime_t te, double tint, int qflg, double *offset,
                   int outtrk, int outpnt, int outalt, int outtime);

/* sbas functions ------------------------------------------------------------*/
EXPORT int  sbsreadmsg (const char *file, int sel, sbs_t *sbs);
EXPORT int  sbsreadmsgt(const char *file, int sel, gtime_t ts, gtime_t te,
                        sbs_t *sbs);
EXPORT void sbsoutmsg(FILE *fp, sbsmsg_t *sbsmsg);
EXPORT int  sbsdecodemsg(gtime_t time, int prn, const unsigned int *words,
                         sbsmsg_t *sbsmsg);
EXPORT int sbsupdatecorr(const sbsmsg_t *msg, nav_t *nav);
EXPORT int sbssatcorr(gtime_t time, int sat, const nav_t *nav, double *rs,
                      double *dts, double *var);
EXPORT int sbsioncorr(gtime_t time, const nav_t *nav, const double *pos,
                      const double *azel, double *delay, double *var);
EXPORT double sbstropcorr(gtime_t time, const double *pos, const double *azel,
                          double *var);

/* options functions ---------------------------------------------------------*/
EXPORT opt_t *searchopt(const char *name, const opt_t *opts);
EXPORT int str2opt(opt_t *opt, const char *str);
EXPORT int opt2str(const opt_t *opt, char *str);
EXPORT int opt2buf(const opt_t *opt, char *buff);
EXPORT int loadopts(const char *file, opt_t *opts);
#ifdef MATLAB
EXPORT int loadmatopts(MATFile *pmat);
#endif
EXPORT int saveopts(const char *file, const char *mode, const char *comment,
                    const opt_t *opts);
EXPORT void resetsysopts(void);
EXPORT void getsysopts(prcopt_t *popt, solopt_t *sopt, filopt_t *fopt);
EXPORT void setsysopts(const prcopt_t *popt, const solopt_t *sopt,
                       const filopt_t *fopt);

/* stream data input and output functions ------------------------------------*/
EXPORT void strinitcom(void);
EXPORT void strinit  (stream_t *stream);
EXPORT void strlock  (stream_t *stream);
EXPORT void strunlock(stream_t *stream);
EXPORT int  stropen  (stream_t *stream, int type, int mode, const char *path);
EXPORT void strclose (stream_t *stream);
EXPORT int  strread  (stream_t *stream, unsigned char *buff, int n);
EXPORT int  strwrite (stream_t *stream, unsigned char *buff, int n);
EXPORT void strsync  (stream_t *stream1, stream_t *stream2);
EXPORT int  strstat  (stream_t *stream, char *msg);
EXPORT int  strstatx (stream_t *stream, char *msg);
EXPORT void strsum   (stream_t *stream, int *inb, int *inr, int *outb, int *outr);
EXPORT int  strgetsel(stream_t *stream, char *sel);
EXPORT int  strsetsel(stream_t *stream, const char *sel);
EXPORT int  strsetsrctbl(stream_t *stream, const char *file);
EXPORT void strsetopt(const int *opt);
EXPORT gtime_t strgettime(stream_t *stream);
EXPORT void strsendnmea(stream_t *stream, const sol_t *sol);
EXPORT void strsendcmd(stream_t *stream, const char *cmd);
EXPORT void strsettimeout(stream_t *stream, int toinact, int tirecon);
EXPORT void strsetdir(const char *dir);
EXPORT void strsetproxy(const char *addr);

/* integer ambiguity resolution ----------------------------------------------*/
EXPORT int lambda(int n, int m, const double *a, const double *Q, double *F,
                  double *s);
EXPORT int lambda_reduction(int n, const double *Q, double *Z);
EXPORT int lambda_search(int n, int m, const double *a, const double *Q,
                         double *F, double *s);


/* observation model */
EXPORT void matchcposb(int type,const obsd_t *obs,const nav_t *nav,int f,double *cbias,double *pbias);
EXPORT void getcorrobs(const prcopt_t *popt,const obsd_t *obs,const nav_t *nav,const int *frq_idxs,
                         const double *dantr,const double *dants, double phw, double *L, double *P,
                         double *Lc, double *Pc,double *freqs,double *dcbs,ssat_t *sat_info);

EXPORT int test_sys(int sys, int m);
EXPORT void ddcov(const int *nb, int n, const double *Ri, const double *Rj, int nv, double *R);
EXPORT void init_prires(const double *v,const int *vflag,int nv,res_t *res);
EXPORT void init_postres(rtk_t *rtk,const double *post_v,res_t *res,const double *R);
EXPORT void freeres(res_t *res);
EXPORT int resqc(gtime_t t,rtk_t *rtk,res_t *res,int *exc,int ppp);
EXPORT int pri_res_check(gtime_t t,rtk_t *rtk,const double *pri_v,const int *vflag,int nv,int *exc);
EXPORT int valins(const prcopt_t *opt,const double *x);
/* standard positioning ------------------------------------------------------*/
EXPORT int pntpos(int iep,const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel,
                  ssat_t *ssat, char *msg,const insopt_t *iopt,kf_t *ins_kf);

/* precise positioning -------------------------------------------------------*/
EXPORT void rtkinit(rtk_t *rtk, const prcopt_t *opt,const imup_t *imup);
EXPORT void rtkfree(rtk_t *rtk);
EXPORT int  rtkpos (rtk_t *rtk, obsd_t *obs, int nobs, const nav_t *nav);
EXPORT int iamb_ppk(const prcopt_t *opt,int sat,int f);
EXPORT int  rtkopenstat(const char *file, int level);
EXPORT void rtkclosestat(void);
EXPORT int  rtkoutstat(rtk_t *rtk, char *buff);
EXPORT int  rtkopenfcbstat(const char *file_wl,const char *file_nl,const char *file_lc);
EXPORT int resamb(rtk_t *rtk,double *bias,double *xa,double *Pa,int gps,int glo,int sbs,int ns);

/* precise point positioning -------------------------------------------------*/
EXPORT int scan_ppp(rtk_t *rtk,const prcopt_t *opt,obsd_t *obs,int n);
EXPORT void clk_repair(obsd_t* obs,int n,rtk_t *rtk);
EXPORT void bd2_multipath(rtk_t *rtk,obsd_t *obs,int n);
EXPORT int pppos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav);
EXPORT int pppnx(const prcopt_t *opt);
EXPORT int pppna(const prcopt_t *opt);
EXPORT int iamb_ppp(const prcopt_t *opt,int sat,int f);
EXPORT int pppoutstat(rtk_t *rtk, char *buff);
EXPORT int manage_pppar(rtk_t *rtk,double *bias,double *xa,double *Pa,int nf,const obsd_t *obs,int ns,const nav_t *nav,int *exc);
EXPORT void holdamb_ppp(rtk_t *rtk,const double *xa);

EXPORT int pppcorr_read(pppcorr_t *corr, const char *file);
EXPORT void pppcorr_free(pppcorr_t *corr);
EXPORT int pppcorr_trop(const pppcorr_t *corr, gtime_t time, const double *pos,
                        double *ztd, double *std);
EXPORT int pppcorr_stec(const pppcorr_t *corr, gtime_t time, const double *pos,
                        double *ion, double *std);

/*semi-tightly coupled*/
EXPORT int ig_stc(const rtk_t *rtk,const double *v,const int nv,const double *H,const int *vflag,const double *R,const int upd_type);

/* post-processing positioning -----------------------------------------------*/
EXPORT int postpos(gtime_t ts, gtime_t te, double ti, double tu,
                   const prcopt_t *popt, const solopt_t *sopt,
                   const filopt_t *fopt, char **infile, int n, char *outfile,
                   const char *rov, const char *base);
EXPORT int couplepos(prcopt_t *popt,const filopt_t *fopt,solopt_t *solopt,stream_t *moni);

/* stream server functions ---------------------------------------------------*/
EXPORT void strsvrinit (strsvr_t *svr, int nout);
EXPORT int  strsvrstart(strsvr_t *svr, int *opts, int *strs, char **paths,
                        char **logs, strconv_t **conv, char **cmds,
                        char **cmds_priodic, const double *nmeapos);
EXPORT void strsvrstop (strsvr_t *svr, char **cmds);
EXPORT void strsvrstop (strsvr_t *svr, char **cmds);
EXPORT void strsvrstat (strsvr_t *svr, int *stat, int *log_stat, int *byte,
                        int *bps, char *msg);
EXPORT strconv_t *strconvnew(int itype, int otype, const char *msgs, int staid,
                             int stasel, const char *opt);
EXPORT void strconvfree(strconv_t *conv);

/* rtk server functions ------------------------------------------------------*/
EXPORT int  rtksvrinit  (rtksvr_t *svr);
EXPORT void rtksvrfree  (rtksvr_t *svr);
EXPORT int  rtksvrstart (rtksvr_t *svr, int cycle, int buffsize, int *strs,
                         char **paths, int *formats, int navsel, char **cmds,
                         char **cmds_periodic, char **rcvopts, int nmeacycle,
                         int nmeareq, const double *nmeapos, prcopt_t *prcopt,
                         solopt_t *solopt, stream_t *moni, char *errmsg);
EXPORT void rtksvrstop  (rtksvr_t *svr, char **cmds);
EXPORT int  rtksvropenstr(rtksvr_t *svr, int index, int str, const char *path,
                          const solopt_t *solopt);
EXPORT void rtksvrclosestr(rtksvr_t *svr, int index);
EXPORT void rtksvrlock  (rtksvr_t *svr);
EXPORT void rtksvrunlock(rtksvr_t *svr);
EXPORT int  rtksvrostat (rtksvr_t *svr, int type, gtime_t *time, int *sat,
                         double *az, double *el, int **snr, int *vsat);
EXPORT void rtksvrsstat (rtksvr_t *svr, int *sstat, char *msg);
EXPORT int  rtksvrmark(rtksvr_t *svr, const char *name, const char *comment);

/* downloader functions ------------------------------------------------------*/
EXPORT int dl_readurls(const char *file, char **types, int ntype, url_t *urls,
                       int nmax);
EXPORT int dl_readstas(const char *file, char **stas, int nmax);
EXPORT int dl_exec(gtime_t ts, gtime_t te, double ti, int seqnos, int seqnoe,
                   const url_t *urls, int nurl, char **stas, int nsta,
                   const char *dir, const char *usr, const char *pwd,
                   const char *proxy, int opts, char *msg, FILE *fp);
EXPORT void dl_test(gtime_t ts, gtime_t te, double ti, const url_t *urls,
                    int nurl, char **stas, int nsta, const char *dir,
                    int ncol, int datefmt, FILE *fp);

/* gis data functions --------------------------------------------------------*/
EXPORT int gis_read(const char *file, gis_t *gis, int layer);
EXPORT void gis_free(gis_t *gis);

/* qzss lex functions --------------------------------------------------------*/
EXPORT int lexupdatecorr(const lexmsg_t *msg, nav_t *nav, gtime_t *tof);
EXPORT int lexreadmsg(const char *file, int sel, lex_t *lex);
EXPORT void lexoutmsg(FILE *fp, const lexmsg_t *msg);
EXPORT int lexconvbin(int type, int format, const char *infile,
                      const char *outfile);
EXPORT int lexeph2pos(gtime_t time, int sat, const nav_t *nav, double *rs,
                      double *dts, double *var);
EXPORT int lexioncorr(gtime_t time, const nav_t *nav, const double *pos,
                      const double *azel, double *delay, double *var);

/* application defined functions ---------------------------------------------*/
extern int showmsg(char *format,...);
extern void settspan(gtime_t ts, gtime_t te);
extern void settime(gtime_t time);

/* ins function ------------------------------------------------------------- */
#define EPS 1E-50
#define MAXLINELEN  512      /**< char numbner limit of line when reading file */
#define MAXIMUOBS   1000000  /**< max memory allocate for imu_t struct */
#define MAXODOBS    1000000  /**< max memory allocate for od_t sturct */
#define MAXPVAOBS   1000000  /**< max memory allocate for pva_t struct */
#define ARCDEG      D2R
#define ARCMIN      ARCDEG/60
#define ARCSEC      ARCMIN/60
#define DPH2RPS     4.84813681109536e-06    /**< deg/h to rad/s */
#define RPS2DPH     206264.806247096        /**< rad/s to deg/h */
#define DPSH2RPSS   2.90888208665722e-4     /**< deg/sqrt(h) to rad/sqrt(s) */
#define RPSS2DPSH   3437.74677078493        /**< rad/sqrt(s) to deg/sqrt(h) */
#define DPS2DPH     3600.0
#define DPHPSHZ2RPSS 4.84813681109536e-06   /**< deg/h/sqrt(Hz) to rad/sqrt(s) */
#define RPSS2DPHPSHZ 206264.806247096       /**< rad/sqrt(s) to deg/h/sqrt(Hz) */
#define G2MPS2      9.7803267714            /**< g0 to m/s^2 */
#define MPS22G      0.101971621297793       /**< m/s^2 to g0 */
#define MG2MPS2     9.7803267714e-3         /**< millo-g0(mg) to m/s^2 */
#define UG2MPS2     9.7803267714e-6
#define MPS22MG     101.971621297793        /**< m/s^2 to milli-g0(mg) */
#define MPS22UG     101971.621297793        /**< m/s^2 to micro-g0(ug) */
#define GAL2MPS2    0.01                    /**< gal to m/s^2 */
#define MPS22GAL    100.0                   /**< m/s^2 to gal */

#define INS_ALIGN_MANUAL        0
#define INS_ALIGN_CORSE         1
#define INS_ALIGN_GNSS_PV       2
#define INS_ALIGN_GNSS_SPP      3
#define INS_ALIGN_GNSS_DGPS     4
#define INS_ALIGN_GNSS_PPK      5
#define INS_ALIGN_GNSS_PPP      6

typedef struct {
    gtime_t time;
    int nx;
    double cCbe[9],pCbe[9],sCbe[9];
    double cre[3],pre[3],sre[3];
    double cve[3],pve[3],sve[3];
    double cba[3],pba[3],sba[3];
    double cbg[3],pbg[3],sbg[3];

    double *Pc,*Pp,*Ps;
    double *F;
}insstate_t;


/* match files  --------------------------------------------------------*/
EXPORT int match_navfile(const prcopt_t *popt,gtime_t ts,gtime_t te,const char prc_dir[],char *navpaths[]);
EXPORT int match_clkfile(const prcopt_t *popt,gtime_t ts,gtime_t te,const char prc_dir[],char *clkpaths[]);
EXPORT int match_sp3file(const prcopt_t* popt,gtime_t ts,gtime_t te,const char prc_dir[],char *sp3paths[]);
EXPORT int match_eopfile(gtime_t ts,gtime_t te,const char dir[],char eoppath[],const prcopt_t *popt);
EXPORT int match_dcbfile(const prcopt_t *popt,gtime_t ts,gtime_t te,const char prc_dir[],char dcbpath[]);
EXPORT int match_ionfile(gtime_t ts, gtime_t te, const char *dir, char *ionpath);
EXPORT int match_fcbfile(const prcopt_t *opt,gtime_t ts,gtime_t te,const char dir[],char fcbpath[]);
EXPORT int match_biafile(const prcopt_t *opt,gtime_t ts,gtime_t te,const char dir[],char biapath[]);
EXPORT int match_updfile(const prcopt_t *opt,gtime_t ts,gtime_t te,const char dir[],char *biapath[]);
EXPORT int match_baseofile(const prcopt_t *opt, gtime_t ts,const char *prcdir,char *basepath);
EXPORT int match_mgexdcbfile(gtime_t ts, gtime_t te, const char dir[], char dcbpath[],int opt);
EXPORT int match_reffile(const prcopt_t *opt,const char *obs_file,gtime_t ts,gtime_t te,const char dir[],char refath[],int dynamic );
EXPORT int match_blqfile(gtime_t ts, gtime_t te, const char dir[], char blqpath[]);
EXPORT int match_atxfile(gtime_t ts, gtime_t te, const char dir[], char atxPath[], const prcopt_t *popt);
EXPORT void matchout(const prcopt_t *popt,const char *prcdir,filopt_t *fopt,const solopt_t *sopt);

/* preprocess function  -------------------------------------------------*/
EXPORT int parsecmd(int arc,char *arv[],prcopt_t *popt,solopt_t *sopt, filopt_t *fopt,int *port);
EXPORT int loadconf(const char *conf_file, prcopt_t *popt, solopt_t *sopt, filopt_t *fopt);
EXPORT int loadprcfiles(const char *dir,const prcopt_t *popt,filopt_t *fopt,sta_t *sta,int* nsta);
EXPORT int freeprcfiles(const prcopt_t *popt, filopt_t *fopt);

EXPORT double sdobs(const obsd_t *obs, int i, int j, int f);
EXPORT int doppler(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, const prcopt_t *opt);
#ifdef __cplusplus
}
#endif
#endif /* RTKLIB_H */
