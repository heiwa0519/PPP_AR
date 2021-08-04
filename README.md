This is demo for multi-GNSS precise point positioning with ambiguity resolution (PPP-AR), which is based on RTKLIB and RTKLIB_demo5.

FEATURES

1. ppp-ar with different type products
2. others

ENVIRONMENT

WIN10 + CLION2019.3 + TDM(can be foun in ./ide folder)

DATA

Example data can be found in /PPP_AR/GNSS_DATA.7z, please unzip.

CMD

fix solution：  YOUR_PATH/PPP_AR/build/Bin/ppp_ar.exe -C YOUR_PATH/conf/PPP/ppp_mgex_wum.conf -S G -M PPP-KINE(or PPP-STATIC) -A 7 -L 0  

float solution: YOUR_PATH/PPP_AR/build/Bin/ppp_ar.exe -C YOUR_PATH/conf/PPP/ppp_mgex_wum.conf -S G -M PPP-KINE(or PPP-STATIC) -A 0 -L 0 

NOTE

Please set 'pos1-prcdir' in configuration file to your local path

Support by 'Assessment of GPS/Galileo/BDS precise point positioning with ambiguity resolution using products from different analysis centers, Remote Sensing (under review)'

SOMETHING SHOULD BE IMPROVED BY YOURSELF.