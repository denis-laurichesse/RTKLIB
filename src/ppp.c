/*------------------------------------------------------------------------------
* ppp.c : precise point positioning
*
*          Copyright (C) 2010-2020 by T.TAKASU, All rights reserved.
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
*           2020/11/30 1.14 use sat2freq() to get carrier frequency
*                           use E1-E5b for Galileo iono-free LC
*           2026/03/16 1.15 full multi-frequency support (receiver DCBs in state vector)
*                           correct bugs in PCOs and PCVs
*           2026/03/18 1.16 add Doppler observations in PPP with receiver
*                           clock-drift states
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define MAX_ITER    8               /* max number of iterations */
#define MAX_STD_FIX 0.15            /* max std-dev (3d) to fix solution */
#define MIN_NSAT_SOL 4              /* min satellite number for solution */
#define THRES_REJECT 4.0            /* reject threshold of posfit-res (sigma) */

#define THRES_MW_JUMP 10.0
#define MAXINNO_DOP 10000.0         /* reject threshold of Doppler prefit residual (m/s) */

#define VAR_POS     SQR(60.0)       /* init variance receiver position (m^2) */
#define VAR_VEL     SQR(10.0)       /* init variance of receiver vel ((m/s)^2) */
#define VAR_ACC     SQR(10.0)       /* init variance of receiver acc ((m/ss)^2) */
#define VAR_CLK     SQR(60.0)       /* init variance receiver clock (m^2) */
#define VAR_CLKD0   SQR(1000.0)     /* initial variance receiver clock drift ((m/s)^2) */
#define PRN_CLKD    10.0            /* receiver clock-drift process noise (m/s/sqrt(s)) */
#define VAR_ZTD     SQR( 0.6)       /* init variance ztd (m^2) */
#define VAR_GRA     SQR(0.01)       /* init variance gradient (m^2) */
#define VAR_DCB     SQR(30.0)       /* init variance dcb (m^2) */
#define VAR_BIAS    SQR(60.0)       /* init variance phase-bias (m^2) */
#define VAR_IONO    SQR(60.0)       /* init variance iono-delay */
#define VAR_GLO_IFB SQR( 0.6)       /* variance of glonass ifb */
#define DOP_NOISE   0.1             /* zenithal Doppler noise (m/s) */

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
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics?9:3)
#define NC(opt)     (NSYS)
#define NCD(opt)    ((opt)->dynamics?1:0)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))
#define NI(opt)     ((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)
#define ND(opt)     ((NF(opt)<=1)?0:(NC(opt)*(NF(opt)-1)))
#define NR(opt)     (NP(opt)+NC(opt)+NCD(opt)+NT(opt)+NI(opt)+ND(opt))
#define NB(opt)     (NF(opt)*MAXSAT)
#define NX(opt)     (NR(opt)+NB(opt))
#define IC(s,opt)   (NP(opt)+(s))
#define ICD(s,opt)  (NP(opt)+NC(opt))
#define IT(opt)     (NP(opt)+NC(opt)+NCD(opt))
#define II(s,opt)   (NP(opt)+NC(opt)+NCD(opt)+NT(opt)+(s)-1)
#define ID(s,f,opt) (NP(opt)+NC(opt)+NCD(opt)+NT(opt)+NI(opt)+((s)*((NF(opt)<=1)?0:(NF(opt)-1)))+((f)-1))
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)

#define AUX_FOR_IONO 1
#define AUX_FOR_BIAS 2
#define AUX_FOR_SLIP_GF 3
#define AUX_FOR_SLIP_MW 4

#define MEAS_PHASE 0
#define MEAS_CODE  1
#define MEAS_DOP   2

/* standard deviation of state -----------------------------------------------*/
static double STD(rtk_t *rtk, int i)
{
    if (rtk->sol.stat==SOLQ_FIX) return SQRT(rtk->Pa[i+i*rtk->nx]);
    return SQRT(rtk->P[i+i*rtk->nx]);
}
/* active constellation to idxx */
static int sys2clkidx(int sys)
{
    int idx=1; /* 0 = groupe de référence (GPS/QZS/SBS/...) */

#ifdef ENAGLO
    if (sys==SYS_GLO) return idx;
    idx++;
#endif
#ifdef ENAGAL
    if (sys==SYS_GAL) return idx;
    idx++;
#endif
#ifdef ENACMP
    if (sys==SYS_CMP) return idx;
    idx++;
#endif
#ifdef ENAIRN
    if (sys==SYS_IRN) return idx;
    idx++;
#endif
    return 0;
}
/* auxiliary frequency index by constellation -------------------------------*/
static int aux_freq_index(const prcopt_t *opt, int sys, const obsd_t *obs,
                          const nav_t *nav, int purpose)
{
    double freq;
    int f,fr=0;
    int nf=MIN(opt->nf,NFREQ);

    for (f=0;f<nf;f++) {
        if (f==fr) continue;
        if (!(freq=sat2freq(obs->sat,obs->code[f],nav))) continue;
        switch (purpose) {
            case AUX_FOR_IONO:
            case AUX_FOR_BIAS:
                if (obs->P[f]==0.0) continue;
                break;
            case AUX_FOR_SLIP_GF:
                if (obs->L[f]==0.0) continue;
                break;
            case AUX_FOR_SLIP_MW:
                if (obs->L[f]==0.0||obs->P[f]==0.0) continue;
                break;
            default:
                break;
        }
        return f;
    }
    return -1;
}
/* code dcb for constellation/frequency -------------------------------------*/
static double code_dcb(const double *x, int sysidx, int f, const prcopt_t *opt)
{
    if (!x||NF(opt)<=1||f==0) return 0.0;
    return x[ID(sysidx,f,opt)];
}
/* write solution status for PPP ---------------------------------------------*/
extern int pppoutstat(rtk_t *rtk, char *buff)
{
    ssat_t *ssat;
    double tow,pos[3],vel[3],acc[3],*x;
    int i,j,k,week;
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
    if (rtk->opt.ionoopt==IONOOPT_EST) {
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
    /* code DCB parameters */
    if (ND(&rtk->opt)>0) {
        for (i=0;i<NC(&rtk->opt);i++) for (j=1;j<NF(&rtk->opt);j++) {
            k=ID(i,j,&rtk->opt);
            if (x[k]==0.0) continue;
            p+=sprintf(p,"$DCB,%d,%.3f,%d,%d,%d,%.4f,%.4f\n",
                       week,tow,rtk->sol.stat,i+1,j+1,x[k],STD(rtk,k));
        }
    }
#ifdef OUTSTAT_AMB
    /* ambiguity parameters */
    for (i=0;i<MAXSAT;i++) for (j=0;j<NF(&rtk->opt);j++) {
        k=IB(i+1,j,&rtk->opt);
        if (rtk->x[k]==0.0) continue;
        satno2id(i+1,id);
        p+=sprintf(p,"$AMB,%d,%.3f,%d,%s,%d,%.4f,%.4f\n",week,tow,
                   rtk->sol.stat,id,j+1,x[k],STD(rtk,k));
    }
#endif
    return (int)(p-buff);
}
/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/
static void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
    double rsun[3],esun[3],r,ang,erpv[5]={0},cosa;
    int i,j;
    const char *type;
    
    trace(3,"testeclipse:\n");
    
    /* unit vector of sun direction (ecef) */
    sunmoonpos(gpst2utc(obs[0].time),erpv,rsun,NULL,NULL);
    normv3(rsun,esun);
    
    for (i=0;i<n;i++) {
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
extern int model_phw(gtime_t time, int sat, const char *type, int opt,
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
static double varerr(int sat, int sys, double el, int idx, int type,
                     const prcopt_t *opt)
{
    double fact=1.0,sinel=sin(el);

    (void)sat;
    if (type==1) fact*=opt->eratio[idx==0?0:1];
    fact*=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);

    if ((sys==SYS_GPS||sys==SYS_QZS)&&idx>=2) {
        fact*=EFACT_GPS_L5; /* apply extra factor to higher GPS/QZS freqs */
    }
    if (opt->ionoopt==IONOOPT_IFLC) fact*=3.0;
    return SQR(fact*opt->err[1])+SQR(fact*opt->err[2]/sinel);
}
/* Doppler measurement error variance ----------------------------------------*/
static double varerr_dop(int sat, int sys, double el, int idx,
                         const prcopt_t *opt)
{
    double fact=1.0,sinel=sin(el);

    (void)sat;
    (void)opt;
    if (sinel<0.1) sinel=0.1;

    fact*=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);

    if ((sys==SYS_GPS||sys==SYS_QZS)&&idx>=2) {
        fact*=EFACT_GPS_L5;
    }
    return SQR(fact*DOP_NOISE/sinel);
}
/* iono-free Doppler measurement error variance -------------------------------*/
static double varerr_dop_iflc(const obsd_t *obs, const nav_t *nav, int sat,
                              int sys, double el, const prcopt_t *opt)
{
    double freq1,freq2,C1,C2,var1,var2;

    if (obs->code[0]==CODE_NONE||obs->code[1]==CODE_NONE) return 0.0;
    if ((freq1=sat2freq(sat,obs->code[0],nav))==0.0||
        (freq2=sat2freq(sat,obs->code[1],nav))==0.0) return 0.0;

    C1=SQR(freq1)/(SQR(freq1)-SQR(freq2));
    C2=-SQR(freq2)/(SQR(freq1)-SQR(freq2));
    var1=varerr_dop(sat,sys,el,0,opt);
    var2=varerr_dop(sat,sys,el,1,opt);
    return SQR(C1)*var1+SQR(C2)*var2;
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
/* geometry-free phase measurement -------------------------------------------*/
static double gfmeas(const obsd_t *obs, const nav_t *nav, int fa, int fb)
{
    double freqa,freqb;

    freqa=sat2freq(obs->sat,obs->code[fa],nav);
    freqb=sat2freq(obs->sat,obs->code[fb],nav);
    if (freqa==0.0||freqb==0.0||obs->L[fa]==0.0||obs->L[fb]==0.0) return 0.0;
    return (obs->L[fa]/freqa-obs->L[fb]/freqb)*CLIGHT;
}
/* Melbourne-Wubbena linear combination --------------------------------------*/
static double mwmeas(const obsd_t *obs, const nav_t *nav, int fa, int fb)
{
    double freqa,freqb;

    freqa=sat2freq(obs->sat,obs->code[fa],nav);
    freqb=sat2freq(obs->sat,obs->code[fb],nav);

    if (freqa==0.0||freqb==0.0||obs->L[fa]==0.0||obs->L[fb]==0.0||
        obs->P[fa]==0.0||obs->P[fb]==0.0) return 0.0;
    return (obs->L[fa]-obs->L[fb])*CLIGHT/(freqa-freqb)-
           (freqa*obs->P[fa]+freqb*obs->P[fb])/(freqa+freqb);
}
/* helper for phase use */
static int use_phase_osb(const prcopt_t *opt, int sys)
{
    if (!opt) return 0;

    /* GLONASS controlled by pos2-gloarmode */
    if (sys==SYS_GLO) {
        return opt->glomodear != 0;
    }
    /* all other constellations controlled by pos2-armode */
    return opt->modear != ARMODE_OFF;
}
/* antenna corrected measurements --------------------------------------------*/
static void corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel,
                      const prcopt_t *opt, const double *dantr,
                      const double *dants, const int *antok,
                      double phw, double *L, double *P,
                      double *Lc, double *Pc)
{
    static const int fallback_code_mask  = SYS_GAL|SYS_CMP;
    static const int fallback_phase_mask = SYS_GPS|SYS_GLO|SYS_GAL|SYS_QZS|SYS_CMP|SYS_IRN;

    double freq[NFREQ]={0},C1,C2;
    double bcode,bphase,Lraw,Praw;
    uint8_t code_used_code,code_used_phase;
    int i,sys=satsys(obs->sat,NULL);
    int apply_phase_bias=use_phase_osb(opt,sys);

    for (i=0;i<NFREQ;i++) {
        L[i]=P[i]=0.0;
        freq[i]=sat2freq(obs->sat,obs->code[i],nav);

        if (antok && !antok[i]) continue;

        /* keep original RTKLIB behaviour:
         * a frequency is processed only if both phase and code exist */
        if (freq[i]==0.0||obs->L[i]==0.0||obs->P[i]==0.0) continue;
        if (testsnr(0,0,azel[1],obs->SNR[i]*SNR_UNIT,&opt->snrmask)) continue;

        Lraw=obs->L[i];
        Praw=obs->P[i];
        bcode=0.0;
        bphase=0.0;
        code_used_code = obs->code[i];
        code_used_phase = obs->code[i];

        if (nav->osb) {
            /* code bias: always applied if available */
            if (!getsinexbias(nav->osb, obs->time, obs->sat, obs->code[i],
                              BIAS_CODE, nav, fallback_code_mask,
                              &bcode, NULL, &code_used_code)) {
                trace(3,"corr_meas: no SINEX code bias sat=%d code=%d\n",
                      obs->sat,(int)obs->code[i]);
                continue;
            }

           Praw -= bcode; /* code bias already in meters */

            if (code_used_code != obs->code[i]) {
                trace(3,"corr_meas: code fallback sat=%d code=%d->%d bias=%.4f m\n",
                      obs->sat,(int)obs->code[i],(int)code_used_code,bcode);
            }

            /* phase bias: apply only if AR is enabled for this constellation */
            if (apply_phase_bias) {
                if (!getsinexbias(nav->osb, obs->time, obs->sat, obs->code[i],
                                  BIAS_PHASE, nav, fallback_phase_mask,
                                  &bphase, NULL, &code_used_phase)) {
                    trace(3,"corr_meas: no SINEX phase bias sat=%d code=%d\n",
                          obs->sat,(int)obs->code[i]);
                    continue;
                }

                Lraw -= bphase * freq[i] / CLIGHT; /* phase bias m -> cycles */

                if (code_used_phase != obs->code[i]) {
                    trace(3,"corr_meas: phase fallback sat=%d code=%d->%d bias=%.4f m\n",
                          obs->sat,(int)obs->code[i],(int)code_used_phase,bphase);
                }
            }
            else {
                trace(3,"corr_meas: phase OSB disabled sat=%d sys=%d code=%d modear=%d glomodear=%d\n",
                      obs->sat,sys,(int)obs->code[i],opt->modear,opt->glomodear);
            }
        }
        else {
            /* legacy RTKLIB behaviour if no SINEX OSB loaded */
            if (sys==SYS_GPS||sys==SYS_GLO) {
                if (obs->code[i]==CODE_L1C) Praw += nav->cbias[obs->sat-1][1];
                if (obs->code[i]==CODE_L2C) Praw += nav->cbias[obs->sat-1][2];
            }
        }

        /* then continue with antenna / wind-up corrections */
        L[i] = Lraw*CLIGHT/freq[i] - dants[i] - dantr[i] - phw*CLIGHT/freq[i];
        P[i] = Praw                - dants[i] - dantr[i];
    }

    /* iono-free LC */
    *Lc=*Pc=0.0;
    if (freq[0]==0.0||freq[1]==0.0) return;
    C1= SQR(freq[0])/(SQR(freq[0])-SQR(freq[1]));
    C2=-SQR(freq[1])/(SQR(freq[0])-SQR(freq[1]));

    if (L[0]!=0.0&&L[1]!=0.0) *Lc=C1*L[0]+C2*L[1];
    if (P[0]!=0.0&&P[1]!=0.0) *Pc=C1*P[0]+C2*P[1];
}
/* corrected Doppler measurements (m/s) --------------------------------------*/
static void corr_dop_meas(const obsd_t *obs, const nav_t *nav, const double *azel,
                          const prcopt_t *opt, double *D, double *Dc)
{
    double freq[NFREQ]={0},C1,C2;
    int i;

    for (i=0;i<NFREQ;i++) {
        D[i]=0.0;
        if (obs->code[i]==CODE_NONE||obs->D[i]==0.0) continue;
        if ((freq[i]=sat2freq(obs->sat,obs->code[i],nav))==0.0) continue;
        if (testsnr(0,0,azel[1],obs->SNR[i]*SNR_UNIT,&opt->snrmask)) continue;

        /* Doppler in RTKLIB observations is in Hz: convert to range-rate (m/s) */
        D[i]=-CLIGHT/freq[i]*obs->D[i];
    }
    *Dc=0.0;
    if (freq[0]==0.0||freq[1]==0.0) return;
    C1= SQR(freq[0])/(SQR(freq[0])-SQR(freq[1]));
    C2=-SQR(freq[1])/(SQR(freq[0])-SQR(freq[1]));
    if (D[0]!=0.0&&D[1]!=0.0) *Dc=C1*D[0]+C2*D[1];
}
/* satellite clock drift for Doppler model (s/s) -----------------------------*/
static int satclkdrift_dop(const obsd_t *obs, const nav_t *nav, double *dclk)
{
    double rsb[6],dtsb[2],varb;
    int svhb;

    satposs(obs->time,obs,1,nav,EPHOPT_BRDC,rsb,dtsb,&varb,&svhb);
    if (dtsb[1]==0.0) return 0;
    *dclk=dtsb[1];
    return 1;
}
/* detect cycle slip by LLI --------------------------------------------------*/
static void detslp_ll(rtk_t *rtk, const obsd_t *obs, int n)
{
    int i,j;
    
    trace(3,"detslp_ll: n=%d\n",n);
    
    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<rtk->opt.nf;j++) {
        if (obs[i].L[j]==0.0||!(obs[i].LLI[j]&3)) continue;
        
        trace(3,"detslp_ll: slip detected sat=%2d f=%d\n",obs[i].sat,j+1);
        
        rtk->ssat[obs[i].sat-1].slip[j]=1;
    }
}
/* detect cycle slip by geometry free phase jump -----------------------------*/
static void detslp_gf(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    double g0,g1;
    int i,j,f,sat,sys,fa;

    trace(3,"detslp_gf: n=%d\n",n);

    for (i=0;i<n&&i<MAXOBS;i++) {
        sat=obs[i].sat;
        if (!(sys=satsys(sat,NULL))) continue;
        fa=0;

        for (f=1; f<NF(&rtk->opt)&&f<NFREQ; f++) {
            if ((g1=gfmeas(obs+i,nav,fa,f))==0.0) {
                rtk->ssat[sat-1].gf[f-1]=0.0;
                continue;
            }
            g0=rtk->ssat[sat-1].gf[f-1];
            rtk->ssat[sat-1].gf[f-1]=g1;

            trace(4,"detslip_gf: sat=%2d f=%d gf0=%8.3f gf1=%8.3f\n",sat,f+1,g0,g1);

            if (g0!=0.0&&fabs(g1-g0)>rtk->opt.thresslip) {
                trace(3,"detslip_gf: slip detected sat=%2d f=%d gf=%8.3f->%8.3f\n",
                      sat,f+1,g0,g1);

                for (j=0;j<rtk->opt.nf;j++) rtk->ssat[sat-1].slip[j]|=1;
            }
        }
    }
}
/* detect slip by Melbourne-Wubbena linear combination jump ------------------*/
static void detslp_mw(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    double w0,w1;
    int i,j,f,sat,sys,fa;

    trace(3,"detslp_mw: n=%d\n",n);

    for (i=0;i<n&&i<MAXOBS;i++) {
        sat=obs[i].sat;
        if (!(sys=satsys(sat,NULL))) continue;
        fa=0;

        for (f=1; f<NF(&rtk->opt)&&f<NFREQ; f++) {
            if ((w1=mwmeas(obs+i,nav,fa,f))==0.0) {
                rtk->ssat[sat-1].mw[f-1]=0.0;
                continue;
            }
            w0=rtk->ssat[sat-1].mw[f-1];
            rtk->ssat[sat-1].mw[f-1]=w1;

            trace(4,"detslip_mw: sat=%2d f=%d mw0=%8.3f mw1=%8.3f\n",sat,f+1,w0,w1);

            if (w0!=0.0&&fabs(w1-w0)>THRES_MW_JUMP) {
                trace(3,"detslip_mw: slip detected sat=%2d f=%d mw=%8.3f->%8.3f\n",
                      sat,f+1,w0,w1);

                for (j=0;j<rtk->opt.nf;j++) rtk->ssat[sat-1].slip[j]|=1;
            }
        }
    }
}
/* satellite antenna phase center variation ----------------------------------*/
extern void satantpcv(const double *rs, const double *rr, const pcv_t *pcv,
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
/* satellite offset and obs code to antenna offset (with fallbacks) */
static double sat2ant(int sat, uint8_t code, double *dant)
{
int sys = sat?satsys(sat, NULL):SYS_GPS;
double x = 0.0;
int idx = code2idx(sys, code);

if (idx != -1) x = dant[idx];

if (sat) {
 sys = satsys(sat, NULL);
 if (((sys == SYS_GPS) || (sys == SYS_GLO)) && !x && (NFREQ >= 2)) x = dant[code2idx(sys, CODE_L2P)]; 
} else {
 if (!x && (NFREQ >= 2)) x = dant[code2idx(sys, CODE_L2P)];
}
return x;
}
/* satellite antenna phase center offset in ECEF for all frequencies ---------*/
static void satantofffreq(gtime_t time, const double *rs, int sat,
                          const nav_t *nav, double *dantx, double *danty,
                          double *dantz)
{
    const pcv_t *pcv=nav->pcvs+sat-1;
    double ex[3],ey[3],ez[3],es[3],r[3],rsun[3],gmst,erpv[5]={0};
    int i,j;

    for (j=0;j<NFREQ;j++) {
        dantx[j]=0.0;
        danty[j]=0.0;
        dantz[j]=0.0;
    }
    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,&gmst);

    for (i=0;i<3;i++) r[i]=-rs[i];
    if (!normv3(r,ez)) return;

    for (i=0;i<3;i++) r[i]=rsun[i]-rs[i];
    if (!normv3(r,es)) return;

    cross3(ez,es,r);
    if (!normv3(r,ey)) return;

    cross3(ey,ez,ex);

    for (j=0;j<NFREQ;j++) {
        dantx[j]=pcv->off[j][0]*ex[0]+pcv->off[j][1]*ey[0]+pcv->off[j][2]*ez[0];
        danty[j]=pcv->off[j][0]*ex[1]+pcv->off[j][1]*ey[1]+pcv->off[j][2]*ez[1];
        dantz[j]=pcv->off[j][0]*ex[2]+pcv->off[j][1]*ey[2]+pcv->off[j][2]*ez[2];
    }
}
/* build receiver/satellite antenna corrections by observation code ----------*/
static void model_antcorr_ppp(const obsd_t *obs, const nav_t *nav,
                              const prcopt_t *opt, const double *rs,
                              const double *rr, const double *e,
                              const double *azel, int sat,
                              double *dantr, double *dants, int *antok)
{
    double dantr0[NFREQ]={0},dants0[NFREQ]={0};
    double dantx[NFREQ]={0},danty[NFREQ]={0},dantz[NFREQ]={0};
    double dantn[NFREQ]={0};
    double dsx,dsy,dsz;
    int k,sys=satsys(sat,NULL);

    for (k=0;k<NFREQ;k++) {
        dantr[k]=0.0;
        dants[k]=0.0;
        antok[k]=0;
    }

    /* stock receiver model: off + var */
    antmodel(opt->pcvr,opt->antdel[0],azel,opt->posopt[1],dantr0);

    /* receiver remapping */
    for (k=2;k<NFREQ;k++) dantr0[k]=dantr0[1];
    if (sys==SYS_CMP&&NFREQ>=2) dantr0[1]=dantr0[0];

    /* satellite PCV + explicit satellite PCO */
    if (opt->posopt[0]) {
        satantpcv(rs,rr,nav->pcvs+sat-1,dants0);
        satantofffreq(obs->time,rs,sat,nav,dantx,danty,dantz);

        for (k=0;k<NFREQ;k++) {
            dantn[k]=sqrt(SQR(dantx[k])+SQR(danty[k])+SQR(dantz[k]));
        }
    }

    for (k=0;k<NFREQ;k++) {
        if (obs->code[k]==CODE_NONE) continue;

        /* receiver correction selected by actual signal code */
        dantr[k]=sat2ant(0,obs->code[k],dantr0);

        if (opt->posopt[0]) {
            /* reject this frequency if satellite PCO for this signal is not found */
            if (sat2ant(sat,obs->code[k],dantn)==0.0) {
                trace(3,
                      "model_antcorr_ppp: missing satellite PCO sat=%2d code=%d f=%d\n",
                      sat,(int)obs->code[k],k+1);
                continue;
            }

            /* satellite correction = stock PCV + projected PCO */
            dsx=sat2ant(sat,obs->code[k],dantx);
            dsy=sat2ant(sat,obs->code[k],danty);
            dsz=sat2ant(sat,obs->code[k],dantz);

            dants[k]=sat2ant(sat,obs->code[k],dants0)
                    +(e[0]*dsx+e[1]*dsy+e[2]*dsz);
        }
        else {
            dants[k]=0.0;
        }

        antok[k]=1;
    }
}
/* temporal update of position -----------------------------------------------*/
static void udpos_ppp(rtk_t *rtk)
{
    double *F,*P,*FP,*x,*xp,pos[3],Q[9]={0},Qv[9];
    int i,j,*ix,nx;
    
    trace(3,"udpos_ppp:\n");
    
    /* fixed mode */
    if (rtk->opt.mode==PMODE_PPP_FIXED) {
        for (i=0;i<3;i++) initx(rtk,rtk->opt.ru[i],1E-8,i);
        return;
    }
    /* initialize position for first epoch */
    if (norm(rtk->x,3)<=0.0) {
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
    /* kinematic mode without dynamics */
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
static void udclk_ppp(rtk_t *rtk)
{
    double dtr;
    int i;
    
    trace(3,"udclk_ppp:\n");
    
    /* initialize every epoch for clock (white noise) */
    for (i=0;i<NSYS;i++) {
        if (rtk->opt.sateph==EPHOPT_PREC) {
            /* time of prec ephemeris is based gpst */
            /* negelect receiver inter-system bias  */
            dtr=rtk->sol.dtr[0];
        }
        else {
            dtr=i==0?rtk->sol.dtr[0]:rtk->sol.dtr[0]+rtk->sol.dtr[i];
        }
        initx(rtk,CLIGHT*dtr,VAR_CLK,IC(i,&rtk->opt));
    }
}
/* temporal update of clock-drift --------------------------------------------*/
static void udclkd_ppp(rtk_t *rtk)
{
    int j;

    trace(3,"udclkd_ppp:\n");

    if (NCD(&rtk->opt)<=0) return;

    j=ICD(0,&rtk->opt);

    if (rtk->P[j+j*rtk->nx]<=0.0) {
        initx(rtk,1E-6,VAR_CLKD0,j);
    }
    else {
        rtk->P[j+j*rtk->nx]+=SQR(PRN_CLKD)*fabs(rtk->tt);
    }
}
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop_ppp(rtk_t *rtk)
{
    double pos[3],azel[]={0.0,PI/2.0},ztd,var;
    int i=IT(&rtk->opt),j;
    
    trace(3,"udtrop_ppp:\n");
    
    if (rtk->x[i]==0.0) {
        ecef2pos(rtk->sol.rr,pos);
        ztd=sbstropcorr(rtk->sol.time,pos,azel,&var);
        initx(rtk,ztd,var,i);
        
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
    double freqr,freqa,ion,sinel,pos[3],*azel;
    double L[NFREQ],P[NFREQ],Lc,Pc,dantr[NFREQ]={0},dants[NFREQ]={0};
    char *p;
    int i,j,gap_resion=GAP_RESION,sys,sysi,fr,fa;

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
        if (!(sys=satsys(obs[i].sat,NULL))) continue;
        fr=0;
        fa=aux_freq_index(&rtk->opt,sys,obs+i,nav,AUX_FOR_IONO);
        if (fa<0) continue;
        j=II(obs[i].sat,&rtk->opt);
        if (rtk->x[j]==0.0) {
            corr_meas(obs+i,nav,rtk->ssat[obs[i].sat-1].azel,&rtk->opt,dantr,dants,NULL,
                0.0,L,P,&Lc,&Pc);
            freqr=sat2freq(obs[i].sat,obs[i].code[fr],nav);
            freqa=sat2freq(obs[i].sat,obs[i].code[fa],nav);
            if (P[fr]==0.0||P[fa]==0.0||freqr==0.0||freqa==0.0||
                fabs(1.0-SQR(freqr/freqa))<1E-12) {
                continue;
            }
            sysi=sys2clkidx(sys);
            ion=((P[fr]-code_dcb(rtk->x,sysi,fr,&rtk->opt))-
                 (P[fa]-code_dcb(rtk->x,sysi,fa,&rtk->opt)))/(1.0-SQR(freqr/freqa));
            azel=rtk->ssat[obs[i].sat-1].azel;
            ion/=ionmapf(pos,azel);
            initx(rtk,ion,VAR_IONO,j);
        }
        else {
            sinel=sin(MAX(rtk->ssat[obs[i].sat-1].azel[1],5.0*D2R));
            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[1]/sinel)*fabs(rtk->tt);
        }
    }
}
/* temporal update of L5-receiver-dcb parameters -----------------------------*/
static void uddcb_ppp(rtk_t *rtk)
{
    int i,sys,f;

    trace(3,"uddcb_ppp:\n");

    for (sys=0;sys<NC(&rtk->opt);sys++) for (f=1;f<NF(&rtk->opt);f++) {
        i=ID(sys,f,&rtk->opt);
        if (rtk->x[i]==0.0) {
            initx(rtk,1E-6,VAR_DCB,i);
        }
    }
}
/* temporal update of phase biases -------------------------------------------*/
static void udbias_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                       const double *rs, const double *dr)
{
    double L[NFREQ],P[NFREQ],Lc,Pc,bias[MAXOBS],offset=0.0;
    double freqr,freqf,freqa,ion,dantr[NFREQ],dants[NFREQ];
    double lam[MAXOBS]={0.0},rr[3],pos[3],e[3],azel[2],r,phw;
    int antok[NFREQ];
    int i,j,k,f,sat,sys,sysi,fr,fa,slip[MAXOBS]={0},clk_jump=0;

    trace(3,"udbias  : n=%d\n",n);

    /* handle day-boundary clock jump */
    if (rtk->opt.posopt[5]) {
        clk_jump=ROUND(time2gpst(obs[0].time,NULL)*10)%864000==0;
    }

    for (i=0;i<MAXSAT;i++) for (j=0;j<rtk->opt.nf;j++) {
        rtk->ssat[i].slip[j]=0;
    }

    /* detect cycle slip by LLI / GF / MW */
    detslp_ll(rtk,obs,n);
    detslp_gf(rtk,obs,n,nav);
    detslp_mw(rtk,obs,n,nav);

    for (i=0;i<3;i++) rr[i]=rtk->x[i]+dr[i];
    ecef2pos(rr,pos);

    for (f=0;f<NF(&rtk->opt);f++) {
        offset=0.0;

        /* reset phase-bias if expire obs outage counter */
        for (i=0;i<MAXSAT;i++) {
            if (++rtk->ssat[i].outc[f]>(uint32_t)rtk->opt.maxout||
                rtk->opt.modear==ARMODE_INST||clk_jump) {
                initx(rtk,0.0,0.0,IB(i+1,f,&rtk->opt));
            }
        }

        for (i=k=0;i<n&&i<MAXOBS;i++) {
            sat=obs[i].sat;
            bias[i]=0.0;
            slip[i]=0;
            lam[i]=0.0;

            if (!(sys=satsys(sat,NULL))) continue;
            sysi=sys2clkidx(sys);
            fr=0;
            j=IB(sat,f,&rtk->opt);

            /* current geometry, same spirit as ppp_res() */
            if ((r=geodist(rs+i*6,rr,e))<=0.0||
                satazel(pos,e,azel)<rtk->opt.elmin) {
                continue;
            }

            model_antcorr_ppp(obs+i,nav,&rtk->opt,rs+i*6,rr,e,azel,sat,
                              dantr,dants,antok);

            phw=rtk->ssat[sat-1].phw;
            if (!model_phw(rtk->sol.time,sat,nav->pcvs[sat-1].type,
                           rtk->opt.posopt[2]?2:0,rs+i*6,rr,&phw)) {
                trace(3,"udbias_ppp: skip sat=%2d f=%d (model_phw failed)\n",
                      sat,f+1);
                continue;
            }

            corr_meas(obs+i,nav,azel,&rtk->opt,dantr,dants,antok,phw,
                      L,P,&Lc,&Pc);

            if (rtk->opt.ionoopt==IONOOPT_IFLC) {
                bias[i]=Lc-Pc;
                for (j=0;j<rtk->opt.nf;j++) slip[i]|=rtk->ssat[sat-1].slip[j];

                if (bias[i]==0.0&&slip[i]) {
                    trace(3,"udbias_ppp: slip sat=%2d f=%d but IFLC unusable after corr_meas\n",
                          sat,f+1);
                }
            }
            else if (L[f]!=0.0&&P[f]!=0.0) {
                slip[i]=rtk->ssat[sat-1].slip[f];

                freqr=sat2freq(sat,obs[i].code[fr],nav);
                freqf=sat2freq(sat,obs[i].code[f],nav);
                if (freqr==0.0||freqf==0.0) continue;

                lam[i]=CLIGHT/freqf; /* lambda of IB(sat,f) */

                if (f!=fr) {
                    if (P[fr]==0.0||fabs(1.0-SQR(freqr/freqf))<1E-12) continue;
                    ion=((P[fr]-code_dcb(rtk->x,sysi,fr,&rtk->opt))-
                         (P[f]-code_dcb(rtk->x,sysi,f,&rtk->opt)))/
                        (1.0-SQR(freqr/freqf));
                }
                else {
                    if (rtk->opt.ionoopt==IONOOPT_EST&&rtk->x[II(sat,&rtk->opt)]!=0.0) {
                        ion=rtk->x[II(sat,&rtk->opt)]*ionmapf(pos,azel);
                    }
                    else {
                        fa=aux_freq_index(&rtk->opt,sys,obs+i,nav,AUX_FOR_BIAS);
                        if (fa<0) continue;
                        freqa=sat2freq(sat,obs[i].code[fa],nav);
                        if (P[fr]==0.0||P[fa]==0.0||freqa==0.0||
                            fabs(1.0-SQR(freqr/freqa))<1E-12) continue;
                        ion=((P[fr]-code_dcb(rtk->x,sysi,fr,&rtk->opt))-
                             (P[fa]-code_dcb(rtk->x,sysi,fa,&rtk->opt)))/
                            (1.0-SQR(freqr/freqa));
                    }
                }
                bias[i]=L[f]-(P[f]-code_dcb(rtk->x,sysi,f,&rtk->opt))+
                        2.0*ion*SQR(freqr/freqf);
            }
            else {
                if (rtk->ssat[sat-1].slip[f]) {
                    trace(3,"udbias_ppp: slip sat=%2d f=%d but phase/code rejected after corr_meas antok=%d\n",
                          sat,f+1,antok[f]);
                }
            }

            if (rtk->x[j]==0.0||slip[i]||bias[i]==0.0) continue;

            if (rtk->opt.ionoopt==IONOOPT_IFLC) {
                offset+=bias[i]-rtk->x[j];
            }
            else {
                if (lam[i]<=0.0) continue;
                offset+=bias[i]-lam[i]*rtk->x[j];
            }
            k++;
        }

        /* correct phase-code jump to ensure phase-code coherency */
        if (k>=2&&fabs(offset/k)>0.0005*CLIGHT) {
            for (i=0;i<n&&i<MAXOBS;i++) {
                sat=obs[i].sat;
                j=IB(sat,f,&rtk->opt);
                if (rtk->x[j]==0.0) continue;

                if (rtk->opt.ionoopt==IONOOPT_IFLC) {
                    rtk->x[j]+=offset/k; /* meters */
                }
                else {
                    double freq=sat2freq(sat,obs[i].code[f],nav);
                    double lami=freq>0.0?CLIGHT/freq:0.0;
                    if (lami<=0.0) continue;
                    rtk->x[j]+=offset/k/lami; /* meters -> cycles */
                }
            }
            trace(2,"phase-code jump corrected: %s n=%2d dt=%12.9fs\n",
                  time_str(rtk->sol.time,0),k,offset/k/CLIGHT);
        }

        for (i=0;i<n&&i<MAXOBS;i++) {
            sat=obs[i].sat;
            j=IB(sat,f,&rtk->opt);

            if (rtk->opt.ionoopt==IONOOPT_IFLC) {
                rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[0])*fabs(rtk->tt);

                if (bias[i]==0.0||(rtk->x[j]!=0.0&&!slip[i])) continue;

                initx(rtk,bias[i],VAR_BIAS,IB(sat,f,&rtk->opt));
            }
            else {
                double freq=sat2freq(sat,obs[i].code[f],nav);
                double lami=freq>0.0?CLIGHT/freq:0.0;
                if (lami<=0.0) continue;

                rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[0]/lami)*fabs(rtk->tt);

                if (bias[i]==0.0||(rtk->x[j]!=0.0&&!slip[i])) continue;

                initx(rtk,bias[i]/lami,SQR(SQRT(VAR_BIAS)/lami),IB(sat,f,&rtk->opt));
            }

            for (k=0;k<MAXSAT;k++) rtk->ambc[sat-1].flags[k]=0;

            trace(5,"udbias_ppp: sat=%2d f=%d bias=%.3f%s slip=%d\n",sat,f+1,
                  rtk->opt.ionoopt==IONOOPT_IFLC?bias[i]:
                  (lam[i]>0.0?bias[i]/lam[i]:0.0),
                  rtk->opt.ionoopt==IONOOPT_IFLC?" m":" cyc",slip[i]);
        }
    }
}/* temporal update of states --------------------------------------------------*/
static void udstate_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    trace(3,"udstate_ppp: n=%d\n",n);
    
    /* temporal update of position */
    udpos_ppp(rtk);
    
    /* temporal update of clock */
    udclk_ppp(rtk);

    /* temporal update of clock drift */
    udclkd_ppp(rtk);
    
    /* temporal update of tropospheric parameters */
    if (rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG) {
        udtrop_ppp(rtk);
    }
    /* temporal update of code-dcb parameters */
    if (ND(&rtk->opt)>0) {
        uddcb_ppp(rtk);
    }
    /* temporal update of ionospheric parameters */
    if (rtk->opt.ionoopt==IONOOPT_EST) {
        udiono_ppp(rtk,obs,n,nav);
    }
}
/* precise tropospheric model ------------------------------------------------*/
static double trop_model_prec(gtime_t time, const double *pos,
                              const double *azel, const double *x, double *dtdx,
                              double *var)
{
    const double zazel[]={0.0,PI/2.0};
    double zhd,m_h,m_w,cotz,grad_n,grad_e;
    
    /* zenith hydrostatic delay */
    zhd=tropmodel(time,pos,zazel,0.0);
    
    /* mapping function */
    m_h=tropmapf(time,pos,azel,&m_w);
    
    if (azel[1]>0.0) {
        
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
    return m_h*zhd+m_w*(x[0]-zhd);
}
/* tropospheric model ---------------------------------------------------------*/
static int model_trop(gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, const double *x, double *dtdx,
                      const nav_t *nav, double *dtrp, double *var)
{
    double trp[3]={0};
    
    if (opt->tropopt==TROPOPT_SAAS) {
        *dtrp=tropmodel(time,pos,azel,REL_HUMI);
        *var=SQR(ERR_SAAS);
        return 1;
    }
    if (opt->tropopt==TROPOPT_SBAS) {
        *dtrp=sbstropcorr(time,pos,azel,var);
        return 1;
    }
    if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
        matcpy(trp,x+IT(opt),opt->tropopt==TROPOPT_EST?1:3,1);
        *dtrp=trop_model_prec(time,pos,azel,trp,dtdx,var);
        return 1;
    }
    return 0;
}
/* ionospheric model ---------------------------------------------------------*/
static int model_iono(gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, int sat, const double *x,
                      const nav_t *nav, double *dion, double *var)
{
    if (opt->ionoopt==IONOOPT_SBAS) {
        return sbsioncorr(time,nav,pos,azel,dion,var);
    }
    if (opt->ionoopt==IONOOPT_TEC) {
        return iontec(time,nav,pos,azel,1,dion,var);
    }
    if (opt->ionoopt==IONOOPT_BRDC) {
        *dion=ionmodel(time,nav->ion_gps,pos,azel);
        *var=SQR(*dion*ERR_BRDCI);
        return 1;
    }
    if (opt->ionoopt==IONOOPT_EST) {
        *dion=x[II(sat,opt)];
        *var=0.0;
        return 1;
    }
    if (opt->ionoopt==IONOOPT_IFLC) {
        *dion=*var=0.0;
        return 1;
    }
    return 0;
}
/* phase, code and Doppler residuals -----------------------------------------*/
static int ppp_res(int post, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *var_rs, const int *svh,
                   const double *dr, int *exc, const nav_t *nav,
                   const double *x, rtk_t *rtk, double *v, double *H, double *R,
                   double *azel)
{
    prcopt_t *opt=&rtk->opt;
    double y,r,cdtr,cdtrd,bias,C=0.0,rr[3],pos[3],e[3],dtdx[3];
    double L[NFREQ],P[NFREQ],D[NFREQ],Lc,Pc,Dc;
    double var[MAXOBS*3*NFREQ],dtrp=0.0,dion=0.0,vart=0.0,vari=0.0,dcb,freq;
    double dantr[NFREQ]={0},dants[NFREQ]={0};
    int antok[NFREQ]={0};
    double ve[MAXOBS*3*NFREQ]={0},vmax=0,dtsd,modeled;
    double ep[6],sod;
    const char *ocode;
    char str[32],satid[16],obsid[16];
    int ne=0,obsi[MAXOBS*3*NFREQ]={0},frqi[MAXOBS*3*NFREQ],typei[MAXOBS*3*NFREQ],maxobs,maxfrq,maxtype,rej;
    int i,j,k,f,sat,sys,sysi,nv=0,nx=rtk->nx,stat=1,ntype=2,type;
    char mlabel[16];

    time2str(obs[0].time,str,2);

    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].vsat[j]=0;

    for (i=0;i<3;i++) rr[i]=x[i]+dr[i];
    ecef2pos(rr,pos);

    for (i=0;i<n&&i<MAXOBS;i++) {
        sat=obs[i].sat;

        if ((r=geodist(rs+i*6,rr,e))<=0.0||
            satazel(pos,e,azel+i*2)<opt->elmin) {
            exc[i]=1;
            continue;
        }
        if (!(sys=satsys(sat,NULL))||!rtk->ssat[sat-1].vs||
            satexclude(obs[i].sat,var_rs[i],svh[i],opt)||exc[i]) {
            exc[i]=1;
            continue;
        }
        /* tropospheric and ionospheric model */
        if (!model_trop(obs[i].time,pos,azel+i*2,opt,x,dtdx,nav,&dtrp,&vart)||
            !model_iono(obs[i].time,pos,azel+i*2,opt,sat,x,nav,&dion,&vari)) {
            continue;
        }

        /* satellite and receiver antenna model */
        model_antcorr_ppp(obs+i,nav,opt,rs+i*6,rr,e,azel+i*2,sat,dantr,dants,antok);

        /* phase windup model */
        if (!model_phw(rtk->sol.time,sat,nav->pcvs[sat-1].type,
                       opt->posopt[2]?2:0,rs+i*6,rr,&rtk->ssat[sat-1].phw)) {
            continue;
        }
        /* corrected phase and code measurements */
        corr_meas(obs+i,nav,azel+i*2,&rtk->opt,dantr,dants,antok,
          rtk->ssat[sat-1].phw,L,P,&Lc,&Pc);

        /* corrected Doppler measurements (m/s) */
        corr_dop_meas(obs+i,nav,azel+i*2,&rtk->opt,D,&Dc);

        sysi=sys2clkidx(sys);

        /* stack phase/code residuals {L1,P1,L2,P2,...} */
        for (j=0;j<ntype*NF(opt);j++) {
            if (nv>=MAXOBS*3*NFREQ) {
                trace(2,"ppp_res: too many residuals nv=%d\n",nv);
                return 0;
            }
            f=j/ntype;
            type=j%ntype;
            dcb=bias=0.0;
            C=0.0;

            if (opt->ionoopt==IONOOPT_IFLC) {
                if ((y=type==MEAS_PHASE?Lc:Pc)==0.0) continue;
            }
            else {
                if ((y=type==MEAS_PHASE?L[f]:P[f])==0.0) continue;

                if ((freq=sat2freq(sat,obs[i].code[f],nav))==0.0) continue;
                C=SQR(FREQ1/freq)*ionmapf(pos,azel+i*2)*
                  (type==MEAS_PHASE?-1.0:1.0);
            }
            for (k=0;k<nx;k++) H[k+nx*nv]=0.0;

            for (k=0;k<3;k++) H[k+nx*nv]=-e[k];

            /* receiver clock */
            cdtr=x[IC(sysi,opt)];
            H[IC(sysi,opt)+nx*nv]=1.0;

            if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                for (k=0;k<(opt->tropopt>=TROPOPT_ESTG?3:1);k++) {
                    H[IT(opt)+k+nx*nv]=dtdx[k];
                }
            }
            if (opt->ionoopt==IONOOPT_EST) {
                if (x[II(sat,opt)]==0.0) continue;
                H[II(sat,opt)+nx*nv]=C;
            }
            if (type==MEAS_CODE&&f!=0) {
                dcb=code_dcb(x,sysi,f,opt);
                H[ID(sysi,f,opt)+nx*nv]=1.0;
            }
            if (type==MEAS_PHASE) {
                if (opt->ionoopt==IONOOPT_IFLC) {
                    /* we keep ambiguity convention in meters for IFLC */
                    if ((bias=x[IB(sat,f,opt)])==0.0) continue;
                    H[IB(sat,f,opt)+nx*nv]=1.0;
                }
                else {
                    double freq=sat2freq(sat,obs[i].code[f],nav);
                    double lam=freq>0.0?CLIGHT/freq:0.0;
                    if (lam<=0.0) continue;
                    if (x[IB(sat,f,opt)]==0.0) continue;

                    bias=lam*x[IB(sat,f,opt)];      /* state in cycle -> model in meters */
                    H[IB(sat,f,opt)+nx*nv]=lam;     /* d(lambda*N)/dN = lambda */
                }
            }
            /* residual */
            v[nv]=y-(r+cdtr-CLIGHT*dts[i*2]+dtrp+C*dion+dcb+bias);

            if (type==MEAS_PHASE) rtk->ssat[sat-1].resc[f]=v[nv];
            else                  rtk->ssat[sat-1].resp[f]=v[nv];

            /* variance */
            var[nv]=varerr(obs[i].sat,sys,azel[1+i*2],f,type==MEAS_CODE,opt)+
                    vart+SQR(C)*vari+var_rs[i];
            if (sys==SYS_GLO&&type==MEAS_CODE) var[nv]+=VAR_GLO_IFB;

            if (type==MEAS_DOP) strcpy(mlabel,"DIF");
            else sprintf(mlabel,"%c%d",type==MEAS_CODE?'P':'L',f+1);
            trace(4,"%s sat=%2d %s res=%9.4f sig=%9.4f el=%4.1f\n",str,sat,
                  mlabel,v[nv],sqrt(var[nv]),
                  azel[1+i*2]*R2D);

            /* reject satellite by pre-fit residuals */
            if (!post&&opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) {
                trace(2,"outlier (%d) rejected %s sat=%2d %s res=%9.4f el=%4.1f\n",
                      post,str,sat,mlabel,v[nv],
                      azel[1+i*2]*R2D);
                exc[i]=1;
                rtk->ssat[sat-1].rejc[type==MEAS_PHASE?0:1]++;
                continue;
            }
            /* record large post-fit residuals */
            if (post&&fabs(v[nv])>sqrt(var[nv])*THRES_REJECT) {
                obsi[ne]=i;
                frqi[ne]=f;
                typei[ne]=type;
                ve[ne]=v[nv];
                ne++;
            }

            /* residuals */
            if (!post) {
                satno2id(sat,satid);
                time2epoch(obs[i].time,ep);
                sod=ep[3]*3600.0+ep[4]*60.0+ep[5];

                if (opt->ionoopt==IONOOPT_IFLC) {
                    strcpy(obsid,type==MEAS_PHASE?"LIF":"CIF");
                }
                else {
                    ocode=code2obs(obs[i].code[f]);
                    sprintf(obsid,"%c%s",type==MEAS_PHASE?'L':'C',
                            (ocode&&*ocode)?ocode:"?");
                }
                trace(2,"res: %s %.3f %s %.6f\n",satid,sod,obsid,v[nv]+bias);  /* residual without ambiguity */
            }

            if (type==MEAS_PHASE) rtk->ssat[sat-1].vsat[f]=1;
            nv++;
        }

        /* stack one iono-free Doppler residual per satellite */
        if (opt->dynamics&&Dc!=0.0) {
            if (nv>=MAXOBS*3*NFREQ) {
                trace(2,"ppp_res: too many residuals nv=%d\n",nv);
                return 0;
            }
            if (!satclkdrift_dop(obs+i,nav,&dtsd)) continue;

            for (k=0;k<nx;k++) H[k+nx*nv]=0.0;

            H[3+nx*nv]=-e[0];
            H[4+nx*nv]=-e[1];
            H[5+nx*nv]=-e[2];

            cdtrd=x[ICD(sysi,opt)];
            H[ICD(sysi,opt)+nx*nv]=1.0;

            modeled=dot(e,rs+i*6+3,3)-dot(e,x+3,3)+cdtrd-CLIGHT*dtsd;
            v[nv]=Dc-modeled;
            var[nv]=varerr_dop_iflc(obs+i,nav,obs[i].sat,sys,azel[1+i*2],opt);

            trace(4,"%s sat=%2d DIF res=%9.4f sig=%9.4f el=%4.1f\n",str,sat,
                  v[nv],sqrt(var[nv]),azel[1+i*2]*R2D);

            if (!post&&fabs(v[nv])>MAXINNO_DOP) {
                trace(2,"outlier (%d) rejected %s sat=%2d DIF res=%9.4f el=%4.1f\n",
                      post,str,sat,v[nv],azel[1+i*2]*R2D);
                continue;
            }

            /* trace des residus Doppler effectivement envoyes au filtre */
            if (!post) {
                satno2id(sat,satid);
                time2epoch(obs[i].time,ep);
                sod=ep[3]*3600.0+ep[4]*60.0+ep[5];
                trace(3,"%s %.3f DIF %.6f\n",satid,sod,v[nv]);
            }

            /* Do not let Doppler trigger whole-satellite post-fit rejection. */
            nv++;
        }
    }

    /* reject satellite with large and max post-fit residual */
    if (post&&ne>0) {
        vmax=ve[0];
        maxobs=obsi[0];
        maxfrq=frqi[0];
        maxtype=typei[0];
        rej=0;
        for (j=1;j<ne;j++) {
            if (fabs(vmax)>=fabs(ve[j])) continue;
            vmax=ve[j];
            maxobs=obsi[j];
            maxfrq=frqi[j];
            maxtype=typei[j];
            rej=j;
        }
        sat=obs[maxobs].sat;
        if (maxtype==MEAS_DOP) strcpy(mlabel,"DIF");
        else sprintf(mlabel,"%c%d",maxtype==MEAS_CODE?'P':'L',maxfrq+1);
        trace(2,"outlier (%d) rejected %s sat=%2d %s res=%9.4f el=%4.1f\n",
              post,str,sat,mlabel,vmax,
              azel[1+maxobs*2]*R2D);
        exc[maxobs]=1;
        rtk->ssat[sat-1].rejc[maxtype==MEAS_PHASE?0:1]++;
        stat=0;
        ve[rej]=0;
    }
    for (i=0;i<nv;i++) for (j=0;j<nv;j++) {
        R[i+j*nv]=i==j?var[i]:0.0;
    }
    return post?stat:nv;
}
/* number of estimated states ------------------------------------------------*/
extern int pppnx(const prcopt_t *opt)
{
    return NX(opt);
}
/* update solution status ----------------------------------------------------*/
static void update_stat(rtk_t *rtk, const obsd_t *obs, int n, int stat)
{
    const prcopt_t *opt=&rtk->opt;
    int i,j;
    
    /* test # of valid satellites */
    rtk->sol.ns=0;
    for (i=0;i<n&&i<MAXOBS;i++) {
        for (j=0;j<opt->nf;j++) {
            if (!rtk->ssat[obs[i].sat-1].vsat[j]) continue;
            rtk->ssat[obs[i].sat-1].lock[j]++;
            rtk->ssat[obs[i].sat-1].outc[j]=0;
            if (j==0) rtk->sol.ns++;
        }
    }
    rtk->sol.stat=rtk->sol.ns<MIN_NSAT_SOL?SOLQ_NONE:stat;
    
    if (rtk->sol.stat==SOLQ_FIX) {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->xa[i];
            rtk->sol.qr[i]=(float)rtk->Pa[i+i*rtk->na];
        }
        rtk->sol.qr[3]=(float)rtk->Pa[1];
        rtk->sol.qr[4]=(float)rtk->Pa[1+2*rtk->na];
        rtk->sol.qr[5]=(float)rtk->Pa[2];
    }
    else {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->x[i];
            rtk->sol.qr[i]=(float)rtk->P[i+i*rtk->nx];
        }
        rtk->sol.qr[3]=(float)rtk->P[1];
        rtk->sol.qr[4]=(float)rtk->P[2+rtk->nx];
        rtk->sol.qr[5]=(float)rtk->P[2];
    }
    rtk->sol.dtr[0]=rtk->x[IC(0,opt)];
    rtk->sol.dtr[1]=rtk->x[IC(1,opt)]-rtk->x[IC(0,opt)];
    
    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<opt->nf;j++) {
        rtk->ssat[obs[i].sat-1].snr[j]=obs[i].SNR[j];
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) {
        if (rtk->ssat[i].slip[j]&3) rtk->ssat[i].slipc[j]++;
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
/* precise point positioning -------------------------------------------------*/
extern void pppos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    const prcopt_t *opt=&rtk->opt;
    double *rs,*dts,*var,*v,*H,*R,*azel,*xp,*Pp,dr[3]={0},std[3];
    char str[32];
    int i,j,nv,info,svh[MAXOBS],exc[MAXOBS]={0},stat=SOLQ_SINGLE;
    
    time2str(obs[0].time,str,2);
    trace(3,"pppos   : time=%s nx=%d n=%d\n",str,rtk->nx,n);
    
    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel=zeros(2,n);
    
    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].fix[j]=0;
    
    /* temporal update of ekf states */
    udstate_ppp(rtk,obs,n,nav);
    
    /* satellite positions and clocks */
    satposs(obs[0].time,obs,n,nav,rtk->opt.sateph,rs,dts,var,svh);
    
    /* exclude measurements of eclipsing satellite (block IIA) */
    if (rtk->opt.posopt[3]) {
        testeclipse(obs,n,nav,rs);
    }
    /* earth tides correction */
    if (opt->tidecorr) {
        tidedisp(gpst2utc(obs[0].time),rtk->x,opt->tidecorr==1?1:7,&nav->erp,
                 opt->odisp[0],dr);
    }

    /* update/reinit ambiguities only after sat geometry is available */
    udbias_ppp(rtk,obs,n,nav,rs,dr);

    nv=n*rtk->opt.nf*3+MAXSAT+3;
    xp=mat(rtk->nx,1); Pp=zeros(rtk->nx,rtk->nx);
    v=mat(nv,1); H=mat(rtk->nx,nv); R=mat(nv,nv);
    
    for (i=0;i<MAX_ITER;i++) {
        
        matcpy(xp,rtk->x,rtk->nx,1);
        matcpy(Pp,rtk->P,rtk->nx,rtk->nx);
        
        /* prefit residuals */
        if (!(nv=ppp_res(0,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel))) {
            trace(2,"%s ppp (%d) no valid obs data\n",str,i+1);
            break;
        }
        /* measurement update of ekf states */
        if ((info=filter(xp,Pp,H,v,R,rtk->nx,nv))) {
            trace(2,"%s ppp (%d) filter error info=%d\n",str,i+1,info);
            break;
        }
        /* postfit residuals */
        if (ppp_res(i+1,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel)) {
            matcpy(rtk->x,xp,rtk->nx,1);
            matcpy(rtk->P,Pp,rtk->nx,rtk->nx);
            stat=SOLQ_PPP;
            break;
        }
    }
    if (i>=MAX_ITER) {
        trace(2,"%s ppp (%d) iteration overflows\n",str,i);
    }
    if (stat==SOLQ_PPP) {
        
        /* ambiguity resolution in ppp */
/*        if (ppp_ar(rtk,obs,n,exc,nav,azel,xp,Pp)&&
            ppp_res(9,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel)) { */
        if (ppp_ar(rtk,obs,n,exc,nav,azel,xp,Pp)) {
            
            matcpy(rtk->xa,xp,rtk->nx,1);
            matcpy(rtk->Pa,Pp,rtk->nx,rtk->nx);
            
            for (i=0;i<3;i++) std[i]=sqrt(Pp[i+i*rtk->nx]);
            if (norm(std,3)<MAX_STD_FIX) stat=SOLQ_FIX;
        }
        else {
            rtk->nfix=0;
        }
        /* update solution status */
        update_stat(rtk,obs,n,stat);
        
        /* hold fixed ambiguities */
        if (stat==SOLQ_FIX&&test_hold_amb(rtk)) {
            matcpy(rtk->x,xp,rtk->nx,1);
            matcpy(rtk->P,Pp,rtk->nx,rtk->nx);
            trace(2,"%s hold ambiguity\n",str);
            rtk->nfix=0;
        }
    }
    free(rs); free(dts); free(var); free(azel);
    free(xp); free(Pp); free(v); free(H); free(R);
}
