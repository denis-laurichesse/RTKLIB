/*------------------------------------------------------------------------------
* ppp_ar.c : ppp ambiguity resolution
*
*    [1] Banville et al., Enabling ambiguity resolution in CSRS-PPP
*         NAVIGATION, 2021
*
* Two-step partial AR using single-differences (SD):
*  1) SD wide-lanes (WL):
*       [(N_fx^i-N_f1^i) - (N_fx^p-N_f1^p)]
*  2) SD ambiguities on the reference frequency:
*       [N_f1^i - N_f1^p]
*
* Main properties:
*  - explicit ambiguity subset build (xa,Qa)
*  - lambda() on transformed ambiguity subsets
*  - accept a transformed ambiguity as valid if all candidates have the same
*    rounded integer value [1]
*  - reinject valid linear constraints into global x,P in sparse batch form
*
* The second step is performed on the state already constrained by the valid
* SD-WL fixes.
*
* The function returns 1 only if at least one constraint has been injected
* (either WL or reference-frequency SD).
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include <stdlib.h>
#include <math.h>

/* same state-index macros as in current ppp.c -------------------------------*/
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

#ifndef ROUND
#define ROUND(x)    (floor((x)+0.5))
#endif

/* lambda settings -----------------------------------------------------------*/
#define AR_LAMBDA_NCAND 5

/* enable(1)/disable(0) SD fixing on reference frequency N1 -----------------*/
#define PPPAR_ENABLE_N1_FIX 0

static int is_valid_amb(const rtk_t *rtk, const double *x, const double *P,
                        int sat, int f)
{
    int j=IB(sat,f,&rtk->opt);

    if (j<0 || j>=rtk->nx) return 0;
    if (x[j]==0.0) return 0;
    if (P[j+j*rtk->nx] <= 0.0) return 0;
    return 1;
}

/* helper to group satellites by constellation ------------------------------*/
static int get_sys_idx(int sys)
{
    int idx=0;

    while ((sys & 1)==0) {
        sys >>= 1;
        idx++;
    }
    return idx;
}

/* select one pivot satellite per constellation -----------------------------*/
static void select_pivot_sat(const int *isat, const int *ifreq, int na,
                             int *pivot_sat)
{
    int i,j,sys,sys_idx;

    for (i=0;i<10;i++) pivot_sat[i]=0;

    for (i=0;i<na;i++) {
        if (ifreq[i]!=1) continue;
        sys=satsys(isat[i],NULL);
        if (!sys) continue;

        sys_idx=get_sys_idx(sys);
        if (pivot_sat[sys_idx]!=0) continue;

        /* Keep the same pivot policy as the WL stage: choose an f1 satellite
         * that also owns at least one additional active ambiguity. */
        for (j=0;j<na;j++) {
            if (isat[j]==isat[i] && ifreq[j]>1) {
                pivot_sat[sys_idx]=isat[i];
                break;
            }
        }
    }
}

/* build explicit active ambiguity subset -----------------------------------*/
static int build_active_amb_set(const rtk_t *rtk,
                                const double *x, const double *P,
                                int *na_out,
                                int **ix_out, int **isat_out, int **ifreq_out,
                                double **xa_out, double **Qa_out)
{
    const prcopt_t *opt=&rtk->opt;
    int na,sat,f,j,i,nx,k;
    int *ix=NULL,*isat=NULL,*ifreq=NULL;
    double *xa=NULL,*Qa=NULL;

    if (!na_out||!ix_out||!isat_out||!ifreq_out||!xa_out||!Qa_out) return 0;

    *na_out=0; *ix_out=NULL; *isat_out=NULL; *ifreq_out=NULL;
    *xa_out=NULL; *Qa_out=NULL;

    nx=rtk->nx;
    na=0;

    for (f=0;f<NF(opt);f++) {
        for (sat=1;sat<=MAXSAT;sat++) {
            if (is_valid_amb(rtk,x,P,sat,f)) na++;
        }
    }

    *na_out=na;
    if (na<=0) return 1;

    ix    = imat(na,1);
    isat  = imat(na,1);
    ifreq = imat(na,1);
    xa    = mat(na,1);
    Qa    = mat(na,na);

    if (!ix || !isat || !ifreq || !xa || !Qa) {
        free(ix); free(isat); free(ifreq); free(xa); free(Qa);
        trace(2,"ppp_ar: memory allocation error in build_active_amb_set\n");
        return 0;
    }

    k=0;
    for (f=0;f<NF(opt);f++) {
        for (sat=1;sat<=MAXSAT;sat++) {
            if (is_valid_amb(rtk,x,P,sat,f)) {
                j=IB(sat,f,opt);
                ix[k]=j;
                isat[k]=sat;
                ifreq[k]=f+1;
                xa[k]=x[j];
                k++;
            }
        }
    }

    for (j=0;j<na;j++) {
        for (i=0;i<na;i++) {
            Qa[i+na*j]=P[ix[i]+nx*ix[j]];
        }
    }

    *ix_out=ix; *isat_out=isat; *ifreq_out=ifreq;
    *xa_out=xa; *Qa_out=Qa;

    return 1;
}

/* build T for Single-Difference wide-lanes ---------------------------------*/
static int build_wl_transform(int na, const int *isat, const int *ifreq,
                              const int *use_sys,
                              int *nwl_out, double **Twl_out,
                              int **wl_sat_out, int **wl_f_out)
{
    int i,j,k,nwl,sys,sys_idx,p_sat,idx_p1,idx_px,row;
    int pivot_sat[10];
    double *Twl=NULL;
    int *wl_sat=NULL,*wl_f=NULL;

    if (!nwl_out||!Twl_out||!wl_sat_out||!wl_f_out) return 0;

    *nwl_out=0; *Twl_out=NULL; *wl_sat_out=NULL; *wl_f_out=NULL;
    if (na<=0) return 1;

    select_pivot_sat(isat,ifreq,na,pivot_sat);

    /* count possible SD-WLs */
    nwl=0;
    for (i=0;i<na;i++) {
        if (ifreq[i]!=1) continue;
        sys=satsys(isat[i],NULL);
        if (!sys) continue;

        sys_idx=get_sys_idx(sys);
        if (use_sys && !use_sys[sys_idx]) continue;
        p_sat=pivot_sat[sys_idx];
        if (p_sat==0 || isat[i]==p_sat) continue;

        for (j=0;j<na;j++) {
            if (isat[j]==isat[i] && ifreq[j]>1) {
                for (k=0;k<na;k++) {
                    if (isat[k]==p_sat && ifreq[k]==ifreq[j]) {
                        nwl++;
                        break;
                    }
                }
            }
        }
    }

    *nwl_out=nwl;
    if (nwl<=0) return 1;

    Twl    = zeros(nwl,na);
    wl_sat = imat(nwl,1);
    wl_f   = imat(nwl,1);

    if (!Twl||!wl_sat||!wl_f) {
        free(Twl); free(wl_sat); free(wl_f);
        return 0;
    }

    /* build the SD-WL transformation matrix */
    row=0;
    for (i=0;i<na;i++) {
        if (ifreq[i]!=1) continue;
        sys=satsys(isat[i],NULL);
        if (!sys) continue;

        sys_idx=get_sys_idx(sys);
        if (use_sys && !use_sys[sys_idx]) continue;
        p_sat=pivot_sat[sys_idx];
        if (p_sat==0 || isat[i]==p_sat) continue;

        for (j=0;j<na;j++) {
            if (isat[j]==isat[i] && ifreq[j]>1) {
                idx_p1=-1;
                idx_px=-1;
                for (k=0;k<na;k++) {
                    if (isat[k]==p_sat) {
                        if (ifreq[k]==1)         idx_p1=k;
                        if (ifreq[k]==ifreq[j]) idx_px=k;
                    }
                }
                if (idx_p1>=0 && idx_px>=0) {
                    /* (N_fx^i - N_f1^i) - (N_fx^p - N_f1^p) */
                    Twl[row+nwl*j]      =  1.0;
                    Twl[row+nwl*i]      = -1.0;
                    Twl[row+nwl*idx_px] = -1.0;
                    Twl[row+nwl*idx_p1] =  1.0;
                    wl_sat[row]=isat[i];
                    wl_f[row]  =ifreq[j];
                    row++;
                }
            }
        }
    }

    *nwl_out=row;
    *Twl_out=Twl;
    *wl_sat_out=wl_sat;
    *wl_f_out=wl_f;
    return 1;
}

/* build T for SD ambiguities on the reference frequency --------------------*/
static int build_ref_transform(int na, const int *isat, const int *ifreq,
                               const int *use_sys,
                               int *nref_out, double **Tref_out,
                               int **ref_sat_out)
{
    int i,j,nref,row,sys,sys_idx,p_sat;
    int pivot_sat[10];
    double *Tref=NULL;
    int *ref_sat=NULL;

    if (!nref_out||!Tref_out||!ref_sat_out) return 0;

    *nref_out=0;
    *Tref_out=NULL;
    *ref_sat_out=NULL;
    if (na<=0) return 1;

    select_pivot_sat(isat,ifreq,na,pivot_sat);

    /* count all possible SD on f1 */
    nref=0;
    for (i=0;i<na;i++) {
        if (ifreq[i]!=1) continue;
        sys=satsys(isat[i],NULL);
        if (!sys) continue;

        sys_idx=get_sys_idx(sys);
        if (use_sys && !use_sys[sys_idx]) continue;
        p_sat=pivot_sat[sys_idx];
        if (p_sat==0 || isat[i]==p_sat) continue;
        nref++;
    }

    *nref_out=nref;
    if (nref<=0) return 1;

    Tref   = zeros(nref,na);
    ref_sat= imat(nref,1);
    if (!Tref || !ref_sat) {
        free(Tref); free(ref_sat);
        return 0;
    }

    row=0;
    for (i=0;i<na;i++) {
        if (ifreq[i]!=1) continue;
        sys=satsys(isat[i],NULL);
        if (!sys) continue;

        sys_idx=get_sys_idx(sys);
        if (use_sys && !use_sys[sys_idx]) continue;
        p_sat=pivot_sat[sys_idx];
        if (p_sat==0 || isat[i]==p_sat) continue;

        for (j=0;j<na;j++) {
            if (isat[j]==p_sat && ifreq[j]==1) {
                /* N_f1^i - N_f1^p */
                Tref[row+nref*i] =  1.0;
                Tref[row+nref*j] = -1.0;
                ref_sat[row]=isat[i];
                row++;
                break;
            }
        }
    }

    *nref_out=row;
    *Tref_out=Tref;
    *ref_sat_out=ref_sat;
    return 1;
}

/* build transformed system y = T*xa, Py = T*Qa*T' --------------------------*/
static int build_linear_system(int na, int ny, const double *T,
                               const double *xa, const double *Qa,
                               double **y_out, double **Py_out)
{
    double *y=NULL,*Py=NULL,*tmp=NULL;
    int i,j,k;

    if (!y_out || !Py_out) return 0;
    *y_out=NULL; *Py_out=NULL;
    if (ny<=0 || na<=0) return 1;

    y  = mat(ny,1);
    Py = mat(ny,ny);
    tmp= mat(ny,na);

    if (!y || !Py || !tmp) {
        free(y); free(Py); free(tmp);
        return 0;
    }

    /* y = T * xa */
    for (i=0;i<ny;i++) {
        y[i]=0.0;
        for (k=0;k<na;k++) y[i]+=T[i+ny*k]*xa[k];
    }

    /* tmp = T * Qa */
    for (j=0;j<na;j++) {
        for (i=0;i<ny;i++) {
            tmp[i+ny*j]=0.0;
            for (k=0;k<na;k++) tmp[i+ny*j]+=T[i+ny*k]*Qa[k+na*j];
        }
    }

    /* Py = tmp * T' */
    for (j=0;j<ny;j++) {
        for (i=0;i<ny;i++) {
            Py[i+ny*j]=0.0;
            for (k=0;k<na;k++) Py[i+ny*j]+=tmp[i+ny*k]*T[j+ny*k];
        }
    }

    free(tmp);
    *y_out=y;
    *Py_out=Py;
    return 1;
}

static int run_linear_lambda(int ny, const double *y, const double *Py,
                             int *lambda_info_out,
                             double **F_out, double **s_out)
{
    double *F=NULL,*s=NULL;
    int info;

    if (!lambda_info_out || !F_out || !s_out) return 0;

    *lambda_info_out=0;
    *F_out=NULL;
    *s_out=NULL;
    if (ny<=0) return 1;

    F=mat(ny,AR_LAMBDA_NCAND);
    s=mat(AR_LAMBDA_NCAND,1);
    if (!F || !s) {
        free(F); free(s);
        return 0;
    }

    info=lambda(ny,AR_LAMBDA_NCAND,y,Py,F,s);
    *lambda_info_out=info;

    if (info) {
        free(F); free(s);
        trace(2,"ppp_ar: lambda failed info=%d\n",info);
        return 1;
    }

    *F_out=F;
    *s_out=s;
    return 1;
}

/* detect valid transformed ambiguities and reinject them into x,P ----------*/
static int apply_valid_linear_constraints(double *x, double *P, int nx,
                                          int na, const int *ix,
                                          int ny, const double *T,
                                          int lambda_info,
                                          const double *F,
                                          int *valid, double *fix)
{
    double *v=NULL,*HP=NULL,*PHt=NULL,*S=NULL,*K=NULL,*KP=NULL;
    int *idx_p1=NULL,*idx_p2=NULL,*idx_n1=NULL,*idx_n2=NULL;
    double ref,hx,val,sym;
    int i,j,k,m,info,nfix,row,p1,p2,n1,n2,c,r;

    if (valid) for (i=0;i<ny;i++) valid[i]=0;
    if (fix)   for (i=0;i<ny;i++) fix[i]=0.0;

    if (lambda_info!=0 || !F || ny<=0 || na<=0) return 0;

    /* valid if all candidates provide the same rounded integer */
    nfix=0;
    for (i=0;i<ny;i++) {
        ref=ROUND(F[i]);
        for (k=1;k<AR_LAMBDA_NCAND;k++) {
            if (ROUND(F[i+ny*k])!=ref) break;
        }
        if (k<AR_LAMBDA_NCAND) continue;
        if (valid) valid[i]=1;
        if (fix)   fix[i]=ref;
        nfix++;
    }
    if (nfix<=0) return 0;

    m      = nfix;
    idx_p1 = imat(m,1); idx_p2 = imat(m,1);
    idx_n1 = imat(m,1); idx_n2 = imat(m,1);
    v      = mat(m,1);
    HP     = mat(m,nx);
    PHt    = mat(nx,m);
    S      = mat(m,m);
    K      = mat(nx,m);
    KP     = mat(nx,nx);

    if (!idx_p1||!idx_p2||!idx_n1||!idx_n2||!v||!HP||!PHt||!S||!K||!KP) {
        free(idx_p1); free(idx_p2); free(idx_n1); free(idx_n2);
        free(v); free(HP); free(PHt); free(S); free(K); free(KP);
        return 0;
    }

    row=0;
    for (i=0;i<ny;i++) {
        if (!valid[i]) continue;

        p1=-1; p2=-1; n1=-1; n2=-1;
        for (j=0;j<na;j++) {
            if (T[i+ny*j]== 1.0) {
                if (p1<0) p1=ix[j]; else p2=ix[j];
            }
            else if (T[i+ny*j]==-1.0) {
                if (n1<0) n1=ix[j]; else n2=ix[j];
            }
        }

        idx_p1[row]=p1; idx_p2[row]=p2;
        idx_n1[row]=n1; idx_n2[row]=n2;

        hx=0.0;
        if (p1>=0) hx+=x[p1];
        if (p2>=0) hx+=x[p2];
        if (n1>=0) hx-=x[n1];
        if (n2>=0) hx-=x[n2];
        v[row]=fix[i]-hx;
        row++;
    }

    /* HP = H*P, PHt = P*H' */
    for (i=0;i<m;i++) {
        p1=idx_p1[i]; p2=idx_p2[i];
        n1=idx_n1[i]; n2=idx_n2[i];
        for (j=0;j<nx;j++) {
            val=0.0;
            if (p1>=0) val+=P[p1+j*nx];
            if (p2>=0) val+=P[p2+j*nx];
            if (n1>=0) val-=P[n1+j*nx];
            if (n2>=0) val-=P[n2+j*nx];
            HP[i+j*m]=val;
            PHt[j+i*nx]=val;
        }
    }

    /* S = H*P*H' */
    for (c=0;c<m;c++) {
        for (r=0;r<m;r++) {
            val=0.0;
            if (idx_p1[r]>=0) val+=PHt[idx_p1[r]+c*nx];
            if (idx_p2[r]>=0) val+=PHt[idx_p2[r]+c*nx];
            if (idx_n1[r]>=0) val-=PHt[idx_n1[r]+c*nx];
            if (idx_n2[r]>=0) val-=PHt[idx_n2[r]+c*nx];
            S[r+c*m]=val;
        }
    }

    if ((info=matinv(S,m))) {
        trace(2,"ppp_ar: apply_valid_linear_constraints matinv error info=%d\n",info);
        free(idx_p1); free(idx_p2); free(idx_n1); free(idx_n2);
        free(v); free(HP); free(PHt); free(S); free(K); free(KP);
        return 0;
    }

    /* K = P*H' * inv(S) */
    matmul("NN",nx,m,m,1.0,PHt,S,0.0,K);

    /* x = x + K*v */
    matmul("NN",nx,1,m,1.0,K,v,1.0,x);

    /* P = P - K*(H*P) */
    matmul("NN",nx,nx,m,1.0,K,HP,0.0,KP);
    for (j=0;j<nx;j++) {
        for (i=0;i<nx;i++) P[i+nx*j]-=KP[i+nx*j];
    }

    /* light symmetrization */
    for (j=0;j<nx;j++) {
        for (i=j+1;i<nx;i++) {
            sym=0.5*(P[i+nx*j]+P[j+nx*i]);
            P[i+nx*j]=sym;
            P[j+nx*i]=sym;
        }
    }

    free(idx_p1); free(idx_p2); free(idx_n1); free(idx_n2);
    free(v); free(HP); free(PHt); free(S); free(K); free(KP);
    return m;
}

/* ambiguity resolution in ppp ----------------------------------------------*/
extern int ppp_ar(rtk_t *rtk, const obsd_t *obs, int n, int *exc,
                  const nav_t *nav, const double *azel, double *x, double *P)
{
    int na=0;
    int *ix=NULL,*isat=NULL,*ifreq=NULL;
    double *xa=NULL,*Qa=NULL;

    int nwl=0;
    double *Twl=NULL;
    int *wl_sat=NULL,*wl_f=NULL;
    double *wl_x=NULL,*wl_P=NULL;
    int wl_lambda_info=0;
    double *wl_F=NULL,*wl_s=NULL;
    int *wl_valid=NULL;
    double *wl_fix=NULL;
    int nwl_fix=0;

#if PPPAR_ENABLE_N1_FIX
    int nref=0;
    double *Tref=NULL;
    int *ref_sat=NULL;
    double *ref_x=NULL,*ref_P=NULL;
    int ref_lambda_info=0;
    double *ref_F=NULL,*ref_s=NULL;
    int *ref_valid=NULL;
    double *ref_fix=NULL;
    int nref_fix=0;

    int ref_sys_mask[10]={0};
#endif

    int ar_ok=0;

    (void)exc;
    (void)nav;
    (void)azel;
    (void)obs;

    if (!rtk || n<=0 || !x || !P || rtk->opt.modear==ARMODE_OFF) return 0;

    /* ------------------------------------------------------------------
     * pass 1: SD wide-lanes
     * ------------------------------------------------------------------*/
    if (!build_active_amb_set(rtk,x,P,&na,&ix,&isat,&ifreq,&xa,&Qa)) goto cleanup;

    if (!build_wl_transform(na,isat,ifreq,NULL,&nwl,&Twl,&wl_sat,&wl_f)) goto cleanup;

    if (!build_linear_system(na,nwl,Twl,xa,Qa,&wl_x,&wl_P)) goto cleanup;

    wl_valid=imat(nwl>0?nwl:1,1);
    wl_fix  =mat(nwl>0?nwl:1,1);
    if (!wl_valid || !wl_fix) goto cleanup;

    if (nwl>0 && rtk->opt.ionoopt!=IONOOPT_IFLC) {
        if (!run_linear_lambda(nwl,wl_x,wl_P,&wl_lambda_info,&wl_F,&wl_s)) goto cleanup;

        nwl_fix=apply_valid_linear_constraints(x,P,rtk->nx,
                                               na,ix,
                                               nwl,Twl,
                                               wl_lambda_info,
                                               wl_F,
                                               wl_valid,wl_fix);
        if (nwl_fix>0) ar_ok=1;
    }

#if PPPAR_ENABLE_N1_FIX
    /* ------------------------------------------------------------------
     * pass 2: SD ambiguities on the reference frequency, after WL update
     * ------------------------------------------------------------------*/
    if (nwl_fix>0) {
        int i,sys,sys_idx;

        for (i=0;i<nwl;i++) {
            if (!wl_valid || !wl_valid[i]) continue;
            sys=satsys(wl_sat[i],NULL);
            if (!sys) continue;
            sys_idx=get_sys_idx(sys);
            if (sys_idx>=0 && sys_idx<10) ref_sys_mask[sys_idx]=1;
        }

        free(ix);    ix=NULL;
        free(isat);  isat=NULL;
        free(ifreq); ifreq=NULL;
        free(xa);    xa=NULL;
        free(Qa);    Qa=NULL;
        na=0;

        if (!build_active_amb_set(rtk,x,P,&na,&ix,&isat,&ifreq,&xa,&Qa)) goto cleanup;

        if (!build_ref_transform(na,isat,ifreq,ref_sys_mask,&nref,&Tref,&ref_sat)) goto cleanup;

        if (!build_linear_system(na,nref,Tref,xa,Qa,&ref_x,&ref_P)) goto cleanup;

        ref_valid=imat(nref>0?nref:1,1);
        ref_fix  =mat(nref>0?nref:1,1);
        if (!ref_valid || !ref_fix) goto cleanup;

        if (nref>0) {
            if (!run_linear_lambda(nref,ref_x,ref_P,
                                   &ref_lambda_info,&ref_F,&ref_s)) goto cleanup;

            nref_fix=apply_valid_linear_constraints(x,P,rtk->nx,
                                                    na,ix,
                                                    nref,Tref,
                                                    ref_lambda_info,
                                                    ref_F,
                                                    ref_valid,ref_fix);
            if (nref_fix>0) ar_ok=1;
        }
    }
#endif

cleanup:
    free(ix);    free(isat);  free(ifreq);
    free(xa);    free(Qa);

    free(Twl);   free(wl_sat); free(wl_f);
    free(wl_x);  free(wl_P);
    free(wl_F);  free(wl_s);
    free(wl_valid); free(wl_fix);

#if PPPAR_ENABLE_N1_FIX
    free(Tref);  free(ref_sat);
    free(ref_x); free(ref_P);
    free(ref_F); free(ref_s);
    free(ref_valid); free(ref_fix);
#endif

    return ar_ok;
}
