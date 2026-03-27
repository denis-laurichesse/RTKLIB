/*-----------------------------------------------------------------------------
 * sinexbias.c : SINEX Bias format (OSB) reader
 *
 * Public API:
 *   - readsinexbias()
 *   - getsinexbias()
 *   - freesinexbias()
 *
 * Internal policy:
 *   - exact lookup is always attempted first
 *   - optional fallback searches another RTKLIB CODE_* on the same carrier
 *     frequency, only if enabled for the satellite constellation
 *---------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rtklib.h"

/* Internal constants -------------------------------------------------------*/
#define OSB_BUFSIZE    256
#define OSB_INITCAP    4096
#define OSB_EMPTY_SLOT (-1)
#define NS2M           (CLIGHT*1E-9)

/* helpers */
static const char *osb_type_str(int type)
{
    return type==BIAS_CODE ? "CODE" :
           type==BIAS_PHASE ? "PHASE" : "UNKNOWN";
}

static const char *osb_sig_str(uint8_t code)
{
    char *obs=code2obs(code);
    return (obs&&*obs) ? obs : "?";
}

/* Pack (sat, code, type) into a single integer key ------------------------*/
static int osb_makekey(int sat, uint8_t code, int type)
{
    return (sat << 8) | ((int)code << 1) | (type & 1);
}
/* integer log2 for powers of two ------------------------------------------*/
static int osb_ilog2(int size)
{
    int k=0;
    while ((1 << k) < size) k++;
    return k;
}
/* multiplicative hash ------------------------------------------------------*/
static int osb_hash(int key, int size)
{
    return (int)(((unsigned int)key * 0x9e3779b9u) >> (32 - osb_ilog2(size)));
}
/* SINEX YYYY:DOY:SOD -> gtime_t -------------------------------------------*/
static gtime_t sinex_str2time(int year, int doy, int sod)
{
    static const int mdays[12]={31,28,31,30,31,30,31,31,30,31,30,31};
    int mon,day,leap;
    double ep[6];

    if (year < 100) year += (year < 80) ? 2000 : 1900;
    leap = ((year%4==0 && year%100!=0) || year%400==0) ? 1 : 0;

    for (mon=0,day=doy; mon<12; mon++) {
        int dm = mdays[mon] + (mon==1 ? leap : 0);
        if (day <= dm) break;
        day -= dm;
    }
    if (mon >= 12) { mon=11; day=mdays[11]; }

    ep[0]=(double)year;
    ep[1]=(double)(mon+1);
    ep[2]=(double)day;
    ep[3]=(double)(sod/3600);
    ep[4]=(double)((sod%3600)/60);
    ep[5]=(double)(sod%60);
    return epoch2time(ep);
}
/* expand data[] ------------------------------------------------------------*/
static int osb_expand(osbdata_t *osb)
{
    int nmax;
    osb_t *tmp;

    if (osb->n < osb->nmax) return 1;

    nmax = (osb->nmax == 0) ? OSB_INITCAP : osb->nmax * 2;
    tmp  = (osb_t *)realloc(osb->data, sizeof(osb_t) * nmax);

    if (!tmp) {
        trace(1,"readsinexbias: out of memory\n");
        return 0;
    }
    osb->data = tmp;
    osb->nmax = nmax;
    return 1;
}
/* qsort comparator ---------------------------------------------------------*/
static int osb_cmp(const void *a, const void *b)
{
    const osb_t *ra=(const osb_t *)a;
    const osb_t *rb=(const osb_t *)b;
    int d;
    double dt;

    if ((d = ra->sat       - rb->sat      ) != 0) return d;
    if ((d = (int)ra->code - (int)rb->code) != 0) return d;
    if ((d = ra->type      - rb->type     ) != 0) return d;

    dt = timediff(ra->ts, rb->ts);
    return (dt < 0.0) ? -1 : (dt > 0.0) ? 1 : 0;
}
/* next power of two --------------------------------------------------------*/
static int next_pow2(int n)
{
    int p=1;
    while (p < n) p <<= 1;
    return p;
}
/* build open-addressing hash table ----------------------------------------*/
static int osb_buildindex(osbdata_t *osb)
{
    int i,n_keys,size,slot,key;
    osbslot_t *ht;

    if (!osb) return 0;

    if (osb->htbl) {
        free(osb->htbl);
        osb->htbl = NULL;
        osb->htbl_size = 0;
    }

    n_keys = 0;
    for (i=0; i<osb->n; ) {
        int k = osb_makekey(osb->data[i].sat,
                            osb->data[i].code,
                            osb->data[i].type);
        do { i++; }
        while (i < osb->n &&
               osb_makekey(osb->data[i].sat,
                           osb->data[i].code,
                           osb->data[i].type) == k);
        n_keys++;
    }

    size = next_pow2(4 * (n_keys > 0 ? n_keys : 1));
    ht   = (osbslot_t *)malloc(sizeof(osbslot_t) * size);

    if (!ht) {
        trace(1,"osb_buildindex: out of memory\n");
        return 0;
    }
    for (i=0; i<size; i++) {
        ht[i].key = OSB_EMPTY_SLOT;
        ht[i].idx = 0;
        ht[i].n   = 0;
    }

    for (i=0; i<osb->n; ) {
        int run_start = i;

        key = osb_makekey(osb->data[i].sat,
                          osb->data[i].code,
                          osb->data[i].type);
        do { i++; }
        while (i < osb->n &&
               osb_makekey(osb->data[i].sat,
                           osb->data[i].code,
                           osb->data[i].type) == key);

        slot = osb_hash(key, size);
        while (ht[slot].key != OSB_EMPTY_SLOT) {
            slot = (slot + 1) & (size - 1);
        }
        ht[slot].key = key;
        ht[slot].idx = run_start;
        ht[slot].n   = i - run_start;
    }

    osb->htbl      = ht;
    osb->htbl_size = size;

    trace(3,"osb_buildindex: %d unique signals, hash table size %d\n",
          n_keys,size);
    return 1;
}
/* same carrier frequency? --------------------------------------------------*/
static int same_obs_freq(int sat, uint8_t code1, uint8_t code2,
                         const nav_t *nav)
{
    double f1,f2,df,fref;

    if (!nav || code1==CODE_NONE || code2==CODE_NONE) return 0;

    f1 = sat2freq(sat, code1, (nav_t *)nav);
    f2 = sat2freq(sat, code2, (nav_t *)nav);

    if (f1 <= 0.0 || f2 <= 0.0) return 0;

    df   = fabs(f1 - f2);
    fref = f1 > f2 ? f1 : f2;

    return df <= 1E-6 * fref;
}
/* exact lookup only --------------------------------------------------------*/
static int osb_lookup_exact(const osbdata_t *osb, gtime_t time, int sat,
                            uint8_t code, int type, double *bias, double *std)
{
    int key,slot,lo,hi,mid;
    char satid[16];

    if (!osb || !osb->htbl || !bias) return 0;

    satno2id(sat,satid);

    key  = osb_makekey(sat, code, type);
    slot = osb_hash(key, osb->htbl_size);

    for (;;) {
        const osbslot_t *s = &osb->htbl[slot];

        if (s->key == OSB_EMPTY_SLOT) {
            trace(3,"osb_lookup_exact: sat=%s sig=%s type=%s not in OSB file\n",
                  satid,osb_sig_str(code),osb_type_str(type));
            return 0;
        }
        if (s->key == key) {
            lo = s->idx;
            hi = s->idx + s->n - 1;

            while (lo <= hi) {
                const osb_t *r;

                mid = lo + (hi - lo) / 2;
                r   = &osb->data[mid];

                if (timediff(time, r->ts) < -0.1) {
                    hi = mid - 1;
                }
                else if (timediff(time, r->te) >= 0.1) {
                    lo = mid + 1;
                }
				else {
					*bias = r->bias;
					if (std) *std = r->std;

					trace(3,"osb_lookup_exact: sat=%s sig=%s type=%s bias=%.4f m std=%.4f m\n",
						satid,osb_sig_str(code),osb_type_str(type),*bias,std?*std:r->std);

					return 1;
				}
			}
            trace(2,"osb_lookup_exact: sat=%s sig=%s type=%s epoch outside OSB coverage\n",
                  satid,osb_sig_str(code),osb_type_str(type));
            return 0;
        }
        slot = (slot + 1) & (osb->htbl_size - 1);
    }
}
/* lookup with optional same-frequency fallback --------------------*/
static int osb_lookup_with_fallback(const osbdata_t *osb, gtime_t time, int sat,
                                    uint8_t code, int type, const nav_t *nav,
                                    int fallback_mask,
                                    double *bias, double *std,
                                    uint8_t *code_used)
{
    uint8_t c;
    double b=0.0,s=0.0;
    int sys=satsys(sat,NULL);
    char satid[16];

    if (bias) *bias = 0.0;
    if (std ) *std  = 0.0;
    if (code_used) *code_used = CODE_NONE;

    if (!osb || !bias) return 0;

    satno2id(sat,satid);

    /* 1) exact signal first */
    if (osb_lookup_exact(osb,time,sat,code,type,&b,&s)) {
        *bias = b;
        if (std) *std = s;
        if (code_used) *code_used = code;
        return 1;
    }

    /* 2) fallback disabled or impossible */
    if (!nav || !(sys & fallback_mask)) return 0;

    /* 3) same-frequency fallback */
    for (c=1; c<=MAXCODE; c++) {
        if (c == code) continue;
        if (!same_obs_freq(sat, code, c, nav)) continue;

        if (osb_lookup_exact(osb,time,sat,c,type,&b,&s)) {
            *bias = b;
            if (std) *std = s;
            if (code_used) *code_used = c;

            trace(3,"osb_lookup_with_fallback: sat=%s type=%s req_sig=%s fallback_sig=%s bias=%.4f m\n",
                  satid,osb_type_str(type),osb_sig_str(code),osb_sig_str(c),b);
            return 1;
        }
    }
    return 0;
}
/* readsinexbias ------------------------------------------------------------*/
int readsinexbias(const char *file, osbdata_t *osb, const nav_t *nav)
{
    FILE *fp;
    char  buf[OSB_BUFSIZE];
    int   in_solution=0,nread=0;

    if (!file || !osb) return -1;

    fp = fopen(file,"r");
    if (!fp) {
        trace(1,"readsinexbias: cannot open %s\n",file);
        return -1;
    }

    while (fgets(buf,sizeof(buf),fp)) {
        int    linelen;
        char   rec_type[4]={0};
        char   prn_str [4]={0};
        char   obs1    [5]={0};
        char   t_start [15]={0};
        char   t_end   [15]={0};
        char   unit    [5]={0};
        double val=0.0,std_val=0.0;
        int    n_val;
        int    sat,meas_type;
        uint8_t code;
        int    y0,d0,s0,y1,d1,s1;
        gtime_t ts,te;
        double bias_m,std_m;

        if      (strncmp(buf,"+BIAS/SOLUTION",14)==0) { in_solution=1; continue; }
        else if (strncmp(buf,"-BIAS/SOLUTION",14)==0) { in_solution=0; continue; }
        if (!in_solution) continue;
        if (buf[0]=='*' || buf[0]=='%') continue;

        linelen = (int)strlen(buf);
        if (linelen < 70) continue;

        memcpy(rec_type, buf+ 1, 3);
        memcpy(prn_str,  buf+11, 3);
        memcpy(obs1,     buf+25, 3);
        memcpy(t_start,  buf+35,14);
        memcpy(t_end,    buf+50,14);
        memcpy(unit,     buf+65, 4);

        n_val = sscanf(buf+69," %lf %lf",&val,&std_val);

        if (strncmp(rec_type,"OSB",3) != 0) continue;
        if (n_val < 1) continue;

        sat = satid2no(prn_str);
        if (sat == 0) continue;

        switch (obs1[0]) {
            case 'C': meas_type = BIAS_CODE;  break;
            case 'L': meas_type = BIAS_PHASE; break;
            default : continue;
        }
        code = obs2code(obs1+1);
        if (code == CODE_NONE) continue;

        if (sscanf(t_start,"%d:%d:%d",&y0,&d0,&s0) != 3) continue;
        if (sscanf(t_end  ,"%d:%d:%d",&y1,&d1,&s1) != 3) continue;

        ts = sinex_str2time(y0,d0,s0);
        te = sinex_str2time(y1,d1,s1);

        if (strncmp(unit,"ns",2) == 0) {
            bias_m = val * NS2M;
            std_m  = (n_val >= 2) ? std_val * NS2M : 0.0;
        }
        else if (strncmp(unit,"cyc",3) == 0) {
            double freq;
            int sys_flag = satsys(sat,NULL);

            if (sys_flag == SYS_GLO && nav == NULL) {
				trace(2,"readsinexbias: skip cyc record sat=%s sig=%s type=%s (GLO needs nav->glo_fcn)\n",
					prn_str,osb_sig_str(code),osb_type_str(meas_type));
                continue;
            }
            freq = sat2freq(sat, code, (nav_t *)nav);
            if (freq <= 0.0) {
				trace(2,"readsinexbias: skip cyc record sat=%s sig=%s type=%s freq=0\n",
					prn_str,osb_sig_str(code),osb_type_str(meas_type));
                continue;
            }
            bias_m = val * CLIGHT / freq;
            std_m  = (n_val >= 2) ? std_val * CLIGHT / freq : 0.0;
        }
        else {
            bias_m = val;
            std_m  = (n_val >= 2) ? std_val : 0.0;
        }

        if (!osb_expand(osb)) {
            fclose(fp);
            return -1;
        }
        osb->data[osb->n].sat  = sat;
        osb->data[osb->n].code = code;
        osb->data[osb->n].type = meas_type;
        osb->data[osb->n].ts   = ts;
        osb->data[osb->n].te   = te;
        osb->data[osb->n].bias = bias_m;
        osb->data[osb->n].std  = std_m;
        osb->n++;
        nread++;
    }
    fclose(fp);

    qsort(osb->data, osb->n, sizeof(osb_t), osb_cmp);
    if (!osb_buildindex(osb)) return -1;

    trace(3,"readsinexbias: %d records read from %s\n",nread,file);
    return nread;
}
/* getsinexbias -------------------------------------------------------------*/
int getsinexbias(const osbdata_t *osb, gtime_t time, int sat,
                 uint8_t code, int type, const nav_t *nav,
                 int fallback_mask,
                 double *bias, double *std, uint8_t *code_used)
{
    return osb_lookup_with_fallback(osb,time,sat,code,type,nav,fallback_mask,
                                    bias,std,code_used);
}
/* freesinexbias ------------------------------------------------------------*/
void freesinexbias(osbdata_t *osb)
{
    if (!osb) return;

    free(osb->data);
    free(osb->htbl);

    osb->data      = NULL;
    osb->htbl      = NULL;
    osb->n         = 0;
    osb->nmax      = 0;
    osb->htbl_size = 0;
}
