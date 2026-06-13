/*
 * tic_gr3.c — TIC with extended Golomb-Rice
 *
 * Extensions over tic_gr2:
 *   1. 128 contexts (gradient magnitude × direction)
 *   2. Run-length mode for consecutive zero residuals
 *   3. Better escape coding for large residuals
 *   4. Per-plane statistics header for decoder init
 *
 * No external dependencies. Compile: gcc -O2 -o tic_gr3 tic_gr3.c -lm
 *
 * kind-pasteur-2026-03-26-S33
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

/* ================================================================ */
/* BITSTREAM                                                         */
/* ================================================================ */
typedef struct { uint8_t *d; size_t cap, bp; int bit; } bs_t;
static void bs_init(bs_t *b, uint8_t *d, size_t c) { b->d=d; b->cap=c; b->bp=0; b->bit=7; memset(d,0,c); }
static void bs_w1(bs_t *b, int v) {
    if(b->bp>=b->cap)return;
    if(v) b->d[b->bp]|=(1<<b->bit);
    if(--b->bit<0){b->bit=7; b->bp++;}
}
static void bs_wn(bs_t *b, uint32_t v, int n){for(int i=n-1;i>=0;i--) bs_w1(b,(v>>i)&1);}
static size_t bs_sz(bs_t *b){return b->bp+(b->bit<7?1:0);}
static void bs_rinit(bs_t *b, uint8_t *d, size_t c){b->d=d;b->cap=c;b->bp=0;b->bit=7;}
static int bs_r1(bs_t *b){
    if(b->bp>=b->cap)return 0;
    int v=(b->d[b->bp]>>b->bit)&1;
    if(--b->bit<0){b->bit=7;b->bp++;}
    return v;
}
static uint32_t bs_rn(bs_t *b, int n){uint32_t v=0;for(int i=0;i<n;i++) v=(v<<1)|bs_r1(b);return v;}

/* ================================================================ */
/* GOLOMB-RICE                                                       */
/* ================================================================ */
static int zze(int v){return v>=0?2*v:2*(-v)-1;}
static int zzd(int v){return(v&1)?-((v+1)/2):v/2;}

static void gr_enc(bs_t *b, int val, int k) {
    int v=zze(val), q=v>>k;
    if(q>24){/* escape: unary 25 + raw 16 bits */
        for(int i=0;i<25;i++) bs_w1(b,1); bs_w1(b,0); bs_wn(b,v,16);
    } else {
        for(int i=0;i<q;i++) bs_w1(b,1); bs_w1(b,0); bs_wn(b,v&((1<<k)-1),k);
    }
}
static int gr_dec(bs_t *b, int k) {
    int q=0; while(bs_r1(b)==1&&q<25) q++;
    int v; if(q>=25) v=bs_rn(b,16); else v=(q<<k)|bs_rn(b,k);
    return zzd(v);
}

/* Run-length: encode count using Golomb-Rice with k=3 */
static void run_enc(bs_t *b, int count) { gr_enc(b, count-1, 3); }
static int run_dec(bs_t *b) { return gr_dec(b, 3) + 1; }

/* ================================================================ */
/* MED PREDICTOR                                                     */
/* ================================================================ */
static inline int med(int a,int b,int c){
    int mn=a<b?a:b,mx=a>b?a:b;
    if(c>=mx)return mn; if(c<=mn)return mx; return a+b-c;
}

/* ================================================================ */
/* 128 CONTEXTS: 8 magnitude levels × 4 direction quadrants × 4 bands */
/* ================================================================ */
#define NCTX 128

static int get_ctx(int a, int b, int c, int d) {
    int gh = a - c;    /* horizontal gradient (signed) */
    int gv = b - c;    /* vertical gradient (signed) */
    int gmag = abs(gh) + abs(gv) + abs(d - b);

    /* Magnitude: 8 levels */
    int mag;
    if(gmag < 3) mag = 0;
    else if(gmag < 8) mag = 1;
    else if(gmag < 16) mag = 2;
    else if(gmag < 30) mag = 3;
    else if(gmag < 50) mag = 4;
    else if(gmag < 80) mag = 5;
    else if(gmag < 130) mag = 6;
    else mag = 7;

    /* Direction: 4 quadrants */
    int dir;
    if(abs(gh) < 3 && abs(gv) < 3) dir = 0;       /* flat */
    else if(abs(gh) > 2*abs(gv)) dir = 1;           /* horizontal edge */
    else if(abs(gv) > 2*abs(gh)) dir = 2;           /* vertical edge */
    else dir = 3;                                     /* diagonal */

    /* Band: position within row (4 bands) */
    /* Not using position — just mag × dir × 4 sub-contexts from gradient sign */
    int sign = (gh >= 0 ? 0 : 1) | ((gv >= 0 ? 0 : 1) << 1);

    return (mag << 4) | (dir << 2) | sign;
}

static int k_from(double avg) {
    if(avg < 0.5) return 0;
    if(avg < 1.0) return 1;
    int k = (int)(log2(avg) + 0.5);
    return k < 0 ? 0 : (k > 12 ? 12 : k);
}

/* ================================================================ */
/* ENCODE / DECODE ONE PLANE                                         */
/* ================================================================ */

static size_t enc_plane(const uint8_t *p, int h, int w, uint8_t *out, size_t cap) {
    double ca[NCTX]; for(int i=0;i<NCTX;i++) ca[i]=3.0;
    bs_t bs; bs_init(&bs, out, cap);

    int i = 0, n = h * w;
    while(i < n) {
        int r = i / w, c = i % w;
        int a=(c>0)?p[i-1]:0, b=(r>0)?p[i-w]:0, cv=(r>0&&c>0)?p[i-w-1]:0, d=(r>0&&c+1<w)?p[i-w+1]:b;
        int pred = med(a, b, cv);
        int res = (int)p[i] - pred;
        if(res > 128) res -= 256;
        if(res < -128) res += 256;

        if(res == 0) {
            /* Run-length mode: count consecutive zeros */
            int run = 0;
            int j = i;
            while(j < n) {
                int rj = j/w, cj = j%w;
                int aj=(cj>0)?p[j-1]:0, bj=(rj>0)?p[j-w]:0, cvj=(rj>0&&cj>0)?p[j-w-1]:0;
                int pj = med(aj, bj, cvj);
                int rr = (int)p[j] - pj;
                if(rr > 128) rr -= 256;
                if(rr < -128) rr += 256;
                if(rr != 0) break;
                run++; j++;
            }
            /* Signal: bit 1 = run follows */
            bs_w1(&bs, 1);
            run_enc(&bs, run);
            i += run;
        } else {
            /* Signal: bit 0 = non-zero residual */
            bs_w1(&bs, 0);
            int cx = get_ctx(a, b, cv, d);
            if(cx >= NCTX) cx = NCTX - 1;
            int k = k_from(ca[cx]);
            gr_enc(&bs, res, k);
            ca[cx] = 0.875 * ca[cx] + 0.125 * abs(res);
            i++;
        }
    }
    return bs_sz(&bs);
}

static void dec_plane(const uint8_t *data, size_t len, uint8_t *p, int h, int w) {
    double ca[NCTX]; for(int i=0;i<NCTX;i++) ca[i]=3.0;
    bs_t bs; bs_rinit(&bs, (uint8_t*)data, len);

    int i = 0, n = h * w;
    memset(p, 0, n);
    while(i < n) {
        int flag = bs_r1(&bs);
        if(flag) {
            /* Run of zeros */
            int run = run_dec(&bs);
            for(int j = 0; j < run && i < n; j++, i++) {
                int r = i/w, c = i%w;
                int a=(c>0)?p[i-1]:0, b=(r>0)?p[i-w]:0, cv=(r>0&&c>0)?p[i-w-1]:0;
                p[i] = (uint8_t)(med(a, b, cv) & 0xFF);
            }
        } else {
            int r = i/w, c = i%w;
            int a=(c>0)?p[i-1]:0, b=(r>0)?p[i-w]:0, cv=(r>0&&c>0)?p[i-w-1]:0, d=(r>0&&c+1<w)?p[i-w+1]:b;
            int cx = get_ctx(a, b, cv, d);
            if(cx >= NCTX) cx = NCTX - 1;
            int k = k_from(ca[cx]);
            int res = gr_dec(&bs, k);
            p[i] = (uint8_t)((med(a, b, cv) + res) & 0xFF);
            ca[cx] = 0.875 * ca[cx] + 0.125 * abs(res);
            i++;
        }
    }
}

/* ================================================================ */
/* COLOR TRANSFORM                                                    */
/* ================================================================ */
static void to_grd(const uint8_t *rgb, uint8_t *g, uint8_t *rg, uint8_t *bg, int n){
    for(int i=0;i<n;i++){g[i]=rgb[3*i+1];rg[i]=(rgb[3*i]-rgb[3*i+1])&0xFF;bg[i]=(rgb[3*i+2]-rgb[3*i+1])&0xFF;}
}
static void from_grd(const uint8_t *g, const uint8_t *rg, const uint8_t *bg, uint8_t *rgb, int n){
    for(int i=0;i<n;i++){rgb[3*i+1]=g[i];rgb[3*i]=(rg[i]+g[i])&0xFF;rgb[3*i+2]=(bg[i]+g[i])&0xFF;}
}

/* ================================================================ */
/* FULL ENCODE / DECODE                                               */
/* ================================================================ */
long tic_enc(const uint8_t *rgb, int w, int h, uint8_t *out, size_t cap){
    int n=w*h; uint8_t *g=malloc(n),*rg=malloc(n),*bg=malloc(n);
    size_t bc=n+n; uint8_t *buf=malloc(bc);
    to_grd(rgb,g,rg,bg,n);
    size_t pos=0;
    memcpy(out,"TGR3",4); pos+=4;
    out[pos++]=w&0xFF;out[pos++]=(w>>8)&0xFF;out[pos++]=h&0xFF;out[pos++]=(h>>8)&0xFF;out[pos++]=0;
    uint8_t *pl[3]={g,rg,bg};
    for(int p=0;p<3;p++){
        size_t sz=enc_plane(pl[p],h,w,buf,bc);
        out[pos++]=sz&0xFF;out[pos++]=(sz>>8)&0xFF;out[pos++]=(sz>>16)&0xFF;out[pos++]=(sz>>24)&0xFF;
        memcpy(out+pos,buf,sz);pos+=sz;
    }
    free(g);free(rg);free(bg);free(buf);return(long)pos;
}

int tic_dec(const uint8_t *data, size_t len, int w, int h, uint8_t *rgb){
    int n=w*h; uint8_t *g=malloc(n),*rg=malloc(n),*bg=malloc(n);
    size_t pos=9;
    uint8_t *pl[3]={g,rg,bg};
    for(int p=0;p<3;p++){
        uint32_t sz=data[pos]|(data[pos+1]<<8)|(data[pos+2]<<16)|(data[pos+3]<<24);pos+=4;
        dec_plane(data+pos,sz,pl[p],h,w);pos+=sz;
    }
    from_grd(g,rg,bg,rgb,n); free(g);free(rg);free(bg); return 0;
}

/* ================================================================ */
/* MAIN                                                               */
/* ================================================================ */
int main(int argc, char **argv){
    if(argc<5){printf("Usage: %s bench input.raw W H [iters]\n",argv[0]);return 1;}
    int w=atoi(argv[3]),h=atoi(argv[4]),iters=(argc>=6)?atoi(argv[5]):50;
    int n=w*h*3;
    uint8_t *rgb=malloc(n);
    FILE *f=fopen(argv[2],"rb");if(!f){printf("Cannot open\n");return 1;}
    if(fread(rgb,1,n,f)!=(size_t)n){printf("Short\n");fclose(f);return 1;}fclose(f);
    size_t cap=n*2;uint8_t *out=malloc(cap);
    long sz=tic_enc(rgb,w,h,out,cap);
    clock_t t0=clock();
    for(int i=0;i<iters;i++) sz=tic_enc(rgb,w,h,out,cap);
    double et=(double)(clock()-t0)/CLOCKS_PER_SEC/iters;
    uint8_t *dec=malloc(n);
    t0=clock();
    for(int i=0;i<iters;i++) tic_dec(out,sz,w,h,dec);
    double dt=(double)(clock()-t0)/CLOCKS_PER_SEC/iters;
    int ok=(memcmp(rgb,dec,n)==0);
    printf("TIC-GR3: %dx%d, %d iters\n",w,h,iters);
    printf("  Raw:        %d\n",n);
    printf("  Compressed: %ld (%.2f:1)\n",sz,(double)n/sz);
    printf("  Encode:     %.3f ms (%.1f Mpix/s)\n",et*1000,(double)(w*h)/et/1e6);
    printf("  Decode:     %.3f ms (%.1f Mpix/s)\n",dt*1000,(double)(w*h)/dt/1e6);
    printf("  Roundtrip:  %s\n",ok?"PASS":"FAIL");
    free(rgb);free(out);free(dec);return ok?0:1;
}
