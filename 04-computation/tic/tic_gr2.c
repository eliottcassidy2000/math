/*
 * tic_gr2.c — TIC with context-adaptive Golomb-Rice (no run-length)
 *
 * Simple, correct, no zlib dependency.
 * Pipeline: G-RG-BG → MED → 64-context adaptive GR
 *
 * kind-pasteur-2026-03-26-S32
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

/* Bitstream */
typedef struct { uint8_t *d; size_t cap, bp; int bit; } bs_t;
static void bs_init(bs_t *b, uint8_t *d, size_t c) { b->d=d; b->cap=c; b->bp=0; b->bit=7; memset(d,0,c); }
static void bs_wbit(bs_t *b, int v) {
    if(b->bp>=b->cap) return;
    if(v) b->d[b->bp]|=(1<<b->bit);
    if(--b->bit<0){b->bit=7; b->bp++;}
}
static void bs_wbits(bs_t *b, uint32_t v, int n) { for(int i=n-1;i>=0;i--) bs_wbit(b,(v>>i)&1); }
static size_t bs_sz(bs_t *b) { return b->bp+(b->bit<7?1:0); }
static int bs_rbit(bs_t *b) {
    if(b->bp>=b->cap) return 0;
    int v=(b->d[b->bp]>>b->bit)&1;
    if(--b->bit<0){b->bit=7; b->bp++;}
    return v;
}
static uint32_t bs_rbits(bs_t *b, int n) { uint32_t v=0; for(int i=0;i<n;i++) v=(v<<1)|bs_rbit(b); return v; }

/* Zigzag */
static int zze(int v){return v>=0?2*v:2*(-v)-1;}
static int zzd(int v){return(v&1)?-((v+1)/2):v/2;}

/* GR encode/decode */
static void gr_enc(bs_t *b, int val, int k) {
    int v=zze(val), q=v>>k, r=v&((1<<k)-1);
    if(q>30){for(int i=0;i<31;i++) bs_wbit(b,1); bs_wbit(b,0); bs_wbits(b,v,16);}
    else{for(int i=0;i<q;i++) bs_wbit(b,1); bs_wbit(b,0); bs_wbits(b,r,k);}
}
static int gr_dec(bs_t *b, int k) {
    int q=0; while(bs_rbit(b)==1&&q<31) q++;
    int v; if(q>=31) v=bs_rbits(b,16); else v=(q<<k)|bs_rbits(b,k);
    return zzd(v);
}

/* MED */
static inline int med(int a,int b,int c){int mn=a<b?a:b,mx=a>b?a:b; if(c>=mx)return mn; if(c<=mn)return mx; return a+b-c;}

/* Context: gradient magnitude → 0..63 */
#define NCTX 64
static int ctx(int a,int b,int c,int d) {
    int g=abs(a-c)+abs(b-c)+abs(d-b);
    if(g<3)return 0; if(g<8)return 1; if(g<15)return 2; if(g<25)return 3;
    if(g<40)return 4+(g-25)/5; if(g<100)return 7+(g-40)/10;
    int r=13+(g-100)/20; return r<NCTX?r:NCTX-1;
}
static int k_from(double a){if(a<0.5)return 0; int k=(int)(log2(a/0.693)+0.5); return k<0?0:(k>12?12:k);}

/* Encode one plane */
static size_t enc_plane(const uint8_t *p, int h, int w, uint8_t *out, size_t cap) {
    double ca[NCTX]; for(int i=0;i<NCTX;i++) ca[i]=4.0;
    bs_t bs; bs_init(&bs,out,cap);
    for(int r=0;r<h;r++) for(int c=0;c<w;c++) {
        int i=r*w+c;
        int a=(c>0)?p[i-1]:0, b=(r>0)?p[i-w]:0, cv=(r>0&&c>0)?p[i-w-1]:0, d=(r>0&&c+1<w)?p[i-w+1]:b;
        int pred=med(a,b,cv), res=(int)p[i]-pred;
        if(res>128)res-=256; if(res<-128)res+=256;
        int cx=ctx(a,b,cv,d); int k=k_from(ca[cx]);
        gr_enc(&bs,res,k);
        ca[cx]=0.9*ca[cx]+0.1*abs(res);
    }
    return bs_sz(&bs);
}

/* Decode one plane */
static void dec_plane(const uint8_t *data, size_t len, uint8_t *p, int h, int w) {
    double ca[NCTX]; for(int i=0;i<NCTX;i++) ca[i]=4.0;
    bs_t bs; bs.d=(uint8_t*)data; bs.cap=len; bs.bp=0; bs.bit=7;
    for(int r=0;r<h;r++) for(int c=0;c<w;c++) {
        int i=r*w+c;
        int a=(c>0)?p[i-1]:0, b=(r>0)?p[i-w]:0, cv=(r>0&&c>0)?p[i-w-1]:0, d=(r>0&&c+1<w)?p[i-w+1]:b;
        int cx=ctx(a,b,cv,d); int k=k_from(ca[cx]);
        int res=gr_dec(&bs,k);
        p[i]=(uint8_t)((med(a,b,cv)+res)&0xFF);
        ca[cx]=0.9*ca[cx]+0.1*abs(res);
    }
}

/* Color transform */
static void to_grd(const uint8_t *rgb, uint8_t *g, uint8_t *rg, uint8_t *bg, int n) {
    for(int i=0;i<n;i++){g[i]=rgb[3*i+1]; rg[i]=(rgb[3*i]-rgb[3*i+1])&0xFF; bg[i]=(rgb[3*i+2]-rgb[3*i+1])&0xFF;}
}
static void from_grd(const uint8_t *g, const uint8_t *rg, const uint8_t *bg, uint8_t *rgb, int n) {
    for(int i=0;i<n;i++){rgb[3*i+1]=g[i]; rgb[3*i]=(rg[i]+g[i])&0xFF; rgb[3*i+2]=(bg[i]+g[i])&0xFF;}
}

/* Full encode */
long tic_enc(const uint8_t *rgb, int w, int h, uint8_t *out, size_t cap) {
    int n=w*h; uint8_t *g=malloc(n),*rg=malloc(n),*bg=malloc(n);
    size_t bc=n+n/2; uint8_t *buf=malloc(bc);
    to_grd(rgb,g,rg,bg,n);
    size_t pos=0;
    memcpy(out,"TGR2",4); pos+=4;
    out[pos++]=w&0xFF; out[pos++]=(w>>8)&0xFF;
    out[pos++]=h&0xFF; out[pos++]=(h>>8)&0xFF;
    out[pos++]=0;
    uint8_t *pl[3]={g,rg,bg};
    for(int p=0;p<3;p++){
        size_t sz=enc_plane(pl[p],h,w,buf,bc);
        out[pos++]=sz&0xFF; out[pos++]=(sz>>8)&0xFF; out[pos++]=(sz>>16)&0xFF; out[pos++]=(sz>>24)&0xFF;
        memcpy(out+pos,buf,sz); pos+=sz;
    }
    free(g);free(rg);free(bg);free(buf);
    return(long)pos;
}

/* Full decode */
int tic_dec(const uint8_t *data, size_t len, int w, int h, uint8_t *rgb) {
    int n=w*h; uint8_t *g=malloc(n),*rg=malloc(n),*bg=malloc(n);
    size_t pos=9;
    uint8_t *pl[3]={g,rg,bg};
    for(int p=0;p<3;p++){
        uint32_t sz=data[pos]|(data[pos+1]<<8)|(data[pos+2]<<16)|(data[pos+3]<<24); pos+=4;
        memset(pl[p],0,n);
        dec_plane(data+pos,sz,pl[p],h,w); pos+=sz;
    }
    from_grd(g,rg,bg,rgb,n);
    free(g);free(rg);free(bg);
    return 0;
}

int main(int argc, char **argv) {
    if(argc<5){printf("Usage: %s bench input.raw W H [iters]\n",argv[0]);return 1;}
    int w=atoi(argv[3]),h=atoi(argv[4]),iters=(argc>=6)?atoi(argv[5]):50;
    int n=w*h*3;
    uint8_t *rgb=malloc(n);
    FILE *f=fopen(argv[2],"rb"); if(!f){printf("Cannot open\n");return 1;}
    if(fread(rgb,1,n,f)!=(size_t)n){printf("Short\n");fclose(f);return 1;} fclose(f);
    size_t cap=n*2; uint8_t *out=malloc(cap);
    long sz=tic_enc(rgb,w,h,out,cap);
    clock_t t0=clock();
    for(int i=0;i<iters;i++) sz=tic_enc(rgb,w,h,out,cap);
    double et=(double)(clock()-t0)/CLOCKS_PER_SEC/iters;
    uint8_t *dec=malloc(n);
    t0=clock();
    for(int i=0;i<iters;i++) tic_dec(out,sz,w,h,dec);
    double dt=(double)(clock()-t0)/CLOCKS_PER_SEC/iters;
    int ok=(memcmp(rgb,dec,n)==0);
    printf("TIC-GR2: %dx%d, %d iters\n",w,h,iters);
    printf("  Raw:        %d\n",n);
    printf("  Compressed: %ld (%.2f:1)\n",sz,(double)n/sz);
    printf("  Encode:     %.3f ms (%.1f Mpix/s)\n",et*1000,(double)(w*h)/et/1e6);
    printf("  Decode:     %.3f ms (%.1f Mpix/s)\n",dt*1000,(double)(w*h)/dt/1e6);
    printf("  Roundtrip:  %s\n",ok?"PASS":"FAIL");
    free(rgb);free(out);free(dec);
    return ok?0:1;
}
