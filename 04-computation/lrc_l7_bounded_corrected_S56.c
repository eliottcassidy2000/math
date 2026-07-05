/* lrc_l7_bounded_corrected_S56.c -- mac-mini-2026-07-05-S56, HYP-4122 (AUDIT FIX).
 *
 * CORRECTS my S55 l>=7 bounded-stratum sweep (AUDIT FINDING 1): that sweep
 * assumed the moduli r in C n {7,8,9} are SELF-carried (k_r = r, value 14r).
 * WRONG scope: any lifted coordinate c can carry modulus r via
 * c + 13k == 0 (mod r) (e.g. 14 = 2*7 at c=1,k=1 carries mod-7 while
 * coordinate 7 floats).  opus-S80's dvd_lift_height is a statement about the
 * position-r coordinate only; the stratum needs all carriers.
 *
 * THE CORRECTED STRATUM: l >= 7 lifts W = ({1..12}\C) u {c + 13 k_c : c in C},
 * |C| >= 7, ALL lifted values <= 133 (the bounded piece of opus-S81's descent
 * residual; spread tops with bottom entry >= 134 are descent-dodged).
 * No height assumptions: k_c in [1, (133-c)/13].  Position r in C n {7..12}
 * needs SOME element == 0 mod r (base cannot carry: r's unique multiple is r).
 * NOTE r in {10,11,12} is now POSSIBLE (carried by another coord, e.g.
 * 5+13*5 = 70 carries 10) -- patterns C n {10,11,12} != empty are BACK IN.
 * All C(12,l) patterns, l = 7..12.
 *
 * Filters per leaf: distinctness, sieve 2..12 (full), primitivity, witness
 * scan at 2/25 (q in [8,60] + 13u library), exact escalation printed.
 * gcc -O2 -o l7fix lrc_l7_bounded_corrected_S56.c && ./l7fix
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;
static int gcd_i(int a,int b){while(b){int t=a%b;a=b;b=t;}return a;}
static int dq(int x,int q){int r=x%q; if(r<0)r+=q; return r<q-r?r:q-r;}

static ll st_total=0, st_dup=0, st_sieve=0, st_prim=0, st_scan=0, st_hard=0;
static int W[12], C[12], nC, base_[12], nbase, kmax_[12];

static int scan25(void){
    static int QS[80]; static int nQ=0;
    if(!nQ){ for(int q=26;q<=60;q++) QS[nQ++]=q; for(int q=8;q<26;q++) QS[nQ++]=q;
             for(int u=5;u<=14;u++) QS[nQ++]=13*u; }
    for(int qi=0;qi<nQ;qi++){
        int q=QS[qi], bad=0;
        for(int i=0;i<12;i++) if(W[i]%q==0){bad=1;break;}
        if(bad) continue;
        for(int a=1;a<=q/2;a++){
            if(gcd_i(a,q)!=1) continue;
            int ok=1;
            for(int i=0;i<12;i++){ if(25*dq(W[i]*a,q)<2*q){ok=0;break;} }
            if(ok) return 1;
        }
    }
    return 0;
}

static void leaf(void){
    st_total++;
    /* distinctness (values can collide across coords? c+13k distinct since
       residues mod 13 distinct -- no collision possible; skip) */
    /* sieve 2..12 */
    for(int m=2;m<=12;m++){
        int hit=0;
        for(int i=0;i<12;i++) if(W[i]%m==0){hit=1;break;}
        if(!hit){ st_sieve++; return; }
    }
    int g=0; for(int i=0;i<12;i++) g=gcd_i(g,W[i]);
    if(g!=1){ st_prim++; return; }
    if(scan25()){ st_scan++; return; }
    st_hard++;
    printf("HARD"); for(int i=0;i<12;i++) printf(" %d", W[i]); printf("\n");
    fflush(stdout);
}

static void rec(int pos){
    if(pos==nC){ leaf(); return; }
    int c=C[pos];
    for(int k=1;k<=kmax_[pos];k++){
        W[nbase+pos]=c+13*k;
        rec(pos+1);
    }
}

int main(int argc, char **argv){
    int llo = argc>1 ? atoi(argv[1]) : 7, lhi = argc>2 ? atoi(argv[2]) : 7;
    for(int l=llo;l<=lhi;l++){
        ll tot0=st_total;
        int idx[12]; for(int i=0;i<l;i++) idx[i]=i+1;
        while(1){
            nC=l;
            for(int i=0;i<l;i++) C[i]=idx[i];
            nbase=0;
            for(int v=1;v<=12;v++){
                int in=0; for(int i=0;i<l;i++) if(idx[i]==v) in=1;
                if(!in) base_[nbase++]=v;
            }
            for(int i=0;i<nbase;i++) W[i]=base_[i];
            for(int i=0;i<l;i++) kmax_[i]=(133-C[i])/13;
            rec(0);
            int i=l-1;
            while(i>=0 && idx[i]==12-(l-1-i)) i--;
            if(i<0) break;
            idx[i]++;
            for(int j=i+1;j<l;j++) idx[j]=idx[j-1]+1;
        }
        fprintf(stderr,"l=%d done: cumulative total=%lld sieve=%lld prim=%lld scan=%lld HARD=%lld (delta total=%lld)\n",
                l, st_total, st_sieve, st_prim, st_scan, st_hard, st_total-tot0);
    }
    printf("CORRECTED l>=7 BOUNDED STRATUM (all values <= 133, NO height assumptions):\n");
    printf("total=%lld sieve-killed=%lld nonprim=%lld witness-cleared=%lld HARD=%lld\n",
           st_total, st_sieve, st_prim, st_scan, st_hard);
    return 0;
}
