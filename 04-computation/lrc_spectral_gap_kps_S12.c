/* lrc_spectral_gap_kps_S12.c -- kind-pasteur-2026-07-05-S12, HYP-4147.
 *
 * The SPECTRAL-GAP structure of the loose branch.  Exhaustively enumerate primitive
 * compressed covering 12-bases in [1,B], compute EXACT M(base) = max_{a/q} min_k
 * ||v_k a/q|| (merge-point scan q<=2*max, integer arithmetic), and report:
 *   - non-AP FLOORS: M == 1/13 but base is not a dilated AP  (dichotomy needs NONE)
 *   - GAP inhabitants: 1/13 < M < 2/25                        (dichotomy needs NONE)
 *   - the low spectrum (smallest M values)                    (what is sigma_2?)
 * Farey: 1/13 and 2/25 are neighbors (1*25-2*13=-1), so any gap value has denom>=38.
 *
 * gcc -O2 -o sg lrc_spectral_gap_kps_S12.c && ./sg 28
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int B=28;
static int gcd_i(int a,int b){while(b){int t=a%b;a=b;b=t;}return a;}
static int V[12];
static long long n_base=0, n_floor_nonAP=0, n_gap=0;
/* low spectrum: track a few smallest M as (num,den) */
#define NS 8
static int lowN[NS], lowD[NS]; static long long lowCnt[NS]; static int nlow=0;

static int dist_q(int x,int q){int r=x%q; if(r<0)r+=q; return r<q-r?r:q-r;}

/* exact M as reduced-ish fraction (bestN/bestD), bestD<=2*max */
static void M_exact(int *bn,int *bd){
    int mx=0; for(int i=0;i<12;i++) if(V[i]>mx) mx=V[i];
    int Q=2*mx; int BN=0,BD=1;
    for(int q=2;q<=Q;q++){
        for(int a=1;a<q;a++){
            if(gcd_i(a,q)!=1) continue;
            int mm=q; /* min dist over runners */
            for(int i=0;i<12;i++){int d=dist_q(V[i]*a,q); if(d<mm)mm=d;}
            /* compare mm/q > BN/BD : mm*BD > BN*q */
            if((long long)mm*BD > (long long)BN*q){ BN=mm; BD=q; }
        }
    }
    *bn=BN; *bd=BD;
}
static int isAPdil(void){
    int s[12]; memcpy(s,V,sizeof s);
    for(int i=0;i<11;i++)for(int j=i+1;j<12;j++) if(s[j]<s[i]){int t=s[i];s[i]=s[j];s[j]=t;}
    int d=s[1]-s[0];
    for(int i=0;i<11;i++) if(s[i+1]-s[i]!=d) return 0;
    return s[0]==d; /* dilated AP: {d,2d,...,12d} */
}
static int covering(void){
    for(int q=2;q<=12;q++){int h=0;for(int i=0;i<12;i++)if(V[i]%q==0){h=1;break;}if(!h)return 0;}
    return 1;
}
static int compressed(void){
    for(int i=0;i<12;i++){int ok=0;for(int j=0;j<12;j++)if(j!=i&&V[i]<=13*V[j]){ok=1;break;}if(!ok)return 0;}
    return 1;
}
static int primitive(void){int g=0;for(int i=0;i<12;i++)g=gcd_i(g,V[i]);return g==1;}

static void record_low(int n,int d){
    /* insert n/d into sorted low list (ascending by value) */
    for(int k=0;k<nlow;k++){ if(n==lowN[k]&&d==lowD[k]){lowCnt[k]++;return;} }
    if(nlow<NS){lowN[nlow]=n;lowD[nlow]=d;lowCnt[nlow]=1;nlow++;}
    else{
        /* replace the largest if this is smaller */
        int worst=0; for(int k=1;k<nlow;k++) if((long long)lowN[k]*lowD[worst]>(long long)lowN[worst]*lowD[k]) worst=k;
        if((long long)n*lowD[worst] < (long long)lowN[worst]*d){lowN[worst]=n;lowD[worst]=d;lowCnt[worst]=1;}
    }
}

static void dfs(int pos,int start){
    if(pos==12){
        if(!covering()||!primitive()||!compressed()) return;
        n_base++;
        int bn,bd; M_exact(&bn,&bd);
        /* M ? 1/13 : bn*13 == bd ; gap: bn*13>bd && bn*25<2*bd */
        int isfloor = ((long long)bn*13==(long long)bd);
        int ingap = ((long long)bn*13>(long long)bd) && ((long long)bn*25<2LL*bd);
        if(isfloor && !isAPdil()){ n_floor_nonAP++; if(n_floor_nonAP<=8){printf("FLOOR-nonAP");for(int i=0;i<12;i++)printf(" %d",V[i]);printf("  M=%d/%d\n",bn,bd);} }
        if(ingap){ n_gap++; if(n_gap<=20){printf("GAP");for(int i=0;i<12;i++)printf(" %d",V[i]);printf("  M=%d/%d=%.6f\n",bn,bd,(double)bn/bd);} }
        /* print sigma_2 attainers: M == 2/25 */
        if((long long)bn*25==2LL*bd){ static long long s2=0; s2++; if(s2<=12){printf("SIGMA2");for(int i=0;i<12;i++)printf(" %d",V[i]);printf("  M=%d/%d\n",bn,bd);} }
        record_low(bn,bd);
        return;
    }
    for(int v=start;v<=B-(11-pos);v++){ V[pos]=v; dfs(pos+1,v+1); }
}
int main(int argc,char**argv){
    if(argc>1) B=atoi(argv[1]);
    dfs(0,1);
    printf("=== SPECTRAL GAP CENSUS [1,%d] ===\n",B);
    printf("primitive compressed covering 12-bases: %lld\n",n_base);
    printf("non-AP FLOORS (M=1/13): %lld\n",n_floor_nonAP);
    printf("GAP inhabitants (1/13<M<2/25): %lld\n",n_gap);
    /* sort low list ascending */
    for(int i=0;i<nlow;i++)for(int j=i+1;j<nlow;j++) if((long long)lowN[j]*lowD[i]<(long long)lowN[i]*lowD[j]){int t;long long c;t=lowN[i];lowN[i]=lowN[j];lowN[j]=t;t=lowD[i];lowD[i]=lowD[j];lowD[j]=t;c=lowCnt[i];lowCnt[i]=lowCnt[j];lowCnt[j]=c;}
    printf("LOW SPECTRUM (smallest M values: value=count):\n");
    for(int i=0;i<nlow;i++) printf("   %d/%d = %.6f  (count %lld)\n",lowN[i],lowD[i],(double)lowN[i]/lowD[i],lowCnt[i]);
    printf("   [1/13=%.6f  2/25=%.6f]\n",1.0/13,2.0/25);
    return 0;
}
