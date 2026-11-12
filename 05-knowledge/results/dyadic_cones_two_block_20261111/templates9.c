/* Blow-up templates. For a tournament K (order n <= 8) and vertices x (and y), N(a) = H(K[TT_a at x]) and
   N(a,b) = H(K[TT_a at x, TT_b at y]). If no Hamiltonian walk visits x (resp. y) three times these are
   g0 + g1 2^a   and   g00 + g10 2^a + g01 2^b + g11 2^(a+b).
   Input (stdin): one tournament per line: n out-masks (n given as argv[1]).
   Output: valid single templates (offset g0, multiplier g1) summary and all valid pair templates with g00 = 1
   (argv[2] = 1: print all pair templates with odd g00 <= argv[3]). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
static int n, K[10];
static ll dp[1<<13][13];
static ll Hblow(int x,int a,int y,int b){
    /* vertices: originals except x,y keep; x replaced by a copies, y by b copies (y<0: none) */
    int map[13], cp[13], m=0;            /* map[i]=original vertex, cp[i]=copy index */
    for(int v=0;v<n;v++){ int c = (v==x)?a:((v==y)?b:1); for(int t=0;t<c;t++){ map[m]=v; cp[m]=t; m++; } }
    int out[13];
    for(int i=0;i<m;i++){ out[i]=0; for(int j=0;j<m;j++){ if(i==j) continue; int u=map[i], w=map[j];
        if(u==w){ if(cp[i]<cp[j]) out[i]|=1<<j; } else if((K[u]>>w)&1) out[i]|=1<<j; } }
    int full=(1<<m)-1;
    for(int S=1;S<=full;S++) memset(dp[S],0,sizeof(ll)*m);
    for(int v=0;v<m;v++) dp[1<<v][v]=1;
    for(int S=1;S<=full;S++){ ll*r=dp[S]; for(int v=0;v<m;v++){ ll c=r[v]; if(!c) continue; int t=out[v]&~S; while(t){int bb=t&-t; dp[S|bb][__builtin_ctz(bb)]+=c; t^=bb;} } }
    ll h=0; for(int v=0;v<m;v++) h+=dp[full][v]; return h;
}
int main(int argc,char**argv){
    n=atoi(argv[1]); int mode=argc>2?atoi(argv[2]):0; ll lim=argc>3?atoll(argv[3]):1;
    long nK=0, nsingle=0, nunit=0, npair=0, nunitpair=0; char line[512];
    while(fgets(line,sizeof(line),stdin)){
        char *p=line; int ok=1; for(int i=0;i<n;i++){ char *e; long v=strtol(p,&e,10); if(e==p){ok=0;break;} K[i]=(int)v; p=e; } if(!ok) continue; nK++;
        ll H1=Hblow(-1,1,-1,1); int valid[10]; ll g0[10], g1[10];
        for(int x=0;x<n;x++){
            ll N2=Hblow(x,2,-1,1), N3=Hblow(x,3,-1,1);
            ll a1=(N2-H1)/2; ll a0=H1-2*a1;            /* N(a) = a0 + a1 2^a */
            valid[x] = ((N2-H1)%2==0) && (a0+8*a1==N3);
            if(valid[x]){ ll N4=Hblow(x,4,-1,1); if(a0+16*a1!=N4) valid[x]=0; }
            g0[x]=a0; g1[x]=a1;
            if(valid[x]){ nsingle++; if(a0==1){ nunit++; } }
        }
        if(mode!=3) for(int x=0;x<n;x++) if(valid[x]) for(int y=x+1;y<n;y++) if(valid[y]){
            ll N22=Hblow(x,2,y,2);
            /* N(a,b)=c00+c10 u+c01 v+c11 uv, u=2^a, v=2^b. N(1,1)=H1, N(2,1)=g0x+4g1x, N(1,2)=g0y+4g1y */
            ll N21=g0[x]+4*g1[x], N12=g0[y]+4*g1[y];
            ll c11n=N22-N21-N12+H1;                  /* = c11*(4-2)*(4-2) = 4 c11 */
            if(c11n%4) continue; ll c11=c11n/4;
            ll c10=g1[x]-2*c11, c01=g1[y]-2*c11; ll c00=H1-2*c10-2*c01-4*c11;
            /* validate */
            ll N33=Hblow(x,3,y,3); if(c00+8*c10+8*c01+64*c11!=N33) continue;
            if(n<=7){ ll N42=Hblow(x,4,y,2); if(c00+16*c10+4*c01+64*c11!=N42) continue; ll N24=Hblow(x,2,y,4); if(c00+4*c10+16*c01+64*c11!=N24) continue; }
            npair++;
            if(c00==1) nunitpair++;
            if((mode==0&&c00==1)||(mode==1&&c00<=lim&&c00>=-lim)){
                printf("T n=%d H=%lld x=%d y=%d  P = %lld + %lld u + %lld v + %lld uv   K:",n,H1,x,y,c00,c10,c01,c11);
                for(int i=0;i<n;i++) printf(" %d",K[i]); printf("\n");
            }
        }
        if(mode==3) for(int x=0;x<n;x++) if(valid[x]&&g1[x]>0&&(g1[x]&(g1[x]-1))==0){ printf("C n=%d H=%lld x=%d  N(a) = %lld + %lld 2^a  K:",n,H1,x,g0[x],g1[x]); for(int i=0;i<n;i++) printf(" %d",K[i]); printf("\n"); }
        for(int x=0;x<n;x++) if(valid[x]&&g0[x]==1&&mode==2){ printf("S n=%d H=%lld x=%d  N(a) = 1 + %lld 2^a  K:",n,H1,x,g1[x]); for(int i=0;i<n;i++) printf(" %d",K[i]); printf("\n"); }
    }
    printf("# n=%d tournaments=%ld valid single templates=%ld (unit offset: %ld) valid pair templates=%ld (g00=1: %ld)\n",n,nK,nsingle,nunit,npair,nunitpair);
    return 0;
}
