/* Literal recount of the cone witnesses. Reads lines "r=R j=J m=M O=mask B: o_0 ... o_{M-1}" (output of
   cone_census), builds L = B + v with N+(v) = O, appends the tail 0^t (t = 1, 2, 3: new vertices z_1 -> z_2 -> ...,
   every vertex of B beats every z, every z beats v) and counts Hamiltonian paths by a subset DP written
   independently of the census (no use of U, W or the offset formula). Checks H(L o 0^t) = r + c 2^(j+t) (c = 1 when the line has no "c=" field). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static ll dp[1<<14][14];
static ll hampaths(int n,const int*out){ int full=(1<<n)-1; for(int S=1;S<=full;S++) memset(dp[S],0,sizeof(ll)*n);
    for(int v=0;v<n;v++) dp[1<<v][v]=1;
    for(int S=1;S<=full;S++) for(int v=0;v<n;v++){ ll c=dp[S][v]; if(!c) continue; int t=out[v]&~S; while(t){int b=t&-t; dp[S|b][__builtin_ctz(b)]+=c; t^=b;} }
    ll h=0; for(int v=0;v<n;v++) h+=dp[full][v]; return h; }
int main(void){ char line[512]; long ok=0,bad=0,lines=0;
    while(fgets(line,sizeof(line),stdin)){ ll r; int j,m,O; int B[10];
        if(line[0]!='r') continue; char*p=strstr(line,"B:"); if(!p) continue;
        ll c=1; if(sscanf(line,"r=%lld j=%d m=%d O=%d",&r,&j,&m,&O)!=4) continue; { char*q=strstr(line,"c="); if(q) c=atoll(q+2); } p+=2; for(int i=0;i<m;i++){ char*e; B[i]=(int)strtol(p,&e,10); p=e; }
        lines++;
        for(int t=1;t<=3;t++){ int n=m+1+t; if(n>14) break; int out[14]; int v=m;
            for(int i=0;i<m;i++){ out[i]=B[i]; if(!((O>>i)&1)) out[i]|=1<<v; for(int z=0;z<t;z++) out[i]|=1<<(m+1+z); }
            out[v]=O;
            for(int z=0;z<t;z++){ out[m+1+z]=1<<v; for(int z2=z+1;z2<t;z2++) out[m+1+z]|=1<<(m+1+z2); }
            ll h=hampaths(n,out); ll want=r+c*((ll)1<<(j+t));
            if(h==want) ok++; else { bad++; if(bad<5) printf("MISMATCH r=%lld j=%d t=%d: counted %lld, formula %lld\n",r,j,t,h,want); } }
    }
    printf("witnesses %ld ; literal checks passed %ld ; failed %ld\n",lines,ok,bad); return bad!=0; }
