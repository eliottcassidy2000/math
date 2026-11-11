// scrace.c - the mod-4 / mod-8 race of the Hamiltonian-path count over tournament classes of order n,
// split by self-converse (T isomorphic to its converse) and by strong connectivity, with labelled weights n!/|Aut|.
// stdin: gentourng output (upper triangle row by row, '1' = i->j for i<j).  Usage: scrace n
// Output: lines "sc strong r classes labelled" for r = H mod 8 in {1,3,5,7}; then totals on stderr.
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64;
#define MAXN 11
static int n; static u64 f[1<<MAXN][MAXN]; static int outm[MAXN], inm[MAXN], sc[MAXN];
static int perm[MAXN], used; static u64 autcount; static int found;
static void autrec(int k){
  if(k==n){ autcount++; return; }
  for(int c=0;c<n;c++){ if(used>>c&1) continue; if(sc[c]!=sc[k]) continue; int ok=1;
    for(int j=0;j<k;j++){ int a=(outm[k]>>j)&1, b=(outm[c]>>perm[j])&1; if(a!=b){ok=0;break;} }
    if(!ok) continue; perm[k]=c; used|=1<<c; autrec(k+1); used&=~(1<<c); }
}
static void antirec(int k){   // anti-automorphism: k->j  <=>  phi(j)->phi(k)
  if(found) return; if(k==n){ found=1; return; }
  for(int c=0;c<n;c++){ if(used>>c&1) continue; if(sc[c]!=n-1-sc[k]) continue; int ok=1;
    for(int j=0;j<k;j++){ int a=(outm[k]>>j)&1, b=(outm[perm[j]]>>c)&1; if(a!=b){ok=0;break;} }
    if(!ok) continue; perm[k]=c; used|=1<<c; antirec(k+1); used&=~(1<<c); if(found) return; }
}
int main(int argc,char**argv){ n=atoi(argv[1]); int full=(1<<n)-1; char line[256];
  u64 fact=1; for(int i=2;i<=n;i++) fact*=i;
  static u64 cls[2][2][8], lab[2][2][8]; u64 nt=0, bad=0;
  while(fgets(line,sizeof line,stdin)){ if(line[0]!='0'&&line[0]!='1')continue; int k=0;
    for(int i=0;i<n;i++){outm[i]=0;inm[i]=0;}
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){ if(line[k++]=='1'){outm[i]|=1<<j;inm[j]|=1<<i;} else {outm[j]|=1<<i;inm[i]|=1<<j;} }
    nt++;
    for(int S=1;S<=full;S++) for(int v=0;v<n;v++){ if(!(S>>v&1)){f[S][v]=0;continue;} int R=S&~(1<<v);
        if(R==0){f[S][v]=1;continue;} u64 a=0; int mm=R&inm[v]; while(mm){int w=__builtin_ctz(mm);mm&=mm-1;a+=f[R][w];} f[S][v]=a; }
    u64 H=0; for(int v=0;v<n;v++) H+=f[full][v];
    for(int i=0;i<n;i++) sc[i]=__builtin_popcount(outm[i]);
    int s[MAXN]; for(int i=0;i<n;i++) s[i]=sc[i];
    for(int i=1;i<n;i++){ int t=s[i], j=i-1; while(j>=0&&s[j]>t){s[j+1]=s[j];j--;} s[j+1]=t; }
    int strong=1, acc=0; for(int kk=1;kk<n;kk++){ acc+=s[kk-1]; if(acc<=kk*(kk-1)/2){strong=0;break;} }
    autcount=0; used=0; autrec(0); if(autcount==0||fact%autcount){bad++;continue;}
    found=0; used=0; antirec(0);
    cls[found][strong][H&7]++; lab[found][strong][H&7]+=fact/autcount;
  }
  for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int r=1;r<8;r+=2) printf("%d %d %d %llu %llu\n", a, b, r, cls[a][b][r], lab[a][b][r]);
  fprintf(stderr,"order %d: classes %llu bad %llu\n", n, nt, bad);
  return 0; }
