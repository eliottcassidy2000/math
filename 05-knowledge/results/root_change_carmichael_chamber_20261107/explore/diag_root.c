/* diag_root.c - kind-pasteur S14, exploratory. Root-tail (K = 0) word multiplicity M_0(N) = #{m >= 0 : H_t(m) = N},
 * H_t(m) = 1 + 2m - Q(m), at the two-run K = 0 chamber values: O^a Z^a gives 2 + (2^(a+1)-1)^2 (Session200's
 * Mersenne diagonal) and Z^a O^a gives (2^a+1)^2. H_t(m) >= m + 2, so all preimages satisfy m < N (exhaustive). */
#include <stdio.h>
typedef long long ll;
static ll T[32]; static int cnt[32]; static int nt; static ll X;
static void dfs(ll m, ll q){ ll h = 1 + 2*m - q; if (m > X) return;
  for (int i=0;i<nt;i++) if (h==T[i]) cnt[i]++;
  if (m==0){ dfs(1,0); return; }
  dfs(2*m,q); dfs(2*m+1,q+m+(m&1)); }
int main(void){ nt=0;
  for(int a=1;a<=13;a++){ ll r=(1LL<<(a+1))-1; T[nt++]=2+r*r; }
  for(int a=1;a<=13;a++){ ll f=(1LL<<a)+1; T[nt++]=f*f; }
  X=0; for(int i=0;i<nt;i++) if(T[i]>X) X=T[i];
  dfs(0,0);
  for(int i=0;i<13;i++) printf("O^aZ^a a=%d N=2+(2^(a+1)-1)^2=%lld M_0=%d\n",i+1,T[i],cnt[i]);
  for(int i=13;i<26;i++) printf("Z^aO^a a=%d N=(2^a+1)^2=%lld M_0=%d\n",i-12,T[i],cnt[i]);
  return 0; }
