/* i_classmass.c - auditor: representation mass of F_j (exactly j root changes, canonical words) by H mod 3 in the
 * window [X/2, X], and the share of that mass carried by final blocks confined by THM-4888 (C) with p = 3
 * (3 | gcd(d, R) of the block's seed).  Same move-only DFS as g_rootchanges.c.  usage: i_classmass X K */
#include <stdio.h>
#include <stdlib.h>
typedef long long ll;
static ll X; static int K; static ll mass[3][3], conf[3][3];
/* conf3: whether the CURRENT block (since the last J, or t) has 3 | gcd(d, U) at its start: tracked via the block seed */
static void dfs(ll D, ll E, ll U, ll d, int j, int lastO, int c3){
  if (!lastO && D <= X && 2 * D > X){ mass[j][D % 3]++; if (c3) conf[j][D % 3]++; }
  ll D2 = D - E + U;
  if (D2 <= X){
    dfs(D2, d, 2 * U - d, d, j, 0, c3);
    if (j < K){ ll R2 = D + U, d2 = D; dfs(D2, D - E, R2, d2, j + 1, 0, (R2 % 3 == 0 && d2 % 3 == 0)); }
  }
  if (D + 2 * U <= X) dfs(D, 0, 2 * U, d, j, 1, c3);
}
int main(int argc, char **argv){
  X = atoll(argv[1]); K = atoi(argv[2]);
  dfs(1, 0, 2, 1, 0, 0, 0);
  for (int j = 0; j <= K; j++){
    ll t = mass[j][0] + mass[j][1] + mass[j][2];
    printf("exactly %d root changes, words with H in [X/2, X], X=%lld: share by H mod 3 = %.4f %.4f %.4f (total %lld); share from 3-confined final blocks = %.4f, their mass by class = %lld %lld %lld\n",
      j, X, (double)mass[j][0] / t, (double)mass[j][1] / t, (double)mass[j][2] / t, t,
      (double)(conf[j][0] + conf[j][1] + conf[j][2]) / t, conf[j][0], conf[j][1], conf[j][2]);
  }
  return 0;
}
