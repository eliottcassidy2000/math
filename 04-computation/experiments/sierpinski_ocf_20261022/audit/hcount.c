/* H (Hamiltonian paths) of each tournament read from stdin in gentourng upper-triangle format; prints H per line.
   Also mode "interval N": the 3-adic interval tournament on {0..N-1} (x->y iff unit part of y-x is 1 mod 3). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n, adj[24][24];
static u64 H(void){
  int N = 1 << n; u64 *f = calloc((size_t)N * n, sizeof(u64));
  for(int v = 0; v < n; v++) f[(size_t)(1 << v) * n + v] = 1;
  for(int S = 1; S < N; S++) for(int v = 0; v < n; v++){ u64 c = f[(size_t)S * n + v]; if(!c) continue;
    for(int w = 0; w < n; w++) if(!(S >> w & 1) && adj[v][w]) f[(size_t)(S | 1 << w) * n + w] += c; }
  u64 h = 0; for(int v = 0; v < n; v++) h += f[(size_t)(N - 1) * n + v]; free(f); return h; }
static int sgn3(int z){ while(z % 3 == 0) z /= 3; int r = z % 3; if(r < 0) r += 3; return r; }
int main(int argc, char **argv){
  if(argc > 2 && !strcmp(argv[1], "interval")){
    int M = atoi(argv[2]);
    for(n = 1; n <= M; n++){ for(int i = 0; i < n; i++) for(int j = 0; j < n; j++) adj[i][j] = (i != j && sgn3(j - i) == 1);
      printf("%llu%s", H(), n < M ? ", " : "\n"); }
    return 0; }
  n = atoi(argv[1]); char line[512];
  while(fgets(line, sizeof line, stdin)){ if(line[0] == '>' || line[0] == '\n') continue; int b = 0;
    for(int i = 0; i < n; i++){ adj[i][i] = 0; for(int j = i + 1; j < n; j++){ int bit = line[b++] - '0'; adj[i][j] = bit; adj[j][i] = 1 - bit; } }
    line[strcspn(line, "\n")] = 0; printf("%s %llu\n", line, H()); }
  return 0; }
