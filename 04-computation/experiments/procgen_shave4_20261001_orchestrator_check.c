/* Orchestrator's independent embedding engine for the audit of the shave4 lane (2026-10-01).
 * Written from the definitions; the lane's code was not read.
 * usage: emb MODE n "a>b a>b ..."  < gentourng lines
 *   MODE = exist : for every tournament, does the oriented graph S embed (injective arc-preserving map)?
 *                  prints "exist n=.. classes=.. fail=.."
 *   MODE = parity: prints the histogram of emb(S,T) mod 2 over the input tournaments and the number with emb = 0
 *   MODE = count : prints emb(S,T) for each input tournament
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
static int n, m, A[64], B[64];
static uint32_t TOUT[32];
static int phi[32], used[32];
static int order[32];
static uint32_t prevmask[32]; /* for vertex order position i: S-arcs to earlier vertices */
static int adjS[32][32];      /* adjS[a][b] = 1 if a->b in S */
static u64 cnt; static int stop_at_one;
static int pos_of[32];

static void dfs(int i) {
  if (stop_at_one && cnt) return;
  if (i == n) { cnt++; return; }
  int x = order[i];
  for (int y = 0; y < n; y++) {
    if (used[y]) continue;
    int okk = 1;
    for (int j = 0; j < i && okk; j++) {
      int w = order[j];
      if (adjS[x][w] && !(TOUT[y] >> phi[w] & 1)) okk = 0;
      if (adjS[w][x] && !(TOUT[phi[w]] >> y & 1)) okk = 0;
    }
    if (!okk) continue;
    used[y] = 1; phi[x] = y;
    dfs(i + 1);
    used[y] = 0;
    if (stop_at_one && cnt) return;
  }
}
int main(int argc, char **argv) {
  const char *mode = argv[1]; n = atoi(argv[2]);
  char *spec = strdup(argv[3]); m = 0;
  memset(adjS, 0, sizeof adjS);
  for (char *tok = strtok(spec, " "); tok; tok = strtok(0, " ")) { int a, b; if (sscanf(tok, "%d>%d", &a, &b) == 2) { A[m] = a; B[m] = b; adjS[a][b] = 1; m++; } }
  /* vertex order: by decreasing degree in S (better pruning) */
  int deg[32] = {0}; for (int k = 0; k < m; k++) { deg[A[k]]++; deg[B[k]]++; }
  for (int i = 0; i < n; i++) order[i] = i;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (deg[order[j]] > deg[order[i]]) { int t = order[i]; order[i] = order[j]; order[j] = t; }
  char line[256]; long ncls = 0, fail = 0, odd = 0, even = 0, zero = 0;
  stop_at_one = !strcmp(mode, "exist");
  while (fgets(line, sizeof line, stdin)) {
    if (line[0] != '0' && line[0] != '1') continue;
    line[strcspn(line, "\r\n")] = 0;
    memset(TOUT, 0, sizeof TOUT);
    int k = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, k++) { if (line[k] == '1') TOUT[i] |= 1u << j; else TOUT[j] |= 1u << i; }
    cnt = 0; memset(used, 0, sizeof used);
    dfs(0);
    ncls++;
    if (!strcmp(mode, "exist")) { if (!cnt) fail++; }
    else if (!strcmp(mode, "count")) printf("%s %llu\n", line, cnt);
    else { if (cnt & 1) odd++; else even++; if (!cnt) zero++; }
  }
  if (!strcmp(mode, "exist")) printf("exist n=%d arcs=%d classes=%ld fail=%ld\n", n, m, ncls, fail);
  else if (!strcmp(mode, "parity")) printf("parity n=%d arcs=%d classes=%ld odd=%ld even=%ld zero=%ld\n", n, m, ncls, odd, even, zero);
  return 0;
}
