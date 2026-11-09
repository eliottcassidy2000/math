/* Independent audit code (auditor) for THM-5020 (1), (5), (6), THM-5021 counts, HYP-9590, PREREGISTRATION Q1-Q2.
 *
 * P(X) = max 3-separated subset of X, by the DP f(i+1) = max(f(i), [i in X](1 + f(i-2))) (not greedy).
 * k = ceil(L/3) = floor(n/3), L = n - 2.  0-junction = factor 11011 (centre 0), 1-junction = 00100.
 *
 * Modes:
 *   junction LMIN LMAX : every word of length L: the junction bound P(A) + P(B) >= ceil((L + 2J)/3);
 *                        counts of tight colourings and of those with a 0-/1-/both junctions (brute force).
 *   central LMAX       : every palindrome of even length L <= LMAX: if P(X) is odd then both central
 *                        positions lie in X; P(A), P(B) never both odd.  Also every palindrome of any length:
 *                        P(A) = P(B) = k never happens at n = 0, 1 (mod 3) (HYP-9590 (C)).
 *   pairsx NMIN NMAX   : exhaustive census of all words of length n - 2: which unordered pairs {P(A), P(B)} are
 *                        realised by (i) admissible non-palindromic words (black), (ii) admissible palindromes
 *                        (blue), (iii) any non-palindromic word, (iv) any palindrome; missing allowed pairs.
 *   pairsdp NMIN NMAX  : the same for (i), (ii) by a DP on (difference vectors of the packing DP, counts, last
 *                        four letters, junction flags) with palindromes enumerated directly.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef uint64_t u64;
#define MAXL 64

static int P_of_mask(u64 m, int L) { /* m: bit i set iff position i in X */
  int f0 = 0, f1 = 0, f2 = 0; /* f(i-2), f(i-1), f(i) */
  for (int i = 0; i < L; i++) {
    int nf = f2;
    if ((m >> i) & 1) { int c = 1 + f0; if (c > nf) nf = c; }
    f0 = f1; f1 = f2; f2 = nf;
  }
  return f2;
}
static int count_factor(u64 w, int L, unsigned pat5) { /* occurrences of a 5-letter factor, pat bit 4 = first letter */
  int s = 0;
  for (int z = 2; z + 2 < L; z++) {
    unsigned f = 0;
    for (int d = -2; d <= 2; d++) f = (f << 1) | ((w >> (z + d)) & 1);
    if (f == pat5) s++;
  }
  return s;
}
static int is_pal(u64 w, int L) {
  for (int i = 0; i < L / 2; i++) if (((w >> i) & 1) != ((w >> (L - 1 - i)) & 1)) return 0;
  return 1;
}
#define PAT11011 0x1B /* 11011 */
#define PAT00100 0x04 /* 00100 */

static void mode_junction(int lmin, int lmax) {
  for (int L = lmin; L <= lmax; L++) {
    int k = (L + 2) / 3;
    u64 full = (L == 64) ? ~0ULL : ((1ULL << L) - 1);
    long long bad = 0, tight = 0, tj0 = 0, tj1 = 0, tboth = 0, tpal = 0, maxJ = 0, sharp = 0;
    for (u64 w = 0; w <= full; w++) {
      int pa = P_of_mask(w, L), pb = P_of_mask(~w & full, L);
      int J0 = count_factor(w, L, PAT11011), J1 = count_factor(w, L, PAT00100), J = J0 + J1;
      if (J > maxJ) maxJ = J;
      if (3 * (pa + pb) < L + 2 * J) bad++;
      if (3 * (pa + pb) - (L + 2 * J) < 3 && J > 0) sharp++; /* bound attained (pa + pb = ceil((L+2J)/3)) */
      if (pa + pb == k) {
        tight++;
        if (J0) tj0++;
        if (J1) tj1++;
        if (J0 && J1) tboth++;
        tpal += is_pal(w, L);
      }
      if (w == full) break;
    }
    printf("L=%d (n=%d, k=%d): words %llu; junction-bound violations %lld (max J %lld; words with J>0 attaining "
           "the bound %lld); tight colourings %lld (palindromic %lld), with 11011 %lld, with 00100 %lld, with both "
           "%lld\n",
           L, L + 2, k, (unsigned long long)full + 1, bad, maxJ, sharp, tight, tpal, tj0, tj1, tboth);
    fflush(stdout);
  }
}

static void mode_central(int lmax) {
  for (int L = 1; L <= lmax; L++) {
    int h = L / 2, k = (L + 2) / 3, n = L + 2;
    long long bad = 0, both = 0, kk = 0, cnt = 0;
    for (u64 half = 0; half < (1ULL << ((L + 1) / 2)); half++) {
      /* positions 0..ceil(L/2)-1 from half, mirrored */
      u64 w = 0;
      for (int i = 0; i < (L + 1) / 2; i++)
        if ((half >> i) & 1) { w |= 1ULL << i; w |= 1ULL << (L - 1 - i); }
      cnt++;
      u64 full = (1ULL << L) - 1;
      int pa = P_of_mask(w, L), pb = P_of_mask(~w & full, L);
      if (L % 2 == 0) {
        for (int col = 0; col < 2; col++) {
          int p = col ? pa : pb;
          int in1 = (((w >> (h - 1)) & 1) == (u64)col), in2 = (((w >> h) & 1) == (u64)col);
          if (p % 2 == 1 && !(in1 && in2)) bad++;
        }
        if (pa % 2 == 1 && pb % 2 == 1) both++;
      }
      if (pa == k && pb == k) kk++;
    }
    printf("L=%d (n=%d, n mod 3 = %d): palindromes %lld; %s; palindromes with P(A) = P(B) = k: %lld\n", L, n, n % 3,
           cnt,
           (L % 2 == 0) ? (bad || both ? "CENTRAL-PAIR LEMMA VIOLATED" : "central-pair lemma holds, no both-odd pair")
                        : "odd length", kk);
    if (L % 2 == 0 && (bad || both)) printf("   violations %lld both-odd %lld\n", bad, both);
  }
}

/* ---- census output ---- */
static void report_pairs(int n, int k, unsigned char blk[MAXL][MAXL], unsigned char blu[MAXL][MAXL],
                         unsigned char anyb[MAXL][MAXL], unsigned char anyp[MAXL][MAXL], int have_any) {
  printf("n=%d k=%d: black missing [", n, k);
  int first = 1;
  for (int i = 1; i <= k; i++)
    for (int j = i; j <= k; j++)
      if (i + j >= k && !blk[i][j]) { printf("%s(%d, %d)", first ? "" : ", ", i, j); first = 0; }
  printf("]; blue (non-tight) missing [");
  first = 1;
  for (int i = 1; i <= k; i++)
    for (int j = i; j <= k; j++)
      if (i + j > k && !blu[i][j]) { printf("%s(%d, %d)", first ? "" : ", ", i, j); first = 0; }
  printf("]");
  if (have_any) {
    printf("; missing even among ALL words: black [");
    first = 1;
    for (int i = 1; i <= k; i++)
      for (int j = i; j <= k; j++)
        if (i + j >= k && !anyb[i][j]) { printf("%s(%d, %d)", first ? "" : ", ", i, j); first = 0; }
    printf("], blue [");
    first = 1;
    for (int i = 1; i <= k; i++)
      for (int j = i; j <= k; j++)
        if (i + j > k && !anyp[i][j]) { printf("%s(%d, %d)", first ? "" : ", ", i, j); first = 0; }
    printf("]");
  }
  printf("\n");
  fflush(stdout);
}

static void mode_pairsx(int nmin, int nmax) {
  for (int n = nmin; n <= nmax; n++) {
    int L = n - 2, k = n / 3;
    u64 full = (1ULL << L) - 1;
    static unsigned char blk[MAXL][MAXL], blu[MAXL][MAXL], anyb[MAXL][MAXL], anyp[MAXL][MAXL];
    memset(blk, 0, sizeof blk); memset(blu, 0, sizeof blu); memset(anyb, 0, sizeof anyb); memset(anyp, 0, sizeof anyp);
    for (u64 w = 0; w <= full; w++) {
      int pa = P_of_mask(w, L), pb = P_of_mask(~w & full, L);
      int lo = pa < pb ? pa : pb, hi = pa < pb ? pb : pa;
      int pal = is_pal(w, L);
      int adm = !count_factor(w, L, PAT11011) || !count_factor(w, L, PAT00100);
      if (pal) { anyp[lo][hi] = 1; if (adm) blu[lo][hi] = 1; }
      else { anyb[lo][hi] = 1; if (adm) blk[lo][hi] = 1; }
    }
    report_pairs(n, k, blk, blu, anyb, anyp, 1);
  }
}

/* DP: per colour, packing DP differences (d1, d2) with d1 = f(i) - f(i-1), d2 = f(i-1) - f(i-2);
 * reachable values (0,0), (1,0), (0,1); on letter of this colour, f increases iff d1 = d2 = 0. */
static int enc(int a1, int a2, int b1, int b2) { return ((a1 * 2 + a2) * 2 + b1) * 2 + b2; }
static void mode_pairsdp(int nmin, int nmax) {
  /* state: dstate (16) x cA (<= 14) x cB (<= 14) x last4 (16) x flags (4) */
#define CMAX 16
  static long long cur[16][CMAX][CMAX][16][4], nxt[16][CMAX][CMAX][16][4];
  for (int n = nmin; n <= nmax; n++) {
    int L = n - 2, k = n / 3;
    memset(cur, 0, sizeof cur);
    cur[enc(0, 0, 0, 0)][0][0][0][0] = 1;
    for (int pos = 0; pos < L; pos++) {
      memset(nxt, 0, sizeof nxt);
      for (int ds = 0; ds < 16; ds++)
        for (int ca = 0; ca < CMAX; ca++)
          for (int cb = 0; cb < CMAX; cb++)
            for (int l4 = 0; l4 < 16; l4++)
              for (int fl = 0; fl < 4; fl++) {
                long long c = cur[ds][ca][cb][l4][fl];
                if (!c) continue;
                int a1 = (ds >> 3) & 1, a2 = (ds >> 2) & 1, b1 = (ds >> 1) & 1, b2 = ds & 1;
                for (int col = 0; col < 2; col++) {
                  int na1, na2, nb1, nb2, nca = ca, ncb = cb;
                  if (col == 1) {
                    int inc = (a1 == 0 && a2 == 0);
                    na1 = inc; na2 = a1; nca += inc;
                    nb1 = 0; nb2 = b1;
                  } else {
                    int inc = (b1 == 0 && b2 == 0);
                    nb1 = inc; nb2 = b1; ncb += inc;
                    na1 = 0; na2 = a1;
                  }
                  /* last five letters: l4 holds the previous four (bit 3 = oldest), only valid after 4 letters */
                  int five = (l4 << 1) | col;
                  int nfl = fl;
                  if (pos >= 4) {
                    if ((five & 31) == PAT11011) nfl |= 1;
                    if ((five & 31) == PAT00100) nfl |= 2;
                  }
                  if (nfl == 3) continue; /* contains both: not admissible */
                  nxt[enc(na1, na2, nb1, nb2)][nca][ncb][five & 15][nfl] += c;
                }
              }
      memcpy(cur, nxt, sizeof cur);
    }
    static long long tot[MAXL][MAXL];
    memset(tot, 0, sizeof tot);
    long long all = 0;
    for (int ds = 0; ds < 16; ds++)
      for (int ca = 0; ca < CMAX; ca++)
        for (int cb = 0; cb < CMAX; cb++)
          for (int l4 = 0; l4 < 16; l4++)
            for (int fl = 0; fl < 4; fl++) {
              long long c = cur[ds][ca][cb][l4][fl];
              if (!c) continue;
              int lo = ca < cb ? ca : cb, hi = ca < cb ? cb : ca;
              tot[lo][hi] += c;
              all += c;
            }
    /* palindromes directly */
    static long long pal[MAXL][MAXL];
    memset(pal, 0, sizeof pal);
    long long npal = 0;
    u64 full = (1ULL << L) - 1;
    for (u64 half = 0; half < (1ULL << ((L + 1) / 2)); half++) {
      u64 w = 0;
      for (int i = 0; i < (L + 1) / 2; i++)
        if ((half >> i) & 1) { w |= 1ULL << i; w |= 1ULL << (L - 1 - i); }
      if (count_factor(w, L, PAT11011) && count_factor(w, L, PAT00100)) continue;
      int pa = P_of_mask(w, L), pb = P_of_mask(~w & full, L);
      int lo = pa < pb ? pa : pb, hi = pa < pb ? pb : pa;
      pal[lo][hi]++;
      npal++;
    }
    static unsigned char blk[MAXL][MAXL], blu[MAXL][MAXL];
    memset(blk, 0, sizeof blk); memset(blu, 0, sizeof blu);
    for (int i = 0; i < MAXL; i++)
      for (int j = 0; j < MAXL; j++) { blk[i][j] = tot[i][j] - pal[i][j] > 0; blu[i][j] = pal[i][j] > 0; }
    printf("[dp] admissible words %lld, admissible palindromes %lld; ", all, npal);
    report_pairs(n, k, blk, blu, NULL, NULL, 0);
  }
}

int main(int argc, char **argv) {
  if (argc < 3) return 1;
  if (!strcmp(argv[1], "junction")) mode_junction(atoi(argv[2]), atoi(argv[3]));
  else if (!strcmp(argv[1], "central")) mode_central(atoi(argv[2]));
  else if (!strcmp(argv[1], "pairsx")) mode_pairsx(atoi(argv[2]), atoi(argv[3]));
  else if (!strcmp(argv[1], "pairsdp")) mode_pairsdp(atoi(argv[2]), atoi(argv[3]));
  return 0;
}
