// complement_moves_bfs_n8_20261027.cpp -- opus-2026-10-27-S63
// Level-synchronous parallel BFS of the native complement move graph on all 2^28 labelled 8-tournaments, from the
// transitive tournament: a pre-registered test of HYP-9420 (2) (order 8 = 0 mod 4: connected, eccentricity >= 8).
// Move: T -> T^op xor E(P) for a Hamiltonian path P of T (every arc off P reversed).  Frontiers are bitsets, distances
// bytes; the only races are benign (two threads writing the same new distance).  Afterwards the score bound
// L(T) = (1/2) sum |s_(i) - i| is computed for every tournament and d >= L is checked (THM-4810, Theorem L).
// Build: g++ -O2 -std=c++17 -fopenmp -o cmb8 complement_moves_bfs_n8_20261027.cpp ; run: ./cmb8
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
using namespace std;
typedef unsigned int u32; typedef unsigned long long u64;
static const int n = 8, M = 28; static const u32 N = 1u << M, fullV = (1u << n) - 1, fullP = N - 1;
static int pidx[8][8]; static int pi_[M], pj_[M];
static inline void decode(u32 code, u32* out) { for (int v = 0; v < n; v++) out[v] = 0; for (int k = 0; k < M; k++) { if (code >> k & 1) out[pi_[k]] |= 1u << pj_[k]; else out[pj_[k]] |= 1u << pi_[k]; } }
static void hpaths(const u32* out, int v, u32 used, u32 pm, vector<u32>& res) { if (used == fullV) { res.push_back(pm); return; }
    u32 ext = out[v] & ~used; while (ext) { int w = __builtin_ctz(ext); ext &= ext - 1; hpaths(out, w, used | (1u << w), pm | (1u << pidx[v][w]), res); } }
int main() {
    int k = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { pidx[i][j] = pidx[j][i] = k; pi_[k] = i; pj_[k] = j; k++; }
    vector<uint8_t> dist(N, 255); vector<u64> fr(N / 64, 0), nx(N / 64, 0);
    u32 trans = fullP; dist[trans] = 0; fr[trans >> 6] |= 1ULL << (trans & 63);
    vector<u64> level; u64 reached = 1;
    for (int d = 0; ; d++) { u64 cnt = 0; for (u64 w : fr) cnt += __builtin_popcountll(w); if (!cnt) break; level.push_back(cnt);
        printf("level %d: %llu states (reached so far %llu)\n", d, cnt, reached); fflush(stdout); u64 added = 0;
        #pragma omp parallel
        {   u32 out[8]; vector<u32> res; u64 loc = 0;
            #pragma omp for schedule(dynamic, 256)
            for (long w = 0; w < (long)(N / 64); w++) { u64 bits = fr[w]; while (bits) { int b = __builtin_ctzll(bits); bits &= bits - 1; u32 code = (u32)w * 64 + b; decode(code, out); res.clear();
                    for (int v = 0; v < n; v++) hpaths(out, v, 1u << v, 0, res);
                    for (u32 pm : res) { u32 nb = code ^ (fullP ^ pm); if (dist[nb] == 255) { dist[nb] = (uint8_t)(d + 1); __atomic_fetch_or(&nx[nb >> 6], 1ULL << (nb & 63), __ATOMIC_RELAXED); loc++; } } } }
            #pragma omp atomic
            added += loc;
        }
        // exact count of the next level from the bitset (loc may double count benign races)
        u64 nc = 0; for (u64 w : nx) nc += __builtin_popcountll(w); reached += nc; fr.swap(nx); fill(nx.begin(), nx.end(), 0ULL); (void)added;
    }
    int ecc = (int)level.size() - 1;
    printf("n = 8: %u labelled tournaments; reached %llu (%s); eccentricity of the transitive tournament %d\nlevel counts:", N, reached, reached == (u64)N ? "connected" : "NOT connected", ecc);
    for (size_t d = 0; d < level.size(); d++) printf(" %zu:%llu", d, level[d]); printf("\n");
    // score bound check
    u64 sl[40] = {0}; long viol = 0; int maxL = 0;
    #pragma omp parallel
    {   u32 out[8]; u64 lsl[40] = {0}; long lv = 0; int lm = 0;
        #pragma omp for schedule(static)
        for (long c = 0; c < (long)N; c++) { if (dist[c] == 255) continue; decode((u32)c, out); int s[8]; for (int v = 0; v < n; v++) s[v] = __builtin_popcount(out[v]); sort(s, s + n); int t = 0; for (int i = 0; i < n; i++) t += abs(s[i] - i); int L = t / 2; if (L > lm) lm = L; int sk = dist[c] - L; if (sk < 0) lv++; else lsl[min(sk, 39)]++; }
        #pragma omp critical
        { for (int i = 0; i < 40; i++) sl[i] += lsl[i]; viol += lv; if (lm > maxL) maxL = lm; }
    }
    printf("score bound: max L = %d; violations of d >= L: %ld; (d - L) histogram:", maxL, viol); for (int i = 0; i < 40; i++) if (sl[i]) printf(" %d:%llu", i, sl[i]); printf("\n");
    u64 par[2] = {0, 0}; for (u64 c = 0; c < N; c++) if (dist[c] != 255) par[__builtin_popcountll(c) & 1]++;
    printf("forward-arc parity of reached states: even %llu, odd %llu (one move toggles 21 pairs: odd, no parity obstruction)\n", par[0], par[1]);
    return 0;
}
