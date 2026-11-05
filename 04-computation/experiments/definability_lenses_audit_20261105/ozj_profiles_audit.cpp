// ozj_profiles_audit.cpp -- independent audit code.  Replays Session 200's physical source grammar from the seed
// t = (1,0,2,1) in raw order (D,E,U,d):  O(s) = (D,0,2U,d), Z(s) = (D-E+U, d, 2U-d, d), J(s) = (D-E+U, D-E, D+U, D).
// U strictly increases, so the set of distinct profiles with U <= B is finite; we count it, and the fibre sizes
// pf(n) = #{distinct profiles with D = n}.  Usage: ozj B
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <unordered_set>
using namespace std;
typedef unsigned long long u64;
static inline u64 pack(u64 D, u64 E, u64 U, u64 d) { return (D << 48) | (E << 32) | (U << 16) | d; }
int main(int argc, char** argv) {
    u64 B = strtoull(argv[1], 0, 10);
    unordered_set<u64> seen; vector<u64> q; u64 t = pack(1, 0, 2, 1); seen.insert(t); q.push_back(t);
    for (size_t i = 0; i < q.size(); i++) { u64 s = q[i]; u64 D = s >> 48, E = (s >> 32) & 0xFFFF, U = (s >> 16) & 0xFFFF, d = s & 0xFFFF;
        u64 nx[3][4] = { {D, 0, 2 * U, d}, {D - E + U, d, 2 * U - d, d}, {D - E + U, D - E, D + U, D} };
        for (int k = 0; k < 3; k++) { if (nx[k][2] > B) continue; if (nx[k][0] >= 65536 || nx[k][1] >= 65536 || nx[k][3] >= 65536) { fprintf(stderr, "overflow of the 16-bit packing\n"); return 1; }
            u64 p = pack(nx[k][0], nx[k][1], nx[k][2], nx[k][3]); if (seen.insert(p).second) q.push_back(p); } }
    // readout F = D - E + U (Session 200's C1 readout: append J); F >= U + 1 so U <= B covers every F <= B + 1
    vector<u64> pf(B + 2, 0); u64 maxD = 0; for (u64 s : seen) { u64 D = s >> 48, E = (s >> 32) & 0xFFFF, U = (s >> 16) & 0xFFFF; u64 F = D - E + U; if (F <= B + 1) pf[F]++; if (D > maxD) maxD = D; }
    printf("B = %llu: %zu distinct profiles with U <= B; max D = %llu; readout F = D - E + U\n", B, seen.size(), maxD);
    printf("odd n <= B+1 with pf(n) = 0:"); int z = 0; for (u64 n = 1; n <= B + 1; n += 2) if (pf[n] == 0) { printf(" %llu", n); z++; } printf("  (%d values)\n", z);
    printf("odd n <= B+1 with pf(n) = 1:"); for (u64 n = 1; n <= B + 1; n += 2) if (pf[n] == 1) printf(" %llu", n); printf("\n");
    printf("pf(n) for odd n <= 101:"); for (u64 n = 1; n <= 101 && n <= B + 1; n += 2) printf(" %llu:%llu", n, pf[n]); printf("\n");
    return 0; }
