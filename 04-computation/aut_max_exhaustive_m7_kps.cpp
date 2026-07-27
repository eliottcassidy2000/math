// aut_max_exhaustive_m7_kps.cpp -- kind-pasteur 2026-07-26
// EXHAUSTIVE check of o(m) = max |Aut(T)| over ALL labeled tournaments,
// m = 3..7 (2^C(m,2) masks; |Aut| counted over odd-order permutations only,
// valid because Aut(T) always has odd order).
// Also reports s(m) = max |Aut| over STRONG tournaments.
#include <cstdio>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <numeric>
using namespace std;

int M;
int pidx[8][8]; // pair index for i<j

static bool odd_order(const vector<int>& p) {
    vector<char> vis(M, 0);
    for (int i = 0; i < M; i++) {
        if (vis[i]) continue;
        int len = 0, j = i;
        while (!vis[j]) { vis[j] = 1; j = p[j]; len++; }
        if (len % 2 == 0) return false;
    }
    return true;
}

static bool strong(uint32_t mask) {
    // adjacency: bit pidx[i][j] set => i->j (for i<j); else j->i
    uint32_t adj[8], radj[8];
    for (int i = 0; i < M; i++) { adj[i] = 0; radj[i] = 0; }
    for (int i = 0; i < M; i++)
        for (int j = i + 1; j < M; j++) {
            if ((mask >> pidx[i][j]) & 1) { adj[i] |= 1u << j; radj[j] |= 1u << i; }
            else                          { adj[j] |= 1u << i; radj[i] |= 1u << j; }
        }
    for (int pass = 0; pass < 2; pass++) {
        uint32_t* g = pass ? radj : adj;
        uint32_t seen = 1, frontier = 1;
        while (frontier) {
            uint32_t nf = 0;
            for (int v = 0; v < M; v++)
                if ((frontier >> v) & 1) nf |= g[v];
            nf &= ~seen;
            seen |= nf;
            frontier = nf;
        }
        if (seen != (1u << M) - 1) return false;
    }
    return true;
}

int main() {
    for (M = 3; M <= 7; M++) {
        int NP = M * (M - 1) / 2;
        int k = 0;
        for (int i = 0; i < M; i++)
            for (int j = i + 1; j < M; j++) pidx[i][j] = k++;
        // odd-order permutations (excluding identity), as pair maps
        vector<int> base(M);
        iota(base.begin(), base.end(), 0);
        vector<vector<int>> permTargets;  // per perm: for each pair slot, target slot
        vector<uint32_t> permFlips;       // per perm: flip bit per pair slot
        vector<int> p = base;
        do {
            if (p == base) continue;
            if (!odd_order(p)) continue;
            vector<int> tgt(NP);
            uint32_t flip = 0;
            for (int i = 0; i < M; i++)
                for (int j = i + 1; j < M; j++) {
                    int a = p[i], b = p[j];
                    int slot = pidx[i][j];
                    if (a < b) tgt[slot] = pidx[a][b];
                    else { tgt[slot] = pidx[b][a]; flip |= 1u << slot; }
                }
            permTargets.push_back(tgt);
            permFlips.push_back(flip);
        } while (next_permutation(p.begin(), p.end()));

        uint32_t NMASK = 1u << NP;
        vector<uint16_t> cnt(NMASK, 0);
        // split-table application: pi(T) = A[T_low] ^ B[T_high] (bit-gather is
        // XOR-linear; orientation flips folded into the tables via image of 0)
        int LOW = NP / 2, HIGH = NP - LOW;
        vector<uint32_t> TA(1u << LOW), TB(1u << HIGH);
        for (size_t pi = 0; pi < permTargets.size(); pi++) {
            const vector<int>& tgt = permTargets[pi];
            uint32_t flip = 0;
            // image of mask 0: all pairs oriented j->i... careful: bit=0 means
            // j beats i; permutation maps orientation, flips recorded:
            // image bit at tgt[slot] = bit(slot) XOR flipbit(slot)
            for (int slot = 0; slot < NP; slot++)
                if ((permFlips[pi] >> slot) & 1) flip |= 1u << tgt[slot];
            for (uint32_t lo = 0; lo < (1u << LOW); lo++) {
                uint32_t im = 0;
                for (int slot = 0; slot < LOW; slot++)
                    if ((lo >> slot) & 1) im |= 1u << tgt[slot];
                TA[lo] = im;
            }
            for (uint32_t hi = 0; hi < (1u << HIGH); hi++) {
                uint32_t im = 0;
                for (int slot = 0; slot < HIGH; slot++)
                    if ((hi >> slot) & 1) im |= 1u << tgt[LOW + slot];
                TB[hi] = im;
            }
            uint32_t lowmask = (1u << LOW) - 1;
            for (uint32_t T = 0; T < NMASK; T++) {
                uint32_t im = TA[T & lowmask] ^ TB[T >> LOW] ^ flip;
                if (im == T) cnt[T]++;
            }
        }
        int best = 0, bestStrong = 0;
        uint32_t argbest = 0, argbestStrong = 0;
        for (uint32_t T = 0; T < NMASK; T++) {
            int a = cnt[T] + 1;
            if (a > best) { best = a; argbest = T; }
            if (a > bestStrong && strong(T)) { bestStrong = a; argbestStrong = T; }
        }
        printf("m=%d: masks=%u, odd-order perms=%zu, o(m)=%d (mask %u), "
               "strong-max s(m)=%d (mask %u)\n",
               M, NMASK, permTargets.size() + 1, best, argbest,
               bestStrong, argbestStrong);
        fflush(stdout);
    }
    return 0;
}
