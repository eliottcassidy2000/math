#!/usr/bin/env python3
"""audit H: greedy exact refinement test for adaptive block certificates (THM-4611 2', 3', and the Z_4 claim of 4).

Given a map, an integer form Q and levels K0 <= K: every hidden state mod d^K0 is a node; a node whose block A_s is zero
or exactly balanced by Q is a leaf; otherwise it is replaced by ALL its lifts mod d^(s+1) (M-lifts = one lift times the
kernel of H_(s+1) -> H_s; e-lifts e + d^s t); at level K every node must be zero or exactly balanced.  For fixed Q this
succeeds iff some refinement-tree partition with leaves at levels K0..K is certified by Q (an unbalanced node must be
refined in every valid partition), so it independently confirms the existence claim without the author's partition.
Checks along the way: kernel sizes (fibres of H_(s+1) -> H_s all equal), leaves cover every state mod d^K exactly once
(counted with multiplicity), last level processed in chunks.
Usage: python3 h_greedy.py NAME [NAME ...]"""
import math, sys, time
import numpy as np
from hdp import Tables, balance_exact, margins_float

MAPS = {
    # name: (d, m, Q, K0, K)
    'Z4_1357': (4, [1, 3, 5, 7], [[48, -13, -21], [-13, 50, -14], [-21, -14, 48]], 1, 6),
    'Z5_1_6_4_31_1': (5, [1, 6, 4, 31, 1], [[256, -31, -111], [-31, 191, -60], [-111, -60, 243]], 3, 5),
    'Z5_1_4_1_11_34': (5, [1, 4, 1, 11, 34], [[19, -4, -5], [-4, 23, -9], [-5, -9, 24]], 3, 5),
    'Z5_1_6_11_11_4': (5, [1, 6, 11, 11, 4], [[33, -15, -13], [-15, 40, -11], [-13, -11, 32]], 3, 5),
    'Z5_1_1_6_39_11': (5, [1, 1, 6, 39, 11], [[19, -4, -9], [-4, 15, -3], [-9, -3, 20]], 3, 5),
    'Z5_1_14_4_1_29': (5, [1, 14, 4, 1, 29], [[46, -22, -9], [-22, 50, -8], [-9, -8, 35]], 3, 5),
    'Z5_1_4_31_1_6': (5, [1, 4, 31, 1, 6], [[30, -12, -5], [-12, 32, -7], [-5, -7, 25]], 3, 5),
    'Z5_1_1_28_11_7': (5, [1, 1, 28, 11, 7], [[19, -3, -8], [-3, 14, -3], [-8, -3, 20]], 3, 5),
    'Z5_1_3_7_24_1': (5, [1, 3, 7, 24, 1], [[16, -3, -6], [-3, 11, -2], [-6, -2, 14]], 3, 5),
    # THM-4611 statement 6 (audit F's maps, constants as in audit F)
    'Z5_1_8_3_7_12': (5, [1, 8, 3, 7, 12], [[43, 0, -11], [0, 50, -16], [-11, -16, 27]], 2, 5, [0, -3, 4, 4, 2]),
    'Z5_1_2_3_7_6': (5, [1, 2, 3, 7, 6], [[40, 5, -9], [5, 34, -11], [-9, -11, 20]], 2, 5, [0, 3, 4, 4, 1]),
    'Z7_1_2_3_5_1_1_1': (7, [1, 2, 3, 5, 1, 1, 1], [[49, -5, -21], [-5, 32, -6], [-21, -6, 50]], 3, 5, [0, 5, 1, 6, 3, 2, 1]),
}

def bal_or_zero(Q, A):
    flat = A.reshape(len(A), -1)
    nz = np.any(flat != 0, axis=1)
    ok = np.ones(len(A), dtype=bool)
    U = None
    if nz.any():
        An = A[nz]
        U, inv = np.unique(An.reshape(len(An), -1), axis=0, return_inverse=True)
        okU = balance_exact(Q, U.reshape(-1, A.shape[1], A.shape[2]))
        ok[np.nonzero(nz)[0]] = okU[inv.ravel()]
    return ok, nz, U

def run(name, chunk=200000):
    d, m, Q, K0, K = MAPS[name][:5]
    r = MAPS[name][5] if len(MAPS[name]) > 5 else None
    t0 = time.time()
    tb = Tables(d, m, r)
    rho = tb.rho
    H = {s: tb.Hs(s) for s in range(1, K + 1)}
    print(f"{name}: d = {d}, m = {m}, r = {tb.r}, rank {rho}, Lambda = {sum(math.log(x / d) for x in m) / d:+.4f}, "
          f"|H_s| = {[len(H[s]) for s in range(1, K + 1)]}, form {Q}, levels {K0}..{K}", flush=True)
    # property (i) directly
    G = sorted(set(int(x) for x in H[1]))
    prop_i = not any(all(m[(a * j + b) % d] == m[j] for j in range(d)) for a in G for b in range(d) if (a, b) != (1, 0))
    # kernels and canonical lifts
    kern, lift = {}, {}
    for s in range(1, K):
        Hn, mod_s = H[s + 1], d ** s
        kern[s] = Hn[Hn % mod_s == 1]
        lift_arr = -np.ones(mod_s, dtype=np.int64)
        red = Hn % mod_s
        lift_arr[red[::-1]] = Hn[::-1]          # some lift of each element of H_s
        assert np.all(lift_arr[H[s]] >= 0)
        lift[s] = lift_arr
        assert len(kern[s]) * len(H[s]) == len(Hn)
    Mv = np.repeat(H[K0], d ** K0)
    Ev = np.tile(np.arange(d ** K0, dtype=np.int64), len(H[K0]))
    cover = 0; refined = []; nleafblocks = 0; keys = []; failK = 0; minmg = np.inf; zero_leaves = 0
    for s in range(K0, K + 1):
        nbad_M, nbad_E = [], []
        fib = (len(H[K]) // len(H[s])) * d ** (K - s)
        for c0 in range(0, len(Mv), chunk):
            Mc, Ec = Mv[c0:c0 + chunk], Ev[c0:c0 + chunk]
            if s < K:
                T = tb.table(s)
                lookup = -np.ones(d ** s, dtype=np.int64); lookup[H[s]] = np.arange(len(H[s]))
                A = T[lookup[Mc], Ec]
            else:
                A = tb.blocks_at(s, Mc, Ec)
            ok, nz, U = bal_or_zero(Q, A)
            if s == K:
                failK += int((~ok).sum())
                leaf = np.ones(len(A), dtype=bool)
            else:
                leaf = ok
            zero_leaves += int((leaf & ~nz).sum())
            cover += int(leaf.sum()) * fib
            if (leaf & nz).any():
                L = np.unique(A[leaf & nz].reshape(-1, rho * rho), axis=0)
                keys.append(L)
                mg, _ = margins_float(Q, L.reshape(-1, rho, rho))
                minmg = min(minmg, float(mg.min()))
            nbad_M.append(Mc[~leaf]); nbad_E.append(Ec[~leaf])
        Mb = np.concatenate(nbad_M); Eb = np.concatenate(nbad_E)
        refined.append(len(Mb))
        if s == K or len(Mb) == 0:
            break
        # all lifts to level s+1
        base = lift[s][Mb]
        kk = kern[s]
        Ml = (base[:, None] * kk[None, :]) % d ** (s + 1)              # (nb, |ker|)
        Ml = np.repeat(Ml[:, :, None], d, axis=2)                       # (nb, |ker|, d)
        El = Eb[:, None, None] + d ** s * np.arange(d, dtype=np.int64)[None, None, :]
        El = np.broadcast_to(El, Ml.shape)
        Mv = Ml.reshape(-1); Ev = np.ascontiguousarray(El).reshape(-1)
        # sanity: every lift reduces to its parent
        assert np.all(Mv % d ** s == np.repeat(Mb, len(kk) * d)) and np.all(Ev % d ** s == np.repeat(Eb, len(kk) * d))
    allU = np.unique(np.concatenate(keys), axis=0) if keys else np.zeros((0, rho * rho))
    full = len(H[K]) * d ** K
    print(f"   property (i) {prop_i}; greedy exact refinement: nodes refined per level {refined[:-1] if refined[-1] == 0 else refined} "
          f"(levels {K0}..{K - 1}); level-{K} nodes not balanced: {failK}; zero-block leaves {zero_leaves}; leaves cover "
          f"{cover} of {full} states mod d^{K}: {cover == full}; {len(allU)} distinct nonzero leaf blocks, ALL exactly balanced: "
          f"{failK == 0}; float margin of the form over these leaves {minmg:+.5f}  [{time.time() - t0:.1f}s]", flush=True)

if __name__ == '__main__':
    for nm in (sys.argv[1:] or list(MAPS)):
        run(nm)
