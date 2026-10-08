#!/usr/bin/env python3
"""audit H: THM-4611 statement 3' (adaptive certificates for the five Z_5 maps with coupling group {+-1}).

Independent check, with the author's integer forms but NOT the author's partition: the greedy exact refinement.
Every hidden state mod 5^3 is a node; a node whose block A_3 is zero or exactly balanced by Q is a leaf; otherwise it is
replaced by all its lifts mod 5^4 (M-lifts in H_4 over M, e-lifts e + 125 t), and so on; at level K = 5 every node must be
zero or exactly balanced.  For a fixed Q this succeeds iff SOME refinement-tree partition with levels 3..5 works (an
unbalanced node must be refined in every valid partition), so success confirms the existence claim of 3'.
Also checks: the lifts partition the states (fibre sizes), coupling group, sticky reflections, property (i), Lambda,
the minimal block rank at fixed k <= 4, and |H_s|.
Usage: python3 e_adaptive_check.py [NAME ...]"""
import math, sys, time
import numpy as np
from hdp import Tables, balance_exact, margins_float, exact_rank_vec
from hcore import ratio_group

MAPS = {
    'Z5_1_4_1_11_34': ([1, 4, 1, 11, 34], [[19, -4, -5], [-4, 23, -9], [-5, -9, 24]]),
    'Z5_1_6_11_11_4': ([1, 6, 11, 11, 4], [[33, -15, -13], [-15, 40, -11], [-13, -11, 32]]),
    'Z5_1_1_6_39_11': ([1, 1, 6, 39, 11], [[19, -4, -9], [-4, 15, -3], [-9, -3, 20]]),
    'Z5_1_14_4_1_29': ([1, 14, 4, 1, 29], [[46, -22, -9], [-22, 50, -8], [-9, -8, 35]]),
    'Z5_1_4_31_1_6': ([1, 4, 31, 1, 6], [[30, -12, -5], [-12, 32, -7], [-5, -7, 25]]),
}

def bal_or_zero(Q, A):
    flat = A.reshape(len(A), -1)
    nz = np.any(flat != 0, axis=1)
    ok = np.ones(len(A), dtype=bool)
    if nz.any():
        An = A[nz]
        U, inv = np.unique(An.reshape(len(An), -1), axis=0, return_inverse=True)
        okU = balance_exact(Q, U.reshape(-1, A.shape[1], A.shape[2]))
        ok[np.nonzero(nz)[0]] = okU[inv.ravel()]
    return ok, nz

def run(name, K0=3, K=5):
    m, Q = MAPS[name]
    d = 5
    t0 = time.time()
    tb = Tables(d, m)
    Lam = sum(math.log(x / d) for x in m) / d
    G = sorted(ratio_group(m, d))
    mb = [x % d for x in m]
    sticky = []
    for b in range(d):
        if (d - 1) in G and all(mb[(b - j) % d] == mb[j] for j in range(d)):
            Dm = tb.Dab[d - 1, b]
            sticky.append((b, int(np.linalg.matrix_rank(Dm.astype(float)))))
    prop_i = not any(all(m[(a * j + b) % d] == m[j] for j in range(d)) for a in G for b in range(d) if (a, b) != (1, 0))
    Hsz = [len(tb.Hs(s)) for s in range(1, K + 1)]
    print(f"{name}: m = {m}, r = {tb.r}, residues {mb}, Lambda = {Lam:+.4f}, rank {tb.rho}, coords {tb.v} (den {tb.den}), "
          f"G = {G}, |H_s| s=1..{K}: {Hsz}, property (i) {prop_i}, sticky reflections (b, rank) {sticky}", flush=True)
    # fixed-length: minimal exact block rank for k <= 4
    for k in range(1, 5):
        T = tb.table(k).reshape(-1, 3, 3)
        nz = np.any(T.reshape(len(T), -1) != 0, axis=1)
        U = np.unique(T[nz].reshape(-1, 9), axis=0).reshape(-1, 3, 3)
        rk = exact_rank_vec(U)
        okQ = balance_exact(Q, U)
        print(f"   fixed k = {k}: {len(U)} distinct nonzero blocks, min exact rank {rk.min()}, zero states {int((~nz).sum())}, "
              f"blocks balanced by the 3' form: {int(okQ.sum())}/{len(U)}", flush=True)
    # greedy exact refinement
    H = {s: tb.Hs(s) for s in range(1, K + 1)}
    pos = {s: {int(x): i for i, x in enumerate(H[s])} for s in range(1, K)}
    # level K0 nodes: all states
    Mv = np.repeat(H[K0], d ** K0); Ev = np.tile(np.arange(d ** K0, dtype=np.int64), len(H[K0]))
    leaves = []
    refined = []
    total_cover = 0     # number of level-K states covered by leaves (must equal |H_K| d^K)
    fail_K = 0
    for s in range(K0, K + 1):
        if s <= K - 1:
            T = tb.table(s)
            ix = np.array([pos[s][int(x)] for x in Mv], dtype=np.int64)
            A = T[ix, Ev]
        else:
            A = tb.blocks_at(s, Mv, Ev)
        ok, nz = bal_or_zero(Q, A)
        leaf = ok if s < K else np.ones(len(A), dtype=bool)
        if s == K and not ok.all():
            fail_K = int((~ok).sum())
            print(f"   FAIL: {fail_K} level-{K} nodes not balanced", flush=True)
        leaves.append(A[leaf & nz])
        fib = (len(H[K]) // len(H[s])) * d ** (K - s)
        total_cover += int(leaf.sum()) * fib
        bad = ~leaf
        refined.append(int(bad.sum()))
        if s == K or not bad.any():
            break
        # lifts of the bad nodes to level s+1
        Hn = H[s + 1]; red = Hn % d ** s
        Mb, Eb = Mv[bad], Ev[bad]
        lifts_M = {}
        for x in np.unique(Mb):
            lifts_M[int(x)] = Hn[red == x]
        fibre_sizes = {len(v) for v in lifts_M.values()}
        assert fibre_sizes == {len(Hn) // len(H[s])}, fibre_sizes
        nM, nE = [], []
        for x, e in zip(Mb, Eb):
            for Ml in lifts_M[int(x)]:
                for t in range(d):
                    nM.append(int(Ml)); nE.append(int(e) + d ** s * t)
        Mv = np.array(nM, dtype=np.int64); Ev = np.array(nE, dtype=np.int64)
    U = np.unique(np.concatenate([x.reshape(-1, 9) for x in leaves]), axis=0).reshape(-1, 3, 3)
    mg, al = margins_float(Q, U)
    full = len(H[K]) * d ** K
    print(f"   greedy exact refinement with the 3' form: refined per level {refined[:K - K0]} (level-{K} nodes not balanced: {fail_K}); "
          f"leaves cover {total_cover} of {full} level-{K} states: {total_cover == full}; {len(U)} distinct nonzero leaf blocks, "
          f"all exactly balanced; float margin of the form over the leaves {mg.min():+.5f}  [{time.time() - t0:.1f}s]", flush=True)

if __name__ == '__main__':
    names = sys.argv[1:] or list(MAPS)
    for nm in names:
        run(nm)
