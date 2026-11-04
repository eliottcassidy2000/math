#!/usr/bin/env python3
"""audit_B_fact_i_circulants.py -- audit B of THM-4871 statement 7, fact (i) "prime circulants".

Claims tested, for r prime:
 (a) every circulant graph Cay(Z_r, C) with C nonempty, non-complete is PRIME (no module S with 1 < |S| < r),
     by brute force over ALL subsets S (every S for r <= 13; for r = 17, 19 all S containing 0 -- translates of
     modules are modules -- and the symmetric sets up to dilation, vectorised with numpy);
 (b) |Aut Cay(Z_r, C)| = r * |M(C)| (M(C) = dilation stabiliser), i.e. Aut = Z_r x| M(C) <= AGL(1, r)  [pynauty];
 (c) no two vertices of such a circulant are twins (false or true);
 (d) K_r is never multiplier-even: on K_r the reversal sign of every dilation equals its permutation sign, and the
     dilation by a primitive root (cycle type (1, r-1)) is odd; the empty graph is multiplier-even.
"""
import sys, time, itertools
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from auditb_common import *

def modules_bruteforce(r, C):
    """all nontrivial modules by brute force over all subsets (r <= 13)."""
    adj = [0] * r
    for i in range(r):
        for c in C:
            adj[i] |= 1 << ((i + c) % r)
    full = (1 << r) - 1
    mods = []
    for S in range(1, full):
        pc = bin(S).count('1')
        if pc <= 1 or pc >= r:
            continue
        ok = True
        for v in range(r):
            if (S >> v) & 1:
                continue
            x = adj[v] & S
            if x != 0 and x != S:
                ok = False; break
        if ok:
            mods.append(S)
    return mods

def modules_numpy(r, C):
    """nontrivial modules containing vertex 0, vectorised (r = 17, 19)."""
    adj = [0] * r
    for i in range(r):
        for c in C:
            adj[i] |= 1 << ((i + c) % r)
    S = np.arange(1 << (r - 1), dtype=np.int64) * 2 + 1   # all subsets containing 0
    ok = np.ones(S.shape, dtype=bool)
    for v in range(1, r):
        inS = ((S >> v) & 1).astype(bool)
        x = S & adj[v]
        ok &= inS | (x == 0) | (x == S)
    pc = np.array([bin(int(s)).count('1') for s in S]) if r <= 13 else None
    # popcount via bit tricks
    cnt = np.zeros(S.shape, dtype=np.int64)
    for v in range(r):
        cnt += (S >> v) & 1
    ok &= (cnt > 1) & (cnt < r)
    return S[ok]

def main():
    t0 = time.time()
    print("audit B / fact (i): prime circulants")
    for r in (3, 5, 7, 11, 13, 17, 19):
        sets = symmetric_sets(r)
        allsets = (r <= 13)
        reps = sets if allsets else dilation_classes(sets, r)
        n_nontrivial = 0; bad_modules = 0; bad_aut = 0; bad_twin = 0
        for C in reps:
            if len(C) == 0 or len(C) == r - 1:
                continue
            n_nontrivial += 1
            # (a) modules
            if allsets:
                mods = modules_bruteforce(r, C)
            else:
                mods = modules_numpy(r, C)
            if len(mods):
                bad_modules += 1
                print("  r=%d C=%s has nontrivial modules: %s" % (r, sorted(C), [bin(int(m)) for m in mods[:5]]))
            # (b) automorphism group
            edges = circulant_edges(C, r)
            cert, gens, order = graph_cert_and_gens(r, edges)
            M = multiplier_stabiliser(C, r)
            if int(round(order)) != r * len(M):
                bad_aut += 1
                print("  r=%d C=%s |Aut|=%s != r*|M(C)|=%d" % (r, sorted(C), order, r * len(M)))
            # (c) twins
            nb = [set() for _ in range(r)]
            for (u, v) in edges:
                nb[u].add(v); nb[v].add(u)
            for u in range(r):
                for v in range(u + 1, r):
                    if nb[u] - {v} == nb[v] - {u}:
                        bad_twin += 1
        print("r=%2d: %d nonempty non-complete symmetric sets%s; nontrivial modules found in %d; |Aut| != r|M(C)| in %d; twin pairs %d  (%.1fs)"
              % (r, n_nontrivial, "" if allsets else " (dilation classes)", bad_modules, bad_aut, bad_twin, time.time() - t0))
        # (d) K_r and the empty graph
        Kr = frozenset(range(1, r))
        edges = circulant_edges(Kr, r)
        signs_match = True; odd_dilations = []
        for a in units(r):
            g = [(a * i) % r for i in range(r)]
            rp = reversal_parity(g, edges)
            if (-1) ** rp != perm_sign(g):
                signs_match = False
            if rp:
                odd_dilations.append(a)
        print("      K_%d: reversal sign == permutation sign for all dilations: %s; odd dilations: %s; multiplier-even: %s; empty graph multiplier-even: %s"
              % (r, signs_match, odd_dilations, multiplier_even(Kr, r), multiplier_even(frozenset(), r)))
    print("done (%.1fs)" % (time.time() - t0))

if __name__ == "__main__":
    main()
