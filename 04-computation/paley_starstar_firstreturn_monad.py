#!/usr/bin/env python3
"""
paley_starstar_firstreturn_monad.py
monad-explorer-2026-06-07 (deep-research, 6th session)

(★★) reduced to the FREE-CUMULANT RECURSION, and the first-return decomposition
analyzed.  Since S_k = Σ_{even-series σ} μ(0̂,σ) = (−1)^k C_k is the free cumulant
κ_{2k+2}/A^{k+1} of the two-point law, it must satisfy the R-transform / first-return
recursion.  Algebra:  S_k=(−1)^k C_k, C_k=Σ_{i+j=k−1}C_iC_j ⟹

      S_k = − Σ_{i+j=k−1, i,j≥0} S_i S_j ,   S_0 = 1.          (REC)

equivalently the GF F(x)=Σ S_k x^k satisfies  x F² + F − 1 = 0
(F = (−1+√(1+4x))/(2x) = Σ (−1)^k C_k x^k).  This script:

 (1) verifies (REC) and the quadratic GF numerically;
 (2) tests the FIRST-RETURN decomposition of even-series patterns: split σ at the
     first r>0 with block(r)=block(0).  A split is CLEAN if the first excursion
     positions {0..r} and the remainder {r..2k} share ONLY the vertex block(0).
     We compute  Σ_clean μ  vs  Σ_crossing μ  graded by (i,j)=(excursion size,
     remainder size), and check whether the CLEAN part realizes −S_i S_j.

The point: the recursion is the proof target; the open content is that the
CROSSING terms (parts sharing >1 vertex) must cancel.  We quantify them here.

NO number theory.
"""
import math, sys
from collections import defaultdict
from fractions import Fraction
import numpy as np


def set_partitions(c):
    c = list(c)
    if len(c) == 1:
        yield [c]; return
    f = c[0]
    for sm in set_partitions(c[1:]):
        for i, s in enumerate(sm):
            yield sm[:i] + [[f] + s] + sm[i + 1:]
        yield [[f]] + sm


def mu_partition(blocks):
    m = 1
    for B in blocks:
        b = len(B)
        m *= ((-1) ** (b - 1)) * math.factorial(b - 1)
    return m


def catalan(k):
    return math.comb(2 * k, k) // (k + 1)


def edge_flow_lines(edges, nb):
    E = len(edges)
    Bm = np.zeros((nb, E))
    for ei, (u, v) in enumerate(edges):
        Bm[v, ei] += 1; Bm[u, ei] -= 1
    u, s, vh = np.linalg.svd(Bm)
    tol = 1e-9; rank = int((s > tol).sum()); m = E - rank
    if m == 0:
        return [tuple()] * E, 0
    ns = vh[rank:]; lines = []
    for e in range(E):
        v = ns[:, e]
        if np.max(np.abs(v)) < 1e-7:
            lines.append(("ZERO",)); continue
        v = v / np.max(np.abs(v))
        for x in v:
            if abs(x) > 1e-7:
                if x < 0:
                    v = -v
                break
        lines.append(tuple(round(float(x), 6) for x in v))
    return lines, m


def is_even_series(edges, nb):
    adj = defaultdict(list)
    for (u, v) in edges:
        adj[u].append(v); adj[v].append(u)
    seen = {0}; stk = [0]
    while stk:
        x = stk.pop()
        for w in adj[x]:
            if w not in seen:
                seen.add(w); stk.append(w)
    if len(seen) != nb:
        return False
    lines, m = edge_flow_lines(edges, nb)
    if m == 0:
        return False
    if any(l == ("ZERO",) for l in lines):
        return False
    groups = defaultdict(int)
    for l in lines:
        groups[l] += 1
    return all(c % 2 == 0 for c in groups.values())


def build_graph(blocks, L):
    pos2blk = {}
    for bi, B in enumerate(blocks):
        for pos in B:
            pos2blk[pos] = bi
    return [(pos2blk[i], pos2blk[i + 1]) for i in range(L)], len(blocks), pos2blk


def compute_S():
    """S_k for k=1..5 by direct even-series enumeration; also return the per-pattern
       list (mu, pos2blk) for the decomposition analysis."""
    S = {0: 1}
    patterns = {}
    for k in range(1, 6):
        L = 2 * k
        tot = 0
        plist = []
        for blocks in set_partitions(range(L + 1)):
            edges, nb, p2b = build_graph(blocks, L)
            if any(u == v for (u, v) in edges):
                continue
            if not is_even_series(edges, nb):
                continue
            mu = mu_partition(blocks)
            tot += mu
            plist.append((mu, p2b))
        S[k] = tot
        patterns[k] = plist
    return S, patterns


def main():
    S, patterns = compute_S()
    print("=" * 72)
    print("(1) S_k = Σ_{even-series} μ   vs   free-cumulant recursion (REC)")
    print(f"    S_0..S_5 = {[S[k] for k in range(6)]}   (target (−1)^k C_k)")
    ok = True
    for k in range(1, 6):
        rec = -sum(S[i] * S[k - 1 - i] for i in range(k))
        tgt = (-1) ** k * catalan(k)
        flag = (S[k] == rec == tgt)
        ok &= flag
        print(f"    k={k}: S_k={S[k]:>5}  −ΣS_iS_j={rec:>5}  (−1)^kC_k={tgt:>5}  {'OK' if flag else 'X'}")
    print(f"    (REC) holds for k≤5: {ok}.  GF: x F² + F − 1 = 0,  F=(−1+√(1+4x))/(2x).")

    print("\n" + "=" * 72)
    print("(2) FIRST-RETURN decomposition.  Split at first r>0 with block(r)=block(0).")
    print("    i = #edges in first excursion (0..r) / 2 's role ... we grade by")
    print("    (a,b) = (positions in excursion interior structure).  CLEAN = excursion")
    print("    {0..r} and remainder {r..2k} share only vertex block(0).")
    print("    Compare CLEAN sum to the recursion −S_i·S_j.")
    for k in range(1, 6):
        L = 2 * k
        clean_by_r = defaultdict(int)
        cross_by_r = defaultdict(int)
        clean_tot = 0
        cross_tot = 0
        for (mu, p2b) in patterns[k]:
            b0 = p2b[0]
            # first return
            r = None
            for t in range(1, L + 1):
                if p2b[t] == b0:
                    r = t
                    break
            exc_pos = set(range(0, r + 1))
            rem_pos = set(range(r, L + 1))
            exc_v = set(p2b[p] for p in exc_pos)
            rem_v = set(p2b[p] for p in rem_pos)
            shared = exc_v & rem_v
            if shared == {b0}:
                clean_by_r[r] += mu
                clean_tot += mu
            else:
                cross_by_r[r] += mu
                cross_tot += mu
        print(f"\n  k={k}:  S_k={S[k]}   CLEAN Σμ={clean_tot}   CROSSING Σμ={cross_tot}")
        print(f"        clean by first-return r: {dict(sorted(clean_by_r.items()))}")
        print(f"        cross by first-return r: {dict(sorted(cross_by_r.items()))}")
    print("=" * 72)
    print("Reading: if CLEAN Σμ factored as −Σ_{i+j=k−1}S_iS_j the first-return would")
    print("prove (★★) directly.  The CROSSING remainder is the precise obstruction —")
    print("it must cancel (a sign-reversing involution on shared-vertex excursions).")


if __name__ == "__main__":
    main()
