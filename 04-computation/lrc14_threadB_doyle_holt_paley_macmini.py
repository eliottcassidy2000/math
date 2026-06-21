#!/usr/bin/env python3
"""
lrc14_threadB_doyle_holt_paley_macmini.py
mac-mini-2026-06-21 -- THREAD B (the user's Doyle-Holt lead, tournament side).

QUESTION: does the Doyle-Holt / half-arc-transitive lens give the PALEY H-max
extremality (the tournament side of the shared theta' wall), or is it an
OBSTRUCTION?

THE TENSION (from kps HYP-2748 + canon THM-126/438/472):
  * Paley QR_7 is the UNIQUE H-maximizer among Z_7 circulants (THM-126), and the
    tournament-Hadamard ceiling det(I+S) is attained by switchings of DOUBLY-
    REGULAR TOURNAMENTS (DRTs) (THM-472).
  * Paley is ARC-TRANSITIVE as a digraph (vertex- AND arc-transitive). It is
    SELF-CONVERSE: inversion i -> -i realizes T <-> T^op.
  * HALF-ARC-TRANSITIVE = vertex- & edge-transitive but TWO arc-orbits = the
    converse Z_2 is NOT realized by Aut = NON-self-converse (NS).
  * So Paley (SC, arc-transitive, 1 arc-orbit) and the half-arc object (NS,
    2 arc-orbits) sit at OPPOSITE ends of the symmetry spectrum.

THREAD-B HYPOTHESIS (HYP-B): the H-max tournament is characterized by the
ARC-TRANSITIVE / maximal-symmetry end (1 arc-orbit), and half-arc-transitivity
is precisely the *defective* sibling that CANNOT be H-max. If true, Doyle-Holt
is the OBSTRUCTION, not the route: the H-max wants the converse Z_2 REALIZED
(SC), the half-arc object has it UNREALIZED (NS).

We test this concretely. For a tournament T (skew matrix S = A - A^T) define:
  - the DIGRAPH automorphism group Aut(T) (permutations preserving all arcs),
  - the number of ARC-ORBITS of Aut(T) on the directed arcs,
  - whether T is self-converse (some perm sends every arc x->y to y->x),
  - H(T) (Hamiltonian-path count) and the det(I+S) ceiling value.

Then we cross-tabulate symmetry-class vs H across:
  (A) all 8 circulants on Z_7 (THM-126 set),
  (B) all DRTs at n=7 (Paley, unique up to iso),
  (C) all DRTs at n=11 (there are 2 non-iso DRTs at n=11 -- Paley QR_11 and
      another -- a sharp test: are BOTH H-max? are BOTH arc-transitive?),
  (D) all DRTs at n=15 if feasible (there are 3 skew-Hadamard classes of order
      16 -> several DRTs; only some are arc-transitive). This is THE decisive
      test: do NON-arc-transitive DRTs still attain the H ceiling?

A DEAD-END verdict (Doyle-Holt does NOT pull the extremality) is just as
valuable as a route -- we want the precise obstruction.
"""
import sys
from itertools import permutations, combinations
import numpy as np
sys.stdout.reconfigure(line_buffering=True)

def banner(t):
    print("\n" + "=" * 74 + "\n  " + t + "\n" + "=" * 74)

# ----------------------------------------------------------------------------
# Tournament basics. A tournament on n vertices = adjacency A (A[i][j]=1 iff arc i->j).
# ----------------------------------------------------------------------------
def circulant_adj(n, S):
    """A[i][j]=1 iff (j-i) mod n in S."""
    S = set(s % n for s in S)
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S:
                A[i][j] = 1
    return A

def is_tournament(A):
    n = len(A)
    for i in range(n):
        if A[i][i] != 0:
            return False
        for j in range(i + 1, n):
            if A[i][j] + A[j][i] != 1:
                return False
    return True

def skew(A):
    n = len(A)
    return np.array([[A[i][j] - A[j][i] for j in range(n)] for i in range(n)], dtype=float)

def det_I_plus_S(A):
    n = len(A)
    S = skew(A)
    return round(float(np.linalg.det(np.eye(n) + S)))

def H_paths(A):
    """Exact Hamiltonian-path count via Held-Karp DP (counts directed Ham paths
    in the tournament; H(T) in canon = number of Hamiltonian paths)."""
    n = len(A)
    # dp[mask][v] = number of directed Ham paths over vertices in mask ending at v
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            if not (mask >> v) & 1:
                continue
            cur = dp[mask][v]
            if cur == 0:
                continue
            for w in range(n):
                if (mask >> w) & 1:
                    continue
                if A[v][w]:  # arc v -> w continues the path
                    dp[mask | (1 << w)][w] += cur
    return sum(dp[full][v] for v in range(n))

# ----------------------------------------------------------------------------
# Automorphisms (exact, brute force for small n; pruned for n<=11).
# An automorphism is a permutation p with A[p[i]][p[j]] = A[i][j] for all i,j.
# ----------------------------------------------------------------------------
def automorphisms(A, cap=None):
    n = len(A)
    auts = []
    # score-preserving pruning: out-degree must be preserved
    outdeg = [sum(A[i]) for i in range(n)]
    by_deg = {}
    for v in range(n):
        by_deg.setdefault(outdeg[v], []).append(v)
    # backtracking
    perm = [-1] * n
    used = [False] * n
    order = sorted(range(n), key=lambda v: len(by_deg[outdeg[v]]))

    def consistent(i_orig, cand):
        # check arcs to already-placed vertices
        for j_orig in range(n):
            pj = perm[j_orig]
            if pj == -1:
                continue
            if A[cand][pj] != A[i_orig][j_orig]:
                return False
            if A[pj][cand] != A[j_orig][i_orig]:
                return False
        return True

    def bt(k):
        if cap is not None and len(auts) >= cap:
            return
        if k == n:
            auts.append(perm[:])
            return
        v = order[k]
        for cand in by_deg[outdeg[v]]:
            if used[cand]:
                continue
            if consistent(v, cand):
                perm[v] = cand
                used[cand] = True
                bt(k + 1)
                used[cand] = False
                perm[v] = -1

    bt(0)
    return auts

def arc_orbits(A, auts):
    n = len(A)
    arcs = [(i, j) for i in range(n) for j in range(n) if A[i][j]]
    arc_idx = {a: k for k, a in enumerate(arcs)}
    parent = list(range(len(arcs)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    for p in auts:
        for (i, j) in arcs:
            img = (p[i], p[j])
            union(arc_idx[(i, j)], arc_idx[img])
    return len({find(k) for k in range(len(arcs))})

def vertex_orbits(A, auts):
    n = len(A)
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for p in auts:
        for v in range(n):
            ra, rb = find(v), find(p[v])
            if ra != rb:
                parent[ra] = rb
    return len({find(v) for v in range(n)})

def is_self_converse(A):
    """True iff some permutation p sends T to T^op: A[p[i]][p[j]] = A[j][i].
    Brute force with degree pruning (out-degree of T^op = in-degree of T)."""
    n = len(A)
    outdeg = [sum(A[i]) for i in range(n)]
    indeg = [sum(A[j][i] for j in range(n)) for i in range(n)]
    # in T^op, out-degree of vertex v = indeg_T(v); we need a perm matching
    # vertices of T (out-degree outdeg) to image whose required out-degree is indeg.
    # p must satisfy A[p[i]][p[j]] = A[j][i]. Then out-degree of i in T^op view:
    # sum_j A[j][i] = indeg[i] must equal outdeg[p[i]].
    cand_by = {}
    for v in range(n):
        cand_by.setdefault(indeg[v], None)
    # build candidate lists: p[i] must have outdeg = indeg[i]
    out_by = {}
    for v in range(n):
        out_by.setdefault(outdeg[v], []).append(v)
    perm = [-1] * n
    used = [False] * n
    order = sorted(range(n), key=lambda i: len(out_by.get(indeg[i], [])))
    def consistent(i_orig, cand):
        for j_orig in range(n):
            pj = perm[j_orig]
            if pj == -1:
                continue
            # need A[cand][pj] == A[j_orig][i_orig]
            if A[cand][pj] != A[j_orig][i_orig]:
                return False
            if A[pj][cand] != A[i_orig][j_orig]:
                return False
        return True
    def bt(k):
        if k == n:
            return True
        i = order[k]
        for cand in out_by.get(indeg[i], []):
            if used[cand]:
                continue
            if consistent(i, cand):
                perm[i] = cand
                used[cand] = True
                if bt(k + 1):
                    return True
                used[cand] = False
                perm[i] = -1
        return False
    return bt(0)

def summarize(name, A, cap_aut=200000):
    auts = automorphisms(A, cap=cap_aut)
    naut = len(auts)
    vo = vertex_orbits(A, auts)
    ao = arc_orbits(A, auts)
    sc = is_self_converse(A)
    H = H_paths(A)
    det = det_I_plus_S(A)
    vt = (vo == 1)
    at = (vt and ao == 1)
    half = (vt and ao == 2 and not sc)  # half-arc analog: VT, 2 arc-orbits, NS
    print(f"  {name:22s} |Aut|={naut:6d}  vO={vo} arcO={ao}  "
          f"VT={'Y' if vt else 'n'} ArcTrans={'Y' if at else 'n'} "
          f"SelfConv={'Y' if sc else 'n'}  H={H}  det(I+S)={det}")
    return dict(name=name, naut=naut, vo=vo, ao=ao, sc=sc, H=H, det=det,
                vt=vt, at=at, half=half)

# ----------------------------------------------------------------------------
# (A) all 8 circulants on Z_7
# ----------------------------------------------------------------------------
banner("(A) all 8 circulant tournaments on Z_7  (THM-126 maximizer set)")
results_A = []
m = 3
halfsets = []
# circulant tournament on Z_7: choose one of each inverse pair {s, 7-s} for s=1,2,3
for bits in range(8):
    S = []
    for k, s in enumerate([1, 2, 3]):
        if (bits >> k) & 1:
            S.append(s)
        else:
            S.append(7 - s)
    A = circulant_adj(7, S)
    assert is_tournament(A)
    r = summarize(f"Z7 S={sorted(S)}", A)
    r['Sset'] = sorted(S)
    results_A.append(r)

Hmax_A = max(r['H'] for r in results_A)
print(f"\n  H-max among Z_7 circulants = {Hmax_A}")
for r in results_A:
    if r['H'] == Hmax_A:
        print(f"    H-MAX attainer: S={r['Sset']}  ArcTrans={r['at']}  SelfConv={r['sc']}")

print("\n  >> Among Z_7 circulants: is ARC-TRANSITIVE <=> H-MAX?")
at_set = set(tuple(r['Sset']) for r in results_A if r['at'])
hm_set = set(tuple(r['Sset']) for r in results_A if r['H'] == Hmax_A)
print(f"     arc-transitive set = {sorted(at_set)}")
print(f"     H-max set          = {sorted(hm_set)}")
print(f"     EQUAL? {at_set == hm_set}")
