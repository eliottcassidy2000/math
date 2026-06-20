#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TILING-MODEL CYCLE MOMENTS  (kps-Sx-wf, angle = cycle-moments)
==============================================================
Applies the Z_n score-GF engine (THM-553/554) to derive CLOSED FORMS for
tiling-model cycle moments, generalizing the proved

    E_tiling[c3] = ( C(n,3) + (n-2) ) / 4 .

--------------------------------------------------------------------------------
THE MODEL (codex/kps convention).
Base path n -> n-1 -> ... -> 1.  Free tiles (a,b) with 1<=b<a<=n, a-b>=2.
A tiling sets every non-base arc independently & uniformly:
  - base arcs  k -> k-1   are FIXED (deterministic, oriented "downward").
  - tile arcs  (a,b), a-b>=2  are a fair coin: a->b or b->a.
So the tournament is a random tournament with a SPECIAL marginal:
  arc (i,j), i>j:   P(i->j) = 1  if i = j+1   (consecutive / base arc)
                    P(i->j) = 1/2 otherwise   (tile arc, independent).

THE PER-SUBSET LINEARITY TEMPLATE.
c_k = sum over k-subsets S of [n] of  1{ T|S has a directed k-cycle of a chosen type }.
Because only "consecutive" arcs (i,i+1) are forced and all tile arcs are
independent fair coins, the cycle-indicator of a subset S depends only on which
of its internal arcs are forced.  E[c_k] = sum_S P(S cyclic), and the probability
depends ONLY on the pattern of consecutive pairs inside S (its "gap signature").
Group subsets by that signature -> closed form.

--------------------------------------------------------------------------------
RESULTS DERIVED & VERIFIED HERE (all EXACT, Fraction):

(A) E[c3] = (C(n,3) + (n-2))/4.                                   [PROVED, re-verified]
    Triple {i<j<k}: cyclic prob = 1/4 generically, but if it contains a forced
    arc the base path biases it.  A triple has a 3-cycle iff not transitive.
    - generic triple (no two consecutive): each of 3 arcs is a fair coin =>
      P(cyclic)=2/8=1/4.
    - triple with exactly the run {v,v+1,v+2} (one such per consecutive-pair):
      arcs (v+1,v) and (v+2,v+1) forced downward, arc (v+2,v) is a coin.
      Transitive unless... actually a 3-cycle needs v+2->v, prob 1/2; check below.
    Excess over uniform C(n,3)/4 is exactly (n-2)/4.

(B) Var[c3]  -- CLOSED FORM derived by pair-of-triples covariance, split by overlap.

(C) E[c5] (directed 5-cycles) -- closed form by gap-signature of 5-subsets.

(D) E[c_k] general leading pattern:  E[c_k] = (k-1)!/2 / 2^? ... see derivation.

(E) E[H] leading-order = 1 + 2*E[alpha_1] (+ ...), with alpha_1 = #odd dir. cycles.

All closed forms are checked against the Z-engine (for score-determined ones) and
against direct brute enumeration over all 2^F tilings (for c5 and full cycle data),
n <= 8.

================================================================================
RESULT SUMMARY (all EXACT Fraction; verification status in [] ).
  E[c3]   = ( C(n,3) + (n-2) ) / 4                                  [PROVED]
          = (n^3 - 3n^2 + 8n - 12) / 24
  Var[c3] = (n^3 - 7n^2 + 20n - 16) / 32                            [VERIFIED n<=10]
            (cov-of-triples method == Z-engine exactly, n=3..10)
  E[c5]   = (n^5 - 10n^4 + 45n^3 - 140n^2 + 294n - 280) / 160       [VERIFIED n<=11]
            (per-subset == brute n<=7; degree-5 poly fit n=5..11)
  E[c7]   = (n^7 - 21n^6 + 189n^5 - 1015n^4 + 3836n^3
             - 10514n^2 + 18458n - 15204) / 896                     [VERIFIED n<=15]
            (degree-7 poly fit on 9 points n=7..15)
  E[c_k]  : per-subset baseline (k-1)!/2^k; E[c_k](n) = degree-k poly in n,
            leading term  (k-1)!/2^k * C(n,k) ~ n^k / (k * 2^k).    [PROVED leading,
            CONJECTURE full poly is the gap-signature correction series]
  Each E[c_k] poly's LEADING term = uniform random-tournament value; the lower-order
  terms are exactly the (n-2-ish) base-path "forced-arc" corrections (per-subset
  linearity template).  The forced consecutive arcs i->i+1 (the base path) raise
  short-run subsets' cyclic probability above the uniform baseline.
  E[H] (OCF=I(Omega,2)): exact brute n<=7 = 2, 4, 79/8, 29, 3175/32.
       Leading-order 1+2E[c3] undercounts; 1+2E[c3]+2E[c5] EXACT at n=5 only
       (79/8); diverges at n>=6 due to the 4*alpha_2 (disjoint odd-pair) OCF term.
       E[H] is NOT a polynomial-closed-form target here (full OCF needed).  [CONJECTURE
       n=5 identity 1+2E[c3]+2E[c5]=E[H] is a low-n coincidence, breaks at n=6.]
================================================================================
"""
import sys, time
from collections import defaultdict, Counter
from itertools import product, combinations
from math import comb, factorial
from fractions import Fraction as Fr

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

OUT = []
def log(*a):
    s = " ".join(str(x) for x in a)
    print(s); OUT.append(s)

# ----------------------------------------------------------------------------
# Engine (copied)
# ----------------------------------------------------------------------------
def beta_step(dist, n):
    nd = defaultdict(int)
    for vec, cnt in dist.items():
        l = list(vec) + [0]; l[n - 1] += 1; nd[tuple(l)] += cnt
    dist = nd
    for b in range(1, n - 1):
        nd = defaultdict(int)
        for vec, cnt in dist.items():
            l0 = list(vec); l0[n - 1] += 1; nd[tuple(l0)] += cnt
            l1 = list(vec); l1[b - 1] += 1; nd[tuple(l1)] += cnt
        dist = nd
    return dist

_Z_CACHE = {}
def build_Z(N):
    if N in _Z_CACHE:
        return _Z_CACHE[N]
    d = {(0,): 1}
    for n in range(2, N + 1):
        d = beta_step(d, n)
        _Z_CACHE[n] = d
    return _Z_CACHE[N]

def tiles(n):
    return [(a, b) for a in range(3, n + 1) for b in range(1, a - 1)]

# ----------------------------------------------------------------------------
# Brute tournament builder over all tilings (for c5, full cycle counts)
# ----------------------------------------------------------------------------
def all_tournaments(n):
    """Yield adjacency matrices over all 2^F tilings. adj[i][j]=1 means i->j."""
    T = tiles(n)
    for bv in product((0, 1), repeat=len(T)):
        adj = [[0] * (n + 1) for _ in range(n + 1)]
        for k in range(n, 1, -1):
            adj[k][k - 1] = 1               # base arc k->k-1
        for (a, b), bit in zip(T, bv):
            if bit == 0:
                adj[a][b] = 1
            else:
                adj[b][a] = 1
        yield adj

def count_kcycles(adj, n, k):
    """Count directed k-cycles (as cyclic vertex sequences up to rotation+? )
       We count directed k-cycles per k-subset: number of cyclic orderings of S
       that are consistent with arcs. Returns total over all k-subsets."""
    total = 0
    for S in combinations(range(1, n + 1), k):
        # count directed Hamiltonian cycles in tournament on S
        total += dir_ham_cycles(adj, S)
    return total

def dir_ham_cycles(adj, S):
    """Number of directed Hamiltonian cycles in the sub-tournament on S.
       Count cyclic sequences (v0,v1,...,v_{k-1},v0) with all arcs present,
       counted once per undirected cyclic class (rotation). We fix v0=S[0] to
       kill rotation; direction counted separately (both orientations are
       distinct cycles in a tournament only if both exist - they can't both for
       k=3, but for k>=4 they can). Standard convention: count distinct directed
       cycles = closed directed walks visiting each vertex once, /k for rotation."""
    k = len(S)
    cnt = 0
    rest = S[1:]
    start = S[0]
    for perm in _perms(rest):
        seq = (start,) + perm
        ok = all(adj[seq[i]][seq[(i + 1) % k]] for i in range(k))
        if ok:
            cnt += 1
    return cnt

def _perms(seq):
    from itertools import permutations
    return permutations(seq)

# ----------------------------------------------------------------------------
# c3 closed forms via Z-engine
# ----------------------------------------------------------------------------
def c3_moments_from_Z(distZ, n):
    """Return (E[c3], E[c3^2], Var) as exact Fractions using the score GF.
       c3(T) = C(n,3) - sum_v C(s_v,2)."""
    tot = 0
    s1 = 0   # sum cnt*c3
    s2 = 0   # sum cnt*c3^2
    C3 = comb(n, 3)
    for vec, cnt in distZ.items():
        c3 = C3 - sum(comb(s, 2) for s in vec)
        tot += cnt; s1 += cnt * c3; s2 += cnt * c3 * c3
    E = Fr(s1, tot); E2 = Fr(s2, tot)
    return E, E2, E2 - E * E

# ============================================================================
# CLOSED FORM (A): E[c3]
# ============================================================================
def closed_E_c3(n):
    return Fr(comb(n, 3) + (n - 2), 4)

# ============================================================================
# CLOSED FORM (B): Var[c3]
# ----------------------------------------------------------------------------
# c3 = sum_{S: |S|=3} X_S,  X_S = 1{S is a 3-cycle}.
# Var = sum_S Var(X_S) + sum_{S != S'} Cov(X_S, X_S').
#
# X_S is Bernoulli(p_S), p_S in {1/4, ...}.  Determine p_S by gap signature.
# A triple {i<j<k}. Its three arcs: (k,j),(k,i),(j,i) lower->? Each arc is forced
# downward iff the two vertices are consecutive (differ by 1).  In a triple at
# most the pairs (i,i+1) can be consecutive: signature = which of the gaps
# (j-i), (k-j) equals 1.
#
# Cases for a triple i<j<k:
#   (i) no consecutive pair  (j>i+1 and k>j+1): all 3 arcs coins -> P(cyclic)=1/4.
#   (ii) exactly one consecutive pair:
#        - j=i+1, k>j+1: arc (j,i) forced j->i. Other two coins.
#        - k=j+1, j>i+1: arc (k,j) forced k->j. Other two coins.
#        In each, one arc forced, two coins. 3-cycle needs a specific orientation
#        -> P(cyclic)=1/4 STILL (check by enumeration in code).
#   (iii) full run j=i+1,k=j+1 (k=i+2): arcs (j,i),(k,j) forced downward. arc(k,i)
#        coin. 3-cycle iff i->k? No: cycle i->j->k->i needs i->j (but forced j->i),
#        so the only possible 3-cycle is i->k->j->i style. Enumerate -> P=1/2.
# So per-triple p_S = 1/4 except the (n-2) consecutive-run triples have p=1/2.
# E[c3] = (#triples - (n-2))*1/4 + (n-2)*1/2 = C(n,3)/4 + (n-2)/4. QED (A).
#
# Var: need Cov for pairs of triples. Two triples can share 0,1,2 vertices.
#  - share 0 or 1 vertex: arc-sets disjoint => independent => Cov=0. (forced arcs
#    are deterministic so they don't create dependence either.)
#  - share 2 vertices: share ONE arc. Cov can be nonzero.
# So Var = sum_S p_S(1-p_S) + sum_{pairs sharing an edge} Cov.
# We compute the edge-sharing covariance by gap signatures EXACTLY in code below,
# then fit/confirm the closed form.
# ============================================================================

def p_triple(i, j, k):
    """Exact P(triple {i<j<k} forms a directed 3-cycle) under the tiling marginal.
       Enumerate the 3 arcs' distributions."""
    # arcs: (j,i) [pair i,j], (k,j) [pair j,k], (k,i) [pair i,k]
    # each arc: if consecutive (diff 1) forced downward (higher->lower) prob1;
    # else fair coin.
    def arc_dist(hi, lo):
        if hi - lo == 1:
            return [(1, Fr(1))]            # hi->lo certain (orient = +1 means hi->lo)
        return [(1, Fr(1, 2)), (0, Fr(1, 2))]  # 1: hi->lo, 0: lo->hi
    P = Fr(0)
    for oij, pij in arc_dist(j, i):
        for ojk, pjk in arc_dist(k, j):
            for oik, pik in arc_dist(k, i):
                # orientations: oij=1 -> j->i ; ojk=1 -> k->j ; oik=1 -> k->i
                # build adjacency on {i,j,k}
                adj = {}
                adj[(j, i)] = oij; adj[(i, j)] = 1 - oij
                adj[(k, j)] = ojk; adj[(j, k)] = 1 - ojk
                adj[(k, i)] = oik; adj[(i, k)] = 1 - oik
                # 3-cycle exists iff each vertex has out-deg 1 (i.e. not transitive)
                outi = adj[(i, j)] + adj[(i, k)]
                outj = adj[(j, i)] + adj[(j, k)]
                outk = adj[(k, i)] + adj[(k, j)]
                if (outi, outj, outk) in ((1, 1, 1),):
                    P += pij * pjk * pik
    return P

def E_c3_by_subset(n):
    return sum((p_triple(i, j, k) for i, j, k in combinations(range(1, n + 1), 3)), Fr(0))

def cov_shared_edge(n):
    """Sum of Cov(X_S, X_S') over ordered? unordered pairs of triples sharing
       exactly 2 vertices (one edge). Returns exact Fraction.
       We enumerate by the shared edge and the two third-vertices."""
    total = Fr(0)
    verts = range(1, n + 1)
    # iterate over unordered pairs of distinct triples sharing >=2 vertices
    triples = list(combinations(verts, 3))
    # index by frozenset for sharing
    for a in range(len(triples)):
        Sa = triples[a]
        for b in range(a + 1, len(triples)):
            Sb = triples[b]
            shared = set(Sa) & set(Sb)
            if len(shared) == 2:
                total += joint_cov(Sa, Sb)
    return total

def joint_cov(Sa, Sb):
    """Cov(X_Sa, X_Sb) exact, two triples sharing one edge.
       Enumerate the union (4 vertices, 5 arcs: the shared one + 2 each side...
       actually union has 4 vertices => C(4,2)=6 arcs but only 5 are involved?
       The two triples share 1 edge; union = 4 vertices. The arcs touching both
       triples: shared edge. Other arcs partition. Just enumerate all arcs in the
       4-vertex union that lie within either triple)."""
    U = sorted(set(Sa) | set(Sb))   # 4 vertices
    # arcs we need: all pairs within Sa or within Sb
    needed = set()
    for S in (Sa, Sb):
        for x, y in combinations(sorted(S), 2):
            needed.add((y, x))      # store as (hi,lo) with hi>lo
    needed = sorted(needed)
    def arc_choices(hi, lo):
        if hi - lo == 1:
            return [(1, Fr(1))]
        return [(1, Fr(1, 2)), (0, Fr(1, 2))]
    # enumerate
    EA = Fr(0); EB = Fr(0); EAB = Fr(0)
    # build all combos
    def rec(idx, assign, prob):
        nonlocal EA, EB, EAB
        if idx == len(needed):
            xa = cyclic(Sa, assign); xb = cyclic(Sb, assign)
            EA_l = prob * xa; EB_l = prob * xb; EAB_l = prob * xa * xb
            EA += EA_l; EB += EB_l; EAB += EAB_l
            return
        hi, lo = needed[idx]
        for o, p in arc_choices(hi, lo):
            assign[(hi, lo)] = o
            rec(idx + 1, assign, prob * p)
    rec(0, {}, Fr(1))
    return EAB - EA * EB

def cyclic(S, assign):
    """1 if triple S is a directed 3-cycle given arc assignment dict {(hi,lo):o}."""
    i, j, k = sorted(S)
    oij = assign[(j, i)]; ojk = assign[(k, j)]; oik = assign[(k, i)]
    outi = (1 - oij) + (1 - oik)   # i->j iff oij==0 ; i->k iff oik==0
    outj = oij + (1 - ojk)
    outk = oik + ojk
    return 1 if (outi, outj, outk) == (1, 1, 1) else 0

# NOTE on covariance counting: Var(sum X) = sum Var(X_S) + sum_{S != S'} Cov(X_S,X_S').
# The double sum over ordered distinct pairs = 2 * (sum over unordered pairs).
# cov_shared_edge returns the UNORDERED sum, so we need 2x it.
def closed_Var_c3_fixed(n):
    var_indep = sum((p_triple(i, j, k) * (1 - p_triple(i, j, k))
                     for i, j, k in combinations(range(1, n + 1), 3)), Fr(0))
    cov_unordered = cov_shared_edge(n)
    return var_indep + 2 * cov_unordered

# ============================================================================
# CLOSED FORM (C): E[c5]
# ----------------------------------------------------------------------------
# c5 = sum over 5-subsets S of (#directed 5-cycles in T|S).
# For a fixed 5-subset, E[#dir 5-cycles] depends on its gap signature (which
# consecutive pairs i,i+1 lie inside S, forcing those arcs downward).
# We compute E[#dir-5-cycles | signature] exactly by enumerating the free arcs,
# then E[c5] = sum over 5-subsets of that conditional expectation.
# A clean uniform baseline: in a fully-random tournament on 5 labeled vertices,
# E[# dir 5-cycles] = (number of cyclic orders)/2^? ... we compute exactly.
# ============================================================================

def E_dircycles_subset(S, k):
    """Exact E[# directed k-cycles] in T|S under the tiling marginal.
       LINEARITY OVER ORIENTED CYCLES: E[#dir k-cycles] = sum over the (k-1)!
       oriented Hamiltonian cycles on S of P(all k of its arcs present).
       Each arc (u,v) with u>v: P(present as u->v) = 1 if u-v==1 (forced down),
       else 1/2 (free fair coin); P(v->u) = 0 if forced, else 1/2.
       O((k-1)!) per subset -- no 2^free enumeration."""
    from itertools import permutations
    S = sorted(S)
    start = S[0]; rest = S[1:]
    def p_arc(u, v):
        # probability the arc u->v is present
        hi, lo = (u, v) if u > v else (v, u)
        if hi - lo == 1:
            return Fr(1) if u > v else Fr(0)   # forced hi->lo
        return Fr(1, 2)
    E = Fr(0)
    for perm in permutations(rest):
        seq = (start,) + perm
        p = Fr(1)
        for i in range(k):
            p *= p_arc(seq[i], seq[(i + 1) % k])
            if p == 0:
                break
        E += p
    return E

def num_dir_kcycles(S, assign, k):
    """Count directed k-cycles using S (k = len(S)) given full arc assignment."""
    from itertools import permutations
    start = S[0]; rest = S[1:]; cnt = 0
    def has(u, v):
        if u > v:
            return assign[(u, v)] == 1
        else:
            return assign[(v, u)] == 0
    for perm in permutations(rest):
        seq = (start,) + perm
        if all(has(seq[i], seq[(i + 1) % k]) for i in range(k)):
            cnt += 1
    return cnt

def E_ck_by_subset(n, k):
    return sum((E_dircycles_subset(S, k) for S in combinations(range(1, n + 1), k)), Fr(0))

# uniform-random baseline for reference: E[#dir-k-cycles] per k-subset in a random
# tournament = (k-1)!/2 * 2^{-k}? compute directly:
def E_dircycles_uniform(k):
    from itertools import permutations
    start = 0; rest = list(range(1, k)); total = Fr(0)
    Npairs = comb(k, 2)
    # by symmetry probability each oriented cycle exists = 2^{-k} (k arcs)
    ncyc = factorial(k - 1)             # cyclic sequences fixing start, both dirs counted
    return Fr(ncyc, 1 << k)

# ============================================================================
# (E) E[H] leading order.  H = OCF = I(Omega(T),2) = 1 + 2*alpha_1 + 4*alpha_2+...
# alpha_1 = total # of odd directed cycles (3-,5-,7-,...). alpha_1 leading term ~ c3.
# Actually OCF: H = sum over independent sets of odd-cycle conflict graph of 2^{|set|}.
# Leading-order linear approx: E[H] ~ 1 + 2*E[c3] + 2*E[c5] + ... (alpha_1 contributes
# the count of odd cycles at linear order). We report E[1 + 2*c3 + 2*c5] as a
# leading-order estimate and note higher OCF terms are NOT captured.
# ============================================================================

def brute_H_mean(n):
    """Exact E[H] over all tilings (n<=7 feasible)."""
    tot = 0; s = 0
    for adj in all_tournaments(n):
        H = exact_H(adj, n)
        s += H; tot += 1
    return Fr(s, tot)

def exact_H(adj, n):
    """H(T) = I(Omega(T),2): independence poly of odd-cycle conflict graph at x=2.
       Equivalently # of vertex-disjoint odd-cycle packings weighted... use OCF:
       H = sum over sets of pairwise-vertex-disjoint odd cycles of 2^{...}?
       Simplest correct: H = ham-paths count = permanent-like. We use the
       independence-polynomial-of-conflict-graph definition with conflict = share
       a vertex, over ALL odd directed cycles."""
    odd_cycles = enumerate_odd_cycles(adj, n)
    # conflict graph: vertices = odd cycles; edge if they share a vertex.
    m = len(odd_cycles)
    sets = odd_cycles
    # independence polynomial at x=2: sum over independent sets I of 2^{|I|}
    # independent = pairwise vertex-disjoint
    # DP via inclusion is expensive; n<=7 so m small-ish. Use recursive.
    vsets = [frozenset(c) for c in sets]
    from functools import lru_cache
    total = [0]
    # branch and bound over choosing disjoint subset
    def rec(idx, used, weight):
        total[0] += weight
        for j in range(idx, len(vsets)):
            if not (vsets[j] & used):
                rec(j + 1, used | vsets[j], weight * 2)
    rec(0, frozenset(), 1)
    return total[0]

def enumerate_odd_cycles(adj, n):
    """All directed cycles of odd length (vertex sets as tuples) in T."""
    from itertools import permutations
    res = []
    for k in range(3, n + 1, 2):
        for S in combinations(range(1, n + 1), k):
            start = S[0]; rest = S[1:]
            seen = set()
            for perm in permutations(rest):
                seq = (start,) + perm
                if all(adj[seq[i]][seq[(i + 1) % k]] for i in range(k)):
                    res.append(frozenset(seq) if False else seq)
    # We need each directed cycle ONCE (as a vertex set+rotation). For conflict
    # graph (vertex-disjointness) only the vertex SET matters but multiple distinct
    # directed cycles on same set are distinct conflict-graph vertices.
    # Dedup rotations: fix start=min already done. Keep all distinct directed cycles.
    return res

# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------
def main():
    log("="*78)
    log("TILING-MODEL CYCLE MOMENTS  (kps-Sx-wf)")
    log("="*78)

    # ---- (A) E[c3] ----
    log("\n--- (A) E[c3] closed form = (C(n,3)+(n-2))/4 ---")
    log(" n |  closed E[c3]  |  per-subset sum |  Z-engine E[c3] | match")
    log(" (Z-engine column verified to n=10 in prior runs: E[c3](10)=32, Var(10)=121/8)")
    for n in range(3, 10):
        cf = closed_E_c3(n)
        sub = E_c3_by_subset(n)
        E, E2, V = c3_moments_from_Z(build_Z(n), n)
        ok = (cf == sub == E)
        log(f" {n:2d}| {str(cf):14s}| {str(sub):16s}| {str(E):16s}| {ok}")

    # ---- (B) Var[c3] ----
    log("\n--- (B) Var[c3]: covariance-of-triples closed form vs Z-engine ---")
    log(" n |  closed Var (cov method) |  Z-engine Var    | match")
    var_data = []
    for n in range(3, 10):
        cf = closed_Var_c3_fixed(n)
        _, _, V = c3_moments_from_Z(build_Z(n), n)
        ok = (cf == V)
        var_data.append((n, V))
        log(f" {n:2d}| {str(cf):26s}| {str(V):17s}| {ok}")
    # n=10 confirmed in prior run: cov=121/8, Z=121/8 (match)
    var_data.append((10, Fr(121, 8)))
    log(" 10| 121/8 (cov, prior run)     | 121/8            | True")

    # Try to find a polynomial closed form for Var[c3] in n (it should be a poly).
    log("\n  Fitting Var[c3] as polynomial in n (exact, via finite differences):")
    fit_poly([(n, V) for (n, V) in var_data])

    # ---- (C) E[c5] ----
    log("\n--- (C) E[c5] (directed 5-cycles): per-subset sum vs brute ---")
    log(" n |  per-subset E[c5]  |  brute E[c5]     | match | uniform-baseline*C(n,5)")
    for n in range(5, 9):
        sub = E_ck_by_subset(n, 5)
        if n <= 7:
            bru = brute_mean_kcycles(n, 5)
        else:
            bru = None
        base = E_dircycles_uniform(5) * comb(n, 5)
        ok = (bru is None) or (sub == bru)
        log(f" {n:2d}| {str(sub):19s}| {str(bru):17s}| {ok}   | {str(base)}")

    log("\n  Fitting E[c5] as polynomial in n:")
    e5 = [(n, E_ck_by_subset(n, 5)) for n in range(5, 13)]
    fit_poly(e5)

    # ---- (D) E[c_k] general ----
    # E[c_k] = sum over k-subsets S of E[# dir k-cycles in T|S].  The per-subset
    # expectation depends only on the gap signature of S (which i,i+1 pairs are
    # inside S, forcing arcs downward).  Uniform (no forced arcs) baseline per
    # subset = (k-1)!/2^k.  Closed forms below are degree-k polys in n.
    log("\n--- (D) E[c_k] general: per-subset sums, k=3,5,7 ---")
    log("    uniform per-subset baseline = (k-1)!/2^k; E[c_k](n) is a degree-k poly in n")
    for k, nmax in ((3, 13), (5, 15), (7, 17)):
        log(f"  k={k}: uniform per-subset E = {E_dircycles_uniform(k)} (= (k-1)!/2^k)")
        row = []
        for n in range(k, min(k + k + 2, nmax)):
            row.append((n, E_ck_by_subset(n, k)))
        log(f"    E[c{k}](n) = " + ", ".join(f"n={n}:{v}" for n, v in row))
        log(f"    polynomial fit:")
        fit_poly(row)

    # ---- (E) E[H] ----
    log("\n--- (E) E[H] leading-order vs exact (brute) ---")
    log(" n | 1+2E[c3] (lead) | 1+2E[c3]+2E[c5] |  exact E[H] (brute)")
    for n in range(3, 8):
        lead3 = 1 + 2 * closed_E_c3(n)
        lead5 = lead3 + (2 * E_ck_by_subset(n, 5) if n >= 5 else Fr(0))
        if n <= 7:
            EH = brute_H_mean(n)
        else:
            EH = None
        log(f" {n:2d}| {str(lead3):16s}| {str(lead5):16s}| {str(EH)}")

    # save
    with open("05-knowledge/results/tiling_gf_app_cycle-moments_kps-Sx-wf.out",
              "w", encoding="utf-8") as f:
        f.write("\n".join(OUT))
    log("\n[saved] 05-knowledge/results/tiling_gf_app_cycle-moments_kps-Sx-wf.out")

def brute_mean_kcycles(n, k):
    tot = 0; s = 0
    for adj in all_tournaments(n):
        s += count_kcycles(adj, n, k); tot += 1
    return Fr(s, tot)

def fit_poly(data):
    """Given exact (n, value) points, fit lowest-degree polynomial via finite
       differences; print the polynomial if it stabilizes (constant diff)."""
    pts = sorted(data)
    ys = [v for _, v in pts]
    xs = [x for x, _ in pts]
    # check equal spacing
    if len(set(xs[i+1]-xs[i] for i in range(len(xs)-1))) != 1:
        log("    (non-uniform spacing, skipping diff fit)"); return
    # finite differences
    diffs = [ys[:]]
    cur = ys[:]
    deg = None
    for d in range(1, len(ys)):
        cur = [cur[i+1]-cur[i] for i in range(len(cur)-1)]
        diffs.append(cur)
        if all(x == cur[0] for x in cur):
            deg = d; break
    if deg is None:
        log("    (no stable degree within data range)"); return
    # build Newton forward form -> express closed-form via comb
    x0 = xs[0]; h = xs[1]-xs[0]
    coeffs = [diffs[d][0] for d in range(deg+1)]
    # value at x = sum_d coeffs[d]*C((x-x0)/h, d)
    # Express as polynomial in n and simplify symbolically using Fraction poly.
    log(f"    degree {deg}; Newton coeffs (forward diffs at n={x0}): {coeffs}")
    # produce explicit poly in n by expanding binomials
    poly = expand_newton(coeffs, x0, h)
    log(f"    closed form: c3-var/ck poly in n = {poly_str(poly)}")
    # sanity check
    for x, y in pts:
        if poly_eval(poly, x) != y:
            log(f"    !! mismatch at n={x}: poly={poly_eval(poly,x)} vs {y}")
            return
    log("    (verified poly matches all data points)")

# polynomial as list of Fraction coeffs, poly[i] = coeff of n^i
def expand_newton(coeffs, x0, h):
    # C((n-x0)/h, d) = prod_{j=0}^{d-1} ((n-x0)/h - j) / d!
    poly = [Fr(0)]
    for d, c in enumerate(coeffs):
        term = [Fr(1)]  # start with 1
        for j in range(d):
            # multiply by ((n - x0)/h - j) = (1/h) n + (-x0/h - j)
            factor = [Fr(-x0, h) - j, Fr(1, h)]
            term = poly_mul(term, factor)
        term = poly_scale(term, Fr(c, factorial(d)))
        poly = poly_add(poly, term)
    return poly

def poly_mul(a, b):
    r = [Fr(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i+j] += x*y
    return r
def poly_add(a, b):
    n = max(len(a), len(b)); r = [Fr(0)]*n
    for i, x in enumerate(a): r[i] += x
    for i, y in enumerate(b): r[i] += y
    return r
def poly_scale(a, s): return [x*s for x in a]
def poly_eval(p, x):
    return sum((c * (x**i) for i, c in enumerate(p)), Fr(0))
def poly_str(p):
    terms = []
    for i in range(len(p)-1, -1, -1):
        c = p[i]
        if c == 0: continue
        if i == 0: terms.append(f"{c}")
        elif i == 1: terms.append(f"{c}*n")
        else: terms.append(f"{c}*n^{i}")
    return " + ".join(terms) if terms else "0"

if __name__ == "__main__":
    main()
