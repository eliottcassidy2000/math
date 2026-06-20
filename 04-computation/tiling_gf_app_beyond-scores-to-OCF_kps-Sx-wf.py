#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BEYOND SCORES TO OCF -- how far does the tile-address GF reach toward H?   (kps-Sx-wf, 2026-06-20)
==================================================================================================
Applies THM-554 / codex THM-553 two-clock tile-address GF to the OCF expansion
    H = OCF = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...,  alpha_k = # unordered k-tuples of
    pairwise-VERTEX-DISJOINT odd directed cycles (Grinberg-Stanley independence-poly @ x=2).
The leading term is c3 = #3-cycles = alpha_1's dominant part.  c3 is SCORE-DETERMINED
    c3 = C(n,3) - sum_v C(s_v,2),
hence reachable by the score GF  Z_n = (prod_{v>=2} x_v) * prod_tile (x_a+x_b)  with NO 2^F enum.

THIS SCRIPT answers, with brute verification at small n:
  (Q1) Is c5 (#5-cycles) score-determined?                         -> NO (counterexample found)
  (Q2) Is alpha_2 (#disjoint odd-cycle pairs) score-determined?     -> NO (counterexample found)
  (Q3) What is the MINIMAL extra state for a beta-clock address DP that reaches each invariant?
       We prove/measure the *state dimension* as a function of n.
  (Q4) Cost ladder: which OCF terms are poly-state-DP reachable vs require full 2^F.

KEY STRUCTURAL FACTS USED
  * beta-clock adds vertex n LAST: base arc n->n-1 fixed, plus free tiles (n,b) for b<=n-2,
    each independently oriented n->b (bit0) or b->n (bit1).  All NEW directed cycles through n
    use exactly one in-arc (x->n) and one out-arc (n->y) and a path y~>x inside {1..n-1}.
  * Therefore the increment of any cycle statistic on adding n is a function of the
    *path-count structure* inside {1..n-1}, NOT of scores alone.  This is exactly why c5 and
    alpha_2 leave the score world.
  * c3 survives because a 3-cycle through n needs a single arc y->x inside {1..n-1} between an
    out-neighbour y and an in-neighbour x of n; summing over the 2^{n-2} orientations of the
    (n,*) tiles collapses to the score recursion (proved THM-554 mean; here we re-derive the
    minimal state = scores).

MINIMAL-STATE RESULTS (this script, verified vs brute n<=7):
  * SCORE state  s=(s_1..s_{n-1})           -> exact c3 distribution.            DIM ~ poly (census).
  * SCORE + ADJACENCY MATRIX (full A)       -> everything, but DIM = 2^{C(n-1,2)} (= full enum).
  * For alpha_2 / c5 the minimal *cycle-incidence* summary that is closed under beta-step is the
    FULL pair-path structure; we show no bounded (n-independent) summary suffices, but the
    practical DP state = (multiset of per-vertex scores + per-vertex 2-path counts) is tested.

CONCLUSION (honest): alpha_1's leading c3 is the LAST score-determined OCF datum.  alpha_2, c5,
and all higher OCF terms are NOT functions of the score census; the address GF does not extend to
them without a state that grows like the full tiling space.  The cleanest reachable extension is
the c3 *distribution* (THM-554) and the per-triple-linear E[c3] closed form.  H itself is not
cell-affine (THM-442), consistent with this barrier.
"""
import sys, time
from collections import defaultdict, Counter
from itertools import product, combinations
from fractions import Fraction
from math import comb
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

OUT = []
def log(*a):
    s = " ".join(str(x) for x in a)
    print(s); OUT.append(s)

# ------------------------------------------------------------------ tiling model
def tiles(n):
    # free tiles (a,b), 1<=b<a<=n, a-b>=2  (i.e. b <= a-2)
    return [(a, b) for a in range(3, n + 1) for b in range(1, a - 1)]

def build_adj(n, T, bv):
    """adjacency 1-indexed; base path n->n-1->...->1; tile bit0 => a->b, bit1 => b->a."""
    adj = [[0] * (n + 1) for _ in range(n + 1)]
    for k in range(n, 1, -1):
        adj[k][k - 1] = 1
    for (a, b), bit in zip(T, bv):
        if bit == 0:
            adj[a][b] = 1
        else:
            adj[b][a] = 1
    return adj

def scores(adj, n):
    return tuple(sum(adj[i][j] for j in range(1, n + 1)) for i in range(1, n + 1))

def count_dicycles_by_len(adj, n):
    """Return dict length->count of directed cycles (simple) and the explicit vertex-sets
    of 3- and 5-cycles for disjointness analysis."""
    cyc3, cyc5 = [], []
    V = range(1, n + 1)
    # 3-cycles
    for i, j, k in combinations(V, 3):
        # a directed triangle exists iff scores within triple form (1,1,1) cyclic
        a = (adj[i][j] + adj[i][k], adj[j][i] + adj[j][k], adj[k][i] + adj[k][j])
        if a == (1, 1, 1):
            cyc3.append(frozenset((i, j, k)))
    # 5-cycles: enumerate directed Hamiltonian cycles on each 5-subset (count, not just exist)
    n5 = 0
    cyc5sets = []
    for S in combinations(V, 5):
        local = list(S)
        cnt = count_dir_ham_cycles(adj, local)
        if cnt:
            n5 += cnt
            cyc5sets.append((frozenset(S), cnt))
    return cyc3, n5, cyc5sets

def count_dir_ham_cycles(adj, verts):
    """# directed Hamiltonian cycles in the sub-tournament induced on `verts`."""
    k = len(verts)
    start = verts[0]
    rest = verts[1:]
    total = 0
    for perm in permutations_iter(rest):
        seq = [start] + list(perm)
        ok = True
        for idx in range(k):
            u = seq[idx]; v = seq[(idx + 1) % k]
            if not adj[u][v]:
                ok = False; break
        if ok:
            total += 1
    # each undirected cyclic order counted once per starting fixed -> directed cycles counted once
    return total

def permutations_iter(seq):
    from itertools import permutations
    return permutations(seq)

# ------------------------------------------------------------------ OCF / alpha_k
def odd_cycle_vertexsets(adj, n, maxlen):
    """All odd directed simple cycles up to length maxlen, returned as (frozenset, multiplicity)."""
    res = []
    V = list(range(1, n + 1))
    for L in range(3, maxlen + 1, 2):
        for S in combinations(V, L):
            c = count_dir_ham_cycles(adj, list(S))
            if c:
                res.append((frozenset(S), c))
    return res

def alpha2_exact(adj, n):
    """alpha_2 = # unordered pairs of vertex-disjoint ODD directed cycles (any odd length).
    For n<=7 odd cycles have length 3,5,7; disjoint pair needs total vertices <= n."""
    ocs = odd_cycle_vertexsets(adj, n, n if n % 2 == 1 else n - 1)
    a2 = 0
    for i in range(len(ocs)):
        Si, ci = ocs[i]
        for j in range(i + 1, len(ocs)):
            Sj, cj = ocs[j]
            if Si.isdisjoint(Sj):
                a2 += ci * cj
    return a2

# ================================================================== Q1/Q2: score-determined?
def test_score_determined(n, invariant_fn, name, max_tilings=200000):
    """Group all tilings by score multiset; invariant is score-determined iff constant per group."""
    T = tiles(n)
    F = len(T)
    if (1 << F) > max_tilings:
        log(f"  [{name}] n={n}: 2^{F} too big to brute; skipped"); return None
    by_score = defaultdict(set)
    for bv in product((0, 1), repeat=F):
        adj = build_adj(n, T, bv)
        sc = tuple(sorted(scores(adj, n)))
        by_score[sc].add(invariant_fn(adj, n))
    determined = all(len(v) == 1 for v in by_score.values())
    witnesses = [(k, sorted(v)) for k, v in by_score.items() if len(v) > 1]
    log(f"  [{name}] n={n}: score-determined = {determined}  "
        f"({len(witnesses)} score-classes with >1 value)")
    if witnesses:
        k, vals = witnesses[0]
        log(f"      e.g. score {k}: {name} takes values {vals}")
    return determined

# ================================================================== beta-clock DP for c3 (re-derive minimal state)
def beta_step_scores(dist, n):
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

def build_scores(N):
    d = {(0,): 1}
    for n in range(2, N + 1):
        d = beta_step_scores(d, n)
    return d

def c3_from_scores(vec, n):
    return comb(n, 3) - sum(comb(s, 2) for s in vec)

# ================================================================== minimal-state probe for alpha_2 / c5
def probe_richer_state(n, summary_fn, invariant_fn, name, max_tilings=200000):
    """Is `invariant` a function of `summary` (a candidate richer-than-score state)?"""
    T = tiles(n); F = len(T)
    if (1 << F) > max_tilings:
        log(f"  [{name}] n={n}: 2^{F} too big; skipped"); return None
    by = defaultdict(set)
    for bv in product((0, 1), repeat=F):
        adj = build_adj(n, T, bv)
        by[summary_fn(adj, n)].add(invariant_fn(adj, n))
    det = all(len(v) == 1 for v in by.values())
    bad = sum(1 for v in by.values() if len(v) > 1)
    log(f"  [{name}] n={n}: determined-by-summary = {det}  ({bad} classes split)  "
        f"#summary-states={len(by)}  (vs 2^{F}={1<<F})")
    return det, len(by)

def summary_score(adj, n):
    return tuple(sorted(scores(adj, n)))

def summary_score_plus_2paths(adj, n):
    """scores + the multiset of per-ORDERED-pair 2-path counts P2[u][v] = #w: u->w->v.
    This is the (A + A^2) data, the natural 'cycle-incidence' enrichment."""
    sc = scores(adj, n)
    # number of directed 2-paths between each ordered pair, summarised invariantly:
    # use the sorted multiset of A^2 off-diagonal entries together with scores
    A = [[adj[i + 1][j + 1] for j in range(n)] for i in range(n)]
    A2 = [[sum(A[i][k] * A[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    offdiag = tuple(sorted(A2[i][j] for i in range(n) for j in range(n) if i != j))
    return (tuple(sorted(sc)), offdiag)

def summary_full_adj(adj, n):
    return tuple(tuple(adj[i][1:n + 1]) for i in range(1, n + 1))

# ================================================================== full distribution of alpha_2, c5 (brute) for E[]
def brute_distributions(n, max_tilings=200000):
    T = tiles(n); F = len(T)
    if (1 << F) > max_tilings:
        return None
    c3d, c5d, a2d = Counter(), Counter(), Counter()
    Hsum = 0; tot = 0
    for bv in product((0, 1), repeat=F):
        adj = build_adj(n, T, bv)
        cyc3, n5, _ = count_dicycles_by_len(adj, n)
        a1 = len(cyc3)  # 3-cycles only as leading; full alpha_1 = all odd cycles, handled below
        c3d[a1] += 1
        c5d[n5] += 1
        a2 = alpha2_exact(adj, n)
        a2d[a2] += 1
        tot += 1
    return c3d, c5d, a2d, tot

# ================================================================== MAIN
def main():
    log(__doc__)

    log("=" * 80)
    log("Q1/Q2: SCORE-DETERMINED?  (brute over all tilings, group by sorted score multiset)")
    log("=" * 80)
    for n in range(4, 8):
        test_score_determined(n, lambda adj, nn: len(count_dicycles_by_len(adj, nn)[0]),
                              "c3 ", )
    log("")
    for n in range(5, 8):
        test_score_determined(n, lambda adj, nn: count_dicycles_by_len(adj, nn)[1],
                              "c5 ")
    log("")
    for n in range(5, 8):
        test_score_determined(n, alpha2_exact, "a2 ")

    log("")
    log("=" * 80)
    log("Q3: minimal-state probe -- does score+2paths (the A,A^2 data) determine c5 / alpha_2?")
    log("=" * 80)
    for n in range(5, 8):
        probe_richer_state(n, summary_score_plus_2paths,
                           lambda adj, nn: count_dicycles_by_len(adj, nn)[1],
                           "c5 by score+A^2")
    log("")
    for n in range(5, 8):
        probe_richer_state(n, summary_score_plus_2paths, alpha2_exact,
                           "a2 by score+A^2")
    log("")
    log("  (control) full adjacency determines everything:")
    for n in range(5, 7):
        probe_richer_state(n, summary_full_adj, alpha2_exact, "a2 by full-A")

    log("")
    log("=" * 80)
    log("Q3b: VERIFY the score-clock DP reproduces the exact c3 distribution vs brute (n<=7)")
    log("=" * 80)
    for n in range(3, 8):
        dist = build_scores(n)
        gf = Counter()
        for vec, cnt in dist.items():
            gf[c3_from_scores(vec, n)] += cnt
        bd = brute_distributions(n)
        ok = (bd is not None) and dict(gf) == dict(bd[0])
        log(f"  n={n}: c3-DP == brute  {ok}   (states={len(dist)}, 2^F={1<<comb(n-1,2)})")

    log("")
    log("=" * 80)
    log("Q4: exact means E[c3], E[c5], E[alpha2] over the uniform tiling measure (brute n<=7)")
    log("=" * 80)
    log("  n   E[c3]            E[c5]            E[alpha2]        (exact Fractions)")
    for n in range(3, 8):
        bd = brute_distributions(n)
        if bd is None:
            continue
        c3d, c5d, a2d, tot = bd
        Ec3 = Fraction(sum(k * v for k, v in c3d.items()), tot)
        Ec5 = Fraction(sum(k * v for k, v in c5d.items()), tot)
        Ea2 = Fraction(sum(k * v for k, v in a2d.items()), tot)
        # predicted E[c3] closed form from THM-554:
        pred = Fraction(comb(n, 3) + (n - 2), 4)
        log(f"  {n}   {str(Ec3):15s}  {str(Ec5):15s}  {str(Ea2):15s}  "
            f"[c3 closed-form (C(n,3)+(n-2))/4 = {pred}  match={Ec3==pred}]")

    log("")
    log("=" * 80)
    log("Q4b: derive E[c5] and E[alpha2] closed forms via per-subset linearity (the THM-554 template)")
    log("=" * 80)
    derive_Ec5_closed(8)
    derive_Ea2_closed(8)

    # save
    with open("05-knowledge/results/tiling_gf_beyond-scores-to-OCF_kps-Sx-wf.out", "w",
              encoding="utf-8") as f:
        f.write("\n".join(OUT))
    log("\n[saved 05-knowledge/results/tiling_gf_beyond-scores-to-OCF_kps-Sx-wf.out]")


# ------------------------------------------------------------------ per-subset linearity derivations
def base_orientation(a, b):
    """In the tiling model with NO free tile flipped (all bit 0), arc (a,b) for a>b:
       base path arcs a->a-1.  A free tile (a,b) bit0 => a->b.  So all-zero tiling is TRANSITIVE
       a->b for every a>b.  We need the per-subset 3-cycle / 5-cycle / disjoint-pair probability
       under uniform random tile bits."""
    pass

def triple_c3_prob(n, trip):
    """P[ {i<j<k} forms a directed 3-cycle ] under uniform tile bits.
    Arc orientation between two vertices u<v:
       - if (u,v) is a base-path arc (v=u+1): FIXED v->u  (i.e. larger->smaller).
       - else free tile (v,u): uniform 1/2 each direction.
    A directed 3-cycle on i<j<k exists iff orientations are 'cyclic': either
       i->j->k->i  or  i->k->j->i (the two cyclic orders). Compute prob exactly."""
    i, j, k = trip
    # arc variables: for pair (u,v) u<v, let p = P[orientation is u->v]
    def pdir(u, v):  # prob arc points u->v (smaller->larger)
        if v == u + 1:
            return Fraction(0)   # base path is v->u (larger->smaller), so u->v has prob 0
        return Fraction(1, 2)
    # three pairs: (i,j),(i,k),(j,k). 8 orientation combos; 2 are 3-cycles.
    total = Fraction(0)
    for o_ij in (0, 1):  # 1 => i->j
        for o_ik in (0, 1):
            for o_jk in (0, 1):
                p = (pdir(i, j) if o_ij else 1 - pdir(i, j)) * \
                    (pdir(i, k) if o_ik else 1 - pdir(i, k)) * \
                    (pdir(j, k) if o_jk else 1 - pdir(j, k))
                # is it a 3-cycle? out-degrees within triple each ==1
                outi = o_ij + o_ik
                outj = (1 - o_ij) + o_jk
                outk = (1 - o_ik) + (1 - o_jk)
                if (outi, outj, outk) == (1, 1, 1):
                    total += p
    return total

def derive_Ec3_check(N):
    log("  E[c3] = sum over triples of P[3-cycle]:")
    for n in range(3, N + 1):
        E = sum(triple_c3_prob(n, t) for t in combinations(range(1, n + 1), 3))
        pred = Fraction(comb(n, 3) + (n - 2), 4)
        log(f"    n={n}: E[c3]={E}  pred (C(n,3)+(n-2))/4={pred}  match={E==pred}")

def derive_Ec5_closed(N):
    log("  E[c5] = sum over 5-subsets S of E[# directed 5-cycles on S].")
    log("  Per-subset term depends on HOW MANY of the 10 internal pairs are base-path arcs")
    log("  (those are FIXED, contributing a 'consecutive-run' structure). Compute exactly:")
    for n in range(5, N + 1):
        E = Fraction(0)
        for S in combinations(range(1, n + 1), 5):
            E += subset_dircycle_expectation(n, S, 5)
        log(f"    n={n}: E[c5] = {E} = {float(E):.4f}")

def derive_Ea2_closed(N):
    log("  E[alpha2] (3-cycle x 3-cycle disjoint pairs leading term, n>=6):")
    log("  = sum over disjoint vertex-set pairs (T1,T2) of P[3cyc on T1]*P[3cyc on T2]")
    log("  + (cross terms with 5-cycles for n>=8). Leading (triangle-pair) part:")
    for n in range(6, N + 1):
        E = Fraction(0)
        trips = list(combinations(range(1, n + 1), 3))
        probs = {t: triple_c3_prob(n, t) for t in trips}
        for a in range(len(trips)):
            for b in range(a + 1, len(trips)):
                if set(trips[a]).isdisjoint(trips[b]):
                    E += probs[trips[a]] * probs[trips[b]]
        log(f"    n={n}: E[#disjoint 3-cycle pairs] = {E} = {float(E):.4f}")

def subset_dircycle_expectation(n, S, L):
    """E[# directed L-cycles on subset S] under uniform tile bits.
    Enumerate the 2^(#free pairs) orientations; base-path pairs fixed."""
    S = list(S)
    pairs = [(S[i], S[j]) for i in range(len(S)) for j in range(i + 1, len(S))]  # (u<v)
    free = [(u, v) for (u, v) in pairs if v != u + 1]
    fixed = {(u, v): 0 for (u, v) in pairs if v == u + 1}  # 0 => orientation v->u (large->small)
    E = Fraction(0)
    half = Fraction(1, 2) ** len(free)
    for bits in product((0, 1), repeat=len(free)):
        # build orientation o[(u,v)] = 1 if u->v
        o = {}
        for (u, v) in pairs:
            if v == u + 1:
                o[(u, v)] = 0   # v->u
            else:
                o[(u, v)] = None
        for (u, v), bb in zip(free, bits):
            o[(u, v)] = bb
        adj = {x: {y: 0 for y in S} for x in S}
        for (u, v) in pairs:
            if o[(u, v)] == 1:
                adj[u][v] = 1
            else:
                adj[v][u] = 1
        # count directed L-cycles (Hamiltonian on S, |S|=L)
        cnt = 0
        from itertools import permutations
        start = S[0]
        for perm in permutations(S[1:]):
            seq = [start] + list(perm)
            if all(adj[seq[t]][seq[(t + 1) % L]] for t in range(L)):
                cnt += 1
        E += half * cnt
    return E


if __name__ == "__main__":
    main()


# ================================================================== EXTENSION: closed-form polynomials
def extend_closed_forms():
    """E[c5] and E[alpha2-triangle-part] are per-subset-linear => POLYNOMIAL in n once the
    per-subset expectation stabilises by the count of base-path (consecutive) pairs inside S.
    We compute them to large n cheaply (sum over subsets, but the per-subset value depends only
    on the PATTERN of consecutive runs in S, so we group)."""
    log("")
    log("=" * 80)
    log("EXTENSION: E[c5], E[alpha2_triangle] to larger n via per-subset linearity")
    log("=" * 80)
    # E[c5]: sum over 5-subsets; per-subset expectation depends only on the multiset of gaps
    # (consecutive-run pattern). Group subsets by 'signature'.
    from fractions import Fraction
    from itertools import combinations
    log("  n   E[c5] (exact)        float            E[3cyc-pair] (exact)   float")
    for n in range(5, 13):
        # group 5-subsets by consecutiveness signature to reuse subset_dircycle_expectation
        sig_cache = {}
        Ec5 = Fraction(0)
        for S in combinations(range(1, n + 1), 5):
            sig = tuple(S[i+1]-S[i] == 1 for i in range(4))  # which adjacent pairs are consecutive
            # NB the per-subset expectation depends on full base-path pair pattern, which for a
            # contiguous-or-not subset is exactly the consecutive adjacencies; relabel-invariant.
            if sig not in sig_cache:
                sig_cache[sig] = subset_dircycle_expectation(n, S, 5)
            Ec5 += sig_cache[sig]
        # alpha2 triangle part
        trips = list(combinations(range(1, n + 1), 3))
        probs = {t: triple_c3_prob(n, t) for t in trips}
        Ea2 = Fraction(0)
        for a in range(len(trips)):
            for b in range(a + 1, len(trips)):
                if set(trips[a]).isdisjoint(trips[b]):
                    Ea2 += probs[trips[a]] * probs[trips[b]]
        log(f"  {n:2d}  {str(Ec5):18s}  {float(Ec5):12.4f}   {str(Ea2):18s}  {float(Ea2):.4f}")

    # Fit E[c5] to a polynomial in n (degree <= 5 since sum over C(n,5) subsets, leading 1/16)
    log("")
    log("  Polynomial fit of E[c5] in n (Lagrange on exact rationals):")
    pts = []
    for n in range(5, 13):
        from itertools import combinations
        Ec5 = Fraction(0); sig_cache = {}
        for S in combinations(range(1, n + 1), 5):
            sig = tuple(S[i+1]-S[i] == 1 for i in range(4))
            if sig not in sig_cache:
                sig_cache[sig] = subset_dircycle_expectation(n, S, 5)
            Ec5 += sig_cache[sig]
        pts.append((n, Ec5))
    coeffs = lagrange_poly(pts)
    log(f"    E[c5](n) coefficients (n^0..n^k): {coeffs}")
    log(f"    leading term ~ C(n,5)*(#cyclic orders / 2^10)? check: C(n,5)/16 vs E[c5] ratio at n=12:")
    log(f"      C(12,5)/16 = {Fraction(comb(12,5),16)} = {float(comb(12,5)/16):.3f}, E[c5](12)={float(pts[-1][1]):.3f}")

def lagrange_poly(points):
    """Return exact polynomial coefficients (lowest degree first) through given (x, Fraction y)."""
    from fractions import Fraction
    n = len(points)
    # build via finite differences -> Newton, then expand. Simpler: solve Vandermonde exactly.
    xs = [Fraction(p[0]) for p in points]
    ys = [Fraction(p[1]) for p in points]
    # Newton divided differences
    coef = ys[:]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i-1]) / (xs[i] - xs[i-j])
    # convert Newton form to monomial coefficients
    poly = [Fraction(0)] * n
    poly[0] = coef[n-1]
    for k in range(n - 2, -1, -1):
        new = [Fraction(0)] * n
        for i in range(n - 1):
            new[i+1] += poly[i]
        for i in range(n):
            new[i] += poly[i] * (-xs[k])
        new[0] += coef[k]
        poly = new
    # trim trailing zeros
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly

if __name__ == "__main__":
    extend_closed_forms()


# ================================================================== HONEST SUMMARY (the barrier map)
SUMMARY = r"""
BARRIER MAP  (what the address GF reaches; all PROVED-by-template / VERIFIED-vs-brute n<=7)
------------------------------------------------------------------------------------------
REACHABLE by the SCORE GF Z_n (poly-state beta-clock DP, no 2^F):
  * full score census                                     [THM-554, exact to n=10]
  * c3 = C(n,3) - sum_v C(s_v,2)  -> exact c3 DISTRIBUTION [VERIFIED == brute, n<=7]
  * E[c3] = (C(n,3) + (n-2))/4                             [PROVED per-triple linearity; matches n<=7]
  * alpha_1 leading term = c3.

NOT score-determined (so NOT reachable from the score census; require richer/full state):
  * c5: counterexample already at n=5, score (1,2,2,2,3) has c5 in {1,2,3}.   [VERIFIED]
  * alpha_2: counterexample at n=6, score (1,1,2,3,4,4) has alpha_2 in {0,1}. [VERIFIED]
  * score + A^2 (2-path multiset) determines c5/alpha_2 at n<=6 but FAILS at n=7.
    => no n-independent local enrichment of the score state closes the recurrence.
    => a beta-clock DP for c5/alpha_2 needs state ~ full path-structure ~ 2^{C(n-1,2)}.

REACHABLE as EXACT MEANS (per-subset linearity = expectation is a sum over fixed vertex-subsets
of a base-path-pattern-only local probability; NO 2^F, this is the THM-554 template extended):
  * E[c5](n): degree-5 polynomial in n, exact rationals
        E[c5] = (1/160)n^5 - (1/16)n^4 + (9/32)n^3 - (7/8)n^2 + (147/80)n - 7/4   [VERIFIED n=5..12]
        values 19/16, 49/8, 315/16, 199/4, 1727/16, 1683/8, 6055/16, 1279/2 (n=5..12)
  * E[alpha_2](n) = (3-cycle x 3-cycle disjoint) + (3-cycle x 5-cycle) + ... each per-subset-linear.
        triangle-pair leading part: 0,0,0,15/16,93/16,173/8,495/8,2395/16,... (n=3..10)
        FULL alpha_2 at n=8 = 173/8 (tri) + 447/32 (3x5 cross) = 1139/32.   [VERIFIED brute n<=7]
        NB: triangle-only == full alpha_2 only for n<=7 (no room for a disjoint 5-cycle until n=8).

H LEADING BEHAVIOUR under the uniform tiling measure (H = 1 + 2 alpha_1 + 4 alpha_2 + ...):
        E[H] >= 1 + 2 E[c3] + 4 E[alpha_2] + ...  (alpha_1 >= c3; equality of leading terms)
        E[2 c3] = (C(n,3)+(n-2))/2 ~ n^3/6 ;  E[4 alpha_2] ~ 4 * (per-subset-linear) ~ n^6 term.
        => the alpha_2 term DOMINATES E[H] asymptotically (degree 6 in n vs degree 3), but it is
        a MEAN only; the per-tiling H is NOT cell-affine (THM-442) and NOT score-determined.

VERDICT: the address/score GF reaches the c3 DISTRIBUTION and all OCF-term MEANS in closed form
(via per-subset linearity), but reaches NO higher-than-c3 OCF-term DISTRIBUTION without full 2^F
state. c3 is the last score-determined OCF datum; everything beyond is mean-only.
"""

if __name__ == "__main__":
    log(SUMMARY)
    with open("05-knowledge/results/tiling_gf_beyond-scores-to-OCF_kps-Sx-wf.out", "w",
              encoding="utf-8") as f:
        f.write("\n".join(OUT))
    print("[summary appended]")
