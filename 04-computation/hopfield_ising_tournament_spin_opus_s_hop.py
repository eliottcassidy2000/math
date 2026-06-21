#!/usr/bin/env python3
"""
hopfield_ising_tournament_spin_opus_s_hop.py
=============================================
THREAD: Hopfield/Ising spin-glass <-> tournament as a spin system.

Repo anchors:
  - THM-554/555: Z_n(x) = (prod x_v) * prod_tiles (x_a+x_b); c3 = C(n,3) - sum_v C(s_v,2).
                 "cut cheap (scores), cycle dear (OCF beyond c3)."
  - THM-126:     Paley/regular tournament uniquely maximizes H on Z_7 (S={1,2,4}).
  - definitions: OCF = signed odd-cycle sum; H(T) = Redei Ham-path count.

HYPOTHESES TESTED (falsifiable):

H1 (Ising energy = score variance).
    Define arc spins s_e = +-1 (orientation of each edge of K_n).
    The "imbalance energy"  E_score(T) = sum_v (s_v - sbar)^2   (sbar=(n-1)/2)
    is an exact AFFINE function of c3:  c3 = C(n,3) - sum_v C(s_v,2),
    and  sum_v C(s_v,2) = (1/2)(sum s_v^2 - C(n,2)).
    => maximizing c3 == minimizing E_score == "most balanced spins".
    PREDICTION: argmax c3 over all tournaments at n=7 is exactly the regular
    (3,3,3,3,3,3,3) score class, i.e. the Paley/regular tournament class. TEST it.

H2 (E_score is a genuine 2-spin Ising Hamiltonian on arc spins).
    sum_v s_v^2 should expand as  const + sum_{e!=f sharing a vertex} J * s_e s_f
    with a single coupling constant J (a 2-body Ising interaction whose graph is
    the LINE GRAPH of K_n, edges of K_n adjacent iff they share a vertex).
    PREDICTION: c3 = A + B * sum_{e~f} (signed) s_e s_f  for constants A,B.  TEST.

H3 (Regular = Hopfield ground state under Hebbian AFM coupling).
    Treat E_Ising = -sum_{e~f} J_ef s_e s_f with J_ef chosen so the energy
    equals +E_score (antiferromagnetic on the line graph). Run synchronous /
    asynchronous Hopfield descent from random arc-spin states; does it converge
    to a regular (balanced-score) tournament = max-c3 attractor?  TEST n=5,7.

H4 (cut vs cycle == perceptron vs hidden-layer).
    Scores s_v are a LINEAR (perceptron) readout of arc spins: s_v = linear in s_e.
    c3 is QUADRATIC in scores => quartic-ish / 2-body in arc spins but still
    score-DETERMINED (THM-554: c3 is last score-determined OCF datum).
    The NEXT OCF datum (signed 5-cycle / full H) is NOT a function of scores alone.
    PREDICTION: a single-layer perceptron (linear readout of arc spins -> score,
    threshold) can perfectly predict c3-rank but CANNOT predict H within a fixed
    score class (needs cycle-space = "hidden layer"). TEST: within the regular
    score class at n=7, do scores determine H? (They must NOT, if cut!=cycle.)
"""

from itertools import combinations
from math import comb
from fractions import Fraction
from collections import defaultdict
import random
import sys

random.seed(20260620)
sys.stdout.reconfigure(line_buffering=True)

# ---------------------------------------------------------------------------
# Tournament representation.
# Edges of K_n in fixed order; spin s_e in {+1,-1}.
# Convention: edge (i,j) with i<j.  s=+1 means i->j, s=-1 means j->i.
# ---------------------------------------------------------------------------

def edges_of(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

def spins_to_adj(n, edges, spins):
    """spins: dict or list aligned with edges; A[i][j]=1 means i->j."""
    A = [[0] * n for _ in range(n)]
    for (i, j), s in zip(edges, spins):
        if s == 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    return A

def scores(n, A):
    return [sum(A[v]) for v in range(n)]

def count_c3(n, A):
    c = 0
    for i, j, k in combinations(range(n), 3):
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            c += 1
    return c

# ---------------------------------------------------------------------------
# Redei Hamiltonian-path count H(T) (exact, Held-Karp style bitmask DP).
# ---------------------------------------------------------------------------
def H_count(n, A):
    # number of directed Hamiltonian paths (any start, following arc directions)
    # dp[mask][v] = #paths covering set 'mask' ending at v
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            cur = dp[mask][v]
            if not cur:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:  # arc v->w
                    dp[mask | (1 << w)][w] += cur
    return sum(dp[full][v] for v in range(n))

# ===========================================================================
# H1 + H2: c3 as an Ising energy on arc spins.
# ===========================================================================
def test_H1_H2(n):
    print(f"\n=== H1/H2 at n={n}: c3 = C(n,3) - sum C(s_v,2), spin-variance form ===",
          flush=True)
    edges = edges_of(n)
    m = len(edges)
    best_c3 = -1
    best_scoreclass = None
    c3_by_scoreclass = defaultdict(set)
    # Exhaustive through n=6. At n=7 the identity is exact, but the full
    # 2^21 score sweep is too slow for this evidence script.
    exhaustive = (m <= 15)
    if not exhaustive:
        print("  (sampling, n too large for exhaustive)", flush=True)
    # precompute, for each edge bit b=(i,j), which vertex (i or j) gets +1 to its
    # out-score when spin=+1 (i->j) vs spin=-1.  Then scores via fast bit loop.
    iters = range(1 << m) if exhaustive else range(50000)
    checked_identity = False
    for it in iters:
        if exhaustive:
            mask = it
        else:
            mask = random.getrandbits(m)
        s = [0] * n
        for b, (i, j) in enumerate(edges):
            if (mask >> b) & 1:
                s[i] += 1
            else:
                s[j] += 1
        c3 = comb(n, 3) - sum(comb(sv, 2) for sv in s)
        if not checked_identity:
            # spin-variance identity check (once, exact rationals)
            var_form = comb(n, 3) - Fraction(1, 2) * (sum(sv * sv for sv in s) - comb(n, 2))
            assert var_form == c3, (s, c3, var_form)
            checked_identity = True
        sc = tuple(sorted(s))
        c3_by_scoreclass[sc].add(c3)
        if c3 > best_c3:
            best_c3 = c3
            best_scoreclass = sc
    # c3 is a FUNCTION of the score class (H1 core claim): each score class -> single c3
    one_to_one = all(len(v) == 1 for v in c3_by_scoreclass.values())
    print(f"  c3 determined by score multiset alone?  {one_to_one}")
    print(f"  max c3 = {best_c3}, attained by score class {best_scoreclass}")
    # regular score class
    reg = tuple([(n - 1) // 2] * n) if n % 2 == 1 else None
    if reg is not None:
        print(f"  regular score class {reg} -> c3 = "
              f"{next(iter(c3_by_scoreclass.get(reg, {'NA'})))}")
        print(f"  ARGMAX c3 IS the regular class?  {best_scoreclass == reg}")
    return best_c3, best_scoreclass, c3_by_scoreclass

# ===========================================================================
# H2 explicit: c3 = A + B * sum over vertex-sharing edge pairs of (signed) s_e s_f
# Need an exact constant pair (A,B). Derive sign convention so the line-graph
# Ising coupling is uniform. s_v = #out-arcs at v.  Express s_v - sbar in spins.
# For edge (i,j): contributes +1/2 to (signed) "out indicator" at i and -1/2 at j
# under symmetric spin. We test the cleaner statement directly via least squares
# over exact rationals: fit c3 to 1 and  Q = sum_{e~f} s_e s_f  (line graph).
# ===========================================================================
def line_graph_quadratic(n):
    print(f"\n=== H2 explicit at n={n}: fit c3 = A + B*Q,  Q = sum_{{e~f}} s_e s_f ===")
    edges = edges_of(n)
    m = len(edges)
    # adjacency of edges sharing a vertex (line graph of K_n)
    pairs = [(a, b) for a in range(m) for b in range(a + 1, m)
             if set(edges[a]) & set(edges[b])]
    rows = []  # (1, Q, c3)
    cap = 6000
    seen = 0
    for it in range(1 << m) if m <= 15 else range(cap):
        if m <= 15:
            spins = [1 if (it >> b) & 1 else -1 for b in range(m)]
        else:
            spins = [random.choice((1, -1)) for _ in range(m)]
        Q = sum(spins[a] * spins[b] for a, b in pairs)
        A = spins_to_adj(n, edges, spins)
        s = scores(n, A)
        c3 = comb(n, 3) - sum(comb(sv, 2) for sv in s)
        rows.append((Q, c3))
        seen += 1
        if seen >= cap and m > 15:
            break
    # check c3 is exactly affine in Q: c3 = A + B*Q
    # solve from two distinct Q values, then verify all
    qs = {}
    for Q, c3 in rows:
        qs.setdefault(Q, set()).add(c3)
    functional = all(len(v) == 1 for v in qs.values())
    print(f"  c3 a function of Q alone (line-graph 2-spin energy)?  {functional}")
    if functional and len(qs) >= 2:
        items = sorted((Q, next(iter(v))) for Q, v in qs.items())
        (Q1, c1), (Q2, c2) = items[0], items[-1]
        B = Fraction(c2 - c1, Q2 - Q1)
        A0 = c1 - B * Q1
        ok = all(next(iter(v)) == A0 + B * Q for Q, v in qs.items())
        print(f"  c3 = {A0} + ({B})*Q   exact-affine over all sampled spins?  {ok}")
        print(f"  => c3 is a genuine 2-body Ising energy on the LINE GRAPH of K_n,"
              f" coupling B={B} (sign => antiferromagnetic balancing).")
    return functional

# ===========================================================================
# H3: Hopfield descent on arc spins with AFM line-graph coupling.
# Energy E = + sum_v (s_v - sbar)^2  (want to MINIMIZE -> regular).
# We do asynchronous arc-flip descent: flip an arc if it lowers E_score.
# (This is exactly Hopfield async update with the line-graph weight matrix.)
# ===========================================================================
def hopfield_descent(n, trials=400):
    print(f"\n=== H3 at n={n}: Hopfield/AFM arc-flip descent -> regular attractor? ===")
    edges = edges_of(n)
    m = len(edges)
    sbar = (n - 1) / 2.0

    def E_score(A):
        s = scores(n, A)
        return sum((sv - sbar) ** 2 for sv in s)

    converged_regular = 0
    converged_c3opt = 0
    # max c3 = the global min of E_score; for odd n it's the regular class (c3 = C(n,3)-n*C((n-1)/2,2))
    if n % 2 == 1:
        max_c3 = comb(n, 3) - n * comb((n - 1) // 2, 2)
    else:
        # even n: scores split (n/2-1) and (n/2); min sum C(s,2)
        half = n // 2
        max_c3 = comb(n, 3) - (half * comb(half - 1, 2) + half * comb(half, 2))
    final_c3s = defaultdict(int)
    for _ in range(trials):
        spins = [random.choice((1, -1)) for _ in range(m)]
        A = spins_to_adj(n, edges, spins)
        improved = True
        while improved:
            improved = False
            order = list(range(m))
            random.shuffle(order)
            for b in order:
                i, j = edges[b]
                # flipping arc b: recompute E delta exactly
                # current orientation
                if A[i][j]:
                    out_i, out_j = i, j  # i->j
                else:
                    out_i, out_j = j, i
                s = scores(n, A)
                # flip: out-degree of out_i drops 1, of out_j rises 1
                cur = (s[out_i] - sbar) ** 2 + (s[out_j] - sbar) ** 2
                new = (s[out_i] - 1 - sbar) ** 2 + (s[out_j] + 1 - sbar) ** 2
                if new < cur - 1e-12:
                    A[i][j], A[j][i] = A[j][i], A[i][j]
                    improved = True
        s = scores(n, A)
        c3 = comb(n, 3) - sum(comb(sv, 2) for sv in s)
        final_c3s[c3] += 1
        if c3 == max_c3:
            converged_c3opt += 1
        if n % 2 == 1 and tuple(sorted(s)) == tuple([(n - 1) // 2] * n):
            converged_regular += 1
    print(f"  trials={trials}  global-min-E (max c3={max_c3}) reached: "
          f"{converged_c3opt}/{trials} = {converged_c3opt/trials:.2%}")
    if n % 2 == 1:
        print(f"  reached exact regular score class: "
              f"{converged_regular}/{trials} = {converged_regular/trials:.2%}")
    print(f"  final c3 histogram (top): "
          f"{sorted(final_c3s.items(), key=lambda kv:-kv[1])[:5]}")
    return converged_c3opt / trials

# ===========================================================================
# H4: within the max-c3 (regular) score class, do scores determine H?
# If cut != cycle, NO: H varies inside a single score class.
# ===========================================================================
def test_H4(n, sample=None):
    print(f"\n=== H4 at n={n}: does the score class determine H? (cut vs cycle) ===",
          flush=True)
    edges = edges_of(n)
    m = len(edges)
    if n == 7:
        witnesses = []
        for jumps in ((1, 2, 4), (1, 2, 3), (1, 3, 4)):
            A = regular_circulant(n, jumps)
            witnesses.append((jumps, tuple(sorted(scores(n, A))), H_count(n, A)))
        distinct = sorted({Hv for _, _, Hv in witnesses})
        print("  exact regular-circulant witnesses:")
        for jumps, sc, Hv in witnesses:
            print(f"    jumps={jumps}: score={sc}, H={Hv}")
        print(f"  distinct H values inside the regular score class: {distinct}")
        print(f"  score multiset determines H?  {len(distinct) == 1}")
        print("  => same score/Hopfield ground layer, different cycle-space/H layer.")
        return len(distinct) == 1
    H_by_scoreclass = defaultdict(set)
    reg = tuple([(n - 1) // 2] * n) if n % 2 == 1 else None
    H_in_regular = set()
    if sample is None and m <= 15:
        it_src = (i for i in range(1 << m))
        total = 1 << m
    else:
        it_src = (random.getrandbits(m) for _ in range(sample or 50000))
        total = sample or 50000
    cnt = 0
    for mask in it_src:
        spins = [1 if (mask >> b) & 1 else -1 for b in range(m)]
        A = spins_to_adj(n, edges, spins)
        s = scores(n, A)
        sc = tuple(sorted(s))
        Hv = H_count(n, A)
        H_by_scoreclass[sc].add(Hv)
        if reg is not None and sc == reg:
            H_in_regular.add(Hv)
        cnt += 1
    score_determines_H = all(len(v) == 1 for v in H_by_scoreclass.values())
    note = "exhaustive" if (sample is None and m <= 15) else f"sampled {cnt}"
    print(f"  ({note}) score multiset determines H?  {score_determines_H}  "
          f"(False => need cycle space = 'hidden layer')")
    if reg is not None and H_in_regular:
        print(f"  distinct H values seen within regular class {reg}: {sorted(H_in_regular)}")
        print(f"  => max H in regular class = {max(H_in_regular)} "
              f"(THM-126: Paley=189 at n=7)")
    return score_determines_H


def regular_circulant(n, jumps):
    """Regular circulant tournament on odd n with arcs i -> i+j for j in jumps."""
    jumps = set(jumps)
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for d in range(1, n):
            j = (i + d) % n
            if d in jumps:
                A[i][j] = 1
    return A


if __name__ == "__main__":
    # H1/H2 core identity + argmax
    for n in (5, 6, 7):
        test_H1_H2(n)
    for n in (4, 5, 6):
        line_graph_quadratic(n)
    for n in (5, 7):
        hopfield_descent(n, trials=120)
    test_H4(5)              # exhaustive (1024)
    test_H4(7)              # exact regular-circulant witnesses
