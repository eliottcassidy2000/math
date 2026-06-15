#!/usr/bin/env python3
"""
moment_cumulant_law_verify_macmini.py
mac-mini-2026-06-15

FRESH, independent verification of the MOMENT-CUMULANT law for tournament
cycle structure (THM-502 territory), written from scratch.

We verify three distinct, precisely-distinguished objects:

(A) The tr(A^k) = sum-over-cycle-configs decomposition (k=3..8), term-by-term,
    by INDEPENDENTLY counting:
      - tr(A^k)  : exact, via integer matrix power, summed diagonal.
      - c_k      : simple directed k-cycles, by brute cyclic enumeration.
      - the named overlap terms p33, TQ, Q44, TF : by direct subgraph search,
        as the count of (cycle, cycle) pairs that SHARE >= 1 vertex.
    and checking the claimed integer identities hold EXACTLY on every sampled
    tournament (random + all small).

(B) The WITT / necklace transform  W_k = (1/k) sum_{d|k} mu(d) tr(A^{k/d})
    and the claim  W_k = c_k + (overlap configs).  This is the "cumulant" side
    in the NECKLACE poset sense (primitive aperiodic closed walks up to rotation).

(C) The FREE-PROBABILITY moment<->cumulant test.  Treat the normalized traces
    m_k = tr(A^k)/n as the moments of the (real, but possibly with multiplicity)
    spectral distribution of A (A is the {0,+-1}/skew-ish tournament matrix; its
    eigenvalues are generally complex, so we treat m_k formally as a moment
    sequence).  Compute:
      - CLASSICAL cumulants  (all set partitions, exponential-formula Mobius)
      - FREE cumulants       (non-crossing partitions, Speicher Mobius)
    and ask honestly: do EITHER of these equal the cycle-structure quantities
    (c_k, W_k, the necklace transform)?  This is the part that decides whether
    the analogy is LITERAL free probability or only a formal Mobius analogy on a
    DIFFERENT poset (necklace vs NC-partition lattice).

All counts are exact integers / exact Fractions.  No float identities are
trusted; every "==" below is an exact integer / Fraction comparison.
"""

import sys
import random
from fractions import Fraction
from itertools import combinations, permutations
from math import comb


# ----------------------------------------------------------------------------
# Tournament representation and basic invariants
# ----------------------------------------------------------------------------

def random_tournament(n, rng):
    """adj[i][j] = 1 iff arc i->j. Exactly one of (i->j),(j->i) for i<j."""
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return adj


def tournament_from_bits(n, bits):
    adj = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj


def mat_pow_trace(adj, n, k):
    """Exact integer trace of A^k for the 0/1 adjacency matrix A=adj."""
    # A^k via repeated integer matmul; we only need the trace.
    def matmul(X, Y):
        Z = [[0] * n for _ in range(n)]
        for i in range(n):
            Xi = X[i]
            Zi = Z[i]
            for t in range(n):
                a = Xi[t]
                if a:
                    Yt = Y[t]
                    for j in range(n):
                        Zi[j] += a * Yt[j]
        return Z
    P = [row[:] for row in adj]
    for _ in range(k - 1):
        P = matmul(P, adj)
    return sum(P[i][i] for i in range(n))


# ----------------------------------------------------------------------------
# Direct cycle enumeration: simple directed cycles of length L
# Returns the list of vertex-FROZENSETS (one per cycle) so overlaps are easy.
# We canonicalize each directed cycle by its minimal rotation to avoid double
# counting the L rotations of the same cycle.
# ----------------------------------------------------------------------------

def simple_cycles_of_length(adj, n, L):
    """Return list of cyclic-vertex-tuples for each distinct simple directed
    L-cycle. Each cycle counted ONCE (canonical minimal rotation)."""
    cycles = []
    for combo in combinations(range(n), L):
        # over all cyclic orderings of these L vertices, check arcs
        # fix the smallest vertex first to kill rotation; both directions allowed
        base = combo[0]
        rest = combo[1:]
        for perm in permutations(rest):
            seq = (base,) + perm
            ok = True
            for i in range(L):
                u = seq[i]
                v = seq[(i + 1) % L]
                if not adj[u][v]:
                    ok = False
                    break
            if ok:
                cycles.append(seq)
    return cycles


def count_simple_cycles(adj, n, L):
    return len(simple_cycles_of_length(adj, n, L))


def overlap_pairs(cyclesA, cyclesB, same):
    """Count UNORDERED pairs (a in cyclesA, b in cyclesB) with a != b that
    SHARE at least one vertex.  If same=True, A and B are the same list and we
    count unordered distinct pairs; else ordered-set (a from A, b from B)."""
    setsA = [frozenset(c) for c in cyclesA]
    if same:
        cnt = 0
        m = len(setsA)
        for i in range(m):
            for j in range(i + 1, m):
                if setsA[i] & setsA[j]:
                    cnt += 1
        return cnt
    else:
        setsB = [frozenset(c) for c in cyclesB]
        cnt = 0
        for sa in setsA:
            for sb in setsB:
                if sa & sb:
                    cnt += 1
        return cnt


# ----------------------------------------------------------------------------
# (A)+(B): verify the explicit ladder and the Witt transform
# ----------------------------------------------------------------------------

MU = {1: 1, 2: -1, 3: -1, 5: -1, 6: 1, 7: -1}  # Mobius mu(d) for d we need
def mobius(d):
    # full mobius up to our needs
    table = {1: 1, 2: -1, 3: -1, 4: 0, 5: -1, 6: 1, 7: -1, 8: 0}
    return table[d]

def divisors(k):
    return [d for d in range(1, k + 1) if k % d == 0]


def witt_transform(traces, k):
    """W_k = (1/k) sum_{d|k} mu(d) tr A^{k/d}.  traces[m] = tr A^m."""
    s = 0
    for d in divisors(k):
        s += mobius(d) * traces[k // d]
    assert s % k == 0, f"W_{k} not integer: {s}"
    return s // k


def analyze_tournament(adj, n, kmax):
    """Return a dict of all the exact integer invariants we need."""
    traces = {m: mat_pow_trace(adj, n, m) for m in range(1, kmax + 1)}

    # simple cycle counts up to length min(n,kmax)
    cyc = {}
    for L in range(3, min(n, kmax) + 1):
        cyc[L] = simple_cycles_of_length(adj, n, L)
    cnt = {L: len(cyc[L]) for L in cyc}

    res = {"traces": traces, "c": cnt}

    # overlap configurations (only those needed for k<=8)
    # p33 = overlapping triangle pairs (two distinct 3-cycles sharing >=1 vertex)
    if 3 in cyc:
        res["p33"] = overlap_pairs(cyc[3], None, same=True)
    # TQ = overlapping (triangle, 4-cycle) pairs sharing >=1 vertex
    if 3 in cyc and 4 in cyc:
        res["TQ"] = overlap_pairs(cyc[3], cyc[4], same=False)
    # Q44 = overlapping 4-cycle pairs
    if 4 in cyc:
        res["Q44"] = overlap_pairs(cyc[4], None, same=True)
    # TF = overlapping (triangle, 5-cycle) pairs
    if 3 in cyc and 5 in cyc:
        res["TF"] = overlap_pairs(cyc[3], cyc[5], same=False)
    return res


def check_ladder(res):
    """Return list of (label, lhs, rhs, ok) for each ladder identity that is
    fully determined by the invariants present."""
    tr = res["traces"]
    c = res["c"]
    checks = []

    # tr A^3 = 3 c3
    if 3 in c:
        checks.append(("trA3 = 3 c3", tr[3], 3 * c[3], tr[3] == 3 * c[3]))
    # tr A^4 = 4 c4
    if 4 in c:
        checks.append(("trA4 = 4 c4", tr[4], 4 * c[4], tr[4] == 4 * c[4]))
    # tr A^5 = 5 c5
    if 5 in c:
        checks.append(("trA5 = 5 c5", tr[5], 5 * c[5], tr[5] == 5 * c[5]))
    # tr A^6 = 6 c6 + 3 c3 + 6 p33
    if 6 in c and 3 in c and "p33" in res:
        rhs = 6 * c[6] + 3 * c[3] + 6 * res["p33"]
        checks.append(("trA6 = 6c6+3c3+6p33", tr[6], rhs, tr[6] == rhs))
    # tr A^7 = 7 c7 + 7 TQ
    if 7 in c and "TQ" in res:
        rhs = 7 * c[7] + 7 * res["TQ"]
        checks.append(("trA7 = 7c7+7TQ", tr[7], rhs, tr[7] == rhs))
    # tr A^8 = 8 c8 + 4 c4 + 8 Q44 + 8 TF
    if 8 in c and 4 in c and "Q44" in res and "TF" in res:
        rhs = 8 * c[8] + 4 * c[4] + 8 * res["Q44"] + 8 * res["TF"]
        checks.append(("trA8 = 8c8+4c4+8Q44+8TF", tr[8], rhs, tr[8] == rhs))
    return checks


def check_witt(res, kmax):
    """Witt transform identities: W_k = c_k + (overlap configs)."""
    tr = res["traces"]
    c = res["c"]
    checks = []
    # W_3 = c3, W_4 = c4, W_5 = c5
    for L in (3, 4, 5):
        if L in c and L <= kmax:
            w = witt_transform(tr, L)
            checks.append((f"W{L} = c{L}", w, c[L], w == c[L]))
    # W_6 = c6 + p33
    if 6 <= kmax and 6 in c and "p33" in res:
        w = witt_transform(tr, 6)
        rhs = c[6] + res["p33"]
        checks.append(("W6 = c6+p33", w, rhs, w == rhs))
    # W_7 = c7 + TQ
    if 7 <= kmax and 7 in c and "TQ" in res:
        w = witt_transform(tr, 7)
        rhs = c[7] + res["TQ"]
        checks.append(("W7 = c7+TQ", w, rhs, w == rhs))
    # W_8 = c8 + Q44 + TF
    if 8 <= kmax and 8 in c and "Q44" in res and "TF" in res:
        w = witt_transform(tr, 8)
        rhs = c[8] + res["Q44"] + res["TF"]
        checks.append(("W8 = c8+Q44+TF", w, rhs, w == rhs))
    return checks


# ----------------------------------------------------------------------------
# (C) Free-probability moment<->cumulant machinery (independent reimplementation)
# ----------------------------------------------------------------------------

def set_partitions(collection):
    """All set partitions of a list (Bell number many)."""
    collection = list(collection)
    if len(collection) == 1:
        yield [collection]
        return
    first = collection[0]
    for smaller in set_partitions(collection[1:]):
        for i, subset in enumerate(smaller):
            yield smaller[:i] + [[first] + subset] + smaller[i + 1:]
        yield [[first]] + smaller


def is_noncrossing(partition, n):
    """A partition of {0..n-1} (blocks=lists) is non-crossing iff no a<b<c<d
    with a,c in one block and b,d in another."""
    block_of = [0] * n
    for bi, B in enumerate(partition):
        for x in B:
            block_of[x] = bi
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                for d in range(c + 1, n):
                    if block_of[a] == block_of[c] and block_of[b] == block_of[d] \
                       and block_of[a] != block_of[b]:
                        return False
    return True


def noncrossing_partitions(n):
    return [p for p in set_partitions(list(range(n))) if is_noncrossing(p, n)]


def classical_cumulants_from_moments(m, N):
    """m[0]=1, m[1..N] moments. Solve m_n = sum_{all partitions pi of [n]}
    prod_B kappa_{|B|}.  Returns kappa[1..N] (Fractions)."""
    kappa = [None] * (N + 1)
    for n in range(1, N + 1):
        s = Fraction(0)
        for part in set_partitions(list(range(n))):
            if len(part) == 1:
                continue
            prod = Fraction(1)
            for B in part:
                prod *= kappa[len(B)]
            s += prod
        kappa[n] = Fraction(m[n]) - s
    return kappa


def free_cumulants_from_moments(m, N):
    """m[0]=1, m[1..N]. Solve m_n = sum_{NC partitions pi of [n]} prod_B k_{|B|}.
    Returns k[1..N] (Fractions)."""
    k = [None] * (N + 1)
    for n in range(1, N + 1):
        s = Fraction(0)
        for part in noncrossing_partitions(n):
            if len(part) == 1:
                continue
            prod = Fraction(1)
            for B in part:
                prod *= k[len(B)]
            s += prod
        k[n] = Fraction(m[n]) - s
    return k


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def run_ladder_witt(label, tournaments, n, kmax):
    fails = []
    total = 0
    ladder_pass = 0
    witt_pass = 0
    ladder_count = 0
    witt_count = 0
    for adj in tournaments:
        res = analyze_tournament(adj, n, kmax)
        lc = check_ladder(res)
        wc = check_witt(res, kmax)
        for (lab, lhs, rhs, ok) in lc:
            ladder_count += 1
            if ok:
                ladder_pass += 1
            else:
                fails.append((lab, lhs, rhs))
        for (lab, lhs, rhs, ok) in wc:
            witt_count += 1
            if ok:
                witt_pass += 1
            else:
                fails.append((lab, lhs, rhs))
        total += 1
    print(f"[{label}] n={n}, kmax={kmax}: {total} tournaments")
    print(f"   ladder identities: {ladder_pass}/{ladder_count} exact")
    print(f"   witt  identities : {witt_pass}/{witt_count} exact")
    if fails:
        print("   FAILURES (first 5):")
        for f in fails[:5]:
            print("     ", f)
    else:
        print("   ALL EXACT")
    return len(fails) == 0


def run_free_prob_test(adj, n, kmax):
    """For one tournament: compute normalized moments m_k = tr(A^k)/n, then
    classical & free cumulants, and compare to cycle quantities."""
    traces = {m: mat_pow_trace(adj, n, m) for m in range(1, kmax + 1)}
    # moments of the empirical spectral distribution: m_k = (1/n) tr A^k
    m = [Fraction(1)] + [Fraction(traces[k], n) for k in range(1, kmax + 1)]
    cl = classical_cumulants_from_moments(m, kmax)
    fr = free_cumulants_from_moments(m, kmax)
    res = analyze_tournament(adj, n, kmax)
    return traces, m, cl, fr, res


def main():
    rng = random.Random(20260615)
    print("=" * 78)
    print("MOMENT-CUMULANT LAW FOR TOURNAMENT CYCLE STRUCTURE -- fresh verification")
    print("=" * 78)

    # ----- (A)+(B) exhaustive on small n, random on n=7,8 -----
    print("\n--- PART (A)+(B): ladder + Witt identities ---\n")
    ok_all = True

    # exhaustive n=4,5
    for n in (4, 5):
        m_edges = n * (n - 1) // 2
        allt = [tournament_from_bits(n, b) for b in range(1 << m_edges)]
        ok_all &= run_ladder_witt(f"exhaustive n={n}", allt, n, kmax=min(n, 8))

    # exhaustive n=6 (2^15 = 32768)
    n = 6
    m_edges = n * (n - 1) // 2
    allt6 = [tournament_from_bits(n, b) for b in range(1 << m_edges)]
    ok_all &= run_ladder_witt("exhaustive n=6", allt6, n, kmax=6)

    # random n=7 (kmax=7) and n=8 (kmax=8)
    n = 7
    rt7 = [random_tournament(n, rng) for _ in range(300)]
    ok_all &= run_ladder_witt("random n=7 (300)", rt7, n, kmax=7)

    n = 8
    rt8 = [random_tournament(n, rng) for _ in range(200)]
    ok_all &= run_ladder_witt("random n=8 (200)", rt8, n, kmax=8)

    print(f"\n>>> PART (A)+(B) overall: {'ALL EXACT' if ok_all else 'FAILURES PRESENT'}")

    # ----- (C) free-probability moment-cumulant comparison -----
    print("\n" + "=" * 78)
    print("--- PART (C): is the necklace/Witt law a FREE-prob moment-cumulant law? ---")
    print("=" * 78)
    print("""
We treat normalized traces m_k = tr(A^k)/n as a moment sequence of the spectral
distribution of A. We compute CLASSICAL cumulants (all partitions) and FREE
cumulants (non-crossing partitions) and compare to the cycle quantities W_k and
c_k. KEY question: does either equal W_k = (1/k) sum mu(d) tr A^{k/d}?
""")
    n = 7
    sample = [random_tournament(n, rng) for _ in range(6)]
    kmax = 7
    free_eq_witt = 0
    classical_eq_witt = 0
    free_eq_witt_times_n = 0
    samples_done = 0
    for adj in sample:
        traces, m, cl, fr, res = run_free_prob_test(adj, n, kmax)
        W = {k: witt_transform(traces, k) for k in range(1, kmax + 1)}
        samples_done += 1
        print(f"\n  tournament #{samples_done} (n={n}):")
        print(f"    traces trA^k        : {[traces[k] for k in range(1, kmax+1)]}")
        print(f"    W_k (necklace)      : {[W[k] for k in range(1, kmax+1)]}")
        print(f"    moments m_k=trA^k/n : {[str(m[k]) for k in range(1, kmax+1)]}")
        print(f"    classical kappa_k   : {[str(cl[k]) for k in range(1, kmax+1)]}")
        print(f"    free      k_k       : {[str(fr[k]) for k in range(1, kmax+1)]}")
        # comparisons (exact Fraction)
        fe = all(Fraction(W[k]) == fr[k] for k in range(1, kmax + 1))
        ce = all(Fraction(W[k]) == cl[k] for k in range(1, kmax + 1))
        # also test the natural normalization W_k vs n*free_k or free_k vs W_k/n etc.
        fen = all(Fraction(W[k], n) == fr[k] for k in range(1, kmax + 1))
        free_eq_witt += fe
        classical_eq_witt += ce
        free_eq_witt_times_n += fen
        print(f"    free_k  == W_k      ? {fe}")
        print(f"    classical_k == W_k  ? {ce}")
        print(f"    free_k  == W_k/n    ? {fen}")

    print(f"\n  Across {samples_done} tournaments:")
    print(f"    free cumulant == W_k        in {free_eq_witt}/{samples_done}")
    print(f"    classical cumulant == W_k   in {classical_eq_witt}/{samples_done}")
    print(f"    free cumulant == W_k/n      in {free_eq_witt_times_n}/{samples_done}")

    # ----- (C2): the poset-size diagnostic -----
    print("\n" + "-" * 78)
    print("POSET DIAGNOSTIC: necklace count vs non-crossing-partition count")
    print("-" * 78)
    print("""
The Witt/necklace Mobius runs over the DIVISOR poset of k (Mobius mu(d), a
1-dimensional necklace inversion). Free probability's Mobius runs over the
NON-CROSSING PARTITION LATTICE NC(k) (Catalan-many elements). If they were the
same inversion these counts would have to match. They do not:
""")
    print(f"    {'k':>3} {'#divisors(k)':>14} {'#NC(k)=Catalan':>16} {'#all-partitions=Bell':>22}")
    def catalan(k):
        return comb(2 * k, k) // (k + 1)
    def bell(k):
        return sum(1 for _ in set_partitions(list(range(k)))) if k > 0 else 1
    for k in range(1, 8):
        ndiv = len(divisors(k))
        print(f"    {k:>3} {ndiv:>14} {catalan(k):>16} {bell(k):>22}")
    print("""
=> The necklace (Witt) inversion is Mobius on the DIVISOR LATTICE of k, size
   #divisors(k) (e.g. 4 for k=6). The free-probability inversion is Mobius on
   NC(k), size Catalan(k) (e.g. 132 for k=6). DIFFERENT posets, different sizes.
   So the THM-502 law is NOT literally Speicher's free moment-cumulant relation;
   it is the NECKLACE/Witt (cyclic) Mobius inversion -- the SAME inversion that
   defines free cumulants ONLY in the degenerate sense that both are Mobius
   inversions, on different posets.
""")


if __name__ == "__main__":
    main()
