#!/usr/bin/env python3
"""death-star-2026-07-16-S17: THM-726 SHAPE AUDIT — multi-killer covering-min over ALL
small-part shapes (not just interval cores), with free outliers, and PROVED union-tail
thresholds replacing the unproven far-element monotonicity.

Setup. S = P ∪ W, P = S ∩ {1..12} (|P| = k), W = outliers (elements ≥ 13, j = 13-k ≥ 2),
S primitive covering 13-set. THM-726 claims M(S) ≥ 1/13. The S70 finite check enumerated
ONLY interval cores P = {1..k}, k=9,10,11, carrier-multiple outliers ≤ 220.

This audit:
  (E) flagship anchors, exact ({1..14}\\{a}; {1..12}\\{a} ∪ {w,182m} free-outlier family;
      self-tests on known values 7/89, 2/23, 14/183).
  (A) k=11 (j=2): ALL 12 shapes, exhaustive below PROVED per-shape tail thresholds.
  (B) k=10 (j=3), k=9 (j=4): witness pre-filter + nested exact decision below thresholds.

UNION-TAIL LEMMA (the monotonicity replacement; 1/13 version of opus-S32 simultaneous-peel):
  Let G ⊆ [0,1) be the closed 1/13-good set of a finite core C (measure m, r components).
  For extra speeds w_1..w_j' all ≥ W:  meas(good(C ∪ {w_i})) ≥ m(13-2j')/13 - (2jr'/13)(1/W)
  [each w chops ≤ (2/13)|I| + 2/(13w) from each component I]. Positive iff
  W > 2 j' r / (m (13 - 2 j')), for j' ≤ 6. Positive measure of a closed set ⟹ some t with
  ALL clearances ≥ 1/13 ⟹ M ≥ 1/13, unconditionally. NO value-monotonicity assumed.
"""
from fractions import Fraction as F
from math import gcd
import sys, time

ONE13 = F(1, 13)

# ---------------- exact interval machinery (closed intervals in [0,1]) ----------------

def good_arcs_of_speed(w, lo=F(0), hi=F(1)):
    """Closed good set {t in [lo,hi] : ||wt|| >= 1/13} as list of (a,b) Fractions, a<=b."""
    out = []
    # bad arcs: (c/w - 1/(13w), c/w + 1/(13w)) open, for integer c
    # good: [c/w + 1/(13w), (c+1)/w - 1/(13w)] for each period
    c0 = int(lo * w) - 1
    c1 = int(hi * w) + 1
    for c in range(c0, c1 + 1):
        a = F(13 * c + 1, 13 * w)
        b = F(13 * c + 12, 13 * w)
        if b < lo or a > hi:
            continue
        out.append((max(a, lo), min(b, hi)))
    return out

def intersect_ivsets(A, B):
    """Intersect two sorted lists of closed intervals."""
    out = []
    i = j = 0
    while i < len(A) and j < len(B):
        a1, a2 = A[i]; b1, b2 = B[j]
        lo, hi = max(a1, b1), min(a2, b2)
        if lo <= hi:
            out.append((lo, hi))
        if a2 < b2:
            i += 1
        else:
            j += 1
    return out

def good_set(speeds):
    """Exact closed 1/13-good set of a list of speeds, as sorted interval list on [0,1]."""
    G = [(F(0), F(1))]
    for w in sorted(speeds, reverse=True):
        G = intersect_ivsets(G, good_arcs_of_speed(w))
        if not G:
            return []
    return G

def ivset_stats(G):
    m = sum(b - a for a, b in G)
    r = len(G)
    longest = max((b - a for a, b in G), default=F(0))
    return m, r, longest

def clears_via_intervals(G, w):
    """Does some t in G (closed interval list) have ||wt|| >= 1/13? Exact.
    Interval [a,b]: image [wa, wb] fails iff contained in ONE open bad arc
    (c - 1/13, c + 1/13)."""
    for a, b in G:
        wa, wb = w * a, w * b
        c = (13 * wa + 1).__floor__()  # candidate: need wa > c - 1/13 i.e. 13wa+1 > 13c
        # fails iff exists integer c with  c - 1/13 < wa <= wb < c + 1/13
        # i.e. 13*wa > 13c - 1  and  13*wb < 13c + 1  -> c must be round(wa)
        lo13, hi13 = 13 * wa, 13 * wb
        ok = True
        # find integer c with 13c in (lo13 - 1 ... ): c in ((lo13-1)/13, (hi13+1)/13)
        cmin = ((lo13 - 1) / 13).__floor__()
        cmax = ((hi13 + 1) / 13).__ceil__()
        for c in range(cmin, cmax + 1):
            if 13 * c - 1 < lo13 and hi13 < 13 * c + 1:
                ok = False
                break
        if ok:
            return True
    return False

# ---------------- exact M via breakpoint grid (for anchors/self-tests) ----------------

def exact_M(S, qmax=None):
    """max_t min_v ||vt|| exactly: max over Farey breakpoints a/q, q <= 2*max(S)."""
    if qmax is None:
        qmax = 2 * max(S) + 1
    best = F(0)
    for q in range(2, qmax + 1):
        for a in range(1, q // 2 + 1):
            if gcd(a, q) != 1:
                continue
            mind = None
            for v in S:
                r = (a * v) % q
                d = min(r, q - r)
                if mind is None or d < mind:
                    mind = d
                    if F(d, q) <= best:
                        break
            val = F(mind, q)
            if val > best:
                best = val
    return best

# ---------------- witness pre-filter ----------------

def witness_list(P, Q0=131):
    """All (a,q), q<=Q0, gcd(a,q)=1, a<=q/2, good for core P: min 13*dist(a*v%q) >= q."""
    wits = []
    for q in range(14, Q0 + 1):
        for a in range(1, q // 2 + 1):
            if gcd(a, q) != 1:
                continue
            ok = True
            for v in P:
                r = (a * v) % q
                d = min(r, q - r)
                if 13 * d < q:
                    ok = False
                    break
            if ok:
                wits.append((a, q))
    return wits

def wit_clears(wits, ws):
    """Does some witness clear all speeds in ws? (integer exact)"""
    for (a, q) in wits:
        ok = True
        for w in ws:
            r = (a * w) % q
            d = min(r, q - r)
            if 13 * d < q:
                ok = False
                break
        if ok:
            return True
    return False

# ---------------- covering helpers ----------------

def covered_moduli(P):
    return {q for q in range(2, 15) if any(v % q == 0 for v in P)}

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def is_primitive(S):
    g = 0
    for v in S:
        g = gcd(g, v)
    return g == 1

# ---------------- Part E: flagship anchors ----------------

def part_E():
    print("=" * 78)
    print("PART E: flagship anchors (exact)")
    print("=" * 78)
    t0 = time.time()
    # self-tests
    st = [
        ([1,2,3,4,5,6,7,8,9,10,11,12,182], F(14,183), "deep well"),
        ([1,2,3,4,5,6,7,8,9,10,11,13,84], F(7,89), "THM-726 min (interval k=11)"),
        ([1,2,3,4,5,7,8,9,10,11,12,13,182], F(2,23), "HYP-6660 near-AP-with-far"),
    ]
    for S, expect, name in st:
        got = exact_M(S)
        tag = "OK" if got == expect else f"MISMATCH expect {expect}"
        print(f"  self-test {name}: M = {got} = {float(got):.6f}  [{tag}]")
    # {1..14}\{a}, a = 1..7 (covering 13-sets, outliers {13,14})
    print("\n  {1..14}\\{a} family (multi-killer by the >=13 definition, j=2):")
    for a in range(1, 8):
        S = [v for v in range(1, 15) if v != a]
        cov = is_covering(S)
        M = exact_M(S)
        rel = "  >= 1/13 OK" if M >= ONE13 else "  *** BELOW 1/13 ***"
        print(f"    a={a}: covering={cov}  M = {M} = {float(M):.6f}{rel}")
    print(f"  [E1 time {time.time()-t0:.1f}s]")
    # free-outlier family {1..12}\{6} u {w, 182}
    t1 = time.time()
    print("\n  free-outlier family S = {1..12}\\{6} + {w, 182}, w = 13..400:")
    P = [v for v in range(1, 13) if v != 6]
    below = []
    minM = None
    wits = witness_list(P)
    G_P182 = good_set(P + [182])
    for w in range(13, 401):
        if w == 182:
            continue
        S = P + [w, 182]
        if not is_covering(S) or not is_primitive(S):
            continue
        if wit_clears(wits, [w, 182]):
            continue
        # exact decision
        if clears_via_intervals(G_P182, w):
            continue
        M = exact_M(S)
        below.append((M, w))
        if minM is None or M < minM[0]:
            minM = (M, w)
    if below:
        print(f"    {len(below)} values of w with M < 1/13:")
        for M, w in sorted(below)[:10]:
            print(f"      w={w}: M = {M} = {float(M):.6f}" + (" *** BELOW 14/183" if M < F(14,183) else ""))
    else:
        print("    ALL w clear 1/13 (witness or exact interval).")
    print(f"  [E2 time {time.time()-t1:.1f}s]")

# ---------------- optimized arc generation near an interval set ----------------

def good_arcs_near(w, G):
    """Good arcs of speed w restricted to the span of interval set G (exact)."""
    out = []
    for lo, hi in G:
        c0 = int(lo * w) - 1
        c1 = int(hi * w) + 1
        for c in range(c0, c1 + 1):
            a = F(13 * c + 1, 13 * w)
            b = F(13 * c + 12, 13 * w)
            if b < lo or a > hi:
                continue
            out.append((max(a, lo), min(b, hi)))
    return out

def refine(G, w):
    """G ∩ good(w), exact, using locality."""
    return intersect_ivsets(G, good_arcs_near(w, G))

def tail_threshold(m, r, j):
    """Smallest integer W with W > 2 j r / (m (13-2j)); all j extra speeds >= W clear."""
    assert 1 <= j <= 6 and m > 0
    bound = F(2 * j * r, 1) / (m * (13 - 2 * j))
    return int(bound) + 1

def divisors_in_range(w, lo=2, hi=14):
    return {q for q in range(lo, hi + 1) if w % q == 0}

def lcm(xs):
    out = 1
    for x in xs:
        out = out * x // gcd(out, x)
    return out

# ---------------- Part A: k = 11 (j = 2), ALL shapes, exhaustive + proved tails ----------------

def part_A(report_every=0):
    print("=" * 78)
    print("PART A: k=11 (j=2) — all 12 shapes, exhaustive below PROVED union-tail thresholds")
    print("=" * 78)
    grand_min = None
    failures = []
    for a in range(1, 13):
        t0 = time.time()
        P = [v for v in range(1, 13) if v != a]
        R = set(range(2, 15)) - covered_moduli(P)
        G_P = good_set(P)
        m, r, longest = ivset_stats(G_P)
        W_all = tail_threshold(m, r, 2)
        wits = witness_list(P)
        n_pairs = n_exact = 0
        shape_min = None
        # enumerate w1 = smaller outlier in [13, W_all); w2 in (w1, W1(P,w1)) covering-compat
        for w1 in range(13, W_all):
            D1 = divisors_in_range(w1) & R
            R2 = R - D1
            G1 = refine(G_P, w1)
            if not G1:
                failures.append((P, w1, "EMPTY G1 — impossible by prime-13 tightness"))
                continue
            m1, r1, _ = ivset_stats(G1)
            W1 = tail_threshold(m1, r1, 1)
            if R2:
                L = lcm(R2)
                cands = range(((w1 + 1) + L - 1) // L * L, W1, L)
            else:
                cands = range(w1 + 1, W1)
            for w2 in cands:
                if w2 == w1:
                    continue
                if R2 and any(w2 % q for q in R2):
                    continue
                n_pairs += 1
                if wit_clears(wits, [w1, w2]):
                    continue
                n_exact += 1
                if clears_via_intervals(G1, w2):
                    continue
                M = exact_M(P + [w1, w2])
                failures.append((P + [w1, w2], M, "M < 1/13" if M < ONE13 else "interval-miss but M>=1/13?"))
                if shape_min is None or M < shape_min[0]:
                    shape_min = (M, w1, w2)
        status = "ALL CLEAR" if not [f for f in failures if f[0][:11] == P[:11]] else "FAILURES"
        print(f"  shape P = {{1..12}}\\{{{a}}}: R={sorted(R)} m={m}≈{float(m):.5f} r={r} "
              f"W_all={W_all} pairs={n_pairs} exact={n_exact}  [{time.time()-t0:.1f}s]")
        sys.stdout.flush()
    print()
    if failures:
        print(f"  *** {len(failures)} FAILURES:")
        for S, M, why in failures[:20]:
            print(f"      {S}: {M} ({why})")
    else:
        print("  k=11 (j=2): EVERY primitive covering multi-killer config with either outlier")
        print("  below its proved tail threshold CLEARS 1/13; all larger outliers clear by the")
        print("  union-tail lemma. SHAPE-COMPLETE at k=11.")
    return failures

# ---------------- Parts B/C: recursive audit for j = 3, 4 (k = 10, 9) ----------------

def audit_rec(P, wits, G, R_rem, j_rem, w_min, prefix, stats):
    """Complete decision on the stratum: outliers w with w > w_min, exactly j_rem more,
    covering-remainder R_rem. G = exact good set of P + prefix. Everything with smallest
    remaining outlier >= tail_threshold(m, r, j_rem) is cleared by the union-tail lemma;
    below, enumerate and recurse. Returns list of failures."""
    fails = []
    m, r, _ = ivset_stats(G)
    if m <= 0:
        fails.append((P + prefix, F(0), "EMPTY good set at internal node"))
        return fails
    W = tail_threshold(m, r, j_rem)
    stats['maxW'] = max(stats['maxW'], W)
    if j_rem == 1:
        if R_rem:
            L = lcm(R_rem)
            start = max(w_min + 1, 13)
            start = (start + L - 1) // L * L
            cands = range(start, W, L)
        else:
            cands = range(max(w_min + 1, 13), W)
        for w in cands:
            stats['leaves'] += 1
            if wit_clears(wits, prefix + [w]):
                continue
            stats['exact'] += 1
            if clears_via_intervals(G, w):
                continue
            M = exact_M(P + prefix + [w])
            fails.append((P + prefix + [w], M,
                          "M < 1/13" if M < ONE13 else "interval-miss but M >= 1/13"))
        return fails
    for w in range(max(w_min + 1, 13), W):
        stats['nodes'] += 1
        G1 = refine(G, w)
        fails += audit_rec(P, wits, G1, R_rem - divisors_in_range(w), j_rem - 1,
                           w, prefix + [w], stats)
    return fails

def part_BC(kk, time_budget=None):
    from itertools import combinations
    j = 13 - kk
    print("=" * 78)
    print(f"PART {'B' if kk==10 else 'C'}: k={kk} (j={j}) — all shapes, recursive exhaust + proved tails")
    print("=" * 78)
    all_fails = []
    t_start = time.time()
    shapes = list(combinations(range(1, 13), 12 - kk))
    # order by good-set measure descending (easy shapes first, honest progress)
    order = []
    for miss in shapes:
        P = [v for v in range(1, 13) if v not in miss]
        m, r, _ = ivset_stats(good_set(P))
        order.append((-m, miss, m, r))
    order.sort()
    done = 0
    for negm, miss, m, r in order:
        if time_budget and time.time() - t_start > time_budget:
            print(f"  [TIME BUDGET reached: {done}/{len(shapes)} shapes complete; "
                  f"remaining shapes NOT certified this session]")
            for _, miss2, m2, r2 in order[done:done+8]:
                print(f"    remaining: miss={miss2} m={float(m2):.5f}")
            break
        t0 = time.time()
        P = [v for v in range(1, 13) if v not in miss]
        R = set(range(2, 15)) - covered_moduli(P)
        G_P = good_set(P)
        wits = witness_list(P)
        st = {'nodes': 0, 'leaves': 0, 'exact': 0, 'maxW': 0}
        fails = audit_rec(P, wits, G_P, R, j, 12, [], st)
        all_fails += fails
        done += 1
        print(f"  miss={miss}: R={sorted(R)} m≈{float(m):.5f} r={r} "
              f"nodes={st['nodes']} leaves={st['leaves']} exact={st['exact']} maxW={st['maxW']} "
              f"fails={len(fails)}  [{time.time()-t0:.1f}s]")
        sys.stdout.flush()
    print()
    if all_fails:
        print(f"  *** {len(all_fails)} FAILURES at k={kk}:")
        for S, M, why in all_fails[:20]:
            print(f"      {S}: {M} ({why})")
    else:
        print(f"  k={kk} (j={j}): ALL ENUMERATED SHAPES CLEAR — shape-complete on the strata run.")
    return all_fails

if __name__ == "__main__":
    part_E()
    print()
    part_A()
    print()
    part_BC(10)
    print()
    part_BC(9, time_budget=7200)
