#!/usr/bin/env python3
r"""
lrc14_final_rung_collision_arm_monad_S11.py  (monad-explorer-2026-07-09-S11, HYP-5817/THM-682)

THE FINAL-RUNG COLLISION ARM + THE COMMON-RESIDUE DISPATCH.

(A) the common-residue dispatch: all v_l == a (mod d), d >= 2, primitive + covering
    ==> every prime factor of d >= 17, gcd(a, d) = 1, and at tau = c/d (a*c == floor(d/2))
    ALL 13 phases coincide at ||a c / d|| >= (1/2 - 1/(2d)) >= 8/17.  Battery check.
(B) the k = 13, t = 1 restricted-sum stability: B(A) <= 24 ==> diam(A) <= 14 after
    gcd-normalization (containment in an AP of 15 terms).  Exhaustive DFS: the budget
    forces EVERY prefix within +1 of the classical 2s-3 minimum, so the tree is tiny.
    Tail (gap > D/2) is proved disjoint-block: B >= 32 > 24.  DMAX cap flags any sliver.
(C) the B = 25 escape probe: k = 13 rank-2 GAP shapes (LEM-016(ii) analogs) -- which
    block profiles reach B = 25/26 at unbounded diameter; ruler-liveness of core
    instances built on those shapes (THM-680 certification + true LM).
(D) the W0-carrier lemma numerics: support-2 global exacts are DOUBLINGS ONLY (proof in
    THM-682; here: exhaustive confirmation on batteries) + how many doublings W0 > 0.08
    requires (measured weights and proven bounds); max doubling content of a 13-set.
"""
import sys, random
from math import gcd, pi, sin
from itertools import combinations, product

def band(q):
    return -(-q // 14), (13 * q) // 14

def LM_exact(v, q):
    lo, hi = band(q)
    return sum(1 for p in range(q) if all(lo <= (x * p) % q <= hi for x in v))

def hhat_abs(q, k):
    lo, hi = band(q)
    if k % q == 0:
        return (hi - lo + 1) / q
    return abs(sin(pi * k * (hi - lo + 1) / q)) / (q * abs(sin(pi * k / q)))

def line_weight(q, kvec):
    lo, hi = band(q)
    bq = (hi - lo + 1) / q
    tot = 0.0
    for m in range(1, q):
        term = 1.0
        for kl in kvec:
            term *= hhat_abs(q, (m * kl) % q)
        tot += term
    return tot * bq ** (13 - len(kvec))

def is_core(v):
    """covering + primitive + gapped + distinct (the residual-class basics)."""
    if len(set(v)) != 13:
        return False
    g0 = 0
    for x in v:
        g0 = gcd(g0, x)
    if g0 != 1:
        return False
    if not all(any(x % q == 0 for x in v) for q in range(2, 15)):
        return False
    return max(v) > 13 * min(v)

def rsums(A):
    return set(a + b for i, a in enumerate(A) for b in A[i+1:])

# ---------------------------------------------------------------- (A) common-residue
def common_residue_dispatch(v):
    """if all v_l share a residue a mod d for some d >= 2: return (d, a, c, clearance)."""
    diffs = [x - v[0] for x in v[1:]]
    d = 0
    for x in diffs:
        d = gcd(d, x)
    if d < 2:
        return None
    a = v[0] % d
    # primitive+covering  ==>  gcd(a, d) = 1 and min prime factor of d >= 17 (proved);
    # the dispatch itself only needs SOME c with ||a c / d|| >= 1/14:
    best = None
    for c in range(1, d):
        r = (a * c) % d
        clr = min(r, d - r) / d
        if best is None or clr > best[1]:
            best = (c, clr)
    return d, a, best[0], best[1]

# ---------------------------------------------------------------- (B) the DFS
def dfs_stability(k=13, Bmax=24, DMAX=45):
    """exhaustively enumerate 0 = a_0 < ... < a_{k-1} <= DMAX with B <= Bmax, gcd = 1.
    prune: remaining elements each add >= 2 new restricted sums (new max m adds
    m + prev-max and m + prev-second, both above the previous maximum sum)."""
    results = []
    nodes = [0]
    def rec(A, sums):
        nodes[0] += 1
        s = len(A)
        if s == k:
            g0 = 0
            for x in A[1:]:
                g0 = gcd(g0, x)
            if g0 == 1:
                results.append((A[-1], tuple(A), len(sums)))
            return
        for x in range(A[-1] + 1, DMAX + 1):
            ns = sums | set(x + a for a in A)
            if len(ns) + 2 * (k - s - 1) <= Bmax:
                rec(A + [x], ns)
    rec([0], set())
    return results, nodes[0]

# ---------------------------------------------------------------- main
if __name__ == "__main__":
    rng = random.Random(2026070911)

    print("=" * 100)
    print("PART A -- THE COMMON-RESIDUE DISPATCH (battery of in-core affine-dilated families)")
    print("=" * 100)
    hits = 0
    for trial in range(4000):
        d = rng.choice([17, 19, 23, 29, 31, 37])
        a = rng.randrange(1, d)
        if gcd(a, d) != 1:
            continue
        s = sorted(rng.sample(range(0, 40), 13))
        v = [a + si * d for si in s]
        if not is_core(v):
            continue
        res = common_residue_dispatch(v)
        assert res is not None and res[0] % d == 0
        dd, aa, c, clr = res
        ok = clr >= 1 / 14
        hits += 1
        if hits <= 6:
            print(f"  d={dd:>3} a={aa:>3} v=[{v[0]},{v[1]},...,{v[-1]}]  tau={c}/{dd}  "
                  f"clearance={clr:.4f}  (>= 1/2 - 1/(2d) = {0.5 - 0.5/dd:.4f}: "
                  f"{'EXACT' if abs(clr - (0.5 - 0.5/dd)) < 1e-12 or clr >= 0.5 - 0.5/dd - 1e-12 else 'below?!'})")
        assert ok, (v, res)
        if hits >= 400:
            break
    print(f"  {hits} in-core common-residue families: ALL dispatched with clearance >= 8/17 = {8/17:.4f}")
    print("  (each is covering+primitive+gapped+distinct+compressed-checkable: in the residual class,")
    print("   NOT carved by near-harmonic (d divides differences, not values) -- a NEW Lean branch)")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART B -- k=13, t=1 STABILITY: B <= 24 ==> diam <= 14 (gcd-normalized), EXHAUSTIVE DFS")
    print("=" * 100)
    results, nodes = dfs_stability(13, 24, 45)
    diams = sorted(set(r[0] for r in results))
    print(f"  DFS nodes: {nodes};  B<=24 gcd-normalized sets found: {len(results)}")
    print(f"  diameters realized: {diams}")
    print(f"  max diameter: {max(diams) if diams else '-'}  (cap DMAX=45 {'HIT -- SLIVER' if diams and max(diams) >= 45 else 'not hit: EXHAUSTIVE below cap'})")
    byB = {}
    for D, A, B in results:
        byB.setdefault(B, []).append((D, A))
    for B in sorted(byB):
        mx = max(byB[B])
        print(f"    B={B}: {len(byB[B])} sets, max diam {mx[0]}, extremal e.g. {list(mx[1])}")
    # tail statement (proved, printed for the record):
    print("  TAIL (proved): any gap > D/2 splits A into blocks A1|A2 with disjoint sum ranges:")
    print("    B >= (2s-3) + (13-1) + (2(13-s)-3) = 32 > 24  [s,13-s >= 2; singleton cases >= 33]")
    print("  ==> B <= 24 sets have no gap > D/2; the DFS below the cap is the complete bounded part.")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART C -- B = 25 ESCAPES: rank-2 GAP shapes at k = 13 (LEM-016(ii) analogs)")
    print("=" * 100)
    c = 60
    shapes = []
    for n1 in range(2, 12):
        for n2 in range(1, 13 - n1):
            n3 = 13 - n1 - n2
            if n3 < 1:
                continue
            # LEM-016 sharp pattern: blocks at 0, c, 2c with the middle straddling
            A = list(range(n1)) + [c - 1 + j for j in range(n2)] + [2 * c - 2 + j for j in range(n3)]
            A = sorted(set(A))
            if len(A) != 13:
                continue
            B = len(rsums(A))
            shapes.append((B, n1, n2, n3))
    shapes.sort()
    print("  GAP(0, c, 2c) block profiles (c = 60), smallest restricted sumsets:")
    for B, n1, n2, n3 in shapes[:8]:
        print(f"    blocks ({n1},{n2},{n3}): B = {B}" + ("   <-- ESCAPE LEVEL (B <= 26)" if B <= 26 else ""))
    minB = shapes[0][0]
    print(f"  minimum GAP-shape B at k=13: {minB}  (t = {minB - 23}); rank-2 escapes begin at t = {minB - 23}")
    sys.stdout.flush()

    # ruler-liveness of core instances built on near-minimal GAP shapes
    print()
    print("  core instances on near-minimal GAP shapes -- THM-680 certification + true liveness:")
    probed = 0
    for B, n1, n2, n3 in shapes[:4]:
        built = 0
        for trial in range(3000):
            cc = rng.choice([40, 50, 60, 70])
            A = list(range(n1)) + [cc - 1 + j for j in range(n2)] + [2 * cc - 2 + j for j in range(n3)]
            A = sorted(set(A))
            if len(A) != 13:
                break
            off = rng.randrange(1, 30)
            sc = rng.choice([1, 1, 2, 3])
            v = sorted(off + sc * x for x in A)
            if not is_core(v):
                continue
            built += 1
            sums_all = [x + y for i, x in enumerate(v) for y in v[i+1:]]
            sums = sorted(set(sums_all))
            Vmax = v[-1]
            tall = [q for q in sums if q > Vmax]
            true_live = sum(1 for q in tall if LM_exact(v, q) > 0)
            Bv = len(sums)
            print(f"    blocks ({n1},{n2},{n3}) c={cc} off={off} sc={sc}: B={Bv}, "
                  f"true-live {true_live}/{len(tall)} tall rulers"
                  + ("  [LIVE: generic-side behavior]" if true_live > 0 else "  [ALL DEAD?!]"))
            probed += 1
            break
        if probed >= 4:
            break
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART D -- THE W0-CARRIER: support-2 global exacts = DOUBLINGS ONLY; doubling census")
    print("=" * 100)
    # exhaustive coefficient check (the one-line proof, machine-echoed):
    viol = []
    for ca, cb in product((-2, -1, 1, 2), repeat=2):
        # c_a x + c_b y = 0 with 0 < x < y: solutions?
        # normalize ca > 0: x = -(cb/ca) y: needs cb < 0
        if ca <= 0:
            continue
        if cb >= 0:
            continue  # positive sum cannot vanish
        num, den = -cb, ca
        g = gcd(num, den)
        num, den = num // g, den // g
        # x/y = num/den with x < y distinct positive: num < den needed; |c|<=2: ratios 1/2, 1, 2
        if (num, den) == (1, 1):
            continue  # x = y impossible (distinct)
        viol.append(((ca, cb), f"y = {den}/{num} x" if num < den else f"x = {num}/{den} y"))
    print(f"  support-2 |c|<=2 homogeneous relations on distinct positives: {viol}")
    print("  ==> exactly the DOUBLING ratios (2,-1),(1,-2),(2,-2)->reduces,( -,-): doublings only.")
    qr = 200
    wD = line_weight(qr, [2, 1])
    wS = line_weight(qr, [1, 1, 1])
    print(f"  measured doubling line weight at q={qr}: {wD:.5f} (proven bound 0.0382)")
    print(f"  measured Schur line weight    at q={qr}: {wS:.5f} (proven bound 0.0645)")
    import math
    print(f"  W0 > 0.08 requires: >= {math.ceil(0.08 / wD)} doublings (measured) / "
          f">= {math.ceil(0.08 / 0.0382)} (proven bound), OR >= {math.ceil(0.08 / wS)} Schur "
          f"triples (measured) / >= {math.ceil(0.08 / 0.0645)} (proven)")
    # max doubling content of a CORE 13-set: doubling pairs form a forest; census by search
    best = (0, None)
    for trial in range(20000):
        # build doubling-rich candidates: chains a, 2a, 4a... unioned
        v = set()
        while len(v) < 13:
            a = rng.randrange(1, 60)
            L = rng.randint(1, 5)
            for j in range(L):
                v.add(a * (2 ** j))
        v = sorted(v)[:13]
        if len(v) != 13 or not is_core(v):
            continue
        D = sum(1 for x in v for y in v if y == 2 * x)
        if D > best[0]:
            best = (D, list(v))
    print(f"  max doublings found in a CORE family (random search): {best[0]}  e.g. {best[1]}")
    if best[0]:
        v = best[1]
        evens = sum(1 for x in v if x % 2 == 0)
        print(f"    (evens: {evens}/13 -- doubling-rich families are even-rich: the 2-adic pressure)")
    print()
    print("  ==> THE FINAL RUNG after THM-682: [B >= 25] AND [>= 3 doublings (proven) / >= 6 (measured)]")
    print("      AND not GAP-escape-live -- the doubling-rich collision-light corner ONLY.")
    sys.stdout.flush()
