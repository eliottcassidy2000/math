"""
lrc14_largest_arc_witness_denom_kps.py  (kind-pasteur-2026-06-27-S31ag)

Tests HYP-3088: for covering 13-tuples, is there a uniform lower bound on the
LARGEST lonely arc length (=> bounded-denominator witness => Conjecture 7.1(13)
of arXiv:2604.23906 => LRC(14))?

For each test tuple S (13 speeds) we compute exactly:
  - M(S) = max_t min_i ||s_i t||           (the lonely depth)
  - witness denominator D* = reduced denom of the argmax t*
  - largest lonely arc length L_arc at threshold 1/14, and D_arc = ceil(1/L_arc)
  - #lonely arcs (arc count)

Key question: does L_arc stay bounded below as the apex speed V grows (NO, per
the fine-hole argument) — confirming the largest-arc floor needs bounded speeds
(the Node-3 peel) as a prerequisite.  And: is the minimal witness denominator D*
unbounded (HYP-2866)?
"""
from fractions import Fraction as Fr
from math import gcd, ceil

def norm(x: Fr) -> Fr:
    # distance to nearest integer, exact
    f = x - (x.numerator // x.denominator)  # fractional part in [0,1)
    return min(f, 1 - f)

def M_and_witness(S):
    """Exact M(S)=max_t min_i ||s_i t|| via pair candidates t=n/(s_a +/- s_b)."""
    S = sorted(set(S))
    denoms = set()
    for i in range(len(S)):
        for j in range(len(S)):
            if i == j:
                continue
            a, b = S[i], S[j]
            denoms.add(a + b)
            if a != b:
                denoms.add(abs(a - b))
        denoms.add(2 * S[i])  # peak of single tent
    best = Fr(-1)
    best_t = None
    for D in denoms:
        if D == 0:
            continue
        for n in range(1, D):
            t = Fr(n, D)
            # min over speeds
            m = min(norm(Fr(s) * t) for s in S)
            if m > best:
                best = m
                best_t = t
    # reduce witness denom
    if best_t is None:
        return Fr(0), 1
    return best, best_t.denominator

def lonely_arcs(S, thr=Fr(1, 14)):
    """Return (largest_arc_length, arc_count, total_measure) of {t: ||s_i t||>=thr all i}.
    Intersection of unions of intervals. thr=1/14."""
    # Each speed s: lonely intervals [ (j+thr)/s , (j+1-thr)/s ] for j=0..s-1, within [0,1).
    # Intersect over all s. Use a sweep over breakpoints.
    # Collect endpoints (as Fractions). Between consecutive breakpoints, lonely indicator constant.
    pts = set([Fr(0), Fr(1)])
    for s in S:
        for j in range(s):
            lo = Fr(j) / s + thr / s
            hi = Fr(j + 1) / s - thr / s
            if lo < hi:
                if 0 <= lo <= 1:
                    pts.add(lo)
                if 0 <= hi <= 1:
                    pts.add(hi)
    pts = sorted(pts)
    def lonely(t):
        return all(norm(Fr(s) * t) >= thr for s in S)
    arcs = []
    cur_start = None
    for a, b in zip(pts, pts[1:]):
        mid = (a + b) / 2
        if lonely(mid):
            if cur_start is None:
                cur_start = a
            cur_end = b
        else:
            if cur_start is not None:
                arcs.append((cur_start, cur_end))
                cur_start = None
    if cur_start is not None:
        arcs.append((cur_start, cur_end))
    if not arcs:
        return Fr(0), 0, Fr(0)
    lengths = [e - s for s, e in arcs]
    return max(lengths), len(arcs), sum(lengths)

def report(name, S):
    S = list(S)
    g = 0
    for s in S:
        g = gcd(g, s)
    M, Dw = M_and_witness(S)
    try:
        Larc, narc, meas = lonely_arcs(S)
    except Exception as e:
        Larc, narc, meas = Fr(0), -1, Fr(0)
    Darc = ceil(1 / Larc) if Larc > 0 else float('inf')
    tight = (M == Fr(1, 14))
    print(f"{name:28s} gcd={g} max={max(S):>6} "
          f"M={float(M):.5f}{'(=1/14 TIGHT)' if tight else ''} "
          f"Dwit={Dw:>5} | Larc={float(Larc):.5f} Darc={Darc} narcs={narc} meas={float(meas):.4f}")
    return dict(name=name, M=M, Dw=Dw, Larc=Larc, Darc=Darc, narc=narc, meas=meas, tight=tight)

if __name__ == "__main__":
    print("=== Dilated APs (tight family, covering) ===")
    for d in [1, 2, 3]:
        report(f"{d}*{{1..13}}", [d * i for i in range(1, 14)])

    print("\n=== Goddyn-Wong {1..11,13,24} and dilations ===")
    report("GW {1..11,13,24}", list(range(1, 12)) + [13, 24])
    report("2*GW", [2 * x for x in list(range(1, 12)) + [13, 24]])

    print("\n=== HYP-2866 family {1..11,13,84m} (unbounded witness?) ===")
    for m in [1, 2, 3, 5, 7, 11]:
        S = list(range(1, 12)) + [13, 84 * m]
        report(f"{{1..11,13,{84*m}}}", S)

    print("\n=== {1..12, 14j} aliasing family (covering, near 1/13) ===")
    for j in [1, 2, 13, 14, 26]:
        report(f"{{1..12,{14*j}}}", list(range(1, 13)) + [14 * j])

    print("\n=== Residual: covering, <=6 mult of 7, bounded core ===")
    import itertools, random
    random.seed(7)
    worst = None
    for trial in range(40):
        # build a 13-set: include 14 (a mult of 14), fill with small coprime-ish speeds
        base = random.sample(range(1, 28), 12)
        S = sorted(set(base) | {14})
        while len(S) < 13:
            x = random.randint(1, 40)
            if x not in S:
                S.append(x)
        S = sorted(S)[:13]
        # require covering (has mult of 14) and non-tight
        if 14 not in S and 28 not in S:
            continue
        r = report(f"rand#{trial}", S)
        if not r['tight'] and r['Larc'] > 0:
            if worst is None or r['Larc'] < worst['Larc']:
                worst = r
    if worst:
        print(f"\nWORST non-tight largest-arc among random: {worst['name']} "
              f"Larc={float(worst['Larc']):.5f} Darc={worst['Darc']}")
