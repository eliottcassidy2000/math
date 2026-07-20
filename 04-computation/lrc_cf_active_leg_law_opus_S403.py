"""
opus-2026-07-19-S403 (H7 / HYP-7970): THE CF ACTIVE-LEG LAW.

Observation (hand-verified on 1/14, 3/41, 4/127 before this script): at a
straddling maximizer t* = a/q (reduced) of M(V), the SMALLER active leg v is a
CONTINUED-FRACTION CONVERGENT DENOMINATOR q_m of a/q, and the represented
determinant D = M * s_pair equals the convergent's remainder
delta_m = |q_m * a - p_m * q| (times s_pair/q when the pair sum is a multiple
of the reduced denominator).

Test on: extremals, gap/ladder families, cross-N discoveries, loose controls.
For each family: exact M and ALL maximizing t* (q <= QMAX, integer arithmetic);
at each t*: all straddling pairs (v_below-integer, v_above-integer) with
v_i t* = a_i + M, v_j t* = a_j - M; s_pair, D = M*s_pair; then check legs
against convergent denominators of t* and D against delta_m * (s_pair/q).
"""
from math import gcd
from fractions import Fraction

def exact_max(V, qmax):
    bg, bq = 0, 1  # best M = bg/bq
    wits = []
    for q in range(2, qmax + 1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            m = q
            for v in V:
                r = (v * a) % q
                r = min(r, q - r)
                if r < m:
                    m = r
                    if m * bq < bg * q: break   # early exit: already worse
            if m * bq > bg * q:
                bg, bq, wits = m, q, [(a, q)]
            elif m * bq == bg * q and bg > 0 and len(wits) < 4:
                wits.append((a, q))
    return bg, bq, wits

def convergents(a, q):
    """convergents p/qd of a/q (0<a<q), return list of (p, qd)."""
    cf = []
    x, y = a, q
    # CF of a/q = [0; c1, c2, ...] -- compute on q/a
    num, den = q, a
    quots = []
    while den:
        quots.append(num // den)
        num, den = den, num % den
    # convergents of q/a are h/k; convergents of a/q are k/h
    h0, h1 = 1, quots[0]
    k0, k1 = 0, 1
    conv = [(k1, h1)]  # (p, qd) for a/q: p/qd = k/h
    for c in quots[1:]:
        h0, h1 = h1, c * h1 + h0
        k0, k1 = k1, c * k1 + k0
        conv.append((k1, h1))
    return [(0,1)] + conv   # include 0/1

def analyze(name, V, qmax=300):
    bg, bq, wits = exact_max(V, qmax)
    M = Fraction(bg, bq)
    print(f"\n== {name}: M = {M} (scan q<={qmax}); maximizers {wits}")
    for (a, q) in wits[:2]:
        conv = convergents(a, q)
        cdens = {qd: (p, qd) for (p, qd) in conv}
        # active speeds at t* = a/q
        above, below = [], []   # above-integer: frac = M; below: frac = 1-M... careful
        for v in V:
            r = (v * a) % q
            if min(r, q - r) * bq == bg * q // gcd(q,q):  # r/q == M as fractions: r*bq == bg*q? M=bg/bq, dist=r/q: equal iff r*bq==bg*q
                pass
        # simpler exact: dist == M  <=>  r == M*q; M*q = bg*q/bq must be integer when q multiple of bq
        for v in V:
            r = (v * a) % q
            d = min(r, q - r)
            if Fraction(d, q) == M:
                (above if r == d else below).append(v)  # r==d: just above integer
        print(f"   t* = {a}/{q}: conv denoms {sorted(cdens)}; active above-int {above}, below-int {below}")
        for vi in above:
            for vj in below:
                s_pair = vi + vj
                D = M * s_pair
                if D.denominator != 1: continue
                D = int(D)
                legs = (vi, vj)
                hit = None
                for leg in legs:
                    if leg in cdens:
                        p, qd = cdens[leg]
                        delta = abs(qd * a - p * q)
                        scale = Fraction(s_pair, q)
                        pred = delta * scale
                        hit = (leg, qd, delta, scale, pred, pred == D)
                print(f"     pair ({vi},{vj}): s={s_pair}, D={D}; "
                      f"conv-leg check: {hit}")

FAMS = [
 ("AP13 {1..13}",            list(range(1,14)), 300),
 ("GW {1..11,13,24}",        list(range(1,12))+[13,24], 300),
 ("{1..12,14} (1/13)",       list(range(1,13))+[14], 300),
 ("{1..11,13,36} (3/41)",    list(range(1,12))+[13,36], 300),
 ("ladder m=4 {..,48}",      list(range(1,12))+[13,48], 300),
 ("ladder m=5 {..,60}",      list(range(1,12))+[13,60], 300),
 ("ladder m=6 {..,72}",      list(range(1,12))+[13,72], 300),
 ("loose {1..12,15}",        list(range(1,13))+[15], 300),
 ("primes13",                [2,3,5,7,11,13,17,19,23,29,31,37,41], 300),
 ("mixed control",           [3,7,10,14,22,25,31,38,44,52,57,63,70], 300),
 ("N19 {1..17,19,54} (3/59)",list(range(1,18))+[19,54], 300),
 ("N31 {1..29,31,120} (4/127)", list(range(1,30))+[31,120], 300),
 ("N61 {1..59,61,240} (4/247)", list(range(1,60))+[61,240], 320),
]
if __name__ == "__main__":
    print("CF ACTIVE-LEG LAW test -- opus-S403")
    for name, V, qm in FAMS:
        analyze(name, V, qm)
