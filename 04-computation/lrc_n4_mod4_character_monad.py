#!/usr/bin/env python3
"""
lrc_n4_mod4_character_monad.py    monad-researcher-2026-06-01-S1

Push the S526 covering+character methodology to n=4 (3 runners, threshold 1/4).

THE n=4 HARMONIC SUBSTRATE.
  Safe indicator g(x) = 1[||x|| >= 1/4] = indicator of the half-circle [1/4, 3/4].
  Fourier coeffs:  g_0 = 1 - 2/n = 1/2,
                   g_k = -sin(2 pi k / 4)/(pi k) = -sin(pi k/2)/(pi k) = -chi4(k)/(pi k),
  where chi4 is the UNIQUE odd Dirichlet character mod 4:
        chi4(odd=1 mod4)=+1, chi4(3 mod4)=-1, chi4(even)=0,   and  sin(pi k/2)=chi4(k).
  So the n=4 character is the mod-4 odd character -- the exact analogue the S526
  handoff asked for (n=3 used the mod-3 Legendre symbol).

DECOMPOSITION.  For primitive speeds (a,b,c), with Lambda = {k in Z^3 : k.(a,b,c)=0},
    |SAFE| = sum_{k in Lambda} prod_i g_{k_i}
           = 1/8                              (k=0, the independence main term (1-2/n)^{n-1})
           + R2                                (support-2 / pairwise resonances)
           + R3 .                              (support-3 genuine 3-term resonances)

CLAIM (proved in reflection, verified here):  pairwise term has a CLOSED FORM
    R2^{(a,b)} = chi4(a')chi4(b') / (4 a' b'),     a'=a/d, b'=b/d, d=gcd(a,b),
  i.e. zero unless both cofactors are odd; this is the mod-4 analogue of the n=3
  Legendre closed form (2/9)chi(a)chi(b)/(ab).

This script:
  (1) verifies sin(pi k/2)=chi4(k) and the g_k values,
  (2) verifies the R2 closed form against EXACT pairwise safe measure,
  (3) computes EXACT |SAFE| for all primitive triples in a range, extracts R3,
  (4) tests LRC(n=4): is |SAFE|>0 except at tight sets? locate tight sets,
  (5) tests the reduction inequality R3 >= -(1/8 + R2) and its equality cases,
  (6) probes whether R3 has sign structure / a bound.
"""
from math import gcd, sin, pi
from fractions import Fraction as F
from itertools import combinations, product

N = 4
THR = F(1, N)            # 1/4

def chi4(k):
    k %= 4
    return 0 if k % 2 == 0 else (1 if k == 1 else -1)

def frac_part_norm(x):
    """||x|| = distance to nearest integer, x a Fraction."""
    r = x - (x.numerator // x.denominator)   # in [0,1)
    return min(r, 1 - r)

# ---------- EXACT safe measure of a speed set at threshold 1/4 ----------
def exact_safe_measure_and_witness(speeds):
    """Return (measure, nonempty_bool, best_margin).
    Breakpoints: t where ||s t|| = 1/4  <=> s t = (4k +- 1)/4 => t=(4k+-1)/(4 s)."""
    # breakpoints t = num/(4s) in (0,1) with num = 1 or 3 mod 4  <=>  ||s t|| = 1/4
    bps = set([F(0), F(1)])
    for s in speeds:
        num = 1
        while num < 4 * s:
            if num % 4 in (1, 3):
                bps.add(F(num, 4 * s))
            num += 1
    pts = sorted(bps)
    measure = F(0)
    best = F(-1)
    nonempty = False
    # measure: test midpoint of each cell; also test breakpoints (closed >=) for witness
    for lo, hi in zip(pts, pts[1:]):
        mid = (lo + hi) / 2
        d = min(frac_part_norm(F(s) * mid) for s in speeds)
        if d >= THR:
            measure += (hi - lo)
        if d > best:
            best = d
    # witnesses at breakpoints (closed safe set uses >=)
    for t in pts:
        if t == 1:
            continue
        d = min(frac_part_norm(F(s) * t) for s in speeds)
        if d > best:
            best = d
    nonempty = (best >= THR)
    return measure, nonempty, best

# ---------- R2 closed form ----------
def R2_pair(a, b):
    d = gcd(a, b)
    ap, bp = a // d, b // d
    return F(chi4(ap) * chi4(bp), 4 * ap * bp)

def R2_total(speeds):
    return sum(R2_pair(x, y) for x, y in combinations(speeds, 2))

def exact_pairwise_safe(a, b):
    """exact |S_a ∩ S_b| = measure of {t: ||at||>=1/4 and ||bt||>=1/4}."""
    m, _, _ = exact_safe_measure_and_witness((a, b))
    return m

def main():
    print("LRC n=4 via covering + mod-4 odd character (monad-researcher-S1)\n")

    # (1) character identity
    print("="*70)
    print("(1) sin(pi k/2) == chi4(k), and g_k = -chi4(k)/(pi k)")
    print("="*70)
    ok = all(abs(sin(pi*k/2) - chi4(k)) < 1e-12 for k in range(-20, 21))
    print(f"  sin(pi k/2)=chi4(k) for |k|<=20: {ok}")
    print(f"  g_0 = 1/2 = (1-2/n).  Even k => g_k=0 (only ODD harmonics survive at n=4).")

    # (2) verify R2 closed form vs exact pairwise safe measure  (|S_a∩S_b| = 1/4 + R2)
    print("\n" + "="*70)
    print("(2) R2 closed form: |S_a∩S_b| =?= 1/4 + chi4(a')chi4(b')/(4 a' b')")
    print("="*70)
    bad = 0; shown = 0
    for a in range(1, 14):
        for b in range(a+1, 15):
            pred = F(1, 4) + R2_pair(a, b)
            exact = exact_pairwise_safe(a, b)
            if pred != exact:
                bad += 1
                print(f"  MISMATCH (a,b)=({a},{b}): pred={pred} exact={exact}")
            elif shown < 8 and gcd(a, b) == 1:
                print(f"  (a,b)=({a},{b}): |S∩S|={exact}=1/4+({R2_pair(a,b)})  ok")
                shown += 1
    print(f"  closed form exact-match failures: {bad}  (0 = formula proved-consistent)")

    # (3)+(4)+(5) full triple decomposition, LRC test, reduction inequality
    print("\n" + "="*70)
    print("(3-5) primitive triples: |SAFE| = 1/8 + R2 + R3; LRC; R3>=-(1/8+R2)?")
    print("="*70)
    MAXS = 16
    tot = 0; lrc_fail = 0; tight = []; ineq_fail = 0
    r3_min = F(10); r3_max = F(-10); r3_pos = 0; r3_neg = 0; r3_zero = 0
    minmargin_examples = []
    for a, b, c in combinations(range(1, MAXS+1), 3):
        if gcd(gcd(a, b), c) != 1:
            continue
        tot += 1
        meas, nonempty, best = exact_safe_measure_and_witness((a, b, c))
        r2 = R2_total((a, b, c))
        r3 = meas - F(1, 8) - r2
        # bookkeeping
        if r3 < r3_min: r3_min = r3
        if r3 > r3_max: r3_max = r3
        if r3 > 0: r3_pos += 1
        elif r3 < 0: r3_neg += 1
        else: r3_zero += 1
        # LRC test
        if not nonempty:
            lrc_fail += 1
            print(f"  *** LRC FAIL (a,b,c)=({a},{b},{c}) best={float(best):.4f}")
        if meas == 0:
            tight.append((a, b, c, nonempty, best))
        # reduction inequality: |SAFE|=1/8+R2+R3 >= 0  <=>  R3 >= -(1/8+R2)
        if r3 < -(F(1, 8) + r2):
            ineq_fail += 1
        # track smallest positive measure
        if meas > 0:
            minmargin_examples.append((meas, a, b, c))
    print(f"  primitive triples tested (speeds<= {MAXS}): {tot}")
    print(f"  LRC failures: {lrc_fail}")
    print(f"  reduction-inequality failures (R3 < -(1/8+R2)): {ineq_fail}")
    print(f"  R3 range: [{float(r3_min):.5f}, {float(r3_max):.5f}]   neg/zero/pos = {r3_neg}/{r3_zero}/{r3_pos}")
    print(f"  tight sets (measure(SAFE)=0): {len(tight)}")
    for (a, b, c, ne, bm) in tight[:20]:
        print(f"     ({a},{b},{c}): nonempty(closed)={ne} best_margin={bm}={float(bm):.4f}")

    # (6) Detailed look at the AP {1,2,3} and the tight structure
    print("\n" + "="*70)
    print("(6) the AP {1,2,3}: exact decomposition")
    print("="*70)
    for trip in [(1,2,3),(1,2,5),(1,3,5),(1,4,5),(2,3,7)]:
        meas, ne, best = exact_safe_measure_and_witness(trip)
        r2 = R2_total(trip)
        r3 = meas - F(1,8) - r2
        print(f"  {trip}: |SAFE|={meas}={float(meas):.4f}  1/8={F(1,8)} R2={r2}={float(r2):+.4f} R3={r3}={float(r3):+.4f}  nonempty={ne} margin={float(best):.4f}")

if __name__ == "__main__":
    main()
