#!/usr/bin/env python3
"""
lrc14_neardilate_adversary_deathstar_S14.py — death-star-2026-07-12-S14

The near-dilate adversary vs the large-diameter atoms (THM-720 / opus-S243 / THM-721).

CONSTRUCTION: V_L = {L, 2L, ..., 12L, 13L+1}, L divisible by 2^3*3^2*5*7*13 (=> divisor-
complete: iL covers 2..10,12,13,14 via L, 11 via 11L; primitive: gcd(L,13L+1)=1).
Diameter 12L+1 grows without bound.

CLAIMS TESTED (all exact, pair-sum ruler enumeration = THM-668 Part 2, pure integers):
  (1) M(V_L) = 1/13 + O(1/L)  — the adversarial large-diameter floor is 1/13, NOT growing
      (THM-720's sampled "min M grows with diameter" is a random-sample artifact).
  (2) V_L escapes the 1D descent legs as stated:
      - mac-mini cont.49 r<=12 leg: needs <=12 distinct lifts at a compressed scale — V_L has
        13 distinct lifts at every atom-applicable scale (census over ALL scales L').
      - opus-S243 Case A trigger (<=6 coprime-to-30030) FIRES, but no small-lift decomposition
        exists — trigger does not imply mechanism.
      - opus-S243 Case B trigger (speed >= lcm(2..14)=360360) fires for L=32760, but the max
        element is NOT far (ratio 13L+1 : 12L ~ 1.083) — far-peel hypothesis fails.
  (3) The 2D atom (THM-721) handles it: at scale L, j=1 impure runner (b=(0,...,0,1), B=1),
      pure lifts {1..12} (12 distinct, LRC(13) floor 1/13), u-escape kills the impure arc:
      reach(V) >= 1/13 - 3B/(2L). Also the 2D PIN: M <= reach2 + 3B/(2L) with reach2 = 1/13 here.
  (4) j-SERIES: perturbing j runners (j=1..7) — M tracks the u-escape floor
      min(M(pure), 1/(2j)) - 3B/(2L); j=7 is the union-bound boundary (exactly 1/14).
  (5) kps blocker [200,...,1856]: compressibility census — incoherent at EVERY scale.
      Unit tests: M({1..13})=1/14, M({1..12})=1/13, M(blocker)=53/227.
"""
from fractions import Fraction
from math import gcd
from functools import reduce
import sys, time

def exact_M(v):
    """Exact M via THM-668: local maxima only at t=p/(v_i+v_j) (i<=j). Scan every pair-sum
    ruler q, p in [1, q//2] (p <-> q-p symmetry), early-exit inner loop with running-best
    integer threshold. Returns (Fraction M, q*, p*). Pure Python, exact integers."""
    vv = tuple(sorted(v))
    n = len(vv)
    rulers = sorted(set(vv[i] + vv[j] for i in range(n) for j in range(i, n)))
    bn, bd = 0, 1          # best fraction bn/bd
    bq, bp = None, None
    for q in rulers:
        thr = (bn * q) // bd          # candidate must have m > thr (integer floor works:
                                      # m/q > bn/bd  <=>  m*bd > bn*q  <=  m > floor(bn*q/bd) when not exact tie;
                                      # ties resolved by exact compare below)
        half = q >> 1
        for p in range(1, half + 1):
            m = q
            for x in vv:
                r = (x * p) % q
                if r > half:
                    r = q - r
                if r < m:
                    m = r
                    if m <= thr:
                        break
            if m > thr and m * bd > bn * q:   # exact improvement
                bn, bd, bq, bp = m, q, q, p
                g = gcd(bn, bd)
                bn //= g; bd //= g
                thr = (bn * q) // bd
    return Fraction(bn, bd), bq, bp

def is_divisor_complete(v):
    return all(any(x % d == 0 for x in v) for d in range(2, 15))

def compress_census(v, max_scale=None):
    """For every scale L', v_i = L'*k_i + b_i (nearest multiple). Report scales where a descent
    leg applies: legA_1d (mac-mini cont.49 / THM-636 analog: all k>=1, <=12 distinct lifts,
    L' > 182*B) or legB_2d (THM-721: j = #{b!=0 or k==0} <= 6, L' > 273*B)."""
    vv = sorted(v)
    vmax = vv[-1]
    out = []
    for L in range(2, (max_scale or vmax) + 1):
        ks, bs = [], []
        for x in vv:
            k = (x + L // 2) // L
            ks.append(k); bs.append(x - L * k)
        B = max(abs(b) for b in bs)
        j = sum(1 for k, b in zip(ks, bs) if b != 0 or k == 0)
        r_all = len(set(k for k in ks if k > 0))
        legA = all(k >= 1 for k in ks) and r_all <= 12 and L > 182 * max(B, 1)
        legB = j <= 6 and L > 273 * max(B, 1)
        if legA or legB:
            pure = set(k for k, b in zip(ks, bs) if b == 0 and k > 0)
            out.append(dict(L=L, B=B, j=j, r=r_all, legA=legA, legB=legB,
                            pure_distinct=len(pure)))
    return out

def coprime_30030_count(v):
    return sum(1 for x in v if gcd(x, 30030) == 1)

def report_family(name, v, expect=None):
    v = sorted(v)
    g = reduce(gcd, v)
    t0 = time.time()
    M, q, p = exact_M(v)
    dt = time.time() - t0
    print(f"\n=== {name} ===")
    print(f"  v = {v}")
    print(f"  primitive={g==1}  DC={is_divisor_complete(v)}  diam={v[-1]-v[0]}  Vmax={v[-1]}"
          f"  #coprime-30030={coprime_30030_count(v)}")
    print(f"  EXACT M = {M} = {float(M):.7f}   witness t = {p}/{q}   [{dt:.1f}s]")
    print(f"  vs 1/14: {'ABOVE' if M > Fraction(1,14) else 'AT/BELOW'}   M/(1/13) = {float(M*13):.5f}")
    if expect is not None:
        ok = (M == expect)
        print(f"  UNIT TEST expect {expect}: {'PASS' if ok else 'FAIL'}")
        if not ok:
            sys.exit(f"UNIT TEST FAILED for {name}")
    sys.stdout.flush()
    return M

print("=" * 78)
print("UNIT TESTS (THM-668 exact-M enumeration)")
print("=" * 78)
report_family("AP {1..13} (the wall)", list(range(1, 14)), expect=Fraction(1, 14))
report_family("AP {1..12} (12-runner tight)", list(range(1, 13)), expect=Fraction(1, 13))
BLOCKER = [200, 496, 540, 656, 851, 921, 935, 1122, 1482, 1680, 1835, 1849, 1856]
# NOTE: kps cont.47 (HYP-6120) reported "true M = 53/227" for the blocker; that is the margin at
# one pair event (t=265/1135, ruler 1135=200+935) — a LOWER bound, not the max. The exact max is
# 406/1669 at t=1147/3338 (ruler 3338=1482+1856), independently found by klein-S264's Parseval
# floor and by this enumeration. Qualitative conclusion (very loose) unchanged/strengthened.
report_family("kps blocker (HYP-6120)", BLOCKER, expect=Fraction(406, 1669))

print("\n" + "=" * 78)
print("(1)+(3) NEAR-DILATE ADVERSARY SERIES: V_L = {L,2L,...,12L,13L+1}")
print("=" * 78)
for L in [1560, 3900, 8190, 32760]:
    v = [i * L for i in range(1, 13)] + [13 * L + 1]
    M = report_family(f"near-dilate L={L}" + (" [DC]" if is_divisor_complete(v) else ""), v)
    lower = Fraction(1, 13) - Fraction(3, 2 * L)
    upper = Fraction(1, 13) + Fraction(3, 2 * L)
    print(f"  THM-721 floor 1/13 - 3/(2L) = {float(lower):.7f}: {'OK' if M >= lower else 'VIOLATION'}"
          f"   | 2D pin M <= 1/13 + 3/(2L): {'OK' if M <= upper else 'VIOLATION'}"
          f"   (M - 1/13 = {float(M - Fraction(1,13)):+.2e})")
    sys.stdout.flush()

print("\n" + "=" * 78)
print("(2) ESCAPE CENSUS on the DC adversary (L=32760): any scale with a 1D leg?")
print("=" * 78)
L = 32760
vD = [i * L for i in range(1, 13)] + [13 * L + 1]
cs = compress_census(vD)
legA_scales = [c for c in cs if c['legA']]
legB_scales = [c for c in cs if c['legB']]
print(f"  scales where legA (1D r<=12 descent) applies: {len(legA_scales)}"
      + (f"  e.g. {legA_scales[:3]}" if legA_scales else "  — NONE (escapes mac-mini cont.49 r<=12 leg)"))
print(f"  scales where legB (2D j<=6 atom, THM-721) applies: {len(legB_scales)}"
      + (f"  e.g. {[dict(L=c['L'],B=c['B'],j=c['j'],pure=c['pure_distinct']) for c in legB_scales[:4]]}" if legB_scales else "  — NONE"))
print(f"  opus-S243 Case A trigger (#coprime-30030 <= 6): count={coprime_30030_count(vD)} — "
      f"{'FIRES' if coprime_30030_count(vD) <= 6 else 'no'}; but legA census above shows NO small-lift scale")
print(f"  opus-S243 Case B trigger (Vmax >= lcm(2..14)=360360): Vmax={max(vD)} — "
      f"{'FIRES' if max(vD) >= 360360 else 'no'}; far-ratio Vmax/second = {max(vD)/sorted(vD)[-2]:.4f} (NOT far)")

print("\n" + "=" * 78)
print("(4) j-SERIES at L=1560: perturb top-j runners by +1..+j — u-escape floor tracking")
print("=" * 78)
L = 1560
for j in range(1, 8):
    v = sorted([i * L for i in range(1, 14 - j)] +
               [(13 - j + 1 + a) * L + (a + 1) for a in range(j)])
    pure = list(range(1, 14 - j))
    M, q, p = exact_M(v)
    Mpure, _, _ = exact_M(pure)
    B = j  # offsets 1..j
    floor = min(Mpure, Fraction(1, 2 * j)) - Fraction(3 * B, 2 * L)
    g = reduce(gcd, v)
    print(f"  j={j}: M = {str(M):>9} = {float(M):.6f} | floor min({Mpure}, 1/{2*j}) - 3B/(2L) = "
          f"{float(floor):.6f}  {'OK' if M >= floor else 'VIOLATION'}  (prim={g==1}, witness {p}/{q})")
    sys.stdout.flush()

print("\n" + "=" * 78)
print("(5) kps BLOCKER compressibility census — expect: NO leg at any scale")
print("=" * 78)
cs = compress_census(BLOCKER)
print(f"  scales with legA or legB: {len(cs)}" + (f" — {cs[:5]}" if cs else "  (incoherent at every scale: pair-sum/coverage domain)"))

print("\nDONE.")
