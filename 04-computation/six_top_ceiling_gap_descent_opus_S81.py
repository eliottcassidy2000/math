#!/usr/bin/env python3
"""
opus-2026-07-05-S81 -- HYP-4116: the six-top ceiling + gap descent (l >= 7 route correction).

Verifies, numerically:
  PART 1 (CEILING): any valid uniform fee (must dominate teeth mass at EVERY window
    placement) is >= its positional mean 2*rho*L per top.  Hence sum over l tops
    >= 2*rho*l*L, and tower_step_12's criterion  sum < L  is UNSATISFIABLE for
    2*rho*l >= 1  (l >= 7 at rho = 1/13).  We compute true sup-fees by direct scan and
    confirm sup >= mean for many (w, L), and that l=7 always busts while l=6 with huge
    tops passes.
  PART 2 (DESCENT): the gap-descent lemma (w*|J| >= 2 ==> J contains a full inter-tooth
    gap of length (1-2rho)/w, pointwise dist >= rho) and the spread tower (consecutive
    ratios >= 2/(1-2rho)): dodges ANY number of spread tops.  Also demonstrates the
    bottom-cluster failure (tops 14..20 kill naive descent) -- the residual is real.
  PART 3 (GRID PROBE): for bottom-cluster l >= 7 lift patterns (visible lifts,
    k mod 13 != 0), does a strict 169-witness  t = a/169  with min_i dist >= 14/169
    exist?  This is the S77-row route one level up; probes the residual's finite path.
"""
from fractions import Fraction as F
import itertools, random

RHO = F(1, 13)

def dist_to_Z(x: F) -> F:
    fl = x - F(int(x) if x >= 0 else int(x) - (0 if x == int(x) else 1))
    return min(fl, 1 - fl)

# ---------------------------------------------------------------- PART 1: CEILING
print("=" * 72)
print("PART 1: THE SIX-TOP CEILING (fee-mean lower bound)")
print("=" * 72)

def sup_fee(w: int, L: F, samples=2000) -> F:
    """True sup over window placements t* of |bad_w cap [t*-L/2, t*+L/2]|,
    bad_w = {t: dist(w t, Z) < rho}.  Teeth: ((m-rho)/w, (m+rho)/w), width 2rho/w,
    period 1/w.  Sup over t* in one period, scanned at breakpoints exactly:
    mass(t*) is piecewise linear with breakpoints where window edges hit tooth edges;
    we scan candidate placements = tooth edges +- window edges."""
    h = 2 * RHO / w          # tooth width
    P = F(1, w)              # period
    # candidate left-endpoints x of window: align window edges with tooth edges
    cands = set()
    m_lo, m_hi = -2, int(L * w) + 3
    for m in range(m_lo, m_hi):
        for edge in ((F(m) - RHO) / w, (F(m) + RHO) / w):
            cands.add(edge)          # window left at tooth edge
            cands.add(edge - L)      # window right at tooth edge
    best = F(0)
    for x in cands:
        # mass of teeth in [x, x+L]
        mass = F(0)
        m0 = int((x * w)) - 2
        m1 = int(((x + L) * w)) + 2
        for m in range(m0, m1 + 1):
            lo, hi = (F(m) - RHO) / w, (F(m) + RHO) / w
            mass += max(F(0), min(hi, x + L) - max(lo, x))
        best = max(best, mass)
    return best

L = F(7, 468)   # the l=7 window: 2*delta = 2*(1/6-1/13)/12
print(f"window L = 2*delta = {L} ~ {float(L):.6f}   (l=7, beta=1/6, B=12)")
print(f"mean fee per top = 2*rho*L = {2*RHO*L} ~ {float(2*RHO*L):.6f}")
print(f"budget check: l * 2*rho = l * 2/13 >= 1  <=>  l >= 6.5  ==> ceiling at l = 6")
print()
print(" w      sup_fee        mean=2*rho*L   sup>=mean?")
viol = 0
for w in [14, 15, 19, 29, 67, 97, 134, 300, 1000, 2500, 10000]:
    sf = sup_fee(w, L)
    ok = sf >= 2 * RHO * L
    if not ok: viol += 1
    print(f"{w:6d}  {float(sf):.8f}   {float(2*RHO*L):.8f}   {ok}")
assert viol == 0, "MEAN BOUND VIOLATED -- theorem wrong!"
print("mean bound holds at every scanned w  ==> sup-fee >= 2*rho*L  ALWAYS")

# l=7: minimal possible fee sum vs budget
def fee_lb(w, L):
    # certified lower bound on sup_fee: exact scan when cheap, else the MEAN 2*rho*L
    return sup_fee(w, L) if w <= 3000 else 2 * RHO * L
for tops in ([14,15,16,17,18,19,98], [134,317,750,1773,4190,9904,23409],
             [10**4,10**6,10**8,10**9,10**10,10**11,10**12]):
    s = sum(fee_lb(w, L) for w in tops)
    print(f"l=7 tops={tops if tops[0]<10**4 else '[1e4..1e12 spread, mean bound]'}: "
          f"sum fee-lower-bounds = {float(s):.6f}  vs budget L = {float(L):.6f}  "
          f"--> {'BUST (>=L)' if s >= L else 'fits?!'}")
    assert s >= L, "l=7 fee sum below budget -- ceiling wrong!"
print("l = 7: fee criterion UNSATISFIABLE at every scanned configuration  [CONFIRMED]")

# l=6: mean parts total 12L/13 < L; slack L/13 ~ 0.00115 must absorb boundary terms.
# sup <= mean + 2*(tooth width) = 2*rho*L + 4*rho/w  ==> need sum 4*rho/w_j < L/13:
# w_j >= 24336/7 ~ 3477 for equal split.  Tops 6000*4^i: sum 4rho/w ~ 0.00068 < 0.00115.
tops6 = [6000 * 4**i for i in range(6)]
s6_ub = sum(2*RHO*L + 4*RHO/w for w in tops6)
print(f"l=6 tops={tops6}: sum fee-UPPER-bounds (mean + 4rho/w) = {float(s6_ub):.6f} "
      f"vs L = {float(L):.6f} --> {'fits (criterion satisfiable at l=6)' if s6_ub < L else 'bust'}")

# ---------------------------------------------------------------- PART 2: DESCENT
print()
print("=" * 72)
print("PART 2: GAP DESCENT (nesting, measure-free)")
print("=" * 72)

def gap_subinterval(w: int, a: F, Ln: F):
    """exists_gap_subinterval: if w*Ln >= 2, return [c,d] subset [a,a+Ln] with
    d-c = (1-2rho)/w and dist(w t, Z) >= rho pointwise on [c,d]."""
    assert w * Ln >= 2
    # k = ceil(w*a + rho)
    x = w * a + RHO
    k = -((-x.numerator) // x.denominator)  # ceil for Fraction
    c, d = (F(k) + RHO) / w, (F(k) + 1 - RHO) / w
    assert a <= c and d <= a + Ln, "gap not inside window!"
    assert d - c == (1 - 2 * RHO) / w
    return c, d

def check_gap_good(w, c, d, ngrid=17):
    for i in range(ngrid + 1):
        t = c + (d - c) * F(i, ngrid)
        assert dist_to_Z(w * t) >= RHO, f"gap point bad: w={w}, t={t}"

random.seed(4116)
ratio = F(2 * 13, 11)  # 2/(1-2rho) = 26/11
print(f"descent ratio 2/(1-2rho) = {ratio} ~ {float(ratio):.4f}; entry w1 >= 2/L = {float(2/L):.1f}")
for trial in range(5):
    l = random.randint(7, 15)
    w = 134 + random.randint(0, 60)
    tops, a, Ln = [], F(random.randint(0, 999), 1000), L
    for _ in range(l):
        tops.append(w)
        w = int(w * float(ratio)) + 1 + random.randint(0, w)
    for wt in tops:
        assert wt * Ln >= 2, f"entry/ratio violated: w={wt}, L={float(Ln)}"
        c, d = gap_subinterval(wt, a, Ln)
        check_gap_good(wt, c, d)
        a, Ln = c, d - c
    # final interval good for ALL tops (nesting)
    tmid = a + Ln / 2
    assert all(dist_to_Z(wt * tmid) >= RHO for wt in tops)
    print(f"  trial {trial}: l={l} spread tops in [{tops[0]},{tops[-1]}] -- "
          f"nested to |J|={float(Ln):.2e}, midpoint good for all {l} tops  [OK]")
print("spread tower: ANY number of spread tops dodged -- no fee, no measure  [CONFIRMED]")

# bottom-cluster failure: tops 14..20 with the l=7 window
print()
print("bottom-cluster reality check (tops 14..20, window L=7/468):")
a, Ln, died = F(0), L, False
for wt in [14, 15, 16, 17, 18, 19, 20]:
    if wt * Ln >= 2:
        c, d = gap_subinterval(wt, a, Ln); a, Ln = c, d - c
    else:
        # slow case: remove up to 3 teeth, keep largest piece: >= (Ln - 6rho/w)/4
        rem = (Ln - 6 * RHO / wt) / 4
        print(f"  top {wt}: SLOW (w*L={float(wt*Ln):.3f} < 2); largest-piece bound "
              f"(L - 6rho/w)/4 = {float(rem):.6f}"
              + ("  <= 0: DESCENT DIES" if rem <= 0 else ""))
        if rem <= 0: died = True; break
        Ln = rem
print(f"  ==> bottom cluster {'kills naive descent (residual is REAL)' if died else 'survived?!'}")

# ---------------------------------------------------------------- PART 3: GRID PROBE
print()
print("=" * 72)
print("PART 3: 169-GRID PROBE for the bottom-cluster residual (l>=7, visible lifts)")
print("=" * 72)

def best_169_witness(v):
    """max over a in [1,168] of min_i dist(v_i a / 169); return (best, argmax a)."""
    best, arga = F(0), None
    for aa in range(1, 169):
        m = min(dist_to_Z(F(vi * aa, 169)) for vi in v)
        if m > best: best, arga = m, aa
    return best, arga

THRESH = F(14, 169)
print(f"strict threshold: 14/169 = {float(THRESH):.6f} > 1/13 = {float(F(1,13)):.6f}")
random.seed(78)
fails = 0
cases = []
# systematic: lift exactly the pattern positions with k=1 (bottom cluster v=r+13)
for lifted in [set(range(1,8)), set(range(6,13)), {1,2,3,4,5,6,12}, {2,4,6,8,10,11,12}]:
    v = [r + 13 if r in lifted else r for r in range(1, 13)]
    cases.append(("k=1 lifts at " + str(sorted(lifted)), v))
# random k mod 13 visible patterns, l=7..9
for _ in range(8):
    l = random.randint(7, 9)
    lifted = set(random.sample(range(1, 13), l))
    v = [r + 13 * random.choice([1,2,3,5,7,12,14,25]) if r in lifted else r
         for r in range(1, 13)]
    if any((vi - r) // 13 % 13 == 0 and vi != r for r, vi in zip(range(1,13), v)):
        pass  # invisible lift present -- keep anyway, informative
    cases.append((f"random l={l}", v))
for name, v in cases:
    b, aa = best_169_witness(v)
    tag = "STRICT WITNESS" if b >= THRESH else "no 169-witness"
    if b < THRESH: fails += 1
    print(f"  {name:34s} v={v}  best={b} at a={aa}  [{tag}]")
print(f"==> {len(cases)-fails}/{len(cases)} bottom-cluster-type families have a strict a/169 witness")
print("    (failures, if any, are the deep-well analogues -- the finite row stratum)")

# spot exact-M floor via merge-grid attainment (klein HYP-4114: M attained at m/(vi+vj))
print()
print("spot exact-M floors (breakpoint grid m/(vi+vj), klein HYP-4114 attainment):")
def exact_M(v, cap=200):
    best = F(0)
    cands = set()
    for vi, vj in itertools.combinations_with_replacement(v, 2):
        s = vi + vj
        for m in range(1, min(s, cap) + 1):
            cands.add(F(m, s))
    for t in cands:
        best = max(best, min(dist_to_Z(vi * t) for vi in v))
    return best
for name, v in cases[:3]:
    M = exact_M(v)
    print(f"  {name:34s} M = {M} ~ {float(M):.4f}  (vs 1/13 ~ 0.0769, 3/19 ~ 0.158)")
print()
print("ALL PARTS DONE.")
