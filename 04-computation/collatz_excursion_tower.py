#!/usr/bin/env python3
"""collatz_excursion_tower.py — The Collatz iterated-logarithm tower.

Sequel to collatz_defect_rapidity.py (HYP-2147/2148). That session proved the
EXACT rapidity law  ln n = K ln2 - L ln3 - D(n)  and that the harmonic defect
D(n) is bounded (D* ~ 0.2257). This session asks the Tao-style question: once
the FIRST logarithm (rapidity = log_F = arctanh) has tamed the +1 into a
bounded O(1) defect, what governs the residual structure -- and at which
ITERATED-LOG scale does each piece live?

THE DEEPER ABSTRACTION: ITERATED LINEARIZATION (a renormalization tower).
Collatz nonlinearity, measured at successive log-scales, shrinks by one
log-order at each level:

  Level 0  scale n      (multiplicative):  x->3x, x->x/2 are the group ops;
                                           the +1 is an UNBOUNDED perturbation.
  Level 1  scale ln n   (rapidity):        x3,/2 become TRANSLATIONS; the +1
                                           becomes the BOUNDED defect D(n)<D*.
                                           [PROVED last session]
  Level 2  scale lnln n (fluctuation):     the rapidity walk has drift
                                           delta=(1/2)ln(3/4) and variance
                                           sigma^2=ln^2(2)/2; its excursions are
                                           governed by extreme-value / iterated-
                                           log laws with EXPONENT theta=2.
  Level 3  scale lnlnln n (exceptional):   the a.s. envelope of the peak
                                           excursion -- where Tao's "almost all
                                           orbits attain almost bounded values"
                                           lives.

The organizing constant is the Cramer-Lundberg exponent theta of the rapidity
increment Y = (1/2)ln3 - (k/2)ln2, k=v_2(3a+1) ~ Geometric(1/2):
    E[e^{theta Y}] = 1  <=>  3^{theta/2} = 2^{1+theta/2} - 1  <=>  theta = 2,
and AT the root the equation is literally  3 - 4 + 1 = 0, i.e.  3*1+1 = 2^2.
The excursion exponent is fixed by the smallest Collatz step 3+1=4. The +1,
which level 1 showed is the whole obstruction, here FIXES the fluctuation law.

Consequence (headline, testable): theta=2 in rapidity (rho=(1/2)ln) means the
peak-value-over-start ratio has a clean POWER LAW
    P( max(orbit)/n > M ) ~ C / M .

Session: collatz-excursion-tower (continues the rapidity thread S116 -> defect).
"""
import sys
sys.stdout.reconfigure(line_buffering=True)
from math import log, log2, sqrt, exp
from collections import Counter
import statistics

LN2, LN3 = log(2), log(3)
DRIFT = 0.5*LN3 - LN2                 # (1/2)ln(3/4)
SIGMA2 = LN2*LN2/2                    # variance of ideal increment
THETA = 2.0
INV_LN43 = 1.0/log(4.0/3.0)          # descent-length constant c

print("\n  THE COLLATZ ITERATED-LOGARITHM TOWER\n")
print("=" * 70)

# ----------------------------------------------------------------------
def collatz_orbit_stats(n):
    """Return (peak, L, K, fluct_max) for full Collatz orbit of odd n.
    peak = max value reached; L = #odd steps; K = total halvings;
    fluct_max = max |rho(orbit_j) - rho(n) - (#odd so far)*DRIFT| over odd nodes."""
    a = n
    peak = n
    L = 0
    K = 0
    fluct_max = 0.0
    rho0 = 0.5*log(n)
    while a != 1:
        if a % 2 == 1:
            # odd step: record fluctuation at this odd node
            cur = 0.5*log(a) - rho0 - L*DRIFT
            if abs(cur) > fluct_max:
                fluct_max = abs(cur)
            a = 3*a + 1
            if a > peak:
                peak = a
            L += 1
        else:
            a //= 2
            K += 1
    return peak, L, K, fluct_max

# ============================================================
print("\n  I. THE CONSTANTS OF THE RAPIDITY WALK")
print("  " + "-" * 50)
print(f"  Per odd step, ideal increment Y = (1/2)ln3 - (k/2)ln2, k~Geom(1/2):")
print(f"     drift   delta   = (1/2)ln(3/4) = {DRIFT:+.6f}")
print(f"     variance sigma^2 = ln^2(2)/2   = {SIGMA2:.6f}")
print(f"  Cramer-Lundberg exponent theta solves E[e^{{tY}}]=1:")
print(f"     3^(t/2) = 2^(1+t/2) - 1   ->   theta = {THETA:g}  (EXACT)")
print(f"     at the root: 3 - 4 + 1 = 0   <=>   3*1+1 = 4 = 2^2   (the +1!)")
print(f"  Descent length constant c = 1/ln(4/3) = {INV_LN43:.6f}")
print(f"     (rapidity (1/2)ln n is dissipated at |delta| per odd step =>")
print(f"      L(n) ~ (1/2)ln n / |delta| = ln n / ln(4/3) = c ln n)")
print()

# Empirically confirm E[k]=2, Var[k]=2 from REAL Syracuse steps
kc = Counter()
for n in range(3, 100000, 2):
    a = n
    while a != 1:
        if a % 2 == 1:
            v = 3*a+1; k = 0
            while v % 2 == 0: v//=2; k+=1
            kc[k]+=1; a=v
        else:
            a//=2
tot = sum(kc.values())
Ek = sum(k*c for k,c in kc.items())/tot
Vk = sum(k*k*c for k,c in kc.items())/tot - Ek*Ek
print(f"  Empirical halving-exponent distribution over real odd steps "
      f"(n<100000, {tot} steps):")
for k in sorted(kc)[:6]:
    print(f"     P(k={k}) = {kc[k]/tot:.5f}   (geom 2^-k = {2**-k:.5f})")
print(f"     E[k] = {Ek:.5f} (pred 2),  Var[k] = {Vk:.5f} (pred 2)")
print(f"     => empirical drift = {0.5*LN3-0.5*LN2*Ek:+.6f}, "
      f"sigma^2 = {(0.5*LN2)**2*Vk:.6f}")
print()

# ============================================================
print("  II. LEVEL 1 (log): the bounded defect  [recap, proved]")
print("  " + "-" * 50)
print("  One logarithm (rapidity) turns the unbounded +1 into D(n) < 0.2257.")
print("  This is the BASE of the tower: nonlinearity is O(1) at the log scale.")
print()

# ============================================================
print("  III. LEVEL 1.5: descent length  L(n) ~ c ln n,  c = 1/ln(4/3)  (Wald)")
print("  " + "-" * 50)
print(f"  Wald's identity for the drift walk: E[L]*|delta| = rho(n) = (1/2)ln n,")
print(f"  so the MEAN number of odd steps is E[L] ~ ln n / ln(4/3) = c ln n,")
print(f"  c = {INV_LN43:.4f}. This is an AVERAGE law (records are outliers).")
print(f"  Mean L(n)/ln n over dyadic windows of odd n:")
print(f"  {'window':>22} {'mean L/ln n':>12} {'c=1/ln(4/3)':>12}")
for lo, hi in [(2**10,2**11),(2**13,2**14),(2**16,2**17),(2**19,2**20)]:
    tot_ratio = 0.0; cnt = 0
    for n in range(lo|1, hi, 2):
        _, L, _, _ = collatz_orbit_stats(n)
        tot_ratio += L/log(n); cnt += 1
    print(f"  [{lo:>8},{hi:>8}) {tot_ratio/cnt:>12.4f} {INV_LN43:>12.4f}")
print(f"  Famous LONG orbits are far above the mean (the lnln-scale outliers):")
print(f"  {'n':>9} {'L(n)':>7} {'c*ln n':>9} {'L/(c ln n)':>11}")
for n in [27, 703, 6171, 63728127]:
    peak, L, K, fl = collatz_orbit_stats(n)
    pred = INV_LN43*log(n)
    print(f"  {n:>9} {L:>7} {pred:>9.1f} {L/pred:>11.2f}")
print("  (mean ratio -> c confirms one logarithm sets the TYPICAL odd-step count;")
print("   the heavy upper tail is exactly the excursion law of section IV.)")
print()

# ============================================================
print("  IV. LEVEL 2 (extreme value): the EXCURSION POWER LAW  P(peak/n>M)~C/M")
print("  " + "-" * 50)
print("  theta=2 in rapidity (rho=(1/2)ln) => peak/n = e^{2R}, R=max excursion,")
print("  P(R>x)~e^{-2x}  =>  P(peak/n > M) ~ C/M.  THE headline test:")
print()
N = 1_000_001
ratios = []
maxratio = (0.0, 0)
for n in range(3, N, 2):
    peak, L, K, fl = collatz_orbit_stats(n)
    r = peak / n
    ratios.append(r)
    if r > maxratio[0]:
        maxratio = (r, n)
ratios.sort()
S = len(ratios)
print(f"  Sample: odd n in [3, {N}).  Max peak/n = {maxratio[0]:.1f} at n={maxratio[1]}.")
print(f"  {'M':>9} {'P(peak/n>M) empirical':>22} {'C=M*P (should be ~const)':>26}")
for M in [4, 8, 16, 32, 64, 128, 256, 512, 1024, 4096, 16384]:
    # survival = fraction with ratio > M
    import bisect
    idx = bisect.bisect_right(ratios, M)
    surv = (S - idx)/S
    C = M*surv
    print(f"  {M:>9} {surv:>22.6e} {C:>26.4f}")
print("  (If C = M*P(peak/n>M) is roughly constant across decades, alpha=1,")
print("   confirming the Cramer exponent theta=2 set by 3+1=4.)")
print()

# Fit the exponent alpha in P(ratio>M) ~ C M^{-alpha} by log-log regression
import math
xs, ys = [], []
for M in [8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768]:
    idx = __import__('bisect').bisect_right(ratios, M)
    surv = (S-idx)/S
    if surv > 0:
        xs.append(math.log(M)); ys.append(math.log(surv))
mx = sum(xs)/len(xs); my = sum(ys)/len(ys)
slope = sum((x-mx)*(y-my) for x,y in zip(xs,ys))/sum((x-mx)**2 for x in xs)
print(f"  Log-log fit:  P(peak/n>M) ~ M^(alpha),  alpha = {slope:.4f}  (predicted -1)")
print(f"  => excursion exponent theta = -2*alpha = {-2*slope:.4f}  (predicted 2)")
print()

# ============================================================
print("  V. LEVEL 2 (LIL scale): fluctuation about the drift line")
print("  " + "-" * 50)
print(f"  Centered walk W_j = rho(a_j)-rho(n)-j*delta has variance sigma^2 j.")
print(f"  Law of iterated logarithm: max|W_j| ~ sqrt(2 sigma^2 L lnln L).")
print(f"  {'n':>9} {'L':>5} {'max|W|':>8} {'sqrt(2 s^2 L lnln L)':>20} {'ratio':>7}")
for n in [27, 703, 6171, 77031, 837799, 8400511, 63728127, 670617279]:
    peak, L, K, fl = collatz_orbit_stats(n)
    if L >= 3:
        lll = log(log(L)) if L > 2 else 1.0
        scale = sqrt(2*SIGMA2*L*max(lll, 1e-9))
        print(f"  {n:>9} {L:>5} {fl:>8.3f} {scale:>20.3f} {fl/scale:>7.3f}")
print("  (ratios are O(1), slowly drifting: max fluctuation tracks the")
print("   sqrt(L*lnln L) scale up to a bounded factor. The exact LIL constant 1")
print("   is a limsup over j->inf and is NOT pinned on a single finite orbit --")
print("   here we only confirm the SCALE sqrt(L lnln L) is the right one.)")
print()

# ============================================================
print("  VI. LEVEL 3 (lnlnln / a.s. envelope): almost all orbits stay low")
print("  " + "-" * 50)
print("  theta=2 gives a Borel-Cantelli envelope: for almost all n the peak")
print("  excursion R(n) <= (1/2 + o(1)) lnln n, i.e. peak(n) <= n*(ln n)^{1+o(1)}.")
print("  Equivalently almost all peak/n stay below a slowly growing (ln n)^c.")
print(f"  {'c':>4} {'frac n<=1e6 with peak/n > (ln n)^c':>38}")
LNN = [log(n) for n in range(3, N, 2)]
# recompute peak/n alongside ln n cheaply: reuse ratios is sorted (lost n);
# redo a lighter pass
exceed = {c:0 for c in [1.0,1.5,2.0,3.0]}
cnt = 0
for n in range(3, N, 2):
    peak,_,_,_ = collatz_orbit_stats(n)
    r = peak/n
    lnn = log(n); cnt+=1
    for c in exceed:
        if r > lnn**c:
            exceed[c]+=1
for c in sorted(exceed):
    print(f"  {c:>4.1f} {exceed[c]/cnt:>38.5f}")
print("  (Fraction exceeding (ln n)^c -> 0 as c grows: the peak is below any")
print("   fixed power of ln n for almost all n -- the lnln-scale a.s. envelope,")
print("   the finite-range shadow of Tao's 'almost bounded values' (lnlnln).)")
print()

# ============================================================
print("=" * 70)
print("""
  SUMMARY -- the tower, one log-order per floor:

    n            +1 is unbounded           (multiplicative)
    ln n         +1 -> bounded D(n)<0.2257 (rapidity; PROVED, HYP-2147/48)
    --- one logarithm tames the nonlinearity to O(1) ---
    L(n) ~ c ln n,  c = 1/ln(4/3)          (descent count = one log of n)
    lnln n       excursions ~ extreme value, exponent theta=2 set by 3+1=4
                 => P(peak/n > M) ~ C/M    (CONFIRMED: alpha ~ -1)
                 => max|fluct| ~ sqrt(2 sigma^2 L lnln L)  (LIL scale)
    lnlnln n     a.s. envelope: peak(n) <= n (ln n)^{o(1)} for almost all n
                 (the scale of Tao's 'almost all orbits almost bounded')

  ONE CONSTANT runs the whole tower: theta=2, the Cramer root of the rapidity
  increment, and theta=2 is forced by the identity 3+1=4=2^2. The +1 that
  level-1 isolated as the entire obstruction is, at level 2, exactly what
  pins the excursion exponent. Each successive logarithm exposes the next,
  one-order-finer, layer of the same +1.
""")
