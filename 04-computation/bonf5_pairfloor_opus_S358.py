# opus-2026-07-17-S358 -- HYP-7410: where the all-regime pair floor (THM-1025)
# actually lands in the LRC(14) budget.
#   uncovered >= B5 = 1 - S1 + S2 - S3 + S4 - S5   (odd Bonferroni truncation)
# S1 = 13*(2lam) = 13/7 is FIXED.  B5 > 0 needs S2, S4 LARGE and S3, S5 SMALL.
# My THM-1025 floor is exactly an S2 LOWER bound, valid in EVERY regime.
from fractions import Fraction as F
from math import gcd, comb
import random, itertools
LAM = F(1,14)

def rho(a, b):
    g = gcd(a,b); tot = F(0); Mb = LAM*(a+b); m = 0
    while True:
        if F(m) >= Mb*g and m > 0: break
        for mm in ([0] if m == 0 else [m,-m]):
            if mm % g: continue
            if F(abs(mm)) >= Mb: continue
            w = min(2*LAM*a, 2*LAM*b, LAM*(a+b)-abs(mm))
            if w > 0: tot += F(w, a*b)*g
        m += g
        if m > int(Mb)+2: break
    return tot

SMALL = {}
for a in range(1,13):
    for b in range(a,13):
        if gcd(a,b)==1 and a*b <= 12: SMALL[(a,b)] = rho(a,b)

def floor_1025(a, b):
    """THM-1025: all-regime pair floor (sawtooth above threshold, table below)."""
    g = gcd(a,b); ap, bp = a//g, b//g
    if ap*bp >= 13: return F(1,49) - F(1, 4*ap*bp)
    return SMALL[(ap,bp)]

print("THE BUDGET.  S1 = 13*(1/7) =", F(13,7), "= %.4f" % float(F(13,7)))
print("  B5 = 1 - S1 + S2 - S3 + S4 - S5 ; need B5 > 0, so S2 must carry the load.")
print("  equidistributed values: S_k = C(13,k)*(1/7)^k")
eq = {k: F(comb(13,k), 7**k) for k in range(1,6)}
for k in range(1,6): print(f"    S{k}^equid = {eq[k]} = {float(eq[k]):.4f}")
b5eq = 1 - eq[1] + eq[2] - eq[3] + eq[4] - eq[5]
print(f"  B5 at equidistribution = {float(b5eq):+.4f}  (positive: the level-5 wall clears)")
print()
print("MY CONTRIBUTION: THM-1025 gives an S2 LOWER bound valid in every regime.")
print("  S2 >= sum over all 78 pairs of floor_1025.  How close to equidistribution?")
random.seed(358)
rows = []
for trial in range(60):
    m0 = random.randint(20, 150)
    V = sorted(random.sample(range(m0, 13*m0), 13))
    s2_floor = sum(floor_1025(a,b) for a,b in itertools.combinations(V,2))
    s2_exact = sum(rho(a,b) for a,b in itertools.combinations(V,2))
    rows.append((float(s2_floor), float(s2_exact)))
    assert s2_floor <= s2_exact + F(1,10**9)
fl = sorted(r[0] for r in rows); ex = sorted(r[1] for r in rows)
print(f"  60 comparable 13-packets:")
print(f"    S2 floor  range [{fl[0]:.4f}, {fl[-1]:.4f}]   median {fl[len(fl)//2]:.4f}")
print(f"    S2 exact  range [{ex[0]:.4f}, {ex[-1]:.4f}]   median {ex[len(ex)//2]:.4f}")
print(f"    S2 equid  = {float(eq[2]):.4f}")
print(f"    floor recovers {100*fl[len(fl)//2]/float(eq[2]):.1f}% of the equidistributed S2")
print()
print("WHAT REMAINS for B5 > 0 (the honest ledger):")
print("  S2: LOWER bound  -- SUPPLIED by THM-1025 (all regimes, kernel-pure for coprime)")
print("  S4: LOWER bound  -- OPEN (the quadruple analogue of THM-1025)")
print("  S3: UPPER bound  -- OPEN (fleet's T2/THM-926 gives the triple deviation)")
print("  S5: UPPER bound  -- OPEN (level-5; fleet's B5 machinery)")
print()
print("  => LRC(14) is NOT closed by the pair floor. The pair floor closes the")
print("     S2 slot of a five-slot ledger. Claiming a completed proof would be false.")
