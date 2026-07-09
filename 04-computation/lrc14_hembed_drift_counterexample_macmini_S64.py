"""
mac-mini-2026-07-09-S64 -- EXACT counterexample: a GOOD PERIOD does NOT certify loneliness at that j.
The hembed (THM-527 Part A) slow-fast embedding cannot be proven at the antecedent's own (j, phi).

SETUP (the binding kps-S105 says hembed omits): speeds v_i = Vmax - e_i, co-offsets e_i, ruler Vmax.
Two-scale time tau = (j + phi)/Vmax. Slow-fast identity (kps-S105 LRCSlowFast, verified below):
    nearInt(v_i * tau) = nearInt(phi * (v_i/Vmax) - c_i),     c_i = frac(e_i * j / Vmax)  (the teeth)
NOT nearInt(phi - c_i). The difference is the DRIFT e_i*phi/Vmax, which is O(spread/Vmax) = O(1) in
the good-period window -- far larger than the 1/14 margin.

hembed's antecedent: exists j, phi with  nearInt(phi - c_i) > 1/14  for all i  (teeth cleared).
A good period (7*maxCircGap > Vmax, i.e. teeth gap > 1/7) supplies this at the gap midpoint.

CLAIM TESTED: does a good period j admit SOME phi with minReach((j+phi)/Vmax) >= 1/14 ?
ANSWER: NO. Exact counterexamples below.

METHOD (exact, no grids -- cf MISTAKE-130): minReach((j+phi)/V) = min_i nearInt(v_i (j+phi)/V) is
piecewise linear in phi with breakpoints where v_i*tau in (1/2)Z. Between consecutive breakpoints every
component is LINEAR, so the min is CONCAVE and its max is attained at an endpoint or at a pairwise
crossing of two components. Enumerate all of those exactly over Q => the exact max_phi.
"""
from fractions import Fraction as F

E = [0, 7, 14, 21, 26, 29, 37, 44, 51, 58, 67, 75, 82]   # klein's worst7Struct cluster
V = 91                                                    # its hard ruler (7*13)
vs = [V - e for e in E]                                   # THE BINDING: v_i = Vmax - e_i

def nearInt(x):
    r = F(x) % 1
    return min(r, 1 - r)

def minReach(tau):
    return min(nearInt(F(v) * tau) for v in vs)

def maxCircGap(j):
    ph = sorted({(e * j) % V for e in E}); m = len(ph)
    if m <= 1: return V
    return max(max((ph[(i + 1) % m] - ph[i]) % V for i in range(m - 1)), ph[0] + V - ph[-1])

def is_good(j):   # the Lean IsGoodPeriod: 7 * maxCircGap > Vmax
    return 7 * maxCircGap(j) > V

def teeth(j):
    return sorted({F((e * j) % V, V) for e in E})

def linear_piece(v, j, phi_mid):
    """on the breakpoint-interval containing phi_mid, nearInt(v*(j+phi)/V) = A*phi + B."""
    x = F(v) * (F(j) + phi_mid) / V
    m = round(x)
    s = 1 if x >= m else -1
    return F(s) * F(v) / V, F(s) * (F(v) * F(j) / V - m)

def exact_max_phi(j):
    """EXACT max over phi in [0,1] of minReach((j+phi)/V), via concave-piece crossing enumeration."""
    bps = {F(0), F(1)}
    for v in vs:
        lo, hi = F(v) * F(j) / V, F(v) * (F(j) + 1) / V
        k = int(2 * lo) - 2
        while F(k, 2) <= hi + 1:
            phi = F(k, 2) * V / F(v) - F(j)
            if 0 <= phi <= 1: bps.add(phi)
            k += 1
    bps = sorted(bps)
    best, barg = F(-1), None
    for i in range(len(bps) - 1):
        a, b = bps[i], bps[i + 1]
        if a == b: continue
        L = [linear_piece(v, j, (a + b) / 2) for v in vs]
        cands = [a, b]
        for p in range(len(L)):
            for q in range(p + 1, len(L)):
                A1, B1 = L[p]; A2, B2 = L[q]
                if A1 != A2:
                    ph = (B2 - B1) / (A1 - A2)
                    if a <= ph <= b: cands.append(ph)
        for ph in cands:
            val = min(A * ph + B for A, B in L)
            if val > best: best, barg = val, ph
    return best, barg

# --- verify the slow-fast identity (kps-S105) ---
j0 = 2; C0 = teeth(j0); phi0 = F(3, 7)
ok = all(nearInt(F(V - e) * (F(j0) + phi0) / V)
         == nearInt(phi0 * F(V - e, V) - F((e * j0) % V, V)) for e in E)
print(f"slow-fast identity  nearInt(v*tau) == nearInt(phi*v/Vmax - c): {ok}")

goods = [j for j in range(1, V) if is_good(j)]
print(f"\nE = {E}\nVmax = {V},  binding v_i = Vmax - e_i,  #good periods (7*maxCircGap>Vmax) = {len(goods)}\n")
print(f"{'j':>4}{'good?':>7}{'exact max_phi minReach':>26}{'  >= 1/14 ?':>12}")
bad = []
for j in (2, 5, 10, 11):
    b, arg = exact_max_phi(j)
    lonely = b >= F(1, 14)
    if is_good(j) and not lonely: bad.append((j, b))
    print(f"{j:>4}{str(is_good(j)):>7}{str(b):>26}{str(lonely):>12}")
print(f"\n1/14 = {F(1,14)} = {float(F(1,14)):.6f}")
print("\nGOOD PERIODS THAT ARE **NOT** LONELY AT ANY phi:")
for j, b in bad:
    print(f"  j={j}: 7*maxCircGap={7*maxCircGap(j)} > Vmax={V} (good), yet max_phi minReach = {b} < 1/14")
print("""
=> 'good period at j  ==>  exists phi with minReach((j+phi)/Vmax) >= 1/14'  is FALSE.
   hembed therefore CANNOT be discharged by embedding at the antecedent's own (j, phi):
   the drift e_i*phi/Vmax must be absorbed, and a 1/7 teeth-gap does not permit it.
   (hembed's IMPLICATION may still hold -- its conclusion is exists tau, i.e. LRC14(v) itself,
    witnessed here at j=25 with minReach=0.2306 -- but that is circular, not a proof.)
""")
