"""
lrc_farey_roof_apmaxgap_opus_S134.py

THE FAREY ROOF: an exact closed-form mechanism for the AP's max-gap function.

CLAIM 1 (origin-gap saturation, sharpens klein-S153's numeric 'gap@0=maxgap a.e.'):
  For the AP config C(x) = {frac(i*x) : 1<=i<=k}, x irrational, the circular gap
  containing the point 0 IS the maximum gap, ALWAYS.
  Proof sketch (three-distance): in the with-0 config {frac(i*x): 0<=i<=k} the two
  gaps flanking 0 have lengths delta1=||q_L x||, delta2=||q_R x|| (the one-sided
  best-approximation denominators q_L, q_R <= k), and by the three-distance theorem
  every gap length lies in {delta1, delta2, delta1+delta2}. Removing the point 0
  merges the two flanking gaps into one of length delta1+delta2 = the maximal value.

CLAIM 2 (the Farey roof): on each Farey-order-k cell (p/q, p'/q') the AP max-gap is
  EXACTLY the linear function
      maxgap(x) = q*(x - p/q) + q'*(p'/q' - x),
  i.e. it interpolates 1/q (left endpoint value) -> 1/q' (right endpoint value).
  [m+(x) = min_i frac(ix) = q*(x-p/q) and m-(x) = min_i (1-frac(ix)) = q'*(p'/q'-x)
   are the classical 'first return' facts on a Farey cell; maxgap = m+ + m- by CLAIM 1.]

CONSEQUENCES (exact, all k):
  E[maxgap(AP_k)]  = sum over consecutive Farey-k pairs  1/(q*q'^2)
                   = (1/2) * sum (q+q')/(q^2 q'^2)          (symmetric form)
  E[min_i frac(ix)] = (1/2) * E[maxgap(AP_k)]               (reflection symmetry)
  mu_theta(AP_k)   = sum over cells of meas{roof > theta}   (per-cell linear tail)

  This turns the load-bearing (A') RHS values mu_{1/7}(AP_k) -- 691/735, 247/294,
  38/49, 1381/2205, 13823/24255, 477/1078 -- into transparent Farey statistics.

Verification below:
  (a) roof formula vs direct maxgap at random x (k=3..13)
  (b) gap@0 == maxgap identity at random x (k=3..13)
  (c) E[maxgap] Farey sum vs known exact values (k=8..13) and vs numeric integral
  (d) mu_{1/7}(AP_k) from the roof vs the canon table (k=8..13); also mu_{2/7}
  (e) the pigeonhole boundary remark: iid mean maxgap H_k/k vs threshold 2/(k+1)
"""
from fractions import Fraction
import random

def farey(k):
    """Ascending list of Fractions in [0,1] with denominator <= k."""
    fr = set()
    for q in range(1, k+1):
        for p in range(0, q+1):
            fr.add(Fraction(p, q))
    return sorted(fr)

def maxgap_direct(E, x):
    """Circular max gap of {frac(e*x): e in E} (floats)."""
    ph = sorted((e*x) % 1.0 for e in E)
    g = [ph[i+1]-ph[i] for i in range(len(ph)-1)] + [ph[0]+1.0-ph[-1]]
    return max(g)

def gap_at_zero(E, x):
    """Length of the circular gap containing the point 0 (= min + 1 - max)."""
    ph = sorted((e*x) % 1.0 for e in E)
    return ph[0] + 1.0 - ph[-1]

def roof_value(k, x, cells):
    """Farey-roof prediction at x (cells = list of (a, b) consecutive Farey fracs)."""
    # binary search for the cell containing x
    lo, hi = 0, len(cells)-1
    while lo < hi:
        mid = (lo+hi)//2
        if float(cells[mid][1]) < x: lo = mid+1
        else: hi = mid
    a, b = cells[lo]
    q, qp = a.denominator, b.denominator
    return q*(x - float(a)) + qp*(float(b) - x)

print("="*84)
print("(a)+(b) Roof formula and origin-gap saturation: direct check at random x")
print("="*84)
rng = random.Random(20260707)
for k in range(3, 14):
    F = farey(k)
    cells = list(zip(F[:-1], F[1:]))
    E = list(range(1, k+1))
    worst_roof, worst_g0 = 0.0, 0.0
    for _ in range(4000):
        x = rng.random()
        mg = maxgap_direct(E, x)
        worst_roof = max(worst_roof, abs(mg - roof_value(k, x, cells)))
        worst_g0 = max(worst_g0, abs(mg - gap_at_zero(E, x)))
    print(f"  k={k:2d}: max|maxgap-roof| = {worst_roof:.3e}   max|maxgap-gap@0| = {worst_g0:.3e}")

print()
print("="*84)
print("(c) E[maxgap(AP_k)] = sum 1/(q q'^2) over consecutive Farey-k pairs (EXACT)")
print("="*84)
known = {8: Fraction(43,140), 9: Fraction(47,168), 10: Fraction(19,72),
         11: Fraction(151,630), 12: Fraction(796,3465), 13: Fraction(93,440)}
for k in range(1, 14):
    F = farey(k)
    S = Fraction(0)
    Ssym = Fraction(0)
    for a, b in zip(F[:-1], F[1:]):
        q, qp = a.denominator, b.denominator
        S += Fraction(1, q*qp*qp)
        Ssym += Fraction(q+qp, 2*q*q*qp*qp)
    # numeric integral for independent confirmation
    E = list(range(1, k+1))
    N = 20000
    num = sum(maxgap_direct(E, (j+0.5)/N) for j in range(N))/N
    tag = ""
    if k in known:
        tag = f"  known={known[k]}  match={'YES' if S==known[k] else '*** NO ***'}"
    sym_ok = "sym=OK" if S == Ssym else "sym=*** MISMATCH ***"
    print(f"  k={k:2d}: Farey sum = {str(S):>16} = {float(S):.6f}   numeric={num:.6f}  {sym_ok}{tag}")
print("  => E[min_i frac(ix)] = half of the above; e.g. k=13: 93/880")

print()
print("="*84)
print("(d) mu_theta(AP_k) from the roof (EXACT per-cell linear tail)")
print("="*84)
def mu_from_roof(k, theta):
    """meas{x: roof_k(x) > theta} exactly (Fraction theta)."""
    F = farey(k)
    tot = Fraction(0)
    for a, b in zip(F[:-1], F[1:]):
        q, qp = a.denominator, b.denominator
        L = b - a                      # cell length = 1/(q q')
        vl, vr = Fraction(1, q), Fraction(1, qp)   # roof endpoint values
        if vl <= theta and vr <= theta:
            continue
        if vl > theta and vr > theta:
            tot += L
            continue
        # linear from vl to vr across length L: crossing point
        # v(t) = vl + (vr-vl) * t/L  > theta
        t_star = (theta - vl) * L / (vr - vl)
        if vl > theta:   # good part is [0, t_star)
            tot += t_star
        else:            # good part is (t_star, L]
            tot += L - t_star
    return tot

known_mu17 = {8: Fraction(691,735), 9: Fraction(247,294), 10: Fraction(38,49),
              11: Fraction(1381,2205), 12: Fraction(13823,24255), 13: Fraction(477,1078)}
th = Fraction(1,7)
for k in range(8, 14):
    m = mu_from_roof(k, th)
    ok = "match=YES" if m == known_mu17[k] else f"*** MISMATCH vs {known_mu17[k]} ***"
    print(f"  k={k:2d}: mu_1/7(AP) = {str(m):>13} = {float(m):.6f}   {ok}")
print("  bonus mu_2/7 (the refuted-floor threshold), exact from the roof:")
for k in (12, 13):
    m = mu_from_roof(k, Fraction(2,7))
    print(f"    k={k}: mu_2/7(AP) = {m} = {float(m):.6f}")

print()
print("="*84)
print("(e) why 1/7 works at every k: threshold 2/(k+1) vs pigeonhole 1/k vs iid H_k/k")
print("="*84)
import math
for k in range(4, 14):
    Hk = sum(1.0/j for j in range(1, k+1))
    ap = None
    F = farey(k); S = Fraction(0)
    for a,b in zip(F[:-1],F[1:]):
        S += Fraction(1, a.denominator*b.denominator**2)
    print(f"  k={k:2d}: 1/k={1.0/k:.4f}  2/(k+1)={2.0/(k+1):.4f}  AP mean={float(S):.4f}"
          f"  iid mean=H_k/k={Hk/k:.4f}  [pigeonhole {'covers' if 1.0/k>=2.0/(k+1) or k<=7 else 'stops'};"
          f" note e^2={math.e**2:.2f} sits at the k=7/8 boundary]")
