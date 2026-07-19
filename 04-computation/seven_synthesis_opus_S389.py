# opus-2026-07-17-S389 -- HOW THE REPO TREATS 7, AND THE PARITY EXTENSION.
#
# ARCHAEOLOGY FIRST.  There are TWO different sevens in this corpus:
#   (A) the TOURNAMENT seven -- Paley tournament on 7 vertices, QR_7, H7, the
#       heptagon (620+ mentions).  Seven is a PRIME = 3 mod 4 there.
#   (B) the LRC seven -- 2*lam = 1/7, S1 = 13/7, hhat zeros at 7|n, the k=7
#       arity ceiling.  Seven is n/2 = 14/2 there.
# They are NOT the same constant, and this run tests whether the LRC seven is a
# parity artifact -- which would explain why LRC(14)'s coincidences are EXACT.
#
# ALSO: THM-1170 (my S383) rediscovered HYP-2059 (my S557) -- the pinch lemma,
# that the optimum sits at a PAIR-SUM time.  Part (1) tests whether my extra
# claim (differences and half-period peaks also win) is real or a co-occurrence
# artifact, which decides whether THM-1170 adds anything to HYP-2059.
from fractions import Fraction as F
from itertools import combinations
import random
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)

print("(1) DOES THM-1170 ADD ANYTHING TO HYP-2059?  Is the winner ALWAYS a pair SUM?")
print("    (HYP-2059/S557: optimum attained at t = m/(v_a+v_b).  My S383 reported")
print("     sum 55 / diff 20 / peak 11 -- but types can share a denominator.)")
random.seed(389)
only_diff=0; only_peak=0; has_sum=0; n=0
for _ in range(25):
    V=sorted(random.sample(range(1,60),13))
    sums={a+b for a,b in combinations(V,2)}
    diffs={abs(a-b) for a,b in combinations(V,2) if a!=b}
    peaks={2*v for v in V}
    cand=sorted((sums|diffs|peaks)-{0})
    best=F(0); bq=None
    for q in cand:
        for p in range(1,q):
            g=gap_at(V,F(p,q))
            if g>best: best,bq=g,q
    n+=1
    if bq in sums: has_sum+=1
    elif bq in diffs: only_diff+=1
    else: only_peak+=1
print(f"    {n} families: winner is a pair SUM in {has_sum};"
      f" difference-only {only_diff}; peak-only {only_peak}")
print("    => if difference/peak-only is 0, HYP-2059 already covers it and my")
print("       'sum 55 / diff 20 / peak 11' was counting CO-OCCURRENCES, not winners.")

print()
print("(2) THE PARITY EXTENSION: the LRC seven is n/2, an INTEGER only for even n")
print("    hhat(m) = sin(2*pi*m/n)/(pi*m) vanishes iff n | 2m, i.e. m = 0 mod n/gcd(2,n)")
print("      n even -> zeros every n/2 steps (DENSE)")
print("      n odd  -> zeros every n steps   (half as dense)")
from math import gcd, sin, pi
for nn in [12,13,14,15,16,17]:
    step = nn//gcd(2,nn)
    half = F(nn,2) if nn%2==0 else None
    print(f"    n={nn:3d}: hhat zeros every {step:3d} steps;  n/2 = "
          f"{'%d (INTEGER)'%(nn//2) if nn%2==0 else '%.1f (not an integer)'%(nn/2)}"
          f"   arity ceiling k < n/2 {'is exact' if nn%2==0 else 'never degenerates exactly'}")

print()
print("(3) THE TENT/INDEPENDENCE COINCIDENCE ACROSS n  (THM-1195 was n=14)")
print("    tent threshold = lam/2 = 1/(2n); independence E[min] = (1/2)/(n_speeds+1)")
print("    with n_speeds = n-1 these agree iff 1/(2n) = 1/(2n) -- ALWAYS.")
for nn in [10,12,13,14,20]:
    tent=F(1,2*nn); indep=F(1,2)/((nn-1)+1)
    print(f"    n={nn:3d}: tent {tent}  independence {indep}   {'EQUAL' if tent==indep else 'differ'}")
print("    => the coincidence THM-1195 found at n=14 is NOT special to 14: it holds")
print("       for every n, and is the statement that LRC's floor sits exactly where")
print("       the geometric and probabilistic estimates cross.")
