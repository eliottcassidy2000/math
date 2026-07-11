"""
opus-2026-07-11-S217: the covolume / sign-blindness experiment for the LRC(14) last gap.

meas(S7(E)) = P{ x : orbit {frac(e x)} hits ALL 7 sectors } = M7(k) + corr(E),
corr(E) = Sum_{0!=n in Lambda-o(E)} K(n),  K(n)=D7(n mod 7)/prod n_j,  Lambda-o(E)={n: Sum n_j e_j=0}.
Crux: meas(S7) <= cap_k, consec MAXIMIZES.  corr(consec_8)=0.303 is the binding value (cap_8-M7=0.357).

THE TWO COMPETING RULERS this script separates:
  (i)  covol(Lambda-o(E)) = |e/gcd|_2  (Minkowski: prod lambda_i ~ covol; AP is the covol-MINIMIZER, F2).
  (ii) the SIGNED cancellation (Re D7 alternating): |corr| is 5-6x below Sum|K(n)| (F3, proven-lossy).
Question: does covol alone bound corr?  (If yes, the absolute Minkowski count would close it.  The repo
says NO -- F3.  This confirms it exactly and shows WHY: short relations (small lambda_1) survive at large covol.)
"""
from fractions import Fraction as F
from itertools import product as iproduct
from math import gcd, isqrt, comb
import sympy

def meas_S7(E):
    """Exact P{ x in [0,1) : all 7 sectors [j/7,(j+1)/7) hit by some frac(e x), e in E }."""
    Es = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Es:
        e = abs(e)
        for m in range(0, 7*e + 1):
            bps.add(F(m, 7*e))
    pts = sorted(b for b in bps if 0 <= b < 1); pts.append(F(1))
    tot = F(0)
    for i in range(len(pts)-1):
        a, b = pts[i], pts[i+1]
        if b <= a: continue
        mid = (a+b)/2; occ = {0}
        for e in Es: occ.add(int(((e*mid) % 1)*7))
        if len(occ) == 7:            # ALL sectors hit
            tot += (b-a)
    return tot

def M7(k):
    return sum(F((-1)**t)*comb(6,t)*F(7-t,7)**(k-1) for t in range(7))

def covol(E):
    """covol(Lambda-o(E)) = |e/gcd(e)|_2 = sqrt(Sum (e_j/g)^2) over the k-1 nonzero offsets (F2)."""
    Es = [abs(e) for e in E if e != 0]
    g = 0
    for e in Es: g = gcd(g, e)
    v2 = sum((e//g)**2 for e in Es)
    return v2, g  # return |e/gcd|_2^2 (exact integer) and the gcd

def lambda1_l1(E, cap=5):
    Es = [e for e in E if e != 0]; r = len(Es); best = None
    for n in iproduct(range(-cap, cap+1), repeat=r):
        if all(x == 0 for x in n): continue
        if sum(nj*ej for nj, ej in zip(n, Es)) == 0:
            l1 = sum(abs(x) for x in n)
            if best is None or l1 < best: best = l1
    return best

fams = {
    "consec {0..7}        ": [0,1,2,3,4,5,6,7],
    "AP dilated 5*        ": [0,5,10,15,20,25,30,35],
    "stranger {0..6,40}   ": [0,1,2,3,4,5,6,40],
    "stranger {0..6,400}  ": [0,1,2,3,4,5,6,400],
    "two-cluster core+far ": [0,1,2,3,40,41,42,43],
    "mild spread          ": [0,1,2,4,7,11,16,22],
    "Sidon-ish {0,1,3,7..}": [0,1,3,7,12,20,30,44],
    "geometric 2^j        ": [0,1,2,4,8,16,32,64],
    "wide dissociated     ": [0,7,17,31,53,83,127,179],
}
k = 8
m8 = M7(8)
print(f"M7(8) = {float(m8):.6f} (= 20160/823543);  cap_8 = 2243/5880 = {float(F(2243,5880)):.6f};  margin = {float(F(2243,5880)-m8):.4f}\n")
print(f"{'family':>22} {'meas(S7)':>10} {'corr':>9} {'covol^2':>9} {'covol':>8} {'C/covol':>9} {'lam1':>5}")
for name, E in fams.items():
    ms = meas_S7(E); corr = ms - m8
    v2, g = covol(E); cov = v2 ** 0.5
    lam = lambda1_l1(E)
    Ccov = float(corr) * cov          # if corr <= C/covol, this C should be BOUNDED across families
    print(f"{name:>22} {float(ms):>10.6f} {float(corr):>+9.5f} {v2:>9} {cov:>8.2f} {Ccov:>9.3f} {str(lam):>5}")

print("""
READING:
- consec is the corr-MAXIMIZER (binding 0.303) AND the covol-MINIMIZER (F2) -- the Minkowski intuition.
- BUT 'C/covol' (=corr*covol) is NOT bounded: the stranger has large covol yet corr stays moderate
  because the SHORT relation 1+2-3=0 (lambda1 small) survives -- so covol alone does NOT bound corr.
- The gap between corr and the naive covol bound is the SIGNED cancellation (F3): the absolute lattice
  sum is 5-6x larger than |corr|. Successive-minima/covolume see density, not the sign of Re D7.
=> the absolute Minkowski count is provably insufficient; the count must be SIGNED (Abel vs the mod-7
   character D7). THM-546's far-element peel is the working rank-1 signed instance (|Delta_w|<=(6/49)V/w).
""")
