#!/usr/bin/env python3
"""
lrc14_BS_sector_majorant_macmini_0618s7.py   (mac-mini-2026-06-18-S7)  ANGLE A

Beurling-Selberg / Vaaler positive-definite majorant of the seven-sector cover S7(E).

GOAL: a RIGOROUS upper bound  meas(S7(E)) <= U_N(E)  that we can evaluate exactly, replacing
the heuristic corr <= C*W.  Want U_N(E) <= cap_k for the dangerous k=8..11 (AP-extremal).

SETUP.  E={0=e_1,e_2,...,e_k}.  Sector S_j=[j/7,(j+1)/7).  g_j(x)=prod_i (1-psi_j(e_i x))
= "sector j empty".  1_{S7}=prod_{j=1..6}(1-g_j) (sector 0 auto-hit since 0 in E).
Inclusion-exclusion over missed-sector set T subset {1..6}:
   1_{S7}(x) = sum_{T}(-1)^{|T|} prod_{i=1}^k 1_{B_T}(e_i x),    B_T = T^c sectors, |B_T|=1-|T|/7.

THE CERTIFICATE.  Integrate term by term over x in [0,1).  For a fixed T, the integrand is a
product of k periodic indicators of the SAME set B_T evaluated at e_i x.  We REPLACE 1_{B_T} by
its degree-N Vaaler majorant/minorant.  Vaaler (1985): for any measurable B that is a finite
union of arcs, deg-N trig polys  V_B^- <= 1_B <= V_B^+  with
   hat V^pm(0) = |B| +- 1/(N+1),     hat V^pm(n) = hat 1_B(n) + eps^pm(n)/(N+1) * (Fejer kernel),
   |hat V^pm(n) - hat 1_B(n)| <= 1/(N+1)   for all n,  support |n|<=N.
We use the SELBERG construction for a union of arcs:  V_B^pm = sum over component arcs of the
single-arc Vaaler functions, but the SHARP single-pass version majorizes the union directly via
B_T = complement of (union of |T| arcs each length 1/7).  We build V^pm for the UNION directly.

KEY ALGEBRAIC FACT (orbit integral collapses).  Because 1_{B_T}(e_i x) is plugged with multiplier
e_i and integrated against dx, and the Vaaler poly has Fourier support |n|<=N:
   int_0^1 prod_i V(e_i x) dx = sum over (n_1,...,n_k), |n_i|<=N, sum n_i e_i = 0  of  prod hat V(n_i).
This is EXACTLY the offset-relation-lattice sum, now with EXPLICIT, BAND-LIMITED coefficients.
The n=0 term = (|B_T|+-1/(N+1))^k -> reproduces M7 + an O(1/N) main-term inflation; the n!=0 terms
are the relation-lattice tail, but each |hat V(n)| is now explicitly bounded.

This script:
 (1) Builds Vaaler majorant/minorant Fourier coeffs for B_T (union of arcs) EXACTLY (Fraction +
     the half-integer-shift cotangent weights as exact algebraic numbers via the Beurling kernel).
     [We use the cleanest exact route: the "Selberg box" coeffs hat V^pm(n) for an arc, summed.]
 (2) Computes the certificate U_N(E) = sum_T (-1)^{|T|} * [ upper or lower bound on the T-term ]
     where for sign(+) terms we use V^+ on each factor (giving an upper bound on a nonneg product)
     and for sign(-) terms we use V^- (lower bound), so the signed sum is an UPPER bound on meas(S7).
 (3) Compares U_N(E) to cap_k for consec and several shapes, over a range of degrees N.

We FIRST validate the relation-lattice orbit-integral identity exactly against the geometric
measS7 (breakpoint) computation, by using the EXACT indicator (N=inf surrogate: full Fourier of
1_{B_T} truncated) -- a sanity gate before trusting the bound.
"""
import sys, itertools, math
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

# ---------- exact geometric meas(S7) (ground truth) ----------
def measS7_geom(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set()
        for e in E:
            y=(e*xm)%1; secs.add(int(y*7))
        if len(secs)==7: total+=x1-x0
    return total

def M7(k):
    s=F(0)
    for t in range(0,7):
        s += F((-1)**t * math.comb(6,t)) * F(7-t,7)**(k-1)
    return s

cap = {8:F(2243,5880), 9:None, 10:None, 11:None, 12:None, 13:F(1,1)}
# cap_k values from canon table (exact where known; floats elsewhere as guide)
cap_float = {8:0.38146,9:0.4943,10:0.6044,11:0.7253,12:0.8571,13:1.0}

# ---------------------------------------------------------------------------
# PART 1.  The exact orbit-integral identity for the TRUE indicator (sanity gate).
# For each T subset {1..6}, B_T is a union of arcs. The Fourier coeff of 1_{B_T} at integer n:
#   hat 1_{B_T}(n) = integral_{B_T} e(-n y) dy.
# For B_T = complement of union_{j in T} [j/7,(j+1)/7):
#   hat 1_{B_T}(0) = 1 - |T|/7
#   hat 1_{B_T}(n) = -sum_{j in T} hat 1_{S_j}(n),  n!=0,   where
#   hat 1_{S_j}(n) = integral_{j/7}^{(j+1)/7} e(-n y) dy = e(-n j/7)*(1-e(-n/7))/(2 pi i n).
# The orbit integral  int prod_i 1_{B_T}(e_i x) dx = sum_{n: sum n_i e_i =0} prod_i hat 1_{B_T}(n_i).
# This is an infinite sum; we truncate at |n_i|<=NT to sanity-check convergence to measS7.
# (Exact arithmetic on the truncated sum uses complex Fractions via (real,imag) rational pairs is
#  heavy; for the SANITY GATE we use high-precision floats. The RIGOROUS bound in Part 2 is the
#  point; this gate just confirms we set up the relation lattice correctly.)
# ---------------------------------------------------------------------------
import cmath
def hat_indicatorBT_float(T, n):
    # Fourier coeff (e(-n y) convention) of 1_{B_T}
    if n==0:
        return 1.0 - len(T)/7.0
    s=0j
    for j in T:
        # hat 1_{S_j}(n)
        if n==0:
            s += 1.0/7.0
        else:
            s += cmath.exp(-2j*math.pi*n*j/7.0)*(1-cmath.exp(-2j*math.pi*n/7.0))/(-2j*math.pi*n)
    return -s

def orbit_integral_trunc(E, T, NT):
    # sum over (n_2,...,n_k) in [-NT,NT], n_1 determined by sum n_i e_i = 0 (e_1=0 -> n_1 free!)
    # NOTE e_1=0 so the constraint is sum_{i>=2} n_i e_i = 0 and n_1 is unconstrained -> divergent
    #   UNLESS hat at n_1: but factor for i=1 is hat 1_{B_T}(n_1) and e_1=0 means that variable's
    #   exponential e(n_1 e_1 x)=1 always, so effectively the i=1 factor integrates to its mean:
    #   the i=1 factor contributes hat 1_{B_T}(0) = |B_T| (only n_1=0 survives the x-integral since
    #   all other factors carry x; precisely n_1 must satisfy n_1*0 + sum_{i>=2} n_i e_i =0, the 0
    #   coeff makes n_1 unconstrained, but the x-integral kills any leftover frequency -> the i=1
    #   factor is just the constant |B_T| because its argument e_1 x=0 is constant: 1_{B_T}(0)=
    #   [0 in B_T] which is 1 (sector 0 in B_T since 0 not in T). So factor_1 = 1, NOT |B_T|.)
    # Correct: 1_{B_T}(e_1 x)=1_{B_T}(0)=1 (0 in B_T always). So the i=1 factor is the constant 1.
    Enz=[e for e in E if e!=0]
    kk=len(Enz)
    tot=0j
    for ns in itertools.product(range(-NT,NT+1), repeat=kk):
        if sum(n*e for n,e in zip(ns,Enz))!=0: continue
        p=1.0+0j
        for n in ns: p*=hat_indicatorBT_float(T,n)
        tot+=p
    return tot  # * factor_1 = *1

def measS7_lattice_float(E, NT):
    s=0.0
    for r in range(0,7):
        for T in itertools.combinations(range(1,7), r):
            s += ((-1)**r) * orbit_integral_trunc(E, T, NT).real
    return s

if __name__=="__main__":
    print("="*78)
    print("PART 1 sanity gate: relation-lattice orbit integral -> measS7 (truncated)")
    print("="*78)
    for name,E in [("consec k=7 {0..6}",list(range(7))), ("k=8 consec {0..7}",list(range(8)))]:
        g=float(measS7_geom(E))
        print(f"\n{name}: geometric measS7 = {g:.6f}, M7={float(M7(len(E))):.6f}")
        for NT in [6,10,14]:
            v=measS7_lattice_float(E,NT)
            print(f"   |n|<={NT:>3}: lattice = {v:.6f}   (err {abs(v-g):.2e})")
    print("\nIf lattice -> geometric as NT grows, the orbit-relation setup is CORRECT.")
