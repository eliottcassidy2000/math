#!/usr/bin/env python3
"""
lrc14_vaaler_kernel_macmini_0618s7.py  (mac-mini-2026-06-18-S7)

Self-contained, COMPUTATIONALLY-VERIFIED Vaaler majorant/minorant of an arc indicator on R/Z.

Vaaler's sawtooth approximation (Vaaler 1985, Montgomery "Ten Lectures" Thm 1.5):
The periodic sawtooth  psi(x) = {x} - 1/2  (x not integer), psi(0)=0, has trig-poly
majorant/minorant of degree N:
   psi^+(x) =  V_N(x) + (1/(2(N+1))) * Fejer_N(x)
   psi^-(x) =  V_N(x) - (1/(2(N+1))) * Fejer_N(x)
where  V_N(x) = sum_{1<=|n|<=N}  (i/(2 pi n)) * (1 - |n|/(N+1)) * (-?) ...
We DO NOT trust the memorized constants; instead we define the standard Vaaler V_N via its
Fourier coefficients and VERIFY psi^- <= psi <= psi^+ on a fine exact grid + breakpoint check.

Then the arc majorant:  1_{[a,b]}(x) = (b-a) + psi(x-a) - psi(x-b)   (a.e.), so
   M^+_{[a,b]}(x) := (b-a) + psi^+(x-a) - psi^-(x-b)  >= 1_{[a,b]}(x),  deg N.
   M^-_{[a,b]}(x) := (b-a) + psi^-(x-a) - psi^+(x-b)  <= 1_{[a,b]}(x),  deg N.
Constant term = (b-a) + 1/(N+1)  (resp - 1/(N+1)).  Coeffs band-limited |n|<=N.

Vaaler's explicit coefficients (the version we verify):
  hat{psi}(n) = -1/(2 pi i n)   (n!=0),  hat{psi}(0)=0.
  hat{psi^+}(n) = hat{psi}(n)*(1-|n|/(N+1)) + (1/(2(N+1)))   for 1<=|n|<=N   [Fejer part real, even]
                  wait: Fejer_N has hat(n)=(1-|n|/(N+1)); times 1/(2(N+1)).
  hat{psi^+}(0) = 1/(2(N+1)) * 1 = 1/(2(N+1)).
  hat{psi^-}(n) = hat{psi}(n)*(1-|n|/(N+1)) - (1/(2(N+1))) * (1-|n|/(N+1)) ...
We just BUILD candidates and VERIFY numerically; if a candidate fails majorization we report it
rather than trusting it.  The ONLY thing we rely on downstream is a VERIFIED majorant.
"""
import sys, cmath, math
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
TPI=2*math.pi

# ---- Vaaler sawtooth majorant/minorant coefficients (complex), degree N ----
# psi(x) = sum_{n!=0} c_psi(n) e(nx),  c_psi(n) = -1/(2 pi i n) = i/(2 pi n).
def c_psi(n):
    if n==0: return 0j
    return 1j/(TPI*n)

# Vaaler's majorant of psi:  psi^+(x).  Coeffs (standard form, Montgomery Ten Lectures eq 1.??):
#   c+(0) = 1/(2(N+1))
#   c+(n) = (1 - |n|/(N+1)) * c_psi(n)  +  (1/(2(N+1))) * (sign adjust) ...
# We implement the WELL-KNOWN closed form:
#   psi^+(x) = -sum_{n=1}^{N} (1/(pi n)) (1-n/(N+1)) sin(2 pi n x)  + (1/(2(N+1))) sum_{n=-N}^{N}(1-|n|/(N+1)) e(2pi i n x)
# i.e. c+(n) = c_psi(n)(1-|n|/(N+1)) + (1/(2(N+1)))(1-|n|/(N+1)).
def c_plus(n, N):
    if abs(n)>N: return 0j
    fej = (1 - abs(n)/(N+1))   # Fejer weight
    return c_psi(n)*fej + (1.0/(2*(N+1)))*fej
def c_minus(n, N):
    if abs(n)>N: return 0j
    fej = (1 - abs(n)/(N+1))
    return c_psi(n)*fej - (1.0/(2*(N+1)))*fej

def eval_trig(coeffs_fn, x, N):
    s=0j
    for n in range(-N,N+1):
        s+=coeffs_fn(n)*cmath.exp(1j*TPI*n*x)
    return s.real

def psi(x):
    fx=x-math.floor(x)
    return fx-0.5

# ---- VERIFY psi^- <= psi <= psi^+ on dense grid (away from the jump at integers) ----
def verify_sawtooth(N, NG=4000):
    worst_up=1e9; worst_lo=1e9
    for j in range(NG):
        x=(j+0.5)/NG
        p=psi(x); pp=eval_trig(lambda n: c_plus(n,N), x, N); pm=eval_trig(lambda n: c_minus(n,N), x, N)
        worst_up=min(worst_up, pp-p)   # want >=0
        worst_lo=min(worst_lo, p-pm)   # want >=0
    return worst_up, worst_lo

if __name__=="__main__":
    print("Verify Vaaler sawtooth majorant/minorant  psi^- <= psi <= psi^+  (min slack over grid):")
    print(f"{'N':>4}{'min(psi+ - psi)':>18}{'min(psi - psi-)':>18}{'c+(0)=int psi+':>16}{'1/(2(N+1))':>12}")
    for N in [4,8,16,32,64]:
        wu,wl=verify_sawtooth(N)
        c0=c_plus(0,N).real
        print(f"{N:>4}{wu:>18.6f}{wl:>18.6f}{c0:>16.6f}{1/(2*(N+1)):>12.6f}")
    print("\nIf both min-slacks are >= 0 (up to grid resolution), the Vaaler kernel is valid.")
    print("Then arc majorant M+_[a,b] = (b-a)+psi+(x-a)-psi-(x-b), int = (b-a)+1/(N+1).")
    # sanity: arc majorant of [0.2,0.5], check >= indicator and integral
    a,b=0.2,0.5; N=16
    def Mplus(x): return (b-a)+eval_trig(lambda n:c_plus(n,N),x-a,N)-eval_trig(lambda n:c_minus(n,N),x-b,N)
    NG=4000; worst=1e9; integ=0.0
    for j in range(NG):
        x=(j+0.5)/NG; ind=1.0 if a<=x<b else 0.0
        worst=min(worst, Mplus(x)-ind); integ+=Mplus(x)/NG
    print(f"\nArc [0.2,0.5] deg {N}: min(M+ - 1_arc)={worst:.6f} (want>=0), int M+ = {integ:.6f}, target (b-a)+1/(N+1)={b-a+1/(N+1):.6f}")
