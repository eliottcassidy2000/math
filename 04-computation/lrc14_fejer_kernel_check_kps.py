import sys
from fractions import Fraction
from math import sin, pi
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
# Verify the CORE Fejer/sine-kernel insight: chat_T(n) for a union of |T| sectors of width 1/7.
# 1_A(y) = sum_{j in T} 1_{[j/7,(j+1)/7)}(y). Fourier: chat_T(n) = sum_{j in T} e^{-2pi i n (j+.5)/7} * sin(pi n/7)/(pi n).
# |single-sector coeff| = |sin(pi n/7)|/(pi|n|). KEY: this VANISHES at 7|n and is a SINE kernel,
# NOT the crude 0.6973/|n| envelope. The crude bound IGNORES the sin(pi n/7) oscillation/cancellation.
print("=== sine-kernel structure of the per-coordinate Fourier coeff (the source of cancellation) ===")
print("  n : |sin(pi n/7)|/(pi n)  vs crude 0.6973/n   ratio (how much the crude OVER-counts)")
tot_sine=0.0; tot_crude=0.0
for n in range(1,50):
    s=abs(sin(pi*n/7))/(pi*n); c=0.6973/n
    tot_sine+=s; tot_crude+=c
    if n<=14 or n%7==0:
        print(f"  {n:3d}: {s:.5f}   {c:.5f}   {s/c:.3f}{'  <-- =0 (7|n)' if n%7==0 else ''}")
print(f"  partial sum to 49: sine-kernel {tot_sine:.4f} vs crude {tot_crude:.4f} (crude/sine ratio {tot_crude/tot_sine:.2f})")
print()
print("=== the bandlimiting insight (sound in principle) ===")
print("  A degree-D trig majorant M_A of 1_A has Fourier support [-D,D]. In")
print("  K(n)=sum_T (-1)^|T| prod_j chat_T(n_j), using M_A's coeffs => K(n)=0 unless ALL |n_j|<=D.")
print("  => the relation n has all coords <= D => FINITELY many lattice relations contribute.")
print("  So bandlimiting CONVERTS the infinite divergent sum to a FINITE sum + Selberg-Beurling error ~1/D.")
print("  Both the sine-kernel cancellation AND bandlimiting truncation are SOUND; the divergence in")
print("  MISTAKE-078 is purely the crude |chat|<=0.6973/|n| envelope discarding the oscillation.")
# Note: BOTH the sine-kernel partial sum AND crude still ~diverge as log (sum 1/n), so the
# pointwise envelope alone is NOT enough -- the CANCELLATION across the lattice (signed) is what
# converges. Confirm: sum |sin(pi n/7)|/(pi n) still ~ (2/pi^2)*ln + O(1) DIVERGES too.
print()
print("  CAVEAT: sum_n |sin(pi n/7)|/(pi n) ALSO diverges (~log), so the per-term sharper envelope is NOT")
print("  by itself convergent -- the win must come from (a) bandlimiting (truncate to |n|<=D) OR")
print("  (b) the SIGNED lattice cancellation (the relation-lattice constraint + alternating T-sum).")
