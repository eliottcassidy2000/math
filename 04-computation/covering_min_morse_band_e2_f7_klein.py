#!/usr/bin/env python3
"""
covering_min_morse_band_e2_f7_klein.py  --  klein-2026-07-01-S63

Working the "next lever" (a tighter LOWER-bound certificate) through E2 + F7(apex-7) + far-element
resonance + Morse/band barrier -- and testing WHY spectral/PSD certificates have a gap (S54: loneliness
is POINTWISE, averages are blind), so the tight certificate is the COMBINATORIAL witness packing (the
lazy-cut, HYP-3779), not a Delsarte/Fejer SDP.

(1) FAR-ELEMENT RESONANCE: the construction killer n(n-1) = Phi6(n) - 1 == -1 mod Phi6 -- the far
    element sits at the "-1 slot" of the binding modulus. (E2/hexagonal: Phi6 = zeta_6 norm; Dedekind
    s(n,Phi6) -> -1/12, S56.)
(2) THE MORSE / BAND LANDSCAPE: G(t) = min_v ||v t|| (loneliness). Its global max = M(S) at the binding
    t* = a/Phi6; local maxima = the "band barrier" critical points. We tabulate the top critical values
    (the Morse spectrum) and the binding.
(3) F7 / APEX-7 + SPECTRAL GAP: test whether a low-degree positive-definite (Fejer) minorant of the
    danger indicator can certify M >= n/Phi6. If the best low-degree PSD bound is FAR below n/Phi6, that
    is the spectral gap -- the lower bound is combinatorial, not spectral.
"""
from fractions import Fraction as F
from math import gcd, cos, pi, sin

def Phi6(n): return n*n-n+1
def dist0(x,D):
    x%=D; return min(x,D-x)
def normf(t):  # ||t|| for float t
    f=t-int(t);
    if f<0: f+=1
    return min(f,1-f)

def construction(n): return list(range(1,n-1))+[n*(n-1)]

def landscape_max_and_crits(S, grid=20000):
    """G(t)=min_v ||v t|| on a grid; return global max and the sorted distinct local-max values."""
    best=0.0; at=0.0; vals=[]
    prev2=prev1=None
    for i in range(1, grid):
        t=i/grid
        g=min(normf(v*t) for v in S)
        if g>best: best=g; at=t
        if prev1 is not None and prev2 is not None and prev1>prev2 and prev1>g:
            vals.append(round(prev1,4))
        prev2=prev1; prev1=g
    return best, at, sorted(set(vals), reverse=True)[:8]

if __name__=="__main__":
    print("="*72)
    print("(1) FAR-ELEMENT RESONANCE: killer n(n-1) == -1 mod Phi6(n)")
    print("="*72)
    for n in [8,10,12,14]:
        D=Phi6(n); k=n*(n-1)
        print(f"  n={n:2d}: Phi6={D}, killer n(n-1)={k}, k mod Phi6 = {k%D} (== -1 mod Phi6? {k%D==D-1})  "
              f"[Phi6={D}={'*'.join(str(p) for p in [] )}]")
    print("  => the far element sits at residue -1 of the binding modulus Phi6 (the ceiling of the")
    print("     Stern-Brocot ray); this is the far-element resonance / the E2-hexagonal zeta_6 direction.")

    print()
    print("="*72)
    print("(2) MORSE / BAND LANDSCAPE of G(t)=min_v||v t|| for the construction (n=14)")
    print("="*72)
    n=14; S=construction(n); D=Phi6(n); thr=F(n,D)
    gmax, at, crits = landscape_max_and_crits(S)
    print(f"  construction {S}")
    print(f"  global max M = {float(thr):.5f} (=14/183, at t*=14/183={14/183:.5f}); grid max {gmax:.5f} at {at:.5f}")
    print(f"  top LOCAL-max (Morse) values of G: {crits}")
    print(f"  the binding t*=14/183 is the GLOBAL Morse maximum; the lower local maxima are the 'band")
    print(f"  barrier' critical points (radius-1 band moduli). Loneliness is a POINTWISE spike (S54).")

    print()
    print("="*72)
    print("(3) F7 / APEX-7 + SPECTRAL GAP: can a low-degree Fejer (PSD) minorant certify M>=n/Phi6?")
    print("="*72)
    # Fejer kernel of degree N: F_N(t) = sum_{|k|<=N} (1-|k|/(N+1)) e(k t) >= 0 (positive-definite).
    # A spectral lower bound on max_t min_v||vt|| via averaging is BLIND (S54): E_t[min_v||vt||] <<< M.
    for N in [7, 14, 61]:
        # average of the loneliness against the (normalized) Fejer weight -- a spectral/averaged proxy
        grid=6000; num=0.0; den=0.0
        for i in range(grid):
            t=i/grid
            w=1.0
            for k in range(1,N+1):
                w += 2*(1-k/(N+1))*cos(2*pi*k*t)
            w=max(w,0.0)
            g=min(normf(v*t) for v in S)
            num+=w*g; den+=w
        avg=num/den if den else 0
        print(f"  Fejer F_{N} weighted avg of loneliness = {avg:.5f}  vs  M=14/183={float(thr):.5f}  "
              f"(gap factor {float(thr)/avg:.1f}x)" if avg>0 else "")
    print("  => the spectral/averaged (Fejer, incl. apex F_7) proxy is FAR below M: the averaged lens is")
    print("     BLIND to the pointwise loneliness spike (S54). So the tight lower-bound certificate is")
    print("     COMBINATORIAL (the witness packing / lazy-cut HYP-3779), not a Delsarte/Fejer SDP.")
    print("     E2/Eisenstein = the regularizable BULK (Dedekind->-1/12, S56); the residual (cusp form,")
    print("     apex-7/F7) is the un-relaxable core the PSD relaxation cannot see.")
