#!/usr/bin/env python3
"""hyperbessel_boundary_zeros_boxeph_S202.py -- boxeph-2026-07-21-S202

The GMC(2) cross-shell boundary (codex THM-2017): at the sharp degree-gap |e-rd|=r, NC2 holds iff
   Phi_{(p0,q0)}(xi) != 0,   Phi_{(p0,q0)}(x) = sum_{k>=0} x^k / ((q0 k)! (p0 k)!),  xi = alpha/(beta^r d^r).
Phi is a HYPER-BESSEL function. codex leaves "the discrete zero loci of the two boundary functions"
open. This brings the CLASSICAL zero theory (Bessel / Laguerre-Polya / total positivity):

 - Symmetric p0=q0=1: Phi_{(1,1)}(x)=sum x^k/(k!)^2 = I_0(2 sqrt(x)). NO positive real zeros (>0 for x>=0);
   zeros exactly at x = -(j_{0,k}/2)^2 (J_0 zeros), REAL NEGATIVE. => boundary NC2-clear off an EXPLICIT
   discrete negative-real locus.
 - General (p0,q0): are the zeros of Phi real (Laguerre-Polya) or complex? Is the coeff sequence
   c_k=1/((q0 k)!(p0 k)!) a Polya-frequency / log-concave (Turan) sequence? This locates codex's
   "inner resonance band".
"""
import numpy as np
from math import factorial, sqrt, lgamma, exp

def phi_coeffs(p0,q0,K):
    # log-space to avoid overflow: c_k = exp(-lgamma(q0 k +1) - lgamma(p0 k +1))
    return [exp(-lgamma(q0*k+1.0)-lgamma(p0*k+1.0)) for k in range(K+1)]

def zeros_of_truncation(coeffs, keep=None):
    # roots of sum_{k=0}^{K} c_k x^k ; return finite roots (the true entire-function zeros are the
    # stable limits as K grows). numpy.roots wants highest-degree first.
    c=list(coeffs)
    # trim to the prefix where coeffs stay well above underflow, so np.roots stays finite
    thr=abs(c[0])*1e-200
    K=len(c)
    for i in range(1,len(c)):
        if abs(c[i])<thr: K=i; break
    c=c[:max(K,4)]
    poly=c[::-1]
    r=np.roots(poly)
    r=r[np.isfinite(r)]
    if keep: r=sorted(r, key=lambda z: abs(z))[:keep]
    return np.array(r)

print("="*86)
print("SYMMETRIC boundary Phi_(1,1)(x) = I_0(2 sqrt(x)): zeros should be REAL NEGATIVE = -(j0k/2)^2")
print("="*86)
from numpy.polynomial import polynomial as P
# J_0 zeros j_{0,k}
j0=[2.404825558,5.520078110,8.653727913,11.79153444,14.93091771]
pred=[-(j/2)**2 for j in j0]
for K in (40,60,80):
    z=zeros_of_truncation(phi_coeffs(1,1,K))
    z=sorted([zz for zz in z if abs(zz.imag)<1e-6], key=lambda t:t.real, reverse=True)[:5]
    print("  K=%d smallest-|.| real zeros:" % K, ["%.4f"%zz.real for zz in z])
print("  predicted -(j0k/2)^2:", ["%.4f"%p for p in pred])
print("  => I_0 boundary zeros are real-negative (J_0), explicit; positive/complex-off-neg-axis xi => Phi!=0.")

print("\n"+"="*86)
print("GENERAL Phi_(p0,q0): reality of zeros (Laguerre-Polya?) + Turan/log-concavity of coeffs")
print("="*86)
def turan_ok(c):
    # Laguerre-Polya necessary: c_k^2 >= c_{k-1} c_{k+1} * (k+1)/k ... use plain Newton: c_k^2>=c_{k-1}c_{k+1}
    return all(c[k]*c[k]>=c[k-1]*c[k+1]-1e-18 for k in range(1,len(c)-1))
for (p0,q0) in [(1,1),(1,2),(1,3),(2,3),(2,2),(1,4),(3,4),(3,5)]:
    K=28
    c=phi_coeffs(p0,q0,K)
    z=zeros_of_truncation(c)
    # keep the "stable" small zeros (truncation artifacts appear at large |z|); look at smallest 8
    zs=sorted(z, key=lambda t:abs(t))[:8]
    max_im=max((abs(zz.imag)/max(1e-9,abs(zz)) for zz in zs), default=0.0)
    allreal = all(abs(zz.imag)<1e-4*max(1,abs(zz)) for zz in zs)
    allneg = all(zz.real<1e-6 for zz in zs if abs(zz.imag)<1e-4*max(1,abs(zz)))
    print("  (p0,q0)=(%d,%d): Turan-logconcave=%-5s ; smallest zeros %s ; all-real=%s all-neg=%s"
          % (p0,q0, turan_ok(c), ["%.2f%+.2fi"%(zz.real,zz.imag) for zz in zs[:4]], allreal, allneg))
print("  => reality/negativity of Phi_(p0,q0) zeros classifies the boundary exceptional locus;")
print("     complex zeros (if any) are the 'inner resonance band' (codex THM-2017 open piece).")

print("\n"+"="*86)
print("CONNECTION: xi = alpha/(beta^r d^r). For REAL positive-definite leading data xi>0 => Phi(xi)>0")
print("(all coeffs positive) => boundary NC2-clear unconditionally (matches THM-1535/1660 real case).")
print("The only boundary danger is xi landing on a (negative-real or complex) Phi-zero -- a codim>=1")
print("coefficient locus, removed one order down by the hyper-Bessel ODE theta^2 Phi = xi Phi (codex).")
print("="*86)
