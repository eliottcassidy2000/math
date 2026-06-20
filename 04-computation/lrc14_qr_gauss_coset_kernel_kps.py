#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) QR / GAUSS-SUM CANCELLATION in the coset kernel D7  (kind-pasteur-2026-06-19)

ANGLE: Euler's even-factor cancellation (1+x^n)=(1-x^{2n})/(1-x^n) <-> the quadratic
character mod 7.  The 6 inner sectors split into QR(7)={1,2,4} and NQR(7)={3,5,6},
the two <2>-orbits (doubling has order 3 mod 7).  The Legendre symbol chi mod 7
(chi=+1 on QR, -1 on NQR) governs the Gauss sum g=sum_r chi(r) zeta_7^r = i*sqrt(7)
(since 7 = 3 mod 4).

EXACT SETUP.  The support-6 Fourier kernel factorizes as
    K(n) = [ prod_{j=1..6}  1/(2 pi i n_j) ]  *  D7( n mod 7 )
where the ARCHIMEDEAN part is real-up-to-i^6 and depends on the integer magnitudes,
and the COSET KERNEL
    D7(c) = sum_{T subset {1..6}} (-1)^|T|  prod_{j=1..6}  [ (1 - zeta^{-c_j}) * sigma_T(c_j) ]
with sigma_T(m) = sum_{t in T} zeta^{-m t},  zeta = zeta_7,  c in (F_7^*)^6.

D7 is an EXACT element of Z[zeta_7] (degree-6 cyclotomic field).  We:
  (1) compute D7(c) exactly in Z[zeta_7] for all c in (F_7^*)^6 (7776 cosets),
  (2) SPLIT each coordinate by the quadratic character chi (QR vs NQR),
  (3) test whether D7 factors / simplifies / telescopes under the QR/NQR grouping,
  (4) test the Gauss-sum collapse: is sum_c D7(c) (or chi-twisted sums) a multiple
      of the Gauss sum g = i sqrt7, i.e. does the apex quadratic structure give the
      missing cancellation for the conditionally-convergent sum K(n)?

All arithmetic exact in Z[zeta_7] represented in the basis {1, z, z^2, z^3, z^4, z^5}
with the relation 1+z+...+z^6 = 0 (so z^6 = -(1+z+...+z^5)).
"""
import sys, itertools
from fractions import Fraction

sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

# ---------------- exact Z[zeta_7] arithmetic, basis 1..z^5, z^6 = -(sum z^k k<6) -------------
# element = tuple of 6 Fractions (coeff of z^0..z^5)
ZERO = (Fraction(0),)*6
def zfromexp(e):
    """zeta^e as a basis vector (e any integer)."""
    e %= 7
    if e <= 5:
        v = [Fraction(0)]*6; v[e] = Fraction(1); return tuple(v)
    # e == 6: z^6 = -(1+z+z^2+z^3+z^4+z^5)
    return (Fraction(-1),)*6
def zadd(a,b): return tuple(x+y for x,y in zip(a,b))
def zsub(a,b): return tuple(x-y for x,y in zip(a,b))
def zscale(a,s): return tuple(x*s for x in a)
def zmul(a,b):
    # multiply two degree<=5 polys, reduce mod z^6=-(...)
    res = [Fraction(0)]*12
    for i in range(6):
        if a[i]==0: continue
        for j in range(6):
            if b[j]==0: continue
            res[i+j]+=a[i]*b[j]
    # reduce powers 6..10 using z^k = z^(k-7)*z^7 ... easier: z^7=1 so z^k = z^(k mod 7)
    out=[Fraction(0)]*7
    for p in range(11):
        out[p%7]+=res[p] if p<12 else 0
    # now out has coeff of z^0..z^6; substitute z^6 = -(z^0+..+z^5)
    c6=out[6]
    final=[out[i]-c6 for i in range(6)]
    return tuple(final)
def zconj(a):
    """complex conjugate: zeta -> zeta^{-1} = zeta^6, i.e. z^k -> z^{-k}."""
    r=ZERO
    for k in range(6):
        if a[k]!=0:
            r=zadd(r, zscale(zfromexp(-k), a[k]))
    return r

# numeric eval for sanity
import cmath, math
Z = cmath.exp(2j*math.pi/7)
def znum(a): return sum(a[k]*(Z**k) for k in range(6))

# ---------------- the coset kernel D7(c) ----------------
def sigma_T(T, m):
    """sigma_T(m) = sum_{t in T} zeta^{-m t}  in Z[zeta_7].  T subset {1..6}, m in F_7^*."""
    r=ZERO
    for t in T:
        r=zadd(r, zfromexp((-m*t)%7))
    return r

def D7(c):
    """c = tuple of 6 residues in F_7^* (the n_j mod 7).
       D7(c) = sum_{T subset {1..6}} (-1)^|T| prod_j (1 - zeta^{-c_j}) sigma_T(c_j).
       Note (1 - zeta^{-c_j}) is a per-coordinate factor independent of T; pull it out."""
    # prefactor: prod_j (1 - zeta^{-c_j})
    pref=zfromexp(0)  # =1
    for cj in c:
        pref=zmul(pref, zsub(zfromexp(0), zfromexp((-cj)%7)))
    # sum over T of (-1)^|T| prod_j sigma_T(c_j)
    acc=ZERO
    for r in range(7):
        for T in itertools.combinations(range(1,7), r):
            prod=zfromexp(0)
            for cj in c:
                prod=zmul(prod, sigma_T(T, cj))
            term=zscale(prod, Fraction((-1)**r))
            acc=zadd(acc, term)
    return zmul(pref, acc)

# Legendre symbol mod 7
QR={1,2,4}; NQR={3,5,6}
def chi(m):
    m%=7
    if m==0: return 0
    return 1 if m in QR else -1

# Gauss sum g = sum_r chi(r) zeta^r  (should be i*sqrt7)
def gauss_sum():
    g=ZERO
    for r in range(1,7):
        g=zadd(g, zscale(zfromexp(r), Fraction(chi(r))))
    return g

if __name__=="__main__":
    print("LRC(14) QR/GAUSS coset kernel D7 in Z[zeta_7]  (kps)")
    g=gauss_sum()
    print(f"  Gauss sum g = sum chi(r) z^r ; numeric = {znum(g):.6f}  (expect i*sqrt7 = {1j*math.sqrt(7):.6f})")
    print(f"    g^2 numeric = {znum(zmul(g,g)):.6f}  (expect -7 since 7=3 mod 4)")

    # ---- (1) basic D7 values; check it's a genuine cyclotomic integer ----
    print("\n(1) Sample D7(c) values (exact basis coeffs + numeric):")
    for c in [(1,1,1,1,1,1),(1,2,4,3,5,6),(1,1,1,1,1,2),(1,2,3,4,5,6)]:
        d=D7(c)
        print(f"    c={c}: coeffs={[str(x) for x in d]}  num={znum(d):.4f}")

    # ---- (2) FULL sum over all cosets, and chi-twisted sums ----
    print("\n(2) Sums over ALL c in (F_7^*)^6  (7776 cosets):")
    S=ZERO; Sabs=0.0
    Schi=ZERO         # sum chi(prod c_j) D7(c)
    Schi1=ZERO        # sum chi(c_1) D7(c)
    cnt=0
    for c in itertools.product(range(1,7),repeat=6):
        d=D7(c)
        S=zadd(S,d); Sabs+=abs(znum(d)); cnt+=1
        prodchi=1
        for cj in c: prodchi*=chi(cj)
        Schi=zadd(Schi, zscale(d, Fraction(prodchi)))
        Schi1=zadd(Schi1, zscale(d, Fraction(chi(c[0]))))
    print(f"    #cosets={cnt}")
    print(f"    raw   sum_c D7(c)         = {[str(x) for x in S]}  num={znum(S):.6f}")
    print(f"    chi-prod twisted sum      = {[str(x) for x in Schi]}  num={znum(Schi):.6f}")
    print(f"    chi(c1)  twisted sum      = {[str(x) for x in Schi1]}  num={znum(Schi1):.6f}")
    print(f"    sum |D7(c)| (numeric)     = {Sabs:.4f}")
    # is the raw sum a rational multiple of g?  test S / g componentwise via numeric
    ng=znum(g)
    print(f"    raw sum / g  (numeric)    = {znum(S)/ng:.6f}   (rational => QR collapse)")
    print(f"    chi-prod sum / g          = {znum(Schi)/ng:.6f}")
