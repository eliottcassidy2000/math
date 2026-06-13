#!/usr/bin/env python3
"""
paley_starstar_handoff_finitesum_monad.py
monad-explorer-2026-06-07 (deep-research, 11th session)

BOTH THM-438 handoffs reduced to FINITE ALTERNATING-BINOMIAL SUMS over the
genuine triangle entries t(k,m) (k=m..2m-1) -- NOT over the R_s transform.

Derivation (hockey-stick on the s-expansion inversion R_s(m,e)=sum_k (-1)^{e-k}C(e-1,k-1)t(k,m)):
   handoff #1 (deg P_m=m-2  <=>  Q_m(-1)=0):
       sum_{k=m}^{2m-1} (-1)^k C(2m-1, k) t(k,m) = 0
       [uses  sum_{e=k}^{2m-1} C(e-1,k-1) = C(2m-1,k) ]
   handoff #2 (lead P_m = 2^m-1 = Q_m'(-1)):
       sum_{k=m}^{2m-1} (-1)^{k+1} k C(2m, k+1) t(k,m) = 2^m - 1
       [uses  e C(e-1,k-1) = k C(e,k),  sum_{e=k}^{2m-1} C(e,k) = C(2m,k+1) ]

These hold for the column entries directly.  Verified here m=2,3,4 (full columns
known) using the EXACT rational columns to read off t(k,m) for k=m..2m-1.
This relocates the proof to a SIGN-REVERSING INVOLUTION on even-series PATTERNS
(clean objects), pairing patterns across k with binomial multiplicity.
"""
import sympy as sp
from math import comb

x = sp.symbols('x')
# exact columns T_m(x); read t(k,m) = [x^k] T_m
cols = {
 1: x/(1-x),
 2: 3*x**2/(1-x)**3,
 3: (13+7*x)*x**3/(1-x)**5,
 4: (69+97*x+15*x**2)*x**4/(1-x)**7,
}
def t(k,m):
    if k < m: return 0
    ser = sp.series(cols[m], x, 0, k+1).removeO()
    return int(sp.Poly(ser, x).coeff_monomial(x**k)) if ser.has(x) or k==0 else int(ser)

print("column entries t(k,m), k=m..2m-1:")
for m in range(2,5):
    row = [t(k,m) for k in range(m, 2*m)]
    print(f"  m={m}: {dict(zip(range(m,2*m), row))}")

print("\nHANDOFF #1:  sum_{k=m}^{2m-1} (-1)^k C(2m-1,k) t(k,m)  (want 0)")
for m in range(2,5):
    S = sum((-1)**k * comb(2*m-1, k) * t(k,m) for k in range(m, 2*m))
    print(f"  m={m}: {S}   terms: " +
          " ".join(f"{(-1)**k*comb(2*m-1,k)}*{t(k,m)}" for k in range(m,2*m)))

print("\nHANDOFF #2:  sum_{k=m}^{2m-1} (-1)^(k+1) k C(2m,k+1) t(k,m)  (want 2^m-1)")
for m in range(2,5):
    S = sum((-1)**(k+1) * k * comb(2*m, k+1) * t(k,m) for k in range(m, 2*m))
    print(f"  m={m}: {S}   (2^m-1 = {2**m-1})   match={S==2**m-1}")

# Also express as a single statement on the m=5 column (predicts relations among t(7,5),t(8,5),t(9,5)):
print("\nm=5 forms (in terms of t(5,5)=421,t(6,5)=4845,t(7,5),t(8,5),t(9,5)):")
import sympy as sp
t75,t85,t95 = sp.symbols('t75 t85 t95')
tt = {5:421,6:4845,7:t75,8:t85,9:t95}
S1 = sum((-1)**k*comb(9,k)*tt[k] for k in range(5,10))
S2 = sum((-1)**(k+1)*k*comb(10,k+1)*tt[k] for k in range(5,10))
print("  handoff#1: ", sp.expand(S1), "= 0")
print("  handoff#2: ", sp.expand(S2), "= 31")
print("  (two linear relations among the three unknown entries t(7,5),t(8,5),t(9,5);")
print("   together with c0,c1 known they are equivalent to c4=0 and c3=31.)")
