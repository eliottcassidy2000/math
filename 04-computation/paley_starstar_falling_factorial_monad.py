#!/usr/bin/env python3
"""
paley_starstar_falling_factorial_monad.py
monad-explorer-2026-06-07 (deep-research, 11th session)

THE CLEAN COLUMN FORM (granting handoff #1):
   t(k,m) = (k)_m * h_m(k)  =  m! C(k,m) h_m(k),   deg h_m = m-2,
where (k)_m = k(k-1)...(k-m+1) is the falling factorial (= #ordered m-subsets of [k]).
Reason: the column Hilbert polynomial p_m(k) (degree 2m-2, exists by the PROVED pole
order 2m-1) vanishes at k=1,...,m-1 (no rank-m patterns there); handoff #1 adds the
zero at k=0, giving m CONSECUTIVE roots 0..m-1, i.e. the factor (k)_m.  The remaining
factor h_m has degree (2m-2)-m = m-2.

Both handoffs become statements about h_m:
   handoff #1  <=>  (k)_m | p_m(k)  <=>  p_m(0)=0  (h_m is an honest poly).
   handoff #2  <=>  h_m(-1) = (2^m-1)/((-1)^m m!)    [ since (-1)_m=(-1)^m m! ].
And the WILD end: h_m(m) = t(m,m)/(m)_m = A088368(m)/m!  -> e  as m->inf (A088368~e*m!).
So h_m(m)->e (tame rational values), h_m(-1)~Mersenne/m! ; h_m is the BRIDGE polynomial.
"""
import sympy as sp
from math import comb, factorial

x, k = sp.symbols('x k')
cols = {
 1: x/(1-x),
 2: 3*x**2/(1-x)**3,
 3: (13+7*x)*x**3/(1-x)**5,
 4: (69+97*x+15*x**2)*x**4/(1-x)**7,
}
A088368 = {1:1,2:3,3:13,4:69,5:421}

def col_entry(m, kk):
    ser = sp.series(cols[m], x, 0, kk+1).removeO()
    return int(sp.Poly(ser, x).coeff_monomial(x**kk)) if kk>0 else 0

print("Verify t(k,m) = (k)_m * h_m(k), deg h_m = m-2:")
for m in range(2,5):
    # interpolate p_m(k) from enough points (k=m..3m), degree 2m-2
    npts = 2*m   # need 2m-1 points for degree 2m-2; take 2m
    ks = list(range(m, m+npts))
    vals = [col_entry(m, kk) for kk in ks]
    p = sp.interpolate(list(zip(ks, vals)), k)
    p = sp.expand(p)
    deg = sp.Poly(p, k).degree()
    # falling factorial (k)_m
    fall = sp.prod([k-i for i in range(m)])
    hq, hr = sp.div(sp.Poly(p, k), sp.Poly(fall, k))
    print(f"\n m={m}: p_m(k) deg={deg} (want {2*m-2})")
    print(f"        p_m(k) = {p}")
    print(f"        (k)_m divides? remainder = {sp.expand(hr.as_expr())}  (0 <=> handoff#1)")
    h = sp.expand(hq.as_expr())
    print(f"        h_m(k) = {h}   deg={sp.Poly(h,k).degree() if h!=0 else 0} (want {m-2})")
    # checks
    print(f"        p_m(0)={p.subs(k,0)} (handoff#1 wants 0);  p_m(-1)={p.subs(k,-1)} (handoff#2 wants 2^m-1={2**m-1})")
    print(f"        h_m(-1)={sp.nsimplify(h.subs(k,-1))}  vs (2^m-1)/((-1)^m m!)={sp.Rational(2**m-1, (-1)**m*factorial(m))}")
    print(f"        h_m(m)={sp.nsimplify(h.subs(k,m))}  =A088368(m)/m! = {sp.Rational(A088368[m], factorial(m))} (->e)")

print("\nWILD-END limit check: A088368(m)/m! for m=1..5:",
      [sp.Rational(A088368[m], factorial(m)) for m in range(1,6)],
      " (-> e = 2.718...)")
print(" floats:", [float(sp.Rational(A088368[m], factorial(m))) for m in range(1,6)])
