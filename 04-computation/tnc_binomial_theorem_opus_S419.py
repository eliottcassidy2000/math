import sympy as sp
from math import gcd, comb
u=sp.symbols('u')
def CT(R,N,mm):
    return sp.Poly(sp.expand(sp.sympify(R)**mm),u).coeff_monomial(u**(N*mm))
print("THEOREM (to verify): TNC HOLDS for BINOMIAL R = r0 + rd u^d (d=M+N).")
print("  Lambda = r0 u^{-N} + rd u^{M}, charges {-N, M} ONLY.  CT(Lambda^m) = sum over")
print("  a+b=m with -Na+Mb=0.  Minimal solution a=M/g, b=N/g (g=gcd(M,N)), m0=(M+N)/g,")
print("  and CT(m0) = C(m0, M/g) r0^{M/g} rd^{N/g} -- a SINGLE term, no cancellation, != 0.")
print("  So a binomial is NEVER a nullcone element => TNC holds at every bidegree for 2-term R.")
print()
print("  N   M   d=M+N   g=gcd   m0=(M+N)/g   CT(m0) predicted    CT(m0) actual")
ok=True
for N in range(1,5):
    for M in range(1,5):
        d=M+N; g=gcd(M,N); m0=(M+N)//g
        R=f"2 + 3*u**{d}"     # r0=2, rd=3
        a=M//g; b=N//g
        pred=comb(m0,a)*2**a*3**b
        act=CT(R,N,m0)
        # also check ALL earlier m vanish
        earlier=all(CT(R,N,m)==0 for m in range(1,m0))
        good=(pred==act) and act!=0 and earlier
        ok=ok and good
        print(f"  {N}   {M}    {d}      {g}       {m0}          {pred:8d}          {act}   "
              f"{'earlier all 0' if earlier else 'EARLIER NONZERO!'}  {'OK' if good else 'FAIL'}")
print(f"\n  binomial theorem holds on all tested (N,M): {ok}")
print()
print("NECESSARY CONDITION r_N = 0 (Lambda has no constant term), since CT(1)=[u^N]R=r_N:")
for R,N in [("1+u+u**2",2),("1+u**3+u**4",2),("5+u+2*u**4",2)]:
    rN=sp.Poly(sp.sympify(R),u).coeff_monomial(u**N)
    print(f"   R={R:16s} N={N}: r_N={rN}, CT(1)={CT(R,N,1)}  (match: {rN==CT(R,N,1)})")
print()
print("WHERE IT STAYS OPEN: >=3 charges allow CANCELLATION among representations of 0 at")
print("the minimal m.  Smallest genuine test: TRINOMIAL with r_N=0.  Find first nonzero CT:")
for R,N in [("1+u**3+u**4",2),("1+u+u**4",2),("1+u**3+u**5",2),("2+u**3+u**6",3),
            ("1+u**4+u**5",3),("1+u**2+u**5",3)]:
    fm=next((m for m in range(1,40) if CT(R,N,m)!=0),None)
    d=sp.Poly(sp.sympify(R),u).degree()
    print(f"   R={R:16s} N={N} (M={d-N}): first nonzero CT at m={fm}, value {CT(R,N,fm) if fm else '--'}")
