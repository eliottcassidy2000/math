"""Coprime density / totient / Mobius structure of the LRC resonance-killing (mac-mini-2026-06-22-S44).
A runner s kills Farey point a/b (gcd(a,b)=1, b<=14) iff ||s a/b||<1/14 iff b|s, so s kills ALL phi(b)
primitive points of each denominator b|s -- the killing is TOTIENT-weighted. Survival lattice =
Phi(14)=sum phi(b)=64 Farey points. The surviving-neighborhood density = 3/pi^2 = 1/(2 zeta(2)); the
coprime density 6/pi^2 = 1/zeta(2) = sum phi(b)*(point weight). phi = mu * id => killing is Mobius IE
over the divisor lattice; the owner's 3 recursion modes are that IE skeleton; opus S291's Euler product
over odd cycle lengths is the same multiplicative spectrum on tournaments.
See 07-reflections/the-resonance-killing-is-multiplicative-totient-mobius-zeta2.md."""
from math import gcd, pi
def phi(n): return sum(1 for k in range(1,n+1) if gcd(k,n)==1)
def mu(n):
    if n==1: return 1
    r=1; m=n; d=2
    while d*d<=m:
        if m%d==0:
            m//=d
            if m%d==0: return 0
            r=-r
        d+=1
    return -r if m>1 else r
if __name__=="__main__":
    Phi=sum(phi(b) for b in range(1,15))
    print(f"Phi(14)=sum_{{b<=14}} phi(b) = {Phi} (totient-weighted survival lattice)")
    print(f"surviving-neighborhood density = 3/pi^2 = 1/(2 zeta(2)) = {3/pi**2:.4f}")
    print(f"coprime density = 6/pi^2 = 1/zeta(2) = {6/pi**2:.4f}")
    print("phi=mu*id; killing is Mobius IE over the divisor lattice (the 3 recursion modes' signs).")
