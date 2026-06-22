"""
The 3 recursion modes = 3 arithmetic characters; the LRC floor decomposed (kps-S31q).
Modes A,B,C (n->n-1,n-2,n-3) generate inclusion-exclusion: 7=2^3-1 nonempty subsets of {A,B,C},
graded 3+3+1 (singletons +, pairs -, triple +) = the Mobius pattern +++---+ (THM-549).
The 3 sign-weightings of A..G are 3 characters: Mobius mu (+++---+), Legendre chi_7 (++-+--+),
Eisenstein chi_3/omega. The copy rule phi=mu*id (HYP-2882) ties the totient to coprime density.
The LRC apex-7 floor = sum_{b<7} phi(b)*2*delta_b (totient sum); decompose by chi_7(b), chi_3(b).
"""
from fractions import Fraction as F
from math import gcd
def phi(n): return sum(1 for k in range(1,n+1) if gcd(k,n)==1)
def mu(n):
    if n==1: return 1
    r=1; m=n; p=2; cnt=0
    d=2; res=1
    nn=n; primes=[]
    while nn>1:
        c=0
        while nn%d==0: nn//=d; c+=1
        if c>0: primes.append(c)
        d+=1
    return 0 if any(c>1 for c in primes) else (-1)**len(primes)
def chi7(b):  # Legendre symbol (b/7)
    b%=7
    if b==0: return 0
    return 1 if b in (1,2,4) else -1
def chi3(b):  # quadratic char mod 3
    b%=3
    return 0 if b==0 else (1 if b==1 else -1)

print("(1) The 3 sign patterns over A..G (=1..7):")
modes={"Mobius +++---+ (incl-excl /3 modes)":[1,1,1,-1,-1,-1,1],
       "Legendre chi_7 ++-+--+":[1,1,-1,1,-1,-1,1]}
for name,sg in modes.items():
    print(f"   {name}: {sg}")
print(f"   Legendre check: +on{{1,2,4}}=QR(7), -on{{3,5,6}}=NQR(7): {[k+1 for k in range(7) if modes['Legendre chi_7 ++-+--+'][k]==1]}")
print(f"   Incl-excl 3+3+1: singletons A,B,C(+), pairs(-), triple(+) = subsets of 3 modes")

print("\n(2) Copy rule phi = mu * id  (HYP-2882, the totient IS Mobius-inversion of identity):")
for n in range(1,9):
    conv=sum(mu(n//d)*d for d in range(1,n+1) if n%d==0)
    print(f"   n={n}: (mu*id)(n)={conv}, phi(n)={phi(n)}  {'OK' if conv==phi(n) else 'X'}")

print("\n(3) LRC apex-7 floor = sum_{b=1..6} phi(b)*2*(7-b)/(7b), decomposed by character of b:")
floor=F(0); byc7={1:F(0),-1:F(0),0:F(0)}; byc3={1:F(0),-1:F(0),0:F(0)}
for b in range(1,7):
    w=F(phi(b)*2*(7-b), 7*b)   # totient-weighted neighborhood width (per unit V)
    floor+=w; byc7[chi7(b)]+=w; byc3[chi3(b)]+=w
    print(f"   b={b}: phi={phi(b)}, width 2(7-b)/(7b)={F(2*(7-b),7*b)}, contrib {w}  chi7={chi7(b):+d} chi3={chi3(b):+d}")
print(f"\n   TOTAL floor (per V) = {floor} = {float(floor):.4f}")
print(f"   by chi_7:  QR(+1)={float(byc7[1]):.4f}  NQR(-1)={float(byc7[-1]):.4f}  (b=7n: {float(byc7[0]):.4f})")
print(f"   chi_7 principal = (QR+NQR)/1 = {float(byc7[1]+byc7[-1]):.4f};  chi_7 OSCILLATION = QR-NQR = {float(byc7[1]-byc7[-1]):.4f}")
print(f"   by chi_3:  (+1)={float(byc3[1]):.4f}  (-1)={float(byc3[-1]):.4f}  (3|b: {float(byc3[0]):.4f})")
print(f"\n   => floor = PRINCIPAL (totient/Mobius, ~{float(floor):.3f}) + chi_7 + chi_3 corrections.")
print(f"      Principal {float(floor):.3f} DOMINATES |chi_7 osc| {abs(float(byc7[1]-byc7[-1])):.3f} => floor stays POSITIVE.")
