#!/usr/bin/env python3
"""
S650 — (A) the fiberable-vs-prime-hard catalogue of LRC(n), and (B) the sqrt(-19)/chi=5 CM handle on 19.

(A) The S640 reduction fibers n over an odd prime divisor p (CRT), base = LRC(p), fiber = mult-of-p
    runners, with the doubling <2> mod p giving the sub-shell. n is 'reducible' when it has a small
    odd prime divisor p with LRC(p) PROVEN (p<=7) and a small fiber. Primes are the hard end.
(B) For prime n=19: witnesses t=j/19 live in Q(zeta_19); its quadratic subfield is Q(sqrt(-19)) (the
    Heegner field, since 19=3 mod 4 -> Gauss sum g^2 = -19). Verify g^2=-19 numerically; map the QR/
    Paley structure and the chi-tower sqrt(-3)->sqrt(-11)->sqrt(-19).
No external libs (complex via cmath).
"""
import cmath, math
from fractions import Fraction

def isprime(m): return m>1 and all(m%d for d in range(2,int(m**0.5)+1))
def primefactors(m):
    f={}; d=2
    while d*d<=m:
        while m%d==0: f[d]=f.get(d,0)+1; m//=d
        d+=1
    if m>1: f[m]=f.get(m,0)+1
    return f
def ordp(a,p):
    x=a%p; k=1
    while x!=1: x=(x*a)%p; k+=1
    return k

print("="*72)
print("(A) the fiberable-vs-prime-hard catalogue of LRC(n)  (LRC proven for <=7 runners)")
print("="*72)
print("  n  | factor | smallest odd p | LRC(p) proven? | ord_p(2) | reducible? | type")
PROVEN=7
for n in range(8,33):
    pf=primefactors(n)
    odd_primes=[p for p in pf if p%2==1]
    if isprime(n):
        kind="PRIME (no fiber, hard)"; red="NO"; p="-"; op="-"; base="-"
    else:
        p=min(odd_primes) if odd_primes else None
        if p is None:  # power of 2
            kind="2^k (fully 2-adic)"; red="(2-adic)"; op="-"; base="-"
        else:
            op=ordp(2,p); base = (p<=PROVEN)
            red = "YES" if base else "weak"
            kind=f"fiber/{p}" + (" (base proven)" if base else " (base open)")
    fac="*".join(f"{q}^{e}" if e>1 else str(q) for q,e in pf.items())
    pstr=str(p) if not isprime(n) and p else "-"
    bstr=("YES" if (not isprime(n) and odd_primes and min(odd_primes)<=PROVEN) else ("-" if isprime(n) else "no"))
    print(f"  {n:2d} | {fac:7s}| {pstr:>14} | {bstr:>14} | {str(op):>8} | {red:>10} | {kind}")
print("  -> n=2p with p<=7 (n=6,10,14) reduce to a PROVEN base (S640). n=2p with p>7 (22,26) reduce to")
print("     an OPEN base. PRIMES (11,13,17,19,23,29,31) have NO fiber = the hard end. 19 is prime-hard.")

print("\n" + "="*72)
print("(B) the sqrt(-19) / chi=5 CM handle on prime n=19")
print("="*72)
p=19
# Legendre symbol
def leg(a,p):
    a%=p
    if a==0: return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1
QR=sorted(a for a in range(1,p) if leg(a,p)==1)
print(f"  19 = 3 mod 4 ({p%4==3}); QR mod 19 = {QR}  (|QR|=(19-1)/2={ (p-1)//2 })")
print(f"  -1 a residue? {leg(-1,p)==1}  (NO -> Gauss sum gives sqrt(-19), Paley-19 exists, S638)")
print(f"  2  a residue? {leg(2,p)==1}  (NO -> 2 primitive root, ord_19(2)={ordp(2,19)}; <2>!=QR)")

# quadratic Gauss sum g = sum_a leg(a,p) zeta^a ;  g^2 should = (-1)^((p-1)/2) p = -19
zeta=cmath.exp(2j*math.pi/p)
g=sum(leg(a,p)*zeta**a for a in range(1,p))
g2=g*g
print(f"\n  quadratic Gauss sum g = sum_a (a|19) zeta_19^a:")
print(f"    g  = {g.real:+.4f}{g.imag:+.4f}i   (should be purely imaginary = i*sqrt(19))")
print(f"    g^2= {g2.real:+.4f}{g2.imag:+.4f}i   -> = -19 ? {abs(g2-(-19))<1e-6}")
print(f"    => sqrt(-19) = the Gauss sum LIVES in Q(zeta_19): the quadratic subfield of the field where")
print(f"       the 19-runner witnesses t=j/19 live. Q(sqrt(-19)) is the index-2 (Paley/QR) level.")

print(f"\n  Q(zeta_19): degree phi(19)=18 over Q; Galois (Z/19)* = Z/18 (cyclic, ABELIAN, solvable).")
print(f"  subfields <-> divisors of 18 (1,2,3,6,9,18); the degree-2 one = Q(sqrt(-19)) (the QR level).")
print(f"  So the witness tower (opus S704) is the cyclotomic/abelian tower; sqrt(-19) is its QUADRATIC")
print(f"  rung = where the Paley/Gauss/QR structure lives.")

print(f"\n  THE HEEGNER / chi-tower:")
heeg=[1,2,3,7,11,19,43,67,163]
for d,chi,who in [(3,3,'Eisenstein lattice'),(11,4,'Moser spindle'),(19,5,'19-runner / hex(2)')]:
    print(f"    sqrt(-{d:>3}): class no. 1 (Heegner {d in heeg}); chi={chi} step; {who}")
print(f"  19=4*5-1 (HYP-2277 rotation field N=5); 19=1+6+12=hex(2); 2n-1=37=hex(3).")

print("\n" + "="*72)
print("HONEST SUMMARY")
print("="*72)
print("  (A) Reducibility is a DIVISOR property: n=2p (p<=7) reduce to proven LRC(p); primes don't")
print("      reduce at all. 19 prime => prime-hard for fibers (verified: 2 primitive root, single cycle).")
print("  (B) 19's leverage is the CM field sqrt(-19) = the Gauss sum = the quadratic subfield of Q(zeta_19)")
print("      = the index-2 Paley/QR level of the abelian witness tower; Heegner (class 1) = the chi=5 rung.")
print("      This ORGANIZES the 19-runner witnesses (their exact field) -- it does NOT prove LRC(19) (open);")
print("      the genuine attack is the cyclotomic depth q* (S704) at q=19,37, with sqrt(-19) the CM handle.")
