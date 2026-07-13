#!/usr/bin/env python3
"""
lrc14_vdc_clean_klein_S276.py
=============================
klein-2026-07-13-S276. Corrected van der Corput test. The prior delta metric was dominated by e'=1
(||w/1||=0 always); e'=1 gives only 7 bounded terms (harmless). The REAL variable is the Diophantine
separation of w from the LARGER offsets. Denjoy-Koksma predicts per-offset contribution
~ gcd(w,e')^2/e' (clean) blowing to ~e' at resonance; total = a RESONANCE SUM, O(1) clean, O(Sigma e')
at w=lcm.

Tests:
 (1) genuinely CLEAN w (prime, max_over-candidates min_{e'>=2}||w/e'||): err*w vs growing k AND Sigma-e'
     -> is it BOUNDED (O(k)/O(1)) at clean w?
 (2) err*w vs the resonance sum R(C,w)=sum_{e'>=2} min(1, 1/(e'*||w/e'||)) -> does err*w <~ C*R?
"""
import math
from math import gcd
from functools import reduce
NG=999983
def lcm(xs):
    r=1
    for x in xs:
        if x:r=r*x//gcd(r,x)
    return r
def isprime(n):
    if n<2:return False
    for p in range(2,int(n**0.5)+1):
        if n%p==0:return False
    return True
def occ_of(E,x):
    o=0
    for e in E:o|=1<<(int((e*x%1.0)*7.0)%7)
    return o
def moments3(E):
    s1=s2=s3=0
    for k in range(1,NG):
        N=7-bin(occ_of(E,k/NG)).count("1");s1+=N;s2+=N*(N-1);s3+=N*(N-1)*(N-2)
    n=NG-1;return s1/n,s2/n,s3/n
def Phi(E):m1,m2,m3=moments3(E);return 1-(2/3)*m1+(47/252)*m2-(5/252)*m3
def Phi_inf(C):m1,m2,m3=moments3(C);return 1-(2/3)*(6/7)*m1+(47/252)*(5/7)*m2-(5/252)*(4/7)*m3
def errw(C,w):return abs(Phi(C+[w])-Phi_inf(C))*w
def mindelta(C,w):  # min_{e'>=2} ||w/e'||
    d=1.0
    for e in C:
        if e>=2:
            f=(w/e)%1.0;d=min(d,f,1-f)
    return d
def ressum(C,w):  # sum_{e'>=2} min(1, 1/(e'*||w/e'||))
    S=0.0
    for e in C:
        if e>=2:
            f=(w/e)%1.0;nm=min(f,1-f)
            S+=min(1.0, 1.0/(e*nm)) if nm>0 else 1.0
    return S
def clean_prime_near(target, C):  # prime near target maximizing min_{e'>=2}||w/e'||
    best=None;bd=-1
    n=target
    for _ in range(400):
        if isprime(n):
            d=mindelta(C,n)
            if d>bd: bd=d;best=n
        n+=1
    return best,bd

print("="*72)
print("(1) genuinely CLEAN w (prime, well-separated from e'>=2): err*w vs k and Sigma-e'")
print("="*72)
print(f"  {'desc':22s} {'k':>2} {'sum':>4} {'w':>6} {'mindelta':>8} {'err*w':>7}")
# fixed k=7, grow diameter
for D in [12,24,48,96,192,384]:
    C=[0,1,2,3,4,5,D]; w,bd=clean_prime_near(17*D,C)
    print(f"  {'k=7 diam='+str(D):22s} {7:2d} {sum(C):4d} {w:6d} {bd:8.3f} {errw(C,w):7.3f}")
# fixed diam=150, grow k
for k in [4,6,8,10,12]:
    C=sorted(set(round(i*150/(k-1)) for i in range(k)))
    w,bd=clean_prime_near(2550,C)
    print(f"  {'diam=150 k='+str(k):22s} {len(C):2d} {sum(C):4d} {w:6d} {bd:8.3f} {errw(C,w):7.3f}")
print("  => if err*w stays bounded (<~0.2) as k,Sigma-e' grow at clean w: O(k) bound CONFIRMED.")

print()
print("="*72)
print("(2) err*w vs resonance sum R(C,w)=sum_{e'>=2} min(1,1/(e'||w/e'||)); clean->resonant w")
print("="*72)
C=[0,1,2,28,29,30,15]; L=lcm(C)
print(f"  C={C} lcm={L}")
print(f"  {'w':>7} {'R(C,w)':>8} {'err*w':>7} {'err*w/R':>8}")
for w in sorted(set([1009,2003,4001,6007, L//6,L//3,L//2,L,2*L])):
    if not(0<w<NG//30):continue
    R=ressum(C,w); ew=errw(C,w)
    print(f"  {w:7d} {R:8.3f} {ew:7.3f} {ew/max(R,1e-9):8.3f}")
print("  => err*w tracks R: small (clean, R~O(1)) to large (resonant, R~Sigma-e'). Bound err*w <~ c*R.")
print("\ndone.")
