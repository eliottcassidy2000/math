from fractions import Fraction as F
from math import sin, pi, gcd, sqrt
import itertools, random
LAM=F(1,14)
def hhat(n):
    if n==0: return float(2*LAM)
    return sin(2*pi*n*float(LAM))/(pi*n)
def delta(A,N):
    tot=0.0
    for n in itertools.product([x for x in range(-N,N+1) if x!=0], repeat=len(A)):
        if sum(ni*ai for ni,ai in zip(n,A))==0:
            p=1.0
            for ni in n: p*=hhat(ni)
            tot+=p
    return tot

print("(3) THE VANISHING AT MULTIPLES OF 7  [7 = 1/(2*lam)]")
print("    hhat(n) = sin(pi n/7)/(pi n) = 0 exactly when 7 | n:")
for n in [1,2,6,7,8,13,14,21]:
    print(f"      hhat({n:2d}) = {hhat(n):+.8f}   {'<-- ZERO' if n%7==0 else ''}")
print("    => every resonance vector with ANY coordinate divisible by 7")
print("       contributes NOTHING.  The modulus 7 is built into the transform.")

print()
print("(4) HOW THE 3-BODY TERM DECAYS WITH SPEED SIZE")
print("    scale    family              d123        (2lam)^3     ratio")
tl=float(2*LAM)
for A in [(2,3,5),(11,13,17),(31,37,41),(101,103,107)]:
    d3=delta(A,45)
    print(f"    {max(A):5d}    {str(A):18s} {d3:+.6f}   {tl**3:.6f}   {d3/tl**3:+8.2f}")

print()
print("(5) THE SHORT-VECTOR BOUND  |delta| <= sum prod min(2lam, 1/(pi|n_i|))")
print("    tested against the k=2 case, where the answer is known (THM-965):")
print("    (a,b)      |delta| exact   shortest-vector term g^2/(pi^2 ab)")
for (a,b) in [(2,3),(3,5),(5,7),(4,6),(6,9),(10,15)]:
    g=gcd(a,b); d=delta((a,b),300)
    sv=g*g/(pi*pi*a*b)
    print(f"    {str((a,b)):10s} {abs(d):.6f}       {sv:.6f}   ratio {abs(d)/sv:.2f}")
print()
print("    => the k=2 deviation IS the shortest-vector term up to a constant,")
print("       which is exactly the sawtooth floor g^2/(4ab) of THM-1025.")
