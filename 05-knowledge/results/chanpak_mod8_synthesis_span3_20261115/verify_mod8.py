#!/usr/bin/env python3
"""Replay for the mod-8 dictionary (SYNTHESIS.md). Pure Python + sympy, about a minute (Paley counts need ./hp2)."""
import math, subprocess, os
from sympy import jacobi_symbol, primerange, isprime, n_order
def perm_sign(f,m):
    seen=[False]*m; s=1
    for x0 in range(m):
        if seen[x0]: continue
        L=0; x=x0
        while not seen[x]: seen[x]=True; x=f(x); L+=1
        if L%2==0: s=-s
    return s
# 1. Frobenius-Zolotarev for the three Collatz-side multipliers and their composite
bad=sum(1 for m in range(3,400,2) for a in (-1,2,-2,3,-3,6) if math.gcd(a,m)==1 and perm_sign(lambda x:(a*x)%m,m)!=jacobi_symbol(a%m,m))
print('1. sign(x -> a x on Z/m) = (a/m) for odd m < 400, a in {-1, 2, -2, 3, -3, 6}: %d violations'%bad)
# 2. the Collatz step and its 3x-1 mirror mod m are conjugate by negation; their sign is (6/m)
bad=0
for m in range(5,400,2):
    if m%3==0: continue
    inv2=pow(2,-1,m)
    Tp=lambda x:((3*x+1)*inv2)%m; Tm=lambda x:((3*x-1)*inv2)%m
    if any(Tm(x)!=(-Tp((-x)%m))%m for x in range(m)): bad+=1
    if perm_sign(Tp,m)!=jacobi_symbol(6,m) or perm_sign(Tm,m)!=jacobi_symbol(6,m): bad+=1
print('2. T-(x) = -T+(-x) and sign(T+/-) = (6/m) on Z/m, odd m < 400 prime to 3: %d violations'%bad)
# 3. Paley tournaments: x -> a x is an automorphism iff (a/p) = 1 and a converse-isomorphism iff (a/p) = -1
bad=0
for p in primerange(3,80):
    if p%4!=3: continue
    QR={(x*x)%p for x in range(1,p)}
    arc=lambda i,j:((j-i)%p) in QR
    for a in (2,-1,-2,3,-3):
        if a%p==0: continue
        chi=jacobi_symbol(a%p,p); ok=True
        for i in range(p):
            for j in range(p):
                if i==j: continue
                same=arc((a*i)%p,(a*j)%p)==arc(i,j)
                if (chi==1)!=same: ok=False
        bad+=not ok
print('3. Paley P_p (p = 3 mod 4, p < 80): multiplication by a in {2,-1,-2,3,-3} preserves arcs iff (a/p) = 1, reverses all iff -1: %d violations'%bad)
# 4. the mod-8 characters and the Collatz multiplier
units=[1,3,5,7]
chi={'-1':lambda n:1 if n%4==1 else -1,'2':lambda n:1 if n%8 in (1,7) else -1,'-2':lambda n:1 if n%8 in (1,3) else -1}
ker={k:[u for u in units if f(u)==1] for k,f in chi.items()}
print('4. kernels in (Z/8)*: chi_-1 %s, chi_2 %s, chi_-2 %s; <3> = %s; <-1,3> = %s'%(ker['-1'],ker['2'],ker['-2'],sorted({pow(3,k,8) for k in range(4)}),sorted({(s*pow(3,k,8))%8 for s in (1,7) for k in range(4)})))
bad=0
for p in primerange(3,5000):
    rep=any(int(math.isqrt(p-2*y*y))**2==p-2*y*y for y in range(0,math.isqrt(p//2)+1))
    if rep!=(p%8 in (1,3)): bad+=1
print('   p = x^2 + 2y^2 <=> p mod 8 in <3> = {1,3}, p < 5000: %d violations'%bad)
# 5. Collatz and 3x-1 valuations by n mod 8, and the 4n+1 / 4n-1 merges
def v2(x):
    v=0
    while x%2==0: x//=2; v+=1
    return v
for sgn,name in ((1,'3n+1'),(-1,'3n-1')):
    tab={r:sorted({min(v2(3*n+sgn),3) for n in range(r,4000,8)}) for r in (1,3,5,7)}
    print('5. v2(%s) by n mod 8 (3 means >= 3): %s'%(name,tab))
ok=all(((3*(4*n+1)+1)==4*(3*n+1)) and ((3*(4*n-1)-1)==4*(3*n-1)) for n in range(1,10000))
print('   merges 3(4n+1)+1 = 4(3n+1) and 3(4n-1)-1 = 4(3n-1) for n < 10^4: %s'%ok)
# 6. Artin density for base 2 in the non-residue classes 3 and 5 mod 8
for c in (3,5):
    ps=[p for p in primerange(3,200000) if p%8==c]
    print('6. p = %d mod 8, p < 2e5: 2 is a primitive root for %.4f of %d primes (2 x Artin constant = 0.7479)'%(c,sum(n_order(2,p)==p-1 for p in ps)/len(ps),len(ps)))
# 7. Paley Hamiltonian paths: |Aut| = p(p-1)/2 divides H, and diagonal joins lie in 3 mod 8
if os.path.exists('./hp2'):
    ps=[3,7,11,19]; lines=[]
    for p in ps:
        QR={(x*x)%p for x in range(1,p)}
        lines.append('%d %s'%(p,' '.join(str(sum(1<<j for j in range(p) if j!=i and ((j-i)%p) in QR)) for i in range(p))))
    out=subprocess.run(['./hp2'],input='\n'.join(lines)+'\n',capture_output=True,text=True).stdout.split('\n')
    for p,r in zip(ps,out):
        h,f=map(int,r.split()); P2=f-2*h; N=(h+P2)**2+2*h*h
        print('7. P_%d: H = %d (H / (p(p-1)/2) = %s), P2 = %d, diagonal join N mod 8 = %d, N prime: %s'%(p,h,h//(p*(p-1)//2) if h%(p*(p-1)//2)==0 else 'NOT DIVISIBLE',P2,N%8,isprime(N)))
