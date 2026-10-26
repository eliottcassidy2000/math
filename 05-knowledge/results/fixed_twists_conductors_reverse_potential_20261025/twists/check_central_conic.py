# Central-line lifts for a fixed twist d over Q: necessary pair conditions -d1 d2, -d3 d4, -d5 d6 squares; then the conic
# w^2 = (d1 d5) (d3 s)^2 + (d3 d5) t^2 ... derived from x = d1 d3 s^2, z = t^2 parametrization.  Decide by Hilbert symbols.
from fractions import Fraction
from math import gcd, isqrt
import itertools, random
def factor(n):
    n=abs(n); f={}; p=2
    while p*p<=n:
        while n%p==0: f[p]=f.get(p,0)+1; n//=p
        p+=1
    if n>1: f[n]=f.get(n,0)+1
    return f
def legendre(a,p):
    a%=p
    if a==0: return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1
def hilbert(a,b,p):
    """Hilbert symbol (a,b)_p for nonzero integers a,b and prime p (standard formula)."""
    def v(n):
        c=0
        while n%p==0: n//=p; c+=1
        return c,n
    va,ua=v(a); vb,ub=v(b)
    if p!=2:
        eps=(-1)**(((p-1)//2)*va*vb)
        return eps*(legendre(ua,p)**vb)*(legendre(ub,p)**va)
    e=lambda u: ((u-1)//2)%2
    w=lambda u: ((u*u-1)//8)%2
    return (-1)**(e(ua)*e(ub)+va*w(ub)+vb*w(ua))
def hilbert_inf(a,b): return -1 if (a<0 and b<0) else 1
def conic_has_rational_point(alpha,beta):
    # w^2 = alpha s^2 + beta t^2 over Q, alpha,beta nonzero integers; Legendre/Hasse-Minkowski: solvable iff Hilbert symbols (alpha,beta)_v = 1 for all v
    if hilbert_inf(alpha,beta)==-1: return False
    primes=set(factor(2*alpha*beta))
    return all(hilbert(alpha,beta,p)==1 for p in primes)
def sqclass(q):
    q=Fraction(q); n=q.numerator*q.denominator; s=-1 if n<0 else 1; n=abs(n); r=1
    for p,e in factor(n).items():
        if e%2: r*=p
    return s*r
def central_lift_exists(d):
    d1,d2,d3,d4,d5,d6 = d[1:]
    pair_ok = all(sqclass(-a*b)==1 for a,b in [(d1,d2),(d3,d4),(d5,d6)])
    if not pair_ok: return False, 'pair condition fails'
    alpha = d1*d5*d3*d3; beta = d3*d5
    return conic_has_rational_point(alpha,beta), f'conic w^2={alpha}s^2+{beta}t^2'
def brute(d,H=40):
    pts=[]
    for x in range(-H,H+1):
        for z in range(1,H+1):
            if gcd(abs(x),z)!=1: continue
            vals=[0,x,-x,z,-z,x+z,-x-z]
            nz=[sqclass(di*v) for di,v in zip(d,vals) if v!=0]
            if len(nz)==6 and len(set(nz))==1: pts.append((x,z))
    return pts
random.seed(1)
tests=[(1,1,1,1,1,1,1),(1,1,-1,1,-1,1,-1),(1,1,-1,2,-2,1,-1),(1,3,-3,1,-1,3,-3),(1,1,-1,1,-1,7,-7),(1,5,-5,1,-1,5,-5),(1,2,-2,3,-3,6,-6),(1,1,-1,3,-3,3,-3)]
for d in tests:
    ex,why=central_lift_exists(d); pts=brute(d)
    print(d, '->', ex, why, '| brute height<=40 lifts on open central line:', len(pts), pts[:3])
# consistency on random twists: Hilbert decision vs brute-force (brute only certifies existence)
agree=0; tot=0; found_without=0
for _ in range(300):
    d=(1,)+tuple(random.choice([1,2,3,5,6,-1,-2,-3,-5,-6]) for _ in range(6))
    ex,why=central_lift_exists(d); pts=brute(d,25); tot+=1
    if pts and not ex: found_without+=1
    agree += (bool(pts)==ex)
print('random twists: Hilbert-decision agrees with brute height<=25 existence in', agree, 'of', tot, '; brute point found but decision False:', found_without)
