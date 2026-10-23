import json, math
from collections import defaultdict

def need(x):
    if not x: raise RuntimeError('failed invariant')
def orbit(x,y):
    return {(s*a,t*b) for a,b in ((x,y),(y,x)) for s in (-1,1) for t in (-1,1)}
def square(x):
    return x>=0 and math.isqrt(x)**2==x

def r2_formula(n):
    k=n; out=4;p=2
    while p*p<=k:
        e=0
        while k%p==0:k//=p;e+=1
        if p%4==3 and e%2:return 0
        if p%4==1:out*=e+1
        p+=1
    if k>1:
        if k%4==3:return 0
        if k%4==1:out*=2
    return out

shell=defaultdict(set)
for x in range(-100,101):
 for y in range(-100,101):
    if x*x+y*y<=10000:shell[x*x+y*y].add((x,y))
N=K=0; checks=0; first=None
for M in range(10001):
    pts=shell[M]; reps={(max(abs(x),abs(y)),min(abs(x),abs(y))) for x,y in pts}
    N+=len(pts); K+=len(reps)
    A=math.isqrt(M);D=math.isqrt(M//2)
    need(8*K==N+4*A+4*D+7);checks+=1
    if M:
        need(len(pts)==r2_formula(M));checks+=1
        need(8*len(reps)==len(pts)+4*square(M)+4*(M%2==0 and square(M//2)));checks+=1
        if first is None and len(reps)>1:first=M
    invariant={}
    for x,y in reps:
        u=x*x+y*y;v=x*x*y*y
        need((u,v) not in invariant);invariant[u,v]=(x,y)
        disc=math.isqrt(u*u-4*v)
        need(disc*disc==u*u-4*v)
        need((u+disc)//2==x*x and (u-disc)//2==y*y)
        need(len(orbit(x,y))==(1 if x==0 else 4 if y==0 or x==y else 8))
        checks+=3
need(first==25)
# Exact full-disk regeneration by primitive octant rays and D4.
for M in (0,1,2,5,25,100,999,10000):
    expected=set().union(*(shell[m] for m in range(M+1)))
    regenerated={(0,0)}
    for a in range(math.isqrt(M)+1):
      for b in range(a+1):
        if math.gcd(a,b)!=1 or a*a+b*b>M:continue
        for d in range(1,math.isqrt(M//(a*a+b*b))+1):regenerated |= orbit(d*a,d*b)
    need(regenerated==expected);checks+=1
# Rational rotation between every point in each small shell, exact numerators.
rotation_checks=0
for m in range(1,201):
 for x,y in shell[m]:
  for a,b in shell[m]:
    re=a*x+b*y;im=b*x-a*y
    need(re*re+im*im==m*m)
    need(re*x-im*y==a*m and im*x+re*y==b*m)
    rotation_checks+=1
# Norm-preserving mod3 golden->Gaussian map, all products.
def golden(x,y,p):
 a,b=x;c,d=y;return ((a*c+b*d)%p,(a*d+b*c+b*d)%p)
def gaussian(x,y,p):
 a,b=x;c,d=y;return ((a*c-b*d)%p,(a*d+b*c)%p)
def bridge(x):return ((x[0]+2*x[1])%3,x[1])
for a in range(3):
 for b in range(3):
  for c in range(3):
   for d in range(3):need(bridge(golden((a,b),(c,d),3))==gaussian(bridge((a,b)),bridge((c,d)),3));checks+=1
need(gaussian((1,1),(1,1),2)==(0,0))
need(all(golden(x,x,2)!=(0,0) for x in ((1,0),(0,1),(1,1))))
print(json.dumps(dict(status='PASS',disk_thresholds=10001,checks=checks,rational_rotation_pairs=rotation_checks,first_multiple_octant_shell=first,scope='D4 orbit reconstruction, radial decoder, quartic separator, local norm algebras; not a Bott or wiring equivalence'),indent=2))
