from itertools import permutations
from math import gcd
p=31
def v(x): return (0,1) if x==p else (1,x)
def step(x):
 a,b=v(x); a,b=(a+b)%p,a
 return p if not a else b*pow(a,-1,p)%p
sq={x*x%p for x in range(1,p)}
def H(xs):
 V=list(map(v,xs)); E=[[((a*d-b*c)%p in sq) for c,d in V] for a,b in V]
 return sum(all(E[i][j] for i,j in zip(w,w[1:])) for w in permutations(range(len(V))))

def periods(n):
    a,b,c,d=1,0,0,1
    proj=0
    for k in range(1,1000):
        a,b,c,d=((a+b)%n,a,(c+d)%n,c)
        if not proj and b==c==0 and a==d and gcd(a,n)==1:
            proj, scalar=k,a
        if (a,b,c,d)==(1,0,0,1):
            return k,proj,scalar
    raise ValueError(n)
for n, expect in {2:(3,3,1),3:(8,4,2),5:(20,5,3),15:(40,20,11),30:(120,60,11),31:(30,30,1),61:(60,15,11)}.items():
    got=periods(n)
    if got!=expect: raise ValueError((n,got))
    print(n,got)
start=(0,1,2,3,4)
xs=start
a=[]
for _ in range(30):
    a.append(H(xs))
    xs=tuple(map(step,xs))
if xs!=start: raise ValueError(xs)
if a != [13,15,5,3,13,9,9,15,5,5,5,15,9,3,5,5,5,5,5,9,5,5,15,3,5,5,9,13,3,5]:
    raise ValueError(a)
if any(all(a[i]==a[(i+k)%30] for i in range(30)) for k in (1,2,3,5,6,10,15)):
    raise ValueError("shorter period")
print("Exact path-count period 30:", a)
