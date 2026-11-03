from itertools import product
from collections import Counter
from math import gcd

def mul(z,w): return(z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0])
def alpha(v):
    z=(1,0)
    for p,e in zip(((2,1),(3,2)),v):
        for _ in range(abs(e)): z=mul(z,(p[0], p[1] if e>0 else -p[1]))
    return z

def allocations(E):
    return list(product(*[((-e,e) if e else (0,)) for e in E]))
def neg(v):return tuple(-x for x in v)
def triple(v):
    x,y=alpha(v)
    return(abs(x*x-y*y),2*abs(x*y),x*x+y*y)
profiles=list(product(range(4),repeat=2))
cases=outputs=0
for E,F in product(profiles,repeat=2):
    if E==F==(0,0): continue
    cases+=1;seen=set();out=Counter()
    for v,w in product(allocations(E),allocations(F)):
        key=min((v,w),(neg(v),neg(w)))
        if key in seen:continue
        seen.add(key)
        out[triple(tuple(x+y for x,y in zip(v,w)))]+=1
    for T,m in out.items():
        z=sum(e==f and e>0 and T[2]%p!=0 for p,e,f in zip((5,13),E,F))
        predicted=2**(z-(T[2]==1))
        if m!=predicted: raise RuntimeError((E,F,T,m,predicted))
        outputs+=1
print('Expanded fibers:',cases,'profile pairs;',outputs,'output multiplicities PASS')
# A direct amalgamation of the two old arc sets, without using the target formula.
def edges(r):
    return {(i,j) for i in range(r+2) for j in range(i+1,r+2) if (i,j)!=(0,r+1)}|{(r+1,0)}
for r,s in product(range(1,7),repeat=2):
    def f(x):return r+s+1 if x==r+1 else x
    def h(x):return 0 if x==0 else r+s+1 if x==s+1 else r+x
    glued={(f(a),f(b)) for a,b in edges(r)}|{(h(a),h(b)) for a,b in edges(s)}|{(a,b) for a in range(1,r+1) for b in range(r+1,r+s+1)}
    if glued!=edges(r+s): raise RuntimeError(('glue',r,s))
print('Independent retained-arc amalgamation: 36 cases PASS')
f=[0,1]
for _ in range(2004):f.append(f[-1]+f[-2])
for n in range(1,1001):
    C,D,K=f[2*n+1],f[2*n+3],f[2*n+2]
    if gcd(C,D)!=1 or D-C!=K: raise RuntimeError(('gcd',n))
    a,b=(f[n+1],f[n]),(f[n+2],-f[n+1])
    if mul(a,b)!=(K,(-1)**(n+1)):raise RuntimeError(('Cassini',n))
print('All n=1..1000 Cassini and norm-coprimality: PASS')
# Guard counterexample computed from the actual Gaussian powers.
z=(1,0);powers={}
for e in range(1,5):
    z=mul(z,(2,1));powers[e]=z
for e in (1,3,4):
    n,m=sorted(map(abs,powers[e]),reverse=True)
    print('guard',e,powers[e],(n+m,n-m))
