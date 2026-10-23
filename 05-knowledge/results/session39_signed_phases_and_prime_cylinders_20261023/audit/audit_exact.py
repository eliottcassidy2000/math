"""Independent Session39 exact checks. Standard library only."""
from itertools import product, combinations, permutations
from fractions import Fraction as F
from collections import Counter
import json

def need(x, why):
    if not x: raise RuntimeError(why)

def compose(f,g):
    a,b=f;c,d=g
    return a*c%3,(a*d+b)%3
def inv(f):
    a,b=f
    return a,(-a*b)%3
def act(f,x):return (f[0]*x+f[1])%3
aff=list(product((1,2),range(3)))
edges=[(0,1),(0,2),(0,3),(1,2),(2,3)]
counts=Counter()
for gs in product(aff,repeat=5):
    charts=[(1,0)]+list(gs[:3])
    hs=[compose(inv(charts[j]),compose(g,charts[i])) for (i,j),g in zip(edges[3:],gs[3:])]
    roots=[x for x in range(3) if all(act(h,x)==x for h in hs)]
    actual=[xs for xs in product(range(3),repeat=4) if all(xs[j]==act(g,xs[i]) for (i,j),g in zip(edges,gs))]
    expected=[tuple(act(t,x) for t in charts) for x in roots]
    need(sorted(actual)==sorted(expected),'fixed-point gluing')
    need((len(roots)==3)==all(h==(1,0) for h in hs),'frame versus points')
    counts['gluing_instances']+=1
    counts['root_count_'+str(len(roots))]+=1

def valuation(n):
    need(n!=0,'valuation of zero')
    return (abs(n)&-abs(n)).bit_length()-1
def iterate(n,s,k):
    w=[]
    for _ in range(k):
        z=3*n+s;a=valuation(z);w.append(a);n=z//(2**a)
    return tuple(w),n
for k in range(1,5):
 for w in product(range(1,4),repeat=k):
    A=sum(w)
    B=sum(3**(k-1-i)*2**sum(w[:i]) for i in range(k))
    M=2**(A+1)
    for s in (-1,1):
        residue=((2**A-s*B)*pow(3**k,-1,M))%M
        # Independent exhaustive residue calculation, no imported code.
        hits=[n for n in range(1,M,2) if iterate(n,s,k)[0]==w]
        need(hits==[residue],'signed cylinder uniqueness')
        D=2**A-3**k
        if B%D==0 and s*D>0:
            n=s*B//D
            need(iterate(n,s,k)==(w,n),'closed word sufficiency')
        counts['signed_words']+=1

K=((3,4),(4,3));Ki=((F(-3,7),F(4,7)),(F(4,7),F(-3,7)))
def dot(x,A,y):return sum(x[i]*A[i][j]*y[j] for i in range(len(x)) for j in range(len(y)))
def q(x,c):return (dot(x,Ki,x)-dot(x,Ki,c))/2%1
for x in product(range(-3,4),repeat=2):
 for z in product(range(-2,3),repeat=2):
    y=tuple(x[i]+sum(K[i][j]*z[j] for j in range(2)) for i in range(2))
    need(q(x,(3,3))==q(y,(3,3)),'characteristic refinement descent')
    counts['characteristic_shifts']+=1
need(q((0,0),(0,0))!=q((3,4),(0,0)),'naive odd K failure')
need(dot((1,1),Ki,(1,1))==F(2,7),'filling')
need(not any(u*u%512==511 for u in range(1,256,2)),'opposite cyclic refinement')
print(json.dumps(dict(counts),indent=2,sort_keys=True))
