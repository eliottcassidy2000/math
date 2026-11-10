"""Small named controls for Session248's new birth identities; stdlib only."""
from itertools import permutations
from collections import defaultdict
from fractions import Fraction
import json

def aut(t):
    n=len(t)
    return [p for p in permutations(range(n))
            if all(t[i][j]==t[p[i]][p[j]] for i in range(n) for j in range(i))]
def key(t):
    n=len(t)
    return min(tuple(t[p[i]][p[j]] for i in range(n) for j in range(i))
               for p in permutations(range(n)))
def delete(t,x):
    return [[t[i][j] for j in range(len(t)) if j!=x]
            for i in range(len(t)) if i!=x]
def birth(t,s):
    n=len(t)
    return [list(t[i])+[int(i not in s)] for i in range(n)]+[[int(i in s) for i in range(n)]+[0]]
def carousel(q):
    return [[int(0<(j-i)%q<=q//2) for j in range(q)] for i in range(q)]

out=[]
for name,t in [('TT2',[[0,1],[0,0]]),('C3',carousel(3)),
               ('punctured_R5',delete(carousel(5),0)),('R5',carousel(5))]:
    n=len(t); g=aut(t); kt=key(t); classes={}; counts=defaultdict(int); ratios=set()
    for mask in range(1<<n):
        s={i for i in range(n) if mask>>i&1}; u=birth(t,s); a=aut(u)
        k=sum({p[i] for i in s}==s for p in g)
        r=len(g)//k; o=len({p[n] for p in a})
        assert len(a)*r==len(g)*o
        assert o%2==r%2==1
        ku=key(u); classes[ku]=(u,len(a)); counts[ku]+=1
        ratios.add(str(Fraction(o,r)))
    for ku,(u,au) in classes.items():
        d=sum(key(delete(u,x))==kt for x in range(n+1))
        assert counts[ku]*au==d*len(g)
        pi=Fraction(__import__('math').factorial(n),(1<<(n*(n-1)//2))*len(g))
        pj=Fraction(__import__('math').factorial(n+1),(1<<(n*(n+1)//2))*au)
        assert pi*Fraction(counts[ku],1<<n)==pj*Fraction(d,n+1)
    out.append(dict(parent=name,automorphisms=len(g),births=1<<n,
                    child_classes=len(classes),symmetry_ratios=sorted(ratios)))
print(json.dumps(dict(status='PASS',controls=out),indent=2))
