from itertools import product
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from profile_tools import hp_all, strong

states=tuple(product(range(3),range(3),range(2)))
def step(D,x,y):
    R={}
    for (a,b,c),z in D.items():
        cells=[(x+1,b,c)]
        if b or x: cells.append((a,y+1,c))
        if c or y: cells.append((a,b,1))
        for w in cells:R[w]=R.get(w,0)+z
    return R

for x in (0,1):
    for start in states:
        D={start:1}
        for k in (0,1,2,3):
            E=step(D,x,x)
            D={w:E.get(w,0)-k*D.get(w,0) for w in states}
        assert not any(D.values())

def runs(S,z):
    P={(0,0,0):1}
    for x,y in S:
        P=step(P,x,y)
    return sum(d for (a,b,c),d in P.items() if a<2 and b<2 and (b or z))

def fast(S,z):
    return runs(S,z)+runs([(b,a) for a,b in S],1-z)

def masks(S,z):
    n=len(S)
    A=[(1<<n)-(1<<(j+1)) for j in range(n)]+[0,0]
    for j,(x,y) in enumerate(S):
        if x:A[n]|=1<<j
        else:A[j]|=1<<n
        if y:A[n+1]|=1<<j
        else:A[j]|=1<<(n+1)
    A[n if z else n+1]|=1<<(n+1 if z else n)
    return tuple(A)

def blocked(p,q):
    return [(1,0),(0,1)]+[(0,0)]*p+[(1,1)]*q+[(0,1),(1,0)]

def formula(p,q):
    return 2*3**(p+q)+6*(2**p*3**q+3**p*2**q)+14*2**(p+q)+2*(3**p+3**q)-1
for p in range(1,4):
    R=[]
    for q in range(1,4):
        S=blocked(p,q);A=masks(S,1)
        assert strong(A)
        h=hp_all(A)[-1]
        assert h==formula(p,q)==fast(S,1)
        R.append(h)
    print('independent DP interpolation row',p,R)
for p,q in product(range(1,21),repeat=2):
    assert fast(blocked(p,q),1)==formula(p,q)

def lower(b):
    return 2*3**b+4*sum(min(3**l*2**(b-1-l),2**l*3**(b-1-l)) for l in range(b))+4*2**b+1
for b in (1,2):
    m=10**20; n=b+4; count=0
    for right in (0,1):
        a,c=n-2+right,n-1-right
        for X in range(1<<n):
            if not X&1 or (X>>a)&1:continue
            for Y in range(1<<n):
                if not (Y&2) or (Y>>c)&1:continue
                S=[((X>>j)&1,(Y>>j)&1) for j in range(n)]
                for z in (0,1):
                    h=hp_all(masks(S,z))[-1]
                    assert h==fast(S,z) and h>=lower(b)
                    m=min(m,h);count+=1
    assert m=={1:41,2:105}[b]
    print('anchored core',b,'literal cases',count,'min',m,'all-order floor',lower(b))
assert (lower(2),lower(3))==(51,143)
print('two-hub paired-frame debt and exact balanced family: PASS')
