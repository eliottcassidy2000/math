from math import factorial
from itertools import permutations, product
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from profile_tools import hp_all, strong

def host(k,p,q,edges):
    beta=p+q;n=2*k+beta;G=[0]*(n+k)
    for i in range(n):
        for j in range(i+1,n):G[i]|=1<<j
        for j in range(k):
            beat=(i==j) if i<k else (i>=k+p if i<k+beta else i-k-beta!=j)
            if beat:G[n+j]|=1<<i
            else:G[i]|=1<<(n+j)
    r=0
    for i in range(k):
        for j in range(i+1,k):
            if edges>>r&1:G[n+i]|=1<<(n+j)
            else:G[n+j]|=1<<(n+i)
            r+=1
    return G

def full_paths(k,b):
    n=2*k+b
    for perm in permutations(range(k)):
        for assign in product(range(k+1),repeat=b):
            R=[[] for _ in range(k+1)]
            for t,j in enumerate(perm):
                R[t].append(k+b+j)
                R[t+1].append(j)
            for j,t in enumerate(assign):R[t].append(k+j)
            H=[]
            for t,r in enumerate(R):
                H.extend(sorted(r))
                if t<k:H.append(n+perm[t])
            yield tuple(H)
for k in (1,2,3):
    for b in range(1,5):
        G=host(k,b//2,b-b//2,0)
        assert strong(G)
        P=set(full_paths(k,b))
        assert len(P)==factorial(k)*(k+1)**b
        assert all(len(set(w))==len(G) and all(G[i]>>j&1 for i,j in zip(w,w[1:])) for w in P)
        L=factorial(k)*k*sum(min((k+1)**l*k**(b-1-l),k**l*(k+1)**(b-1-l)) for l in range(b))
        h=hp_all(G)[-1]
        assert h>=len(P)+L
        print(k,b,h,len(P),L)
print('many-hub k factorial chamber controls: PASS')
