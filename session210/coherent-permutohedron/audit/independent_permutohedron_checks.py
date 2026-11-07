#!/usr/bin/env python3
from independent_carrier_checks import tournament, hp, require
from itertools import permutations, product, combinations
from collections import Counter
from math import comb, factorial
import json
from pathlib import Path
ROOT=Path(__file__).parent

def faces(n):
    for k in range(1,n+1):
        for assign in product(range(k),repeat=n):
            if set(assign)==set(range(k)):
                yield tuple(sum(1<<v for v,b in enumerate(assign) if b==i) for i in range(k))

def energy(a,p):return sum(not(a[p[i]]>>p[j]&1) for i in range(len(p)) for j in range(i+1,len(p)))

def main():
    out={'orders':[]};first_bad=None
    for n in range(1,6):
        ps=list(permutations(range(n)));fs=list(faces(n));masks=[(1<<n)-1]
        for bits in range(1<<(n*(n-1)//2)):
            a=tournament(n,bits);energies={p:energy(a,p) for p in ps};sink=[];rhs=[0]*n
            for p in ps:
                f=sum(a[p[i]]>>p[i+1]&1 for i in range(n-1))
                for d in range(f+1):rhs[d]+=comb(f,d)
                is_sink=True
                for i in range(n-1):
                    q=p[:i]+(p[i+1],p[i])+p[i+2:]
                    require(abs(energies[q]-energies[p])==1,('gradient',n,bits,p,q))
                    if energies[q]<energies[p]:is_sink=False
                require(is_sink==(f==n-1),('sink bijection',n,bits,p))
                if is_sink:sink.append(p)
            h=[hp(a,m) for m in range(1<<n)];lhs=[0]*n
            for face in fs:
                prod=1
                for b in face:prod*=h[b]
                lhs[n-len(face)]+=prod
            require(lhs==rhs,('face incidence',n,bits,lhs,rhs))
            require(len(sink)==h[-1],('sink count',n,bits))
            if first_bad is None:
                low=min(energies.values())
                for p in sink:
                    if energies[p]>low:
                        first_bad={'n':n,'mask':bits,'rows':a,'hamilton_path':p,'path_energy':energies[p],'global_minimum':low};break
        out['orders'].append({'n':n,'labelled_tournaments':1<<(n*(n-1)//2),'ordered_partition_faces':len(fs)})
    out['first_nonglobal_sink']=first_bad
    rows=[]
    for n in range(3,9):
        m=n-2;a=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
        a[0]^=1<<(n-1);a[n-1]|=1
        actual=Counter(energy(a,p) for p in permutations(range(n)) if all(a[u]>>v&1 for u,v in zip(p,p[1:])))
        expected=Counter({1:1})
        for mask in range(1<<m):
            A=[i+1 for i in range(m) if mask>>i&1];B=[i+1 for i in range(m) if not mask>>i&1]
            inv=sum(x>y for x in A for y in B);p=tuple(A+[n-1,0]+B)
            require(energy(a,p)==m+inv,('extreme energy',n,A))
            expected[m+inv]+=1
        require(actual==expected,('extreme poly',n,actual,expected))
        rows.append({'n':n,'sink_energy_histogram':dict(sorted(actual.items()))})
    out['extreme_reversal_controls']=rows
    # RREF free slots, uniqueness as actual finite-field subspaces, and Gaussian count.
    gs=[]
    for q in (2,3):
        for m in range(1,5):
            seen=set();dims=Counter()
            for k in range(m+1):
                for A in combinations(range(1,m+1),k):
                    B=[b for b in range(1,m+1) if b not in A];piv=sorted(m-a for a in A)
                    slots=[(piv.index(m-a),m-b) for a in A for b in B if a>b]
                    for values in product(range(q),repeat=len(slots)):
                        mat=[[0]*m for _ in piv]
                        for r,p in enumerate(piv):mat[r][p]=1
                        for (r,c),v in zip(slots,values):mat[r][c]=v
                        span=frozenset(tuple(sum(coef[r]*mat[r][c] for r in range(k))%q for c in range(m)) for coef in product(range(q),repeat=k))
                        require(span not in seen,('duplicate RREF span',q,m,A,values));seen.add(span);dims[k]+=1
            gaussian=[[1]]
            for n in range(1,m+1):
                prev=gaussian[-1];gaussian.append([1]+[prev[k-1]+q**k*prev[k] for k in range(1,n)]+[1])
            require([dims[k] for k in range(m+1)]==gaussian[m],('Gaussian',q,m,dims,gaussian[m]))
            gs.append({'q':q,'m':m,'subspaces_by_dimension':[dims[k] for k in range(m+1)]})
    out['rref_controls']=gs;out['status']='PASS'
    (ROOT/'independent_permutohedron_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
