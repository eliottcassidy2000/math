#!/usr/bin/env python3
"""Third-party verifier: generic determinants; arbitrary finite direct-product groups."""
from itertools import combinations,permutations,product
from math import gcd,prod
from functools import reduce
from collections import Counter
import json
from pathlib import Path

def need(a,b,where):
    if a!=b: raise RuntimeError((where,a,b))

def det(M):
    n=len(M)
    return sum((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))*prod(M[i][p[i]] for i in range(n)) for p in permutations(range(n)))

def invfac(rows,width):
    D=[1]
    for k in range(1,min(len(rows),width)+1):
        delta=reduce(gcd,(abs(det([[rows[i][j] for j in cs] for i in rs])) for rs in combinations(range(len(rows)),k) for cs in combinations(range(width),k)),0)
        if not delta:break
        D.append(delta)
    return tuple(D[i]//D[i-1] for i in range(1,len(D)))

def arrangement():
    rows=list(product([1],range(-1,2),range(-1,2)))
    tally=Counter()
    cert=json.loads((Path(__file__).resolve().parent.parent/'graph-magic-degeneration'/'branch-certificate.json').read_text())['snf']
    need(cert['vectors'],[list(x) for x in rows],'input rows')
    for m in range(1<<len(rows)):
        sub=[r for i,r in enumerate(rows) if (m>>i)&1]
        ds=invfac(sub,3)
        need(ds,tuple(cert['subsets'][m]['invariants']),('record',m))
        tally[ds]+=(-1)**len(sub)
    expected={():1,(1,):-9,(1,1):20,(1,2):8,(1,1,1):1,(1,1,2):-10,(1,1,3):-4,(1,1,4):-4,(1,2,2):-3}
    need(dict(tally),expected,'SNF signed table')
    return {str(k):v for k,v in sorted(tally.items())}

def formula(N,h):
    return N**3-9*N*N+(20+8*h[2])*N+1-10*h[2]-4*h[3]-4*h[4]-3*h[2]**2

def direct_groups():
    groups=[(),(2,),(3,),(4,),(5,),(7,),(8,),(9,),(12,),(15,),(2,2),(2,3),(2,4),(2,6),(3,3),(3,9),(4,4),(5,5),(2,2,4),(2,4,4),(2,3,4),(2,2,2,2)]
    result=[]
    for orders in groups:
        els=list(product(*(range(k) for k in orders)))
        N=len(els);zero=(0,)*len(orders)
        h={d:sum(all((d*x)%n==0 for x,n in zip(a,orders)) for a in els) for d in [1,2,3,4,5,6]}
        count=0;injective_pairs=0
        for x,z in product(els,repeat=2):
            offsets={tuple((i*a+j*b)%n for a,b,n in zip(x,z,orders)) for i,j in product(range(-1,2),repeat=2)}
            injective_pairs+=len(offsets)==9
            for e in els:
                # Translate each distinct offset, independently enumerate every center.
                count+=all(any((a+b)%n for a,b,n in zip(e,v,orders)) for v in offsets)
        need(count,formula(N,h),('all-nonzero',orders))
        distinct_formula=N*N-4*(1+h[2])*N-2-4*h[2]+4*h[4]+2*h[5]+4*h[6]+2*h[2]**2+h[2]*h[4]
        need(injective_pairs,distinct_formula,('injective pair',orders))
        result.append({'cyclic_factors':orders,'N':N,'h':h,'all_nonidentity':count,'injective_pairs':injective_pairs})
    return result

def actual_trace(z):
    n=2**len(z)
    def mul(a,b):
        return a^b,prod(z[i] for i in range(len(z)) if (a&b)>>i&1)
    def tr(a):
        # Compute trace of the actual basis multiplication map, not a diagonal ansatz.
        return sum(c for b in range(n) for target,c in [mul(a,b)] if target==b)
    matrix=[]
    for a in range(n):
        matrix.append([c*tr(target) for b in range(n) for target,c in [mul(a,b)]])
    # Bareiss exact determinant.
    M=[r[:] for r in matrix];prev=1;sign=1
    for i in range(n-1):
        j=next((j for j in range(i,n) if M[j][i]),None)
        if j is None:return 0
        if j!=i:M[j],M[i]=M[i],M[j];sign=-sign
        pivot=M[i][i]
        for j in range(i+1,n):
            for k in range(i+1,n):M[j][k]=(M[j][k]*pivot-M[j][i]*M[i][k])//prev
            M[j][i]=0
        prev=pivot
    return sign*M[-1][-1]

def trace_tests():
    inputs=[(0,),(2,),(-2,3),(0,2,3),(2,3,5),(-2,3,4,7),(2,3,5,7,11)]
    out=[]
    for z in inputs:
        r=len(z);actual=actual_trace(z);expected=2**(r*2**r)*prod(z)**(2**(r-1))
        need(actual,expected,('actual multiplication traces',z))
        out.append({'z':z,'degree':2**r,'determinant':actual})
    return out

if __name__=='__main__':
    result={'signed_SNF':arrangement(),'groups':direct_groups(),'trace_checks':trace_tests(),'status':'PASS'}
    print(json.dumps(result,sort_keys=True,indent=2))
