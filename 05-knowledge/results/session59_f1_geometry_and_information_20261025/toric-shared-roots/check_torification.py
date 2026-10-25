#!/usr/bin/env python3
"""Exact characteristic-zero line-arrangement and root-surface checks."""
from collections import Counter
from functools import reduce
from itertools import combinations
from math import comb, gcd
from pathlib import Path
import json
import sympy as sp

NAMES=['E','E+x','E-x','E+z','E-z','E+x+z','E-x-z','E+x-z','E-x+z']
LINES=[(1,0,0),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),
       (1,1,1),(1,-1,-1),(1,1,-1),(1,-1,1)]


def normalize(point):
    g=reduce(gcd,point)
    assert g
    out=tuple(x//g for x in point)
    return out if next(x for x in out if x)>0 else tuple(-x for x in out)


def cross(a,b):
    return normalize((a[1]*b[2]-a[2]*b[1],
                      a[2]*b[0]-a[0]*b[2],
                      a[0]*b[1]-a[1]*b[0]))


def incidences(indices):
    points={cross(LINES[i],LINES[j]) for i,j in combinations(indices,2)}
    return {p:tuple(i for i in indices if sum(a*b for a,b in zip(LINES[i],p))==0)
            for p in sorted(points)}


def check_arrangement(indices):
    n=len(indices)
    lines=sp.Matrix([LINES[j] for j in indices])
    assert lines.rank()==3
    relations=sp.Matrix.hstack(*lines.T.nullspace()).T
    assert relations.shape==(n-3,n) and relations.rank()==n-3
    assert relations*lines==sp.zeros(n-3,3)
    data=incidences(indices)
    histogram=Counter(len(v) for v in data.values())
    assert all(m in [2,3] for m in histogram)
    assert sum(comb(m,2)*v for m,v in histogram.items())==comb(n,2)
    nodes=[]
    for p,incident in data.items():
        if len(incident)!=3:
            continue
        i,j,k=incident
        coefficients=list(sp.linsolve((sp.Matrix.hstack(sp.Matrix(LINES[i]),
                                                        sp.Matrix(LINES[j])),
                                      sp.Matrix(LINES[k]))))[0]
        alpha,beta=coefficients
        assert alpha and beta
        local_quadratic=sp.diag(-alpha,-beta,1)
        assert local_quadratic.rank()==3
        # One exact geometric root point; independent signs rescale nonzero
        # Jacobian columns and hence preserve its rank for all root lifts.
        values=[sum(a*b for a,b in zip(LINES[j],p)) for j in indices]
        roots=[sp.sqrt(value) for value in values]
        assert sum(v==0 for v in values)==3
        jacobian=2*relations*sp.diag(*roots)
        assert jacobian.rank()==n-4
        nodes.append(dict(point=p,incident=[NAMES[j] for j in incident],
                          local_alpha=str(alpha),local_beta=str(beta),
                          jacobian_rank=n-4,geometric_root_lifts=2**(n-4)))
    return dict(n=n,indices=indices,labels=[NAMES[j] for j in indices],
                multiplicities=dict(sorted(histogram.items())),
                triplepoints=nodes,total_nodes=histogram[3]*2**(n-4))


def census():
    output={}
    expected={5:{0:28,1:76,2:22},6:{0:2,1:20,2:46,3:16},
              7:{2:4,3:18,4:12,5:2},8:{4:1,5:4,6:4},9:{8:1}}
    for n in range(5,10):
        bytriple=Counter()
        for indices in combinations(range(9),n):
            multiplicities=Counter(len(v) for v in incidences(indices).values())
            assert max(multiplicities)<=3
            assert sp.Matrix([LINES[j] for j in indices]).rank()==3
            bytriple[multiplicities[3]]+=1
        assert sum(bytriple.values())==comb(9,n)
        assert dict(bytriple)==expected[n]
        output[n]=dict(sorted(bytriple.items()))
    return output


def invariants():
    z=sp.Symbol('z')
    output=[]
    for n in range(4,10):
        degree=2**(n-3)
        m=n-6
        # Independent Hilbert coefficient via alternating Koszul dimensions.
        pg=0 if m<0 else sum((-1)**j*comb(n-3,j)*comb(n-1+m-2*j,n-1)
                            for j in range(min(n-3,m//2)+1))
        if m>=0:
            series=sp.series((1-z*z)**(n-3)/(1-z)**n,z,0,m+1).removeO()
            assert int(sp.expand(series).coeff(z,m))==pg
        euler=2**(n-4)*(n*n-9*n+24)
        assert euler==12*(1+pg)-m*m*degree
        for t3 in range(comb(n,2)//3+1):
            t2=comb(n,2)-3*t3
            # Formal stratum calculation, independent of realizability of
            # these t2,t3 values by a line arrangement.
            singular=(2**(n-1)*(3-2*n+t2+2*t3)
                      +2**(n-2)*(2*n-2*t2-3*t3)
                      +2**(n-3)*t2+2**(n-4)*t3)
            assert singular+2**(n-4)*t3==euler
        output.append(dict(n=n,degree=degree,canonical_multiple=m,
                           canonical_square=m*m*degree,geometric_genus=pg,
                           resolution_euler=euler))
    assert [row['geometric_genus'] for row in output if row['n']>=5]==[0,1,7,31,111]
    return output


def main():
    chosen=[check_arrangement(tuple(range(n))) for n in [5,6,7,9]]
    assert [(r['n'],r['total_nodes']) for r in chosen]==[(5,4),(6,8),(7,24),(9,256)]
    # H- is the coordinate permutation z -> -z of H+.
    minus=check_arrangement((0,1,2,3,4,7,8))
    assert minus['total_nodes']==24
    full=incidences(tuple(range(9)))
    triplepoints={p for p,s in full.items() if len(s)==3}
    assert sum(p[0]==0 for p in triplepoints)==4
    assert {p for p in triplepoints if p[0]}=={(1,1,0),(1,-1,0),(1,0,1),(1,0,-1)}
    data=dict(status='PASS',scope='Characteristic-zero exact incidence and local Jacobian checks',
              selected_surfaces=chosen,hourglass_minus=minus,
              all_subset_census=census(),invariants=invariants(),
              full_arrangement=[dict(point=p,incident=[NAMES[j] for j in s])
                                for p,s in full.items()])
    Path(__file__).with_name('torification_results.json').write_text(json.dumps(data,indent=2)+'\n')
    print('PASS: exact 9-line arrangement; all 256 subsets of sizes 5--9;')
    print('5/6/7/9-root node counts 4/8/24/256, H- nodes24; local A1 and Jacobian ranks;')
    print('Hilbert/Koszul genus values, canonical squares, independent Euler/Noether identities')


if __name__=='__main__':
    main()
