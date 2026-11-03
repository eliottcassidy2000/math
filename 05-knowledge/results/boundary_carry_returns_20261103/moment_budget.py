#!/usr/bin/env python3
"""Physical triangular controls for the matching-cap/path-cover moment budget."""
from pathlib import Path
from math import comb,factorial
import json
from scout import endpoint,clone,require


def witness(d):
    n=2*d+1;rows=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if i==0:forward=j<=d
            elif j<=d or i>d:forward=True
            else:forward=j==d+i
            a,b=(i,j)if forward else(j,i);rows[a]|=1<<b
    return tuple(rows)


def covers(t,k):
    h=[sum(r)for r in endpoint(t)];full=len(h)-1
    previous=[0]*(full+1);previous[0]=1
    for _ in range(k):
        current=[0]*(full+1)
        for mask in range(1,full+1):
            sub=mask
            while sub:
                current[mask]+=h[sub]*previous[mask^sub];sub=(sub-1)&mask
        previous=current
    return previous[full]


def basis(r,m):
    return sum((-1)**(r-j)*comb(r,j)*j**m for j in range(r+1))


if __name__=='__main__':
    checks=[]
    for d,K in ((0,3),(1,3),(2,3)):
        dimension=d+K;triangular=[]
        for j in range(1,dimension+1):
            root_matching=min(d,j-1);k=j-root_matching;t=witness(root_matching)
            values=[covers(clone(t,0,m),k)for m in range(1,dimension+2)]
            cs=[]
            for m in range(1,dimension+1):
                residual=values[m-1]-sum(c*basis(r,m)for r,c in enumerate(cs,1))
                require(residual%factorial(m)==0,'integral skeleton coefficients')
                cs.append(residual//factorial(m))
            require(cs[j-1]>0 and not any(cs[j:]),'positive triangular realization')
            require(sum(c*basis(r,dimension+1)for r,c in enumerate(cs,1))==values[-1],'extra cover size')
            triangular.append({'matching':root_matching,'components':k,'rows':t,'coefficients':cs,'values':values})
        determinant=1
        for i,row in enumerate(triangular):determinant*=row['coefficients'][i]
        checks.append({'matching_cap':d,'component_cap':K,'dimension':dimension,
                       'coefficient_determinant':determinant,'witnesses':triangular})
    Path(__file__).with_name('moment_budget.json').write_text(json.dumps(checks,indent=2)+'\n')
    for r in checks:print(r['matching_cap'],r['component_cap'],r['dimension'],r['coefficient_determinant'])
