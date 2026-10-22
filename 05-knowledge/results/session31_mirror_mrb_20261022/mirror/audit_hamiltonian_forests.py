#!/usr/bin/env python3
"""Independent audit: parent-function forest enumeration; permutation H counts.
No determinant, matrix-forest theorem, subset DP, sympy, or shared worker code.
"""
import itertools, math, json
from fractions import Fraction as F
from pathlib import Path

rows=[]
for n in range(2,6):
    pairs=list(itertools.combinations(range(n),2)); size=1<<len(pairs)
    fs=[0]*n;hfs=[0]*n;h_sum=0;energy_sum=F(0);henergy_sum=F(0)
    for mask in range(size):
        out=[set() for _ in range(n)]
        for bit,(i,j) in enumerate(pairs):
            if mask>>bit&1:out[i].add(j)
            else:out[j].add(i)
        h=sum(all(b in out[a] for a,b in zip(p,p[1:])) for p in itertools.permutations(range(n)))
        forest=[0]*n
        for parents in itertools.product(*[[-1]+sorted(v) for v in out]):
            valid=True
            for start in range(n):
                seen=set();v=start
                while v!=-1:
                    if v in seen:valid=False;break
                    seen.add(v);v=parents[v]
                if not valid:break
            if valid:forest[sum(v!=-1 for v in parents)]+=1
        h_sum+=h
        energy=sum((F(len(v))-F(n-1,2))**2 for v in out)
        energy_sum+=energy;henergy_sum+=h*energy
        for j in range(n):fs[j]+=forest[j];hfs[j]+=h*forest[j]
    mean=[F(x,size) for x in fs]
    biased=[F(x*2**(n-1),size*math.factorial(n)) for x in hfs]
    predicted=[F(0)]*n
    # sum_s C(n-1-s,s)(t^2/4)^s(1+nt/2)^(n-1-2s)
    for s in range((n-1)//2+1):
        for j in range(n-2*s):
            predicted[j+2*s]+=F(math.comb(n-1-s,s),4**s)*math.comb(n-1-2*s,j)*F(n,2)**j
    assert mean==[math.comb(n-1,j)*F(n,2)**j for j in range(n)]
    assert biased==predicted
    delta=[b-a for a,b in zip(mean,biased)]
    assert delta[:2]==[0,0] and all(v>0 for v in delta[2:])
    assert F(h_sum*2**(n-1),size*math.factorial(n))==1
    cov_energy=henergy_sum*F(2**(n-1),size*math.factorial(n))-energy_sum/size
    assert cov_energy==-F(n-2,2)
    rows.append({'n':n,'tournaments':size,'mean_coefficients':[str(v) for v in mean],'biased_coefficients':[str(v) for v in biased],'covariance_coefficients':[str(v) for v in delta],'score_energy_covariance':str(cov_energy)})
result={'status':'PASS: independent exact parent-function forests and permutation Hamiltonian paths, all labelled tournaments n=2..5','rows':rows}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
