#!/usr/bin/env python3
"""Full exact two-step and spectral audit of the105 matching boundary graph."""
from pathlib import Path
from fractions import Fraction as F
import sys,json
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031j'))
from boundary_matchings import matchings,components,need

def mm(A,B):
    columns=list(zip(*B))
    return [[sum(x*y for x,y in zip(row,col)) for col in columns] for row in A]
def eye(n):return [[int(i==j) for j in range(n)] for i in range(n)]
def shift(A,z):return [[v-z*int(i==j) for j,v in enumerate(row)] for i,row in enumerate(A)]
def main():
    states=list(matchings(tuple(range(8))));n=len(states)
    A=[[int(components(P,Q,8)==1) for Q in states] for P in states]
    need(n==105 and all(sum(r)==48 for r in A),'complete105-state48-regular graph')
    need(all(A[i][j]==A[j][i] for i in range(n) for j in range(n)),'symmetric kernel')
    A2=mm(A,A);values=sorted({x for r in A2 for x in r})
    need(values==[20,24,48] and all(A2[i][i]==48 for i in range(n)),'all-pairs exact two-step floors')
    # Verify the entire annihilating matrix, not just a five-cell quotient.
    roots=[48,4,-2,-6,-8];product=eye(n)
    for z in roots:product=mm(product,shift(A,z))
    need(all(v==0 for row in product for v in row),'full integer annihilator')
    powers=[eye(n)]
    for k in range(4):powers.append(mm(powers[-1],A))
    traces=[sum(P[i][i] for i in range(n)) for P in powers]
    multiplicities={}
    for z in roots:
        polynomial=[1];denominator=1
        for w in roots:
            if z==w:continue
            out=[0]*(len(polynomial)+1)
            for k,c in enumerate(polynomial):out[k]-=w*c;out[k+1]+=c
            polynomial=out;denominator*=z-w
        m=F(sum(c*t for c,t in zip(polynomial,traces)),denominator)
        need(m.denominator==1 and m>0,'integer spectral projector trace');multiplicities[z]=int(m)
    need(sum(multiplicities.values())==105 and multiplicities[48]==1,'simple stationary eigenspace')
    need(multiplicities=={48:1,4:56,-2:14,-6:14,-8:20},'exact spectral multiplicities')
    # Every permutation of the original labels induces a graph automorphism;
    # additionally check all rows have the same actual A2 multiplicity profile.
    hist=sorted((v,A2[0].count(v)) for v in values)
    need(all(sorted((v,row.count(v)) for v in values)==hist for row in A2),'source-specific two-step histogram')
    print(json.dumps({'status':'passed','states':n,'degree':48,'two_step_values':values,
          'two_step_row_histogram':hist,'annihilator_roots':roots,'multiplicities':multiplicities,'power_traces':traces,
          'normalized_two_step_floor':'5/576','uniform_minorization':'175/192',
          'l2_contraction_per_step':'1/6','scope':'Exact105-state graph on marked matchings, not a geometric selector or tournament-H coverage proof. Every row and every two-step entry is checked; the full integer matrix satisfies the displayed annihilator.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
