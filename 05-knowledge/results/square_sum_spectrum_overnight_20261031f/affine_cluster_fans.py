#!/usr/bin/env python3
"""Polynomial coefficients of fixed affine apex clusters; symbolic hostile controls."""
import itertools,json,sys
from fractions import Fraction as F
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031e'))
from centered_fans import debt,count

def two(k):return F(1<<k) if k>=0 else F(1,1<<(-k))
def polynomial(groups):
    M=max(groups);out=[F(0)]*(M+1);out[0]=1-sum(debt(s) for s in groups.values())
    for j,S in groups.items():out[j]+=sum(two(b-1) for b in S)
    for u,S in groups.items():
        for v,T in groups.items():
            if u<v:out[v-u]-=sum(two(t-s-1) for s in S for t in T)
    return out

def adjacency(groups,a):
    losing={j*a+b for j,S in groups.items() for b in S};K=max(losing)
    if min(losing)<2:raise ValueError('apex must beat first chain vertex')
    out=[0]*(K+1)
    for i in range(1,K+1):
        for j in range(i+1,K+1):out[i]|=1<<j
        if i in losing:out[i]|=1
        else:out[0]|=1<<i
    return out

def main():
    alphabet=[]
    for mask in range(8):alphabet.append([b for j,b in enumerate([-1,0,1]) if mask>>j&1])
    tested=accepted=physical=0;by_slope={};examples=[]
    for M in range(1,6):
        for D in (0,1):
            for seq in itertools.product(alphabet,repeat=M-1):
                groups={j+1:S for j,S in enumerate(seq) if S};groups[M]=[D]
                coeff=polynomial(groups);valid=all(v==0 for v in coeff[1:M])
                support=set(groups)
                predicted=(support=={M} or (M%2==0 and support=={M//2,M} and set(groups[M//2])=={D-b for b in groups[M//2]}))
                if valid!=predicted:raise RuntimeError('affine-cluster rigidity classification')
                tested+=1
                if valid:
                    accepted+=1;by_slope[M]=by_slope.get(M,0)+1
                    for a in (4,5):
                        rows=adjacency(groups,a);value=sum(v*two(a*j) for j,v in enumerate(coeff))
                        if count(rows,[0])!=value:raise RuntimeError('accepted cluster adjacency')
                        physical+=1
                    examples.append({'groups':groups,'deficit':str(-coeff[0]),'exponent_slope':M,'exponent_offset':D-1})
    # Three nonempty slope classes necessarily leave an uncancelled coefficient.
    hostile={1:[0],2:[0],3:[0]};c=polynomial(hostile)
    if c!=[1,F(-1,2),F(0),F(1,2)]:raise RuntimeError('minimal three-slope hostile')
    for a in (3,4,5):
        if count(adjacency(hostile,a),[0])!=sum(v*two(a*j) for j,v in enumerate(c)):
            raise RuntimeError('hostile literal adjacency')
    print(json.dumps({'status':'passed','complete_profiles':tested,'accepted_profiles':accepted,'accepted_by_slope':by_slope,
                      'actual_adjacency_controls':physical+3,'accepted_examples':examples,
                      'three_slope_hostile':{'groups':hostile,'polynomial_in_2_to_a':[str(v) for v in c]},
                      'scope':'Finite coefficient census verifies the all-size support-difference and reflection proof; it is not a census of all one-apex tournaments.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
