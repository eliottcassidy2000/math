#!/usr/bin/env python3
"""Exact midpoint-fan deficits and an outer-pair controller, not universal coverage."""
import json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031c'))
from audit_port_words import count
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261030'))
from audit_feedback import subset_dp

def weight(g):return 1<<max(g-1,1)
def debt(s):return sum(weight(t-u) for j,t in enumerate(s) for u in s[:j])
def mask(s):return sum(1<<u for u in s)
def patterns(d):
    if d==0:yield [0];return
    orbits=[]
    for a in range(1,d//2+1):orbits.append(sorted({a,d-a}))
    for bits in range(1<<len(orbits)):
        yield sorted({0,d}|{x for j,o in enumerate(orbits) if bits>>j&1 for x in o})

def fan(s,K):
    d=max(s)
    if (K-d)%2 or K<d+4:raise ValueError('translated interior domain')
    losing={x+(K-d)//2 for x in s}|{K}
    rows=[0]*(K+1)
    for i in range(1,K+1):
        for j in range(i+1,K+1):rows[i]|=1<<j
        if i in losing:rows[i]|=1
        else:rows[0]|=1<<i
    return rows,sorted(losing)

def main():
    limit=1000000;D=(limit+1).bit_length();catalog={};byparity=[set(),set()]
    first_collision=None;checked=0;recurrences=0;carry_checks=0
    for d in range(1,D+1):
        local={}
        for s in patterns(d):
            z=debt(s);A=mask(s);c=z-1
            if c<=limit:
                byparity[d%2].add(c);catalog.setdefault(c,[d,s])
            if z in local and first_collision is None:
                prior=local[z]
                first_collision={'diameter':d,'debt':z,'patterns':[prior,s],'masks':[mask(prior),A]}
            local.setdefault(z,s)
            for h in [1,2,5]:
                bigger=sorted([0,d+2*h]+[x+h for x in s])
                expected=z+(1<<h)*A+(1<<(d+2*h-1))+(2 if h==1 else 0)
                if debt(bigger)!=expected or mask(bigger)!=1+(1<<h)*A+(1<<(d+2*h)):
                    raise RuntimeError('outer pair recurrence')
                recurrences+=1
            if d<=9:
                for extra in [4,6,10]:
                    K=d+extra;rows,losing=fan(s,K);expected=(1<<(K-1))+1-z
                    if count(rows,[0])!=expected:raise RuntimeError('native adjacency fan count')
                    if K<=13 and subset_dp(rows)!=expected:raise RuntimeError('independent subset DP fan count')
                    checked+=1
        # Independent adjacent-address audit of the exact palindrome carry law.
        ordered=sorted((mask(s),debt(s),s) for s in patterns(d))
        for before,after in zip(ordered,ordered[1:]):
            old,new=set(before[2]),set(after[2]);a=min(new-old);w=d-2*a
            outside=old&new
            if w==0:
                increment=sum(weight(abs(x-a)) for x in outside)
            elif w==1:
                increment=2+sum(weight(abs(x-a))+weight(abs(x-(d-a))) for x in outside)
            else:
                increment=2
                for x in outside:
                    b=a-x if x<a else x-(d-a)
                    if b<1:raise RuntimeError('carry outside domain')
                    increment+=3*(1<<(b-1))+int(b==1)
            if after[1]-before[1]!=increment or increment<=0:
                raise RuntimeError('strict palindrome carry law')
            carry_checks+=1
    # A translated non-symmetric pattern carries an unavoidable exponential term.
    hostile=[]
    for s in [[0,1,3],[0,2,3],[0,1,4]]:
        d=max(s);A=mask(s);R=sum(1<<(d-x) for x in s)
        vals=[]
        for a in [2,3,4]:
            K=2*a+d;losing={x+a for x in s}|{K};rows=[0]*(K+1)
            for i in range(1,K+1):
                for j in range(i+1,K+1):rows[i]|=1<<j
                if i in losing:rows[i]|=1
                else:rows[0]|=1<<i
            H=count(rows,[0]);expected=(1<<(K-1))+1+(1<<(a-1))*(A-R)-debt(s)
            if H!=expected:raise RuntimeError('unbalanced cluster term')
            vals.append([K,H,H-(1<<(K-1))])
        hostile.append({'pattern':s,'A_minus_reverse':A-R,'samples':[v for v in vals]})
    if 9 in catalog:raise RuntimeError('false claimed first absent positive odd deficit')
    firstmissing=[c for c in range(1,100,2) if c not in catalog]
    print(json.dumps({'status':'passed','complete_deficit_limit':limit,'complete_diameter_limit':D,
         'distinct_positive_deficits':len(catalog),'by_K_parity':[len(x) for x in byparity],
         'first_missing_below100':firstmissing,'first_patterns':[[c,*catalog[c]] for c in sorted(catalog)[:30]],
         'actual_adjacency_controls':checked,'outer_pair_recurrence_controls':recurrences,
         'strict_palindrome_carry_controls':carry_checks,
         'diameter_erasure_hostile':{'same_debt':2,'patterns':[[0,1],[0,2]],'after_outer_gap_one_debts':[14,22]},
         'first_same_debt_diameter_collision':first_collision,'unbalanced_hostiles':hostile,
         'scope':'Complete fixed-deficit midpoint-cluster grammar, not complete tournaments; exponent k=K-1 has opposite parity to the normalized diameter.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
