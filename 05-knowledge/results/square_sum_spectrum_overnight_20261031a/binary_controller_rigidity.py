#!/usr/bin/env python3
"""All six two-letter prefix languages: exact invariant and mod-four audit.

The integral realization is an inherited, independently verified certificate.
No field-rank argument is reduced modulo four: the full finite row orbit is
enumerated in the actual integral coordinates.
"""
import json
from itertools import combinations, product
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
PREVIOUS=HERE.parent/'square_sum_spectrum_overnight_20261030'
sys.path.insert(0,str(PREVIOUS))
from audit_feedback import require, gaps, subset_dp, tournament
from terminal_rigidity import RESTRICTED

data=json.loads((PREVIOUS/'integral_feedback.json').read_text())
matrices=data['integer_matrices']
initial=data['initial'][0]
final=[r[0] for r in data['final']]


def step(row,a,modulus=None):
    result=[sum(row[i]*matrices[a][i][j] for i in range(17)) for j in range(17)]
    return tuple(x%modulus if modulus else x for x in result)


def combine(first,second,sign=-1):
    result={w:c for w,c in first}
    for w,c in second:result[w]=result.get(w,0)+sign*c
    return [(w,c) for w,c in sorted(result.items()) if c]


def signed_row(terms):
    result=[0]*17
    for word,c in terms:
        row=initial
        for a in word:row=step(row,a)
        result=[x+c*y for x,y in zip(result,row)]
    return tuple(result)


def dot(row):return sum(x*y for x,y in zip(row,final))


def main():
    p13=[((3,),2),((1,3),-3),((1,1,3),1)]
    q13=[((1,),-2),((3,),2),((1,1),1),((1,3),-1),((3,1),2),
         ((3,3),-2),((1,3,1),-1),((1,3,3),1)]
    p23=[((3,),2),((2,3),-3),((2,2,3),1)]
    q23=[((2,),-2),((3,),2),((2,2),1),((2,3),-1),((3,2),2),
         ((3,3),-2),((2,3,2),-1),((2,3,3),1)]
    banks={'01':[((),1),((1,),-1),((1,0),-1),((1,1),1)],
           '02':[((),1),((2,),-1),((2,0),-1),((2,2),1)],
           '13':combine(q13,p13), '23':combine(p23,q23)}
    constants={}
    for alphabet,constant,terms in RESTRICTED:
        banks[alphabet]=terms;constants[alphabet]=constant
    receipts=[]
    for alphabet,terms in sorted(banks.items()):
        require(sum(c for _,c in terms)==0,'coefficient sum')
        require(all(set(w)<=set(map(int,alphabet)) for w,_ in terms),'prefix domain')
        D=signed_row(terms)
        require(all(step(D,a)==D for a in map(int,alphabet)),'prefix invariance')
        if alphabet in constants:
            require(all(step(D,a)==D for a in range(4)),'exact constant invariance')
            require(dot(D)==constants[alphabet]!=0,'nonzero exact constant')
            orbit=[]
        else:
            orbit=[tuple(x%4 for x in D)];seen=set(orbit);i=0
            while i<len(orbit):
                require(dot(orbit[i])%4==2,'nonvanishing terminal residue')
                for a in range(4):
                    new=step(orbit[i],a,4)
                    if new not in seen:seen.add(new);orbit.append(new)
                i+=1
                require(len(orbit)<=100000,'guard against incomplete finite search')
        for n in range(3):
            for A in product(range(4),repeat=n):
                literal=sum(c*subset_dp(tournament(w+A)) for w,c in terms)
                row=D
                for a in A:row=step(row,a)
                require(literal==dot(row),'independent signed prefix evaluation')
                require(literal==constants[alphabet] if alphabet in constants else literal%4==2,
                        'independent nonzero invariant')
        receipts.append({'alphabet':alphabet,'terms':terms,'signed_row':D,
                         'nonzero_constant':constants.get(alphabet),
                         'complete_mod4_row_orbit':orbit,'mod4_orbit_size':len(orbit),
                         'independent_suffix_controls':21})
    # Finite controls of the full neutral/non-neutral classification. The proof
    # of positive exponential debt is the odd-cycle injection in README.md.
    debt_controls=0
    suffixes=[w for n in range(3) for w in product(range(4),repeat=n)]
    for pair in combinations(range(4),2):
        alphabet=''.join(map(str,pair))
        for length in (1,2):
            for B in product(pair,repeat=length):
                for A in suffixes:
                    neutral=(alphabet=='01' and set(B)<={1} and set(A)<={1,3}) or (
                        alphabet in ('03','13','23') and set(B+A)<={3})
                    if alphabet=='01':letter=1
                    elif alphabet=='02':letter=2
                    elif alphabet=='12':letter=2 if 1 in B else 1
                    else:letter=3
                    for m in range(1,5):
                        difference=gaps((letter,)*m+B+A)-gaps((letter,)*m+A)
                        require(difference==0 if neutral else difference>=2**m,'insertion debt')
                        debt_controls+=1
                    if neutral:
                        for n in range(4):
                            for w in product(pair,repeat=n):
                                require(gaps(w+B+A)==gaps(w+A),'neutral arbitrary-word control')
    cross_chart_controls=0
    for n in range(5):
        for w in product((0,1),repeat=n):
            require(subset_dp(tournament(w+(3,1)))==3*subset_dp(tournament(w+(1,))),
                    'cross-chart hostile to dropping B membership')
            cross_chart_controls+=1
    print(json.dumps({'status':'passed','alphabets':receipts,
                      'debt_controls':debt_controls,
                      'independent_cross_chart_controls':cross_chart_controls,
                      'scope':'All suffixes by exact invariance or complete finite mod-four orbit. Inserted B must lie in the declared two-letter source alphabet.'},
                     indent=2,sort_keys=True))


if __name__=='__main__':main()
