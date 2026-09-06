#!/usr/bin/env python3
"""Formal cycle defects before geometric averaging; exact degree-nine audit.

No geometry census is repeated. Frozen triple-event multiplicities are
independently weighted by formal logarithms of cyclic-run subset polynomials.
"""
from collections import Counter,defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import json
import sys

sys.stdout.reconfigure(newline='\n')
ROOT=Path(__file__).resolve().parents[1]
LIMIT=9
GATES=0


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def weight(monomial):
    return sum(x//3 for x in monomial)


def clean(P):
    return {m:Q(c) for m,c in P.items() if c}


def add(P,Qpoly,scale=1):
    result=defaultdict(Q,P)
    for m,c in Qpoly.items():
        result[m]+=scale*c
    return clean(result)


def times(P,Qpoly):
    result=defaultdict(Q)
    for a,x in P.items():
        for b,y in Qpoly.items():
            if weight(a)+weight(b)<=LIMIT:
                result[tuple(sorted(a+b))]+=x*y
    return clean(result)


def scaled(P,c):
    return clean({m:c*x for m,x in P.items()})


def homogeneous(P,k):
    return {m:c for m,c in P.items() if weight(m)==k}


def logarithm(P):
    need(P.get((),0)==1,'formal logarithm has constant one')
    U=add(P,{():Q(1)},-1)
    power={():Q(1)}
    result={}
    for k in range(1,LIMIT+1):
        power=times(power,U)
        result=add(result,power,Q((-1)**(k+1),k))
    return result


def cycle_bank(L):
    need(L>=4 and L%2==0,'simple bipartite cycle')
    bank=Counter()
    for mask in range(1<<L):
        if mask.bit_count()>LIMIT:
            continue
        if not mask:
            sig=()
        elif mask==(1<<L)-1:
            sig=(3*L,)
        else:
            tokens=[]
            for i in range(L):
                if mask>>i&1 and not(mask>>((i-1)%L)&1):
                    m=1
                    while mask>>((i+m)%L)&1:
                        m+=1
                    tokens.append(3*m+(0 if m%2 else 2 if i%2==0 else 1))
            sig=tuple(sorted(tokens))
        need(weight(sig)==mask.bit_count(),'cyclic-run component grading')
        bank[sig]+=1
    return clean(bank)


def pretty(P):
    def name(token):
        m,r=divmod(token,3)
        return 'C'+str(m) if m%2==0 and r==0 else 'p'+str(m)+('L' if r==1 else 'R' if r==2 else '')
    return [(('*'.join(name(x) for x in m)) or '1',str(c)) for m,c in sorted(P.items())]


def main():
    banks={L:cycle_bank(L) for L in (4,6,8,10,12,14,16)}
    logs={L:logarithm(P) for L,P in banks.items()}
    bulk=scaled(logs[10],Q(1,10))
    for L in (12,14,16):
        need(scaled(logs[L],Q(1,L))==bulk,('formal long-cycle bulk control through weight9',L))
    need(homogeneous(bulk,1)=={(3,):Q(1)},'one-edge bulk coefficient')
    defects={L:add(logs[L],bulk,-L) for L in (4,6,8)}
    for L,D in defects.items():
        need(all(weight(m)>=L for m in D),('cycle-defect weight onset',L))
    D4=homogeneous(defects[4],4)
    D5=homogeneous(defects[4],5)
    expected={(12,):1,(3,3,3,3):1,(3,3,7):-2,(3,3,8):-2,
              (3,9):4,(7,7):1,(8,8):1,(13,):-2,(14,):-2}
    need(D4==expected,'explicit signed four-cycle defect')
    D8=homogeneous(defects[8],8)
    D9=homogeneous(defects[8],9)
    for name,P in (('D4',D4),('D5',D5),('D8',D8),('D9',D9)):
        need(sum(P.values())==0,('forgetting components and retaining only edge count kills cycle defects',name))
    p1={(3,):Q(1)}
    T5=scaled(add(D5,scaled(times(p1,D4),8)),Q(1,4))
    expected_T={(3,3,9):1,(3,7,8):-1,(3,12):1,(3,13):-1,(3,14):-1,
                (7,9):1,(8,9):1,(15,):-1}
    need(T5==expected_T,'simplified next-order signed cycle defect')
    n=8
    four_squared=times(D4,D4)
    quadratic=add(scaled(four_squared,Q(1,2)),
                  add(scaled(times(p1,four_squared),n),times(D4,D5)))
    need(quadratic==add(scaled(four_squared,Q(1,2)),scaled(times(D4,T5),4)),
         'n8 cancellation in the quadratic copy coefficient')
    eight=add(D8,add(D9,scaled(times(p1,D8),2*n)))
    direct_eight=scaled(add(times(banks[8],banks[8]),banks[16],-1),Q(1,2))
    direct_four=scaled(add(add(times(times(banks[4],banks[4]),times(banks[4],banks[4])),
                                 times(banks[4],banks[12]),-4),banks[16],3),Q(1,12))
    need(eight==direct_eight,'whole polynomial eight-cycle coefficient versus contrast')
    need(quadratic==direct_four,'whole polynomial quadratic four-cycle coefficient versus contrast')
    data_path=ROOT/'05-knowledge/results/overnight2_20260906_no3line_third_certificates.json'
    data=json.loads(data_path.read_text())
    totals=[Q(0),Q(0)]
    grouped=defaultdict(lambda:[0,0,Q(0),Q(0)])
    sign_counts=Counter()
    geometry_weights={}
    for row in data['profiles']:
        s=row['signature']
        sig=tuple(sorted(3*m+(1 if r>c else 2 if c>r else 0) for m,r,c,cycle in s))
        need(eight.get(sig,0)==Q(row['profile_coefficients_A_B_D_E_F'][3]),'full eight-cycle profile coefficient')
        need(quadratic.get(sig,0)==Q(row['profile_coefficients_A_B_D_E_F'][4]),'full quadratic profile coefficient')
        geometry=Q(6*row['unordered_event_triples'],row['grid_copies'])
        geometry_weights[sig]=geometry/6
        contribution=[geometry*eight.get(sig,0),geometry*quadratic.get(sig,0)]
        need(contribution==list(map(Q,row['ordered_weighted_contributions_E_F'])),'geometric weighting normalization')
        totals=[a+b for a,b in zip(totals,contribution)]
        key=(weight(sig),tuple(sorted(m for m,r,c,cycle in s if cycle)))
        grouped[key][0]+=1
        grouped[key][1]+=row['unordered_event_triples']
        for i,x in enumerate(contribution):
            grouped[key][i+2]+=x
            sign_counts[(i,(x>0)-(x<0))]+=1
    need(totals==[Q(172483,529200),Q(11881,50400)],'full inherited third-moment coefficients')
    need(quadratic[(3,3,3,3,3,3,7)]<0 and quadratic[(12,12)]>0,
         'signed square has positive and negative monomial coefficients')
    # Ordinary/factorial cumulant subtractions contain only M1 and M2,
    # which have cycle weights at most six; they cannot alter either result.
    need(totals[0]>0 and totals[1]>0,'both coefficients survive in third cumulants')
    def W(P):
        return sum((c*geometry_weights.get(m,Q(0)) for m,c in P.items()),Q(0))
    need(3*W(four_squared)==Q(456371,2116800),'degree8 geometry-weighted square')
    need(24*W(times(D4,T5))==Q(42631,2116800),'degree9 geometry-weighted cross term')
    hostile={(3,3,7):Q(1),(7,7):Q(-1)}
    hostile_square=times(hostile,hostile)
    need(W(hostile_square)==Q(-397,529200),'geometry averaging is not positive on formal squares')
    needed_weights=[geometry_weights[(3,3,3,3,7,7)],geometry_weights[(3,3,7,7,7)],
                    geometry_weights[(7,7,7,7)]]
    need(needed_weights==[Q(16747,1058400),Q(2749,235200),Q(1,147)],'three forest-type hostile values')
    need(8*totals[1]==Q(11881,6300),'component-additivity hostile contrast for third cumulants')
    print('FORMAL_BULK long-cycle lengths10,12,14,16 agree through weight9')
    print('FOUR_CYCLE_DEFECT_DEGREE4',pretty(D4))
    print('FOUR_CYCLE_DEFECT_DEGREE5',pretty(D5))
    print('SIMPLIFIED_NEXT_DEFECT_T5',pretty(T5))
    print('POLYNOMIAL_COUNTS bulk',len(bulk),'D4',len(D4),'D5',len(D5),
          'D8',len(D8),'D9',len(D9),'eight',len(eight),'quadratic',len(quadratic))
    print('MOMENT_AND_THIRD_CUMULANT_COEFFICIENTS',list(map(str,totals)))
    print('QUADRATIC_WEIGHTED_SQUARE_AND_CROSS',str(3*W(four_squared)),str(24*W(times(D4,T5))))
    print('NEGATIVE_FORMAL_SQUARE Q=p2L*(p1^2-p2L)',str(W(hostile_square)),
          'three_nonnegative_geometry_weights',list(map(str,needed_weights)))
    print('CUMULANT_ADDITIVITY_CONTRAST 2C8 - 2*(2C4+C8) + 4C4',str(8*totals[1]))
    print('SIGNED_PROFILE_COUNTS',sorted(sign_counts.items()))
    for key,values in sorted(grouped.items()):
        print('ACTUAL_CYCLE_GROUP',key,values[:2],list(map(str,values[2:])))
    print('GEOMETRY_CERTIFICATE_SHA256',sha256(data_path.read_bytes()).hexdigest())
    digest=sha256(repr((sorted(D4.items()),sorted(D5.items()),sorted(D8.items()),sorted(D9.items()),
                        sorted(eight.items()),sorted(quadratic.items()))).encode()).hexdigest()
    print('FORMAL_SEMANTIC_SHA256',digest)
    print('PASS_OPTIMIZATION_LIVE_GATES',GATES)


if __name__=='__main__':
    main()
