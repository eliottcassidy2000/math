#!/usr/bin/env python3
"""Separate exact controls for TOP_SUFFIX.md; no large graph census."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def values(bits):
    while bits:
        low=bits & -bits
        yield low.bit_length()-1
        bits ^= low


def step(bits):
    result=bits
    for a in values(bits):
        result |= (bits & ((1<<a)-1)) << a
    return result


def numeric_closure(seed,cap):
    present=[False]*(cap+1)
    bits=0
    for n in range(1,cap+1):
        present[n]=n in seed or any(present[a] and present[n-a]
                                    for a in range(1,(n+1)//2))
        if present[n]:
            bits |= 1<<n
    return bits,present


def run():
    interval_checks=0
    for B in range(4,17):
        L,U=B-2,B-1
        for p in range(1,9):
            for q in range(1,9):
                left=range(p*L-1,p*U+2)
                right=range(q*L-1,q*U+2)
                actual={x+y for x in left for y in right if x != y}
                target=set(range((p+q)*L-1,(p+q)*U+2))
                require(target <= actual,'balanced interval inclusion')
                if p == q:
                    require(target == actual,'equal interval exact strict sum')
                else:
                    require(actual == set(range((p+q)*L-2,(p+q)*U+3)),
                            'distinct interval full sum including overlap')
                interval_checks += 1

    cap=16384
    seeds=[set(x) for x in itertools.combinations(range(1,9),2)]
    seeds += [{1,3,8},{2,3,4},{4,7,11},{6,10,15},{6,9,20},
              {4,5,6,7},{3,4,5},{5,6,7},{6,8,10},{4,14,22}]
    histories=[]
    for seed in seeds:
        d=math.gcd(*seed)
        normalized={x//d for x in seed}
        closure,present=numeric_closure(normalized,cap)
        holes=[n for n in range(1,cap+1) if not present[n]]
        require(all(present[-256:]),'cap extends beyond all normalized holes in controls')
        genus=len(holes)
        conductor=max(holes,default=0)+1
        bits=sum(1<<x for x in normalized)
        records=[]
        first_certificate=None
        for t in range(15):
            M=bits.bit_length()-1
            require(M <= cap,'explicit synchronous target cap')
            saturated=bits == closure & ((1<<(M+1))-1)
            top3=all(bits>>x & 1 for x in [M-2,M-1,M]) if M >= 3 else False
            top4=all(bits>>x & 1 for x in [M-3,M-2,M-1,M]) if M >= 4 else False
            records.append({'t':t,'maximum':M,'size':bits.bit_count(),
                            'saturated':saturated,'top3':top3,'top4':top4})
            if saturated and top3 and M >= 2*conductor and first_certificate is None:
                first_certificate=t
            if first_certificate is not None:
                require(saturated,'saturation persists after full state beyond twice conductor')
                require(bits.bit_count() == M-genus,'eventual exact cardinality')
                if t > first_certificate:
                    require(M == 2*records[-2]['maximum']-1,'closed maximum recurrence')
                    require(bits.bit_count() == 2*records[-2]['size']+genus-1,
                            'closed cardinality recurrence')
                if t == first_certificate+2:
                    break
            bits=step(bits)
        require(first_certificate is not None,'selected seeds reach saturation within scope')
        histories.append({'original_seeds':sorted(seed),'gcd':d,'normalized_holes':holes,
                          'normalized_conductor':conductor,
                          'first_saturation_top3_and_twice_conductor':first_certificate,
                          'generations':records})

    hostile=next(x for x in histories if x['original_seeds'] == [2,3])
    require([r['maximum'] for r in hostile['generations']] == [3,5,8,15,28,55,109],
            'transient saturation maximum control')
    require([r['saturated'] for r in hostile['generations']] ==
            [True,True,True,False,True,True,True],
            'early finite saturation need not persist')
    stronger=next(x for x in histories if x['original_seeds'] == [5,6,7])
    require([r['saturated'] for r in stronger['generations'][:3]] == [True,True,False]
            and all(r['top3'] for r in stronger['generations']),
            'early finite saturation plus stable top3 still need not persist')
    s567=sum(1<<x for x in [5,6,7])
    require(not(step(step(s567))>>21 & 1) and (step(step(step(s567)))>>21 & 1),
            'literal missing-21 witness and third-generation repair')

    # A fixed exact top window, rather than the exponentially growing state.
    # From generation 1 onward the window lower edge exceeds the prior max,
    # so retained old vertices cannot enter. All needed sum parents remain.
    long_controls=[]
    for A in [2,3,7,25,100]:
        B=A+2
        initial={A,A+1,A+2}
        first=initial | {a+b for a,b in itertools.combinations(initial,2)}
        M=2*B-1
        W=B+1
        D={M-x for x in first if M-x <= W}
        first_top4=1 if 3 in D else None
        for t in range(1,A):
            require({0,1,2} <= D,'persistent top3 in exact deficit window')
            if first_top4 is None and 3 in D:
                first_top4=t
            if t == A-1:
                break
            # At A=2 this loop has already terminated before this condition
            # is needed. For A>=3 it holds from generation 1 onward.
            require(M-1 > W,'no retained old value enters fixed top window')
            D={u+v-1 for u,v in itertools.combinations(sorted(D),2)
               if 0 <= u+v-1 <= W}
            M=2*M-1
        require(first_top4 == A-1,'three-consecutive-seed first top-four time')
        long_controls.append({'seeds':sorted(initial),'fixed_window_width':W,
                              'first_top4_generation':first_top4,
                              'maximum_at_top4':M,'final_deficits':sorted(D)})

    return {'status':'FINITE-EXACT controls for the all-seed saturation proof',
            'universes':{'intervals':'B=4..16, p,q=1..8',
                         'seed_sets':'all 28 pairs from1..8 and10 named multi-seed controls',
                         'literal_target_cap':cap,'maximum_generation_cap':14,
                         'long_transient':'five named consecutive-three seed sets; exact fixed top window'},
            'interval_addition_checks':interval_checks,'seed_histories':histories,
            'long_top_four_transients':long_controls,
            'first_finite_saturation_implies_persistence':False}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    data=json.dumps(run(),indent=2,sort_keys=True)+'\n'
    target=Path(__file__).with_name('top_suffix_receipt.json')
    if args.write:
        target.write_text(data)
    else:
        require(target.read_text() == data,'receipt drift')
    print(json.dumps({'passed':True,'receipt_sha256':hashlib.sha256(data.encode()).hexdigest()},sort_keys=True))
