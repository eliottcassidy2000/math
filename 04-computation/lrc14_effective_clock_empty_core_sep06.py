#!/usr/bin/env python3
"""Exact controls for the gcd-aware clock repair and ten-subset gcd cap.

Standard library only. The all-height result is proved in the companion
note; the bounded clock and arc controls below are independently typed.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import gcd, lcm
from hashlib import sha256

CHECKS = 0


def need(value, label):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(label)


def distance(x):
    r = x % 1
    return min(r, 1-r)


def safe_union_measure(speeds):
    intervals = []
    for v in speeds:
        for k in range(v+1):
            a, b = max(Q(0), Q(14*k-1,14*v)), min(Q(1), Q(14*k+1,14*v))
            if a < b:
                intervals.append((a,b))
    merged = []
    for a,b in sorted(intervals):
        if not merged or a > merged[-1][1]:
            merged.append([a,b])
        else:
            merged[-1][1] = max(merged[-1][1],b)
    return 1-sum((b-a for a,b in merged),Q(0))


def safe_wall_measure(speeds):
    walls = {Q(0),Q(1)}
    for v in speeds:
        for k in range(v+1):
            for s in (-1,1):
                p=Q(14*k+s,14*v)
                if 0 <= p <= 1:
                    walls.add(p)
    points=sorted(walls)
    return sum((b-a for a,b in zip(points,points[1:])
                if all(distance(v*(a+b)/2) >= Q(1,14) for v in speeds)),Q(0))


def beta(q):
    return Q((q+6)//7,q)


def arc_max(q):
    events={Q(0),Q(1)}
    for j in range(q):
        for sign in (-1,1):
            events.add((Q(sign,14)-Q(j,q)) % 1)
    points=sorted(events)
    probes=points+[(a+b)/2 for a,b in zip(points,points[1:])]
    return max(sum(distance(t+Q(j,q)) < Q(1,14) for j in range(q)) for t in probes)


def main():
    r=tuple(range(1,13))
    primitive=tuple(sorted([2*x for x in r]+[13]))
    hostile=tuple(sorted([4*x for x in r]+[26]))
    values=[]
    for name,speeds in [('core',r),('primitive',primitive),('hostile',hostile)]:
        a,b=safe_union_measure(speeds),safe_wall_measure(speeds)
        need(a==b,('independent interval engines',name))
        values.append(a)
        print('measure',name,a)
    body,base,actual=values
    need(body==Q(6617,194040) and base==actual==body/2,'recovered exact measures')
    false_floor=body*Q(3,4)
    corrected_floor=body*(1-beta(2))
    need(actual<false_floor and actual==corrected_floor,'gcd-aware repair is sharp')
    need(gcd(26,4)==2 and 26%4!=0,'actual nonmultiple detuning')
    print('false_coprime_floor',false_floor,'corrected_floor',corrected_floor)
    # Complete arc-event probes include strict endpoints and every open cell.
    for q in range(1,71):
        need(arc_max(q)==(q+6)//7,('sharp open arc count',q))
        if q>=3:
            need(beta(q)<=Q(1,3) and (beta(q)==Q(1,3))==(q==3),('order3 extremal',q))
    print('arc_universe q=1..70, endpoints and all event cells, PASS')
    words=[]
    for c in range(2,421):
        orders=[q for q in range(1,c+1) if c%q==0 and c//q<=2]
        for word in combinations_with_replacement(orders,3):
            pairs=[lcm(word[i],word[j]) for i,j in [(0,1),(0,2),(1,2)]]
            if pairs==[c,c,c] and sum(map(beta,word),Q(0))>=1:
                words.append((c,word))
    expected=[(2,(1,2,2)),(2,(2,2,2)),(3,(3,3,3)),(4,(2,4,4))]
    need(words==expected,'complete bounded scalar/hereditary controls')
    print('order_universe c=2..420, all sorted divisor words, hereditary pair-lcm and eleven-gcd filters')
    print('surviving_words',words)
    # A phase-only hostile for the c=2 residual: all chosen body speeds
    # have the same safe residue, but the two named odd tails kill both lifts.
    y=Q(3,8)
    body_control=tuple(1+8*j for j in range(10))
    tails=(3,5,7)
    need(all(distance(v*y)>=Q(1,14) for v in body_control),'body-safe phase control')
    masks=[sum((1<<j) for j in range(2) if distance(t*(y+j)/2)<Q(1,14)) for t in tails]
    need(masks[0]|masks[1]|masks[2]==3,'safe body phase with both sheets killed')
    print('phase_only_hostile body=(1+8j:j=0..9), y=3/8, tails=(3,5,7), masks',masks)
    print('SCOPE: clocks/phase budgets only; no claimed unsafe speed row')
    print('semantic_sha256',sha256(repr((values,words,masks)).encode()).hexdigest())
    print('explicit_checks',CHECKS)


if __name__=='__main__':
    main()
