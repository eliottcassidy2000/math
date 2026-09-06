#!/usr/bin/env python3
"""Exact nine-subset gcd-profile classifier; no body-height census.

The companion proves c<=9 analytically. This source independently classifies
all c=2..12 by gcd words and by effective-order words, then audits divisor
absorption, inherited exact ten-pack signatures, and explicit safe realizations.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
from math import gcd, lcm
import json

CHECKS = 0
EXPECTED = {
    2: [(1,1,2,2),(1,2,2,2),(2,2,2,2)],
    3: [(1,3,3,3),(3,3,3,3)],
    4: [(1,2,4,4),(1,4,4,4),(2,2,4,4),(2,4,4,4),(4,4,4,4)],
    6: [(2,3,3,6),(2,3,6,6),(2,6,6,6),(3,3,6,6)],
    8: [(2,4,8,8),(2,8,8,8),(4,4,8,8),(4,8,8,8),(8,8,8,8)],
    9: [(3,9,9,9)]}
TEN = {2: {(1,2,2),(2,2,2)}, 3: {(3,3,3)}, 4: {(2,4,4)}}


def need(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def beta(q):
    return F((q+6)//7,q)


def distancesafe(v,x):
    z = v*x % 1
    return min(z,1-z) >= F(1,14)


def gcd_route(c):
    answer = []
    for gs in combinations_with_replacement([g for g in range(1,5) if c%g==0],4):
        if any(gcd(*pair)>2 for pair in combinations(gs,2)):
            continue
        if any(gcd(*triple)!=1 for triple in combinations(gs,3)):
            continue
        qs = tuple(sorted(c//g for g in gs))
        if sum(map(beta,qs),F(0)) >= 1:
            answer.append(qs)
    return sorted(answer)


def order_route(c):
    answer = []
    for qs in combinations_with_replacement([q for q in range(1,c+1) if c%q==0],4):
        if any(c//q>4 for q in qs):
            continue
        if any(c//lcm(*pair)>2 for pair in combinations(qs,2)):
            continue
        if any(lcm(*triple)!=c for triple in combinations(qs,3)):
            continue
        if sum(((q+6)//7)*(c//q) for q in qs) >= c:
            answer.append(qs)
    return answer


def divisor_exits(c,qs):
    gs = [c//q for q in qs]
    exits = []
    for d in range(2,c+1):
        if c%d:
            continue
        residual = [g for g in gs if g%d]
        cap = sum(gcd(d,g)*((d//gcd(d,g)+6)//7) for g in residual)
        if cap < d:
            exits.append((d,cap))
    return exits


def ten_exit(c,qs):
    gs = [c//q for q in qs]
    for i,g in enumerate(gs):
        if g==1:
            continue
        orders = tuple(sorted(g//gcd(g,other) for j,other in enumerate(gs) if j!=i))
        if orders not in TEN[g]:
            return True
    return False


def main():
    found = {}
    for c in range(2,13):
        a,b = gcd_route(c),order_route(c)
        need(a==b,('independent full gcd/order enumeration',c))
        need(a==EXPECTED.get(c,[]),('literal complete profile bank',c))
        if a:
            found[c]=a
    need(sum(map(len,found.values()))==20,'twenty inherited-only profiles')
    # Positive beta boundary controls, with the all-q tail proved in the note.
    need([q for q in range(4,9) if beta(q)==F(1,4)] == [4,8], 'quarter equality controls')
    need(all(beta(q)<=F(1,4) for q in range(4,9)), 'quarter bound finite base')
    need(F(9+6,7*9)<F(1,4), 'monotone beta tail begins at9')
    removed = []
    realizations = []
    for c,rows in found.items():
        for qs in rows:
            exits = divisor_exits(c,qs)
            need(bool(exits)==ten_exit(c,qs),'all-divisor exit agrees with exact ten-pack signature')
            if exits:
                removed.append((c,qs,exits))
            gs = [c//q for q in qs]
            body = [c*j for j in range(1,10)]
            tails = [g*(1+c*(1+14*(i+1))) for i,g in enumerate(gs)]
            physical = body+tails
            need(len(set(physical))==13 and min(physical)>0,'thirteen distinct positive realization speeds')
            need(gcd(*physical)==1,'primitive realization')
            need([gcd(c,w) for w in tails]==gs,'exact tail gcd coordinates')
            need(all(gcd(*(body+[tails[i] for i in ids]))<=2 for ids in combinations(range(4),2)),
                 'actual eleven-subset gcd filters')
            need(all(gcd(*(body+[tails[i] for i in ids]))==1 for ids in combinations(range(4),3)),
                 'actual twelve-subset gcd filters')
            witness = F(1,14*c)
            need(all(distancesafe(v,witness) for v in physical),'explicit safe full-row witness')
            need([w*witness % 1 for w in tails] == [F(g,14)+F(1,14*q) for g,q in zip(gs,qs)],
                 'literal tail phase formula')
            realizations.append([c,qs,tails,str(witness),bool(exits)])
            print('profile',c,qs,'budget',sum(map(beta,qs),F(0)),
                  'divisor_exits',exits,'safe_realization_tails',tails,'witness',witness)
    need(removed==[(4,(1,4,4,4),[(4,3)]),(8,(2,8,8,8),[(4,3)])], 'two strict absorption exits')
    need(sum(map(len,found.values()))-len(removed)==18,'eighteen profiles after absorption')
    # Boundary exclusions used to reduce the analytic c<=12 estimate to c<=9.
    need(F(4,5)<1 and F(8,11)<1 and F(1,3)+F(1,4)+2*F(1,6)<1,
         'clocks10,11,12 excluded with strict margins')
    # Actual phase-cover boundary at the maximal residual clock. The complete
    # speed row is safe elsewhere, so this is not an unsafe-row claim.
    core = tuple(range(1,18,2))
    tails = (1,5,6,7)
    killed = [tuple(j for j in range(9) if not distancesafe(t,(F(1,2)+j)/9))
              for t in tails]
    need(killed == [(0,8),(3,5),(1,4,7),(2,6)],'literal clock9 cap-tight partition')
    need(all(distancesafe(v,F(1,2)) for v in core),'strictly safe chosen body phase')
    need(all(distancesafe(v,F(1,14)) for v in tuple(9*v for v in core)+tails),
         'same full-row has an explicit safe witness elsewhere')
    print('clock9_phase_hostile core=(1,3,...,17), tails=(1,5,6,7), y=1/2, killed',killed)
    print('clock9_hostile full-row safety witness x=1/14; no unsafe row claimed')
    print('PROVED bound from companion: every nine-subset gcd belongs to {1,2,3,4,6,8,9}')
    print('FINITE universe: c=2..12; all sorted inherited gcd/order quadruples')
    print('inherited_profiles',20,'after_all_divisor_absorption',18)
    print('all20 gcd profiles have explicit safe primitive realizations; no unsafe row claimed')
    print('semantic_sha256',sha256(json.dumps([realizations,killed],separators=(',',':')).encode()).hexdigest())
    print('explicit_checks',CHECKS)


if __name__=='__main__':
    main()
