#!/usr/bin/env python3
"""Exact body/event interface controls; no repository mathematical imports."""
from fractions import Fraction as F
from math import gcd
from hashlib import sha256
import sys

CHECKS=0


def need(value,label):
    global CHECKS
    CHECKS+=1
    if not value:
        raise RuntimeError(label)


def distance(x):
    r=x%1
    return min(r,1-r)


def blocked(w,c):
    D=14*w
    return {k for k in range(w)
            if min((c*(14*k+3))%D,(-c*(14*k+3))%D)<w}


def blocker_count(w,c):
    g=gcd(w,c)
    q,h=w//g,c//g
    return g*((q-1-3*h)//14+(q-1+3*h)//14+1)


def count_controls():
    rows=0
    digest=sha256()
    for w in range(1,62):
        if w%3==0:
            continue
        for c in range(1,121):
            B=blocked(w,c)
            g=gcd(w,c)
            need(len(B)==blocker_count(w,c),('literal residue/count identity',w,c))
            need((len(B)==w)==(c%(14*w)==0),('single-owner event cover iff',w,c))
            need(len(B)<=g*((w//g+6)//7),('gcd-weighted event count ceiling',w,c))
            digest.update(repr((w,c,tuple(sorted(B)))).encode())
            rows+=1
    return rows,digest.hexdigest()


def phase_loss_controls():
    need(blocked(5,1)==blocked(5,71)=={0},'same full affine residue')
    need(blocked(5,29)=={1},'same gcd/count/quotient-residue but shifted event')
    need(all(gcd(c,5)==1 and c%14==1 for c in (1,29,71)),'identical proposed scalar signatures')
    need(len(blocked(5,1)|blocked(5,71))==1,'coincident blocker union')
    need(len(blocked(5,1)|blocked(5,29))==2,'separated blocker union')


def bohr_packet_controls():
    C=tuple(range(1,11));T=(1,5,11)
    delta=F(1,11)-F(1,14)
    rho=delta/10
    need((delta,rho)==(F(3,154),F(3,1540)),'exact full return packet radius')
    need(10*delta<F(1,2),'first-speed restriction prevents max-speed wrap')
    for y0 in (F(1,11),F(2,11)):
        need(min(distance(c*y0) for c in C)==F(1,11),'both anchors are beta deep')
    left,right=F(2,11)-rho,F(2,11)+rho
    need((left,right)==(F(277,1540),F(283,1540)),'spoiled translated packet')
    need(F(5,28)<left<right<F(13,70),'whole packet lies in inherited spoiled component')
    owners=[]
    for w,n in ((1,0),(5,1),(11,2)):
        need(n-F(3,14)<w*left<w*right<n+F(3,14),'fixed effective tooth throughout packet')
        owners.append((-n*pow(w,-1,3))%3)
    need(owners==[0,1,2],'all three owners remain a permutation')
    left,right=F(1,11)-rho,F(1,11)+rho
    need(F(3,14)<5*left<5*right<1-F(3,14),'same packet at another deep anchor has inactive tail')
    return delta,rho,owners


def all_event_hostile():
    C=(1,2,3,4,6,7,8,14,70,154);T=(1,5,11)
    for w in T:
        need(14*w in C and blocked(w,14*w)==set(range(w)),'a body owner blocks the full event coset')
    S=tuple(3*c for c in C)+T
    need(min(distance(F(v,13)) for v in S)==F(1,13),'event-only hostile is physically safe')
    return C,T,F(1,13)


def divisor_event_control():
    for w in range(37,44):
        need(6*((w+6)//7)<w,'all residue classes at sharp threshold37')
    need(6*((36+6)//7)==36,'threshold36 is a count-statistic equality hostile')
    for w in range(31,37):
        if w%3:
            need(6*((w+6)//7)<w,'ternary-unit filter improves uniform threshold to31')
    need(6*((29+6)//7)>29,'last ternary-unit failure of this count statistic')
    C=(8,9,10,11,13,14,43,86,129,172);T=(1,5,43);w=43
    need(sum(c%w==0 for c in C)==4 and all(c%(14*w) for c in C),'four harmless divisible body speeds')
    need(all(c%w==0 or gcd(c,w)==1 for c in C),'six residual coprime speeds')
    unions=set()
    for c in C:
        unions|=blocked(w,c)
    survivors=sorted(set(range(w))-unions)
    need(survivors==[1,2,7,8,10,11,20,22,25,26,28,29,37,40,41,42],
         'complete positive control event list')
    S=tuple(3*c for c in C)+T
    need(all(any(v%q==0 for v in S) for q in range(2,15)),'positive control survives every small-clock divisor sieve')
    y=F(17,602);x=F(619,1806)
    need(3*x%1==y,'raw body phase and physical sheet map')
    need(all(distance(c*y)>=F(1,14) for c in C),'event belongs to actual closed body set')
    need(distance(w*y)==F(3,14),'tail effective endpoint is weak-safe')
    need(min(distance(v*x) for v in S)>=F(1,14),'all thirteen physical speeds checked')
    return C,T,survivors,y,x,sum(blocker_count(w,c) for c in C)


def norm_five_fixed_roof_hostile():
    w=(10,11,16);v=(1,2,-2);q=F(3,7*w[2])
    terms=[]
    for k in (1,2):
        terms.append(tuple(min(q,F(3*(sum(w)-w[i])-14*k*abs(v[i]),14*w[(i+1)%3]*w[(i+2)%3]))
                           for i in range(3)))
    E=tuple(2*sum(row[i] for row in terms) for i in range(3))
    mass=2*sum(min(row) for row in terms)
    need(E==(F(17,176),F(9,140),F(3,55)),'all exact norm-five projections')
    need(mass==F(331,6160) and min(E)-mass==F(1,1232),'fixed norm-five roof identity is false')
    need(terms[0].index(min(terms[0]))==1 and terms[1].index(min(terms[1]))==2,'minimum roof switches with multiplier')
    return w,E,mass


def main():
    sys.stdout.reconfigure(newline='\n')
    counts=count_controls()
    phase_loss_controls()
    packet=bohr_packet_controls()
    hostile=all_event_hostile()
    positive=divisor_event_control()
    norm5=norm_five_fixed_roof_hostile()
    print('scope=elementary event-interface corollary and exact sidecars; arbitrary entry remains OPEN')
    print('count_universe=1<=w<=61,3 does not divide w,1<=c<=120; rows,digest='+str(counts))
    print('single_owner_event_cover=iff14w divides c')
    print('phase_loss=w5:B1=B71={0},B29={1}; samegcd and h mod14 do not preserve union')
    print('full_Bohr_packet_hostile='+str(packet))
    print('all_event_cosets_blocked_but_safe='+str(hostile))
    print('divisor_event_positive_control='+str(positive))
    print('norm5_fixed_projection_hostile='+str(norm5))
    print('optimization_live_gates='+str(CHECKS))


if __name__=='__main__':
    main()
