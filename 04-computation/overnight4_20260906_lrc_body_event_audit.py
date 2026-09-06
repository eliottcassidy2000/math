#!/usr/bin/env python3
"""Independent rational-danger-interval audit of the body/event interface.

Sorted event phases are located in open danger intervals by binary search;
the producer's modular residue membership engine is not reused or imported.
"""
from bisect import bisect_left,bisect_right
from fractions import Fraction as Q
from hashlib import sha256
from math import gcd
from pathlib import Path
import csv
import sys

sys.stdout.reconfigure(newline='\n')
ROOT=Path(__file__).resolve().parents[1]
GATES=0


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def event_points(w):
    return [Q(14*k+3,14*w) for k in range(w)]


def interval_blockers(w,c):
    points=event_points(w)
    bad=set()
    for m in range(c+1):
        left,right=Q(14*m-1,14*c),Q(14*m+1,14*c)
        # Strict left and right endpoints: equality is body-safe.
        first=bisect_right(points,left)
        last=bisect_left(points,right)
        bad.update(range(first,last))
    return bad


def count_formula(w,c):
    g=gcd(w,c)
    q,h=w//g,c//g
    return g*((q-1-3*h)//14+(q-1+3*h)//14+1)


def safe(x,c):
    m=(c*x).numerator//(c*x).denominator
    return not any(Q(14*j-1,14*c)<x<Q(14*j+1,14*c) for j in (m,m+1))


def clearance(x,c):
    m=(c*x).numerator//(c*x).denominator
    return c*min(abs(x-Q(m,c)),abs(x-Q(m+1,c)))


def full_interval_bad(left,right,w):
    return [m for m in range(w+1)
            if Q(14*m-1,14*w)<left<right<Q(14*m+1,14*w)]


def disjoint_from_danger(left,right,w):
    return all(max(left,Q(14*m-1,14*w))>=min(right,Q(14*m+1,14*w)) for m in range(w+1))


def main():
    digest=sha256()
    rows=0
    for w in range(1,62):
        if w%3==0:
            continue
        for c in range(1,121):
            bad=interval_blockers(w,c)
            need(len(bad)==count_formula(w,c),('rational danger intervals versus floor count',w,c))
            need((len(bad)==w)==(c%(14*w)==0),('complete one-owner cover iff',w,c))
            g=gcd(w,c)
            need(len(bad)<=g*((w//g+6)//7),('gcd-weighted marginal bound',w,c))
            digest.update(repr((w,c,tuple(sorted(bad)))).encode())
            rows+=1
    need(rows==4920,'complete independent marginal universe')
    need(digest.hexdigest()=='bbeeb45e66237a65b66f0837833787a2bb2275abd5ae5f1eb1c7bc8a306d3671',
         'every full affine event word agrees with independent residue certificate')
    need(interval_blockers(5,1)==interval_blockers(5,71)=={0} and interval_blockers(5,29)=={1},
         'gcd and quotient residue lose overlap phase')
    C=(8,9,10,11,13,14,43,86,129,172)
    T=(1,5,43)
    w=43
    words={c:interval_blockers(w,c) for c in C}
    union=set().union(*words.values())
    survivors=sorted(set(range(w))-union)
    need(survivors==[1,2,7,8,10,11,20,22,25,26,28,29,37,40,41,42],
         'full positive event-survivor list from geometric intervals')
    need(sum(map(len,words.values()))==37 and len(union)==27,'marginal versus union loss')
    y=event_points(w)[survivors[0]]
    need(y==Q(17,602) and all(safe(y,c) for c in C),'actual body-safe event')
    physical=tuple(3*c for c in C)+T
    safe_labels=[j for j in range(3) if all(safe((y+j)/3,v) for v in physical)]
    need(safe_labels==[1,2],'all valid labels reconstructed from literal danger intervals')
    x=(y+1)/3
    need(x==Q(619,1806) and all(safe(x,v) for v in physical),'positive physical witness')
    need(min(clearance(x,v) for v in physical)==Q(1,7),'positive witness exact clearance')
    need(all(any(v%d==0 for v in physical) for d in range(2,15)),'all small-clock divisors are present')
    for w in range(37,44):
        need(w-6*((w+6)//7)>0 and (w+7)-6*((w+13)//7)==w-6*((w+6)//7)+1,
             'seven residue proof of uniform count threshold')
    need(36==6*((36+6)//7),'all-integer count equality boundary')
    for w in range(31,37):
        if w%3:
            need(6*((w+6)//7)<w,'ternary-unit threshold31 corner')
    need(6*((29+6)//7)>29,'last ternary-unit numerical failure')
    H=(1,2,3,4,6,7,8,14,70,154)
    tails=(1,5,11)
    for w in tails:
        need(interval_blockers(w,14*w)==set(range(w)),'full positive event coset blocked by literal intervals')
        need(all(not safe(1-y,14*w) for y in event_points(w)),'negative event coset also blocked')
    safe_row=tuple(3*c for c in H)+tails
    need(all(safe(Q(1,13),v) for v in safe_row),'event incidence is not necessary')
    need(min(clearance(Q(1,13),v) for v in safe_row)==Q(1,13),'safe hostile clearance')
    beta=Q(1,11)
    delta=beta-Q(1,14)
    rho=delta/10
    need((delta,rho)==(Q(3,154),Q(3,1540)) and 10*delta<Q(1,2),'whole Bohr set has one centered interval')
    body=tuple(range(1,11))
    for center in (Q(1,11),Q(2,11)):
        need(min(clearance(center,c) for c in body)==beta,'both anchors have the same deep margin')
    left,right=Q(2,11)-rho,Q(2,11)+rho
    need((left,right)==(Q(277,1540),Q(283,1540)),'complete spoiled translated packet')
    center=Q(2,11)
    safe_left=[]
    safe_right=[]
    for c in body:
        m=(c*center).numerator//(c*center).denominator
        safe_left.append(Q(14*m+1,14*c))
        safe_right.append(Q(14*m+13,14*c))
    component=(max(safe_left),min(safe_right))
    need(component==(Q(5,28),Q(13,70)),'spoiled body component reconstructed from safe intervals')
    need(component[0]<left<right<component[1],'strict packet containment in actual body component')
    tooth_certificate=[]
    for j in range(3):
        intervals=[(w,m) for w in tails for m in full_interval_bad((left+j)/3,(right+j)/3,w)]
        need(bool(intervals),('entire physical packet lift is spoiled',j))
        tooth_certificate.append(intervals)
    need(tooth_certificate==[[(1,0)],[(5,2)],[(11,8)]],'independent physical tooth ownership')
    good_left,good_right=Q(1,11)-rho,Q(1,11)+rho
    need(all(disjoint_from_danger((good_left+j)/3,(good_right+j)/3,5) for j in range(3)),
         'same entire return packet has inactive tail5 at the other deep anchor')
    w=(10,11,16)
    terms=[]
    joint=[]
    for k in (1,2):
        nearest=(2*k,2*k,3*k)
        intervals=[(Q(14*n-3,14*s),Q(14*n+3,14*s)) for n,s in zip(nearest,w)]
        cap=min(b-a for a,b in intervals)
        row=[]
        for omitted in range(3):
            others=[intervals[j] for j in range(3) if j!=omitted]
            row.append(min(cap,min(b for a,b in others)-max(a for a,b in others)))
        terms.append(tuple(row))
        joint.append(min(b for a,b in intervals)-max(a for a,b in intervals))
    E=tuple(2*sum(row[i] for row in terms) for i in range(3))
    mass=2*sum(joint)
    need(E==(Q(17,176),Q(9,140),Q(3,55)) and mass==Q(331,6160),'norm5 independent interval widths')
    need(min(E)-mass==Q(1,1232) and [r.index(min(r)) for r in terms]==[1,2],
         'physical and selected projections differ through a switching roof')
    tsv=ROOT/'05-knowledge/results/overnight4_20260906_lrc_parityfree_probe_head63.tsv'
    with tsv.open(newline='',encoding='utf-8') as f:
        row=next(r for r in csv.DictReader(f,delimiter='\t') if tuple(int(r[z]) for z in ('a','b','c'))==w)
    denominator=int(row['denominator'])
    need(tuple(Q(int(row['E'+str(i)+'_numerator']),denominator) for i in range(3))==E
         and Q(int(row['mass_numerator']),denominator)==mass and int(row['raw_carriers'])==4,
         'norm5 matches every native-audited frozen head column')
    print('COUNT_UNIVERSE',rows,'method=rational danger endpoints and binary search')
    print('FULL_EVENT_WORD_SEMANTIC_SHA256',digest.hexdigest())
    print('POSITIVE_CONTROL survivors',survivors,'marginals',37,'union',27,'y',y,'safe_labels',safe_labels,'x',x,'clearance',Q(1,7))
    print('COUNT_THRESHOLDS all_positive37 ternary_unit31; no event-sharpness inference')
    print('ALL_EVENT_COSETS_BLOCKED_BUT_SAFE x=1/13 clearance=1/13')
    print('FULL_PACKET',left,right,'body_component',component,'physical_danger_teeth',tooth_certificate)
    print('NORM5_INTERVAL_CONTROL',E,mass,'selected_minus_physical',min(E)-mass,'terms',terms)
    print('NATIVE_AUDITED_HEAD_SHA256',sha256(tsv.read_bytes()).hexdigest())
    print('PASS_OPTIMIZATION_LIVE_GATES',GATES)


if __name__=='__main__':
    main()
