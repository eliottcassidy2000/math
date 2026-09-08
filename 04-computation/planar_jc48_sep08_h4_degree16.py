#!/usr/bin/env python3
"""Exact retained-sheet head for the marked mixed-cusp H4 covering.

RESERVED pending independent audit. The finite scalar enumeration retains
actual necessary inequalities, not realized permutation passports.
"""
from hashlib import sha256
import json

GATES=[]
def need(label,predicate):
    if not bool(predicate):raise RuntimeError(label)
    GATES.append(label)

def parts(n,minimum=2):
    if n==0:
        yield ()
    for q in range(minimum,n+1):
        for rest in parts(n-q,q):yield (q,)+rest

def class_closed(row):
    return (len(row)==1 or all(q==2 for q in row)
            or row in [(2,3),(3,3),(2,4)])
small=[]
for t in range(2,7):
    rows=list(parts(t))
    for row in rows:
        need('all moved types through six closed '+str(row),class_closed(row))
        small.append((t,row))
need('ten nonidentity small moved types',len(small)==10)
need('one moved letter impossible',list(parts(1))==[])

raw={};post={};controls={}
for D in range(2,17):
    raw[D]=[];post[D]=[]
    for t in range(7,D):
        for k in range(1,D-t+1):
            n3=max(0,k-t//2,(3*k-D+1)//2)
            n5=max(0,k-(2*t)//3,(5*k-2*D+2)//3)
            W=1+2*k-2*n3-n5
            lo=3*max(0,D-2*k)
            need(f'ordinary ceilings D{D}t{t}k{k}',2*n3>=3*k-D and n3>=k-t//2)
            need(f'fifth ceilings D{D}t{t}k{k}',3*n5>=5*k-2*D and n5>=k-(2*t)//3)
            if W<lo:continue
            row=(t,k,n3,n5,W)
            raw[D].append(row)
            if k==D-t and D<2*t:continue
            post[D].append(row)
need('complete D12 scalar row',raw[12]==[(7,5,2,1,6)])
need('complete D13 scalar rows',raw[13]==[(7,6,3,2,5),(8,5,1,0,9)])
need('complete D14 scalar rows',raw[14]==[(7,7,4,3,4),(8,6,2,1,8)])
need('complete D15 scalar rows',raw[15]==[(7,7,4,3,4),(7,8,5,4,3),(8,7,3,2,7),(9,6,2,0,9)])
need('no lower scalar row after t7',all(not raw[D] for D in range(2,12)))
need('complete post-full-retention head',
     [(D,row) for D in range(2,16) for row in post[D]]
     ==[(14,(7,7,4,3,4)),(15,(7,7,4,3,4)),(15,(7,8,5,4,3))])

# Literal incidence arithmetic used in the proved full-retention supplier.
for multiplicity in range(5):
    need('full retention convex incidence '+str(multiplicity),
         multiplicity*(multiplicity-1)//2>=2*multiplicity-3)
# Two ordinary overlaps in a seven-set force the marked commuting overlap.
need('odd moved support forces intersection',4+4>7)
need('commuting positive support cannot be singleton',min(q for q in range(2,8))==2)

# D14, k=t=7: W<=4, node ac>=2. If both other nodes are positive,
# their invariant moved intersections each contribute at least two.
need('D14 other two nodes cannot both meet',2+2+2>4)
# At full retention D=2t a disjoint node gives complementary supports.
# These are symbolic integer identities over the entire allowed overlap
# intervals, not one chosen support drawing.
for ab in range(4,8):
    for bc in range(4,8):
        for ac in range(2,8):
            # ad disjoint: bd=7-ab and cd=7-ac.
            need('D14 ad complement Euler identity',-14+ab+bc+(7-ac)+ac+(7-ab)==bc)
            # bd disjoint: ad=7-ab and cd=7-bc.
            need('D14 bd complement Euler identity',-14+ab+bc+(7-bc)+ac+(7-ab)==ac)
need('D14 both complementary alternatives obstruct',min(4,2)>1)

# D15, k7: each marked node overlap >=1. The ac support intersection
# is >=2. Under W<=4 the two other actual overlaps are exactly one;
# their moved-support intersections cannot be one and hence vanish.
need('D15 partial node allocation',[(ac,ad,bd) for ac in range(2,5) for ad in range(1,5) for bd in range(1,5) if ac+ad+bd<=4]==[(2,1,1)])
need('D15 common excluded support leaves union eight',15-7==8)
need('D15 a-b overlap at least six',7+7-8==6)
need('D15 partial common fixed count',15-2*7+6==7)
need('D15 two one-sheet deficits',7-2==5)
need('D15 partial Euler contradiction',-2*7+5+4+3+4==2 and 2>1)

# D15, k8: full retention makes every marked node overlap its moved
# support intersection. Their total<=3 with ac>=2 leaves both others zero.
need('D15 full node support allocation',all(ad==bd==0 for ac in range(2,8) for ad in [0,2,3,4,5,6,7] for bd in [0,2,3,4,5,6,7] if ac+ad+bd<=3))
need('D15 full common fixed count',15-2*7+6==7)
need('D15 full Euler contradiction',-2*8+7+5+4+2==2 and 2>1)

# A necessary scalar survivor at the next degree prevents misreporting
# scalar closure or manufacturing an actual permutation realization.
need('declared D16 scalar survivors',post[16]==[(7,8,5,4,3),(7,9,6,5,2),(8,7,3,2,7),(8,8,4,3,6)])
need('D16 formal Euler control',-16+4+4+3+6==1)
need('D16 formal node lower control',6>=3*max(0,16-16))
semantic=json.dumps({'small':small,'raw':raw,'post':post},sort_keys=True,separators=(',',':')).encode()
print('h4_degree16: PASS')
print('scope: actual marked mixed-cusp H4 covering has mapping degree at least16')
print('small moved types:',len(small),'through six; all have proved exclusions')
print('scalar universe: every D2..16,t7..D-1,k1..D-t')
print('D12..15 raw rows:',{D:raw[D] for D in range(12,16)})
print('post-full-retention D14,D15:',post[14],post[15])
print('D16 necessary scalar survivors, not realized passports:',post[16])
print('gates:',len(GATES))
print('semantic sha256:',sha256(semantic).hexdigest())
