#!/usr/bin/env python3
"""Independent literal-predicate audit of the scale-three Haar consumers.

No producer imports and no permutation-of-three interval intersection.
All breakpoints and midpoint predicates use exact integer rational grids.
Body boundary points are tested separately, retaining isolated safe points.
The adjacent-family proof uses exact residue-class polynomial certificates.
"""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from math import gcd,lcm
import json
import sys

sys.stdout.reconfigure(newline="\n")
GATES=0


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def literal_tail(speeds):
    denominator=42*lcm(*speeds)
    # F(x+1/3)=F(x), merely permuting the three literal sheet predicates.
    end=denominator//3
    points={0,end}
    for w in speeds:
        unit=denominator//(14*w)
        for k in range(w):
            for sign in (-1,1):
                for j in range(3):
                    point=((14*k+sign)*unit-j*end)%denominator
                    if point<=end:
                        points.add(point)
    points=sorted(points)
    twice=2*denominator
    length=0
    for left,right in zip(points,points[1:]):
        x=left+right
        blocked=True
        for j in range(3):
            one=False
            for w in speeds:
                residue=w*(x+2*j*end)%twice
                if 14*min(residue,twice-residue)<twice:
                    one=True
                    break
            if not one:
                blocked=False
                break
        if blocked:
            length+=right-left
    return Q(3*length,denominator)


def body_predicate(speeds,numerator,denominator):
    return all(14*min((w*numerator)%denominator,
                      denominator-(w*numerator)%denominator)>=denominator for w in speeds)


def literal_body(speeds):
    denominator=14*lcm(*speeds)
    points={0,denominator}
    for w in speeds:
        unit=denominator//(14*w)
        for k in range(w):
            points.add((14*k+1)*unit)
            points.add((14*k+13)*unit)
    points=sorted(points)
    inside=[body_predicate(speeds,a+b,2*denominator) for a,b in zip(points,points[1:])]
    endpoint=[body_predicate(speeds,a,denominator) for a in points]
    intervals=[]
    for i,ok in enumerate(inside):
        if not ok:
            continue
        need(endpoint[i] and endpoint[i+1],"closed body keeps both interval endpoints")
        a,b=points[i],points[i+1]
        if intervals and intervals[-1][1]==a:
            intervals[-1]=(intervals[-1][0],b)
        else:
            intervals.append((a,b))
    isolated=[a for i,a in enumerate(points) if endpoint[i]
              and not(i and inside[i-1]) and not(i<len(inside) and inside[i])]
    mass=Q(sum(b-a for a,b in intervals),denominator)
    return mass,tuple((Q(a,denominator),Q(b,denominator)) for a,b in intervals),tuple(Q(a,denominator) for a in isolated)


def nonmultiples(k):
    multiples=k//3
    return k-multiples,k*(k+1)//2-3*multiples*(multiples+1)//2


def adjacent_numerator(b):
    bulk_count,bulk_sum=nonmultiples(3*b//14)
    all_count,all_sum=nonmultiples((3*b+2)//14)
    return (3*(2*b+1)*bulk_count-14*bulk_sum
            +(b+1)*(3*(b+1)*(all_count-bulk_count)-14*(all_sum-bulk_sum)))


def shift(poly,q0):
    A,B,C=poly
    return A,2*A*q0+B,A*q0*q0+B*q0+C


def adjacent_proof():
    table=[]
    for r in range(1,42,3):
        values=[adjacent_numerator(42*q+r) for q in range(3)]
        C=values[0]
        twiceA=values[2]-2*values[1]+C
        need(twiceA%2==0,"quadratic integer coefficient")
        A=twiceA//2
        B=values[1]-C-A
        # Degree <=2 follows algebraically from the count/sum formulas:
        # both cutoffs are 9q+constant on each residue class.
        need(A==1134,"universal asymptotic numerator coefficient")
        for q in range(6):
            need(adjacent_numerator(42*q+r)==A*q*q+B*q+C,"count formula versus residue polynomial")
        lower=(11*A-6*42**2,11*B-6*42*(2*r+1),11*C-6*r*(r+1))
        lower=shift(lower,int(r<7))
        need(all(x>0 for x in lower),"all-height strict 6/77 violation for b>=7")
        upper=(42*42**2-55*A,42*42*(2*r+1)-55*B,42*r*(r+1)-55*C)
        upper=shift(upper,int(r<4))
        need(upper[0]>0 and upper[1]>0 and upper[2]>=0,"all-height 6/55 upper bound")
        need((upper[2]==0)==(r==10),"unique equality at b10")
        table.append((r,A,B,C,lower,upper))
        for q in (0,1,2):
            b=42*q+r
            if b>=4:
                need(literal_tail((1,b,b+1))==Q(A*q*q+B*q+C,7*b*(b+1)),
                     "literal physical predicate verifies residue formula")
    need(Q(1134,7*42**2)==Q(9,98),"exact limiting mass")
    print("ADJACENT_RESIDUE_TABLE",table)
    return table


def main():
    anchor_tails=((1,5,7),(1,5,11),(1,2,4),(1,2,5),(1,4,5),(2,5,7),(1,10,11))
    print("LITERAL_TAIL_CONTROLS",[(w,str(literal_tail(w))) for w in anchor_tails])
    need(literal_tail((2,5,7))==Q(22,245)>Q(6,77),"mixed-parity hostile")
    need(literal_tail((1,10,11))==Q(6,55),"adjacent maximum witness")
    units=[w for w in range(1,61) if w%3]
    tails={w:literal_tail(w) for w in combinations(units,3) if gcd(*w)==1}
    violations=[(w,value) for w,value in tails.items() if value>Q(6,77)]
    first=min(violations,key=lambda item:(item[0][-1],item[0]))
    maximum=max(tails.values())
    maximizers=[w for w,value in tails.items() if value==maximum]
    need(len(tails)==8664,"complete primitive ternary-unit head")
    need(len(violations)==136,"independent literal violation count")
    need(first==((2,5,7),Q(22,245)),"smallest maximum-speed hostile")
    need(maximum==Q(6,55) and maximizers==[(1,10,11)],"independent head maximum and uniqueness")
    nonadditive=[(w,value) for w,value in violations if w[0]+w[1]!=w[2]]
    need(nonadditive==[((2,11,20),Q(11,140))],"sole nonadditive head violation")
    print("LITERAL_HEAD60",len(tails),"VIOLATIONS",len(violations),"FIRST",first,"MAXIMUM",maximum,"AT",maximizers)
    print("NONADDITIVE_HEAD_VIOLATIONS",nonadditive)
    C=(1,2,3,5,7,8,9,11,12,13)
    mass,components,isolated=literal_body(C)
    need(mass==Q(14249,252252)<Q(6,77),"recovered ten-body hostile")
    need(isolated==tuple(Q(r,14) for r in (1,3,5,9,11,13)),"six isolated weakly-safe body points")
    lower=((Q(15,154),Q(13,126)),(Q(29,182),Q(9,56)),(Q(29,168),Q(27,154)),
           (Q(43,182),Q(27,112)),(Q(29,112),Q(41,154)),(Q(29,98),Q(55,182)))
    need(components==lower+tuple((1-b,1-a) for a,b in reversed(lower)),"complete recovered component addresses")
    print("RECOVERED_BODY",C,"MASS",mass,"INTERVALS",len(components),"ISOLATED",list(map(str,isolated)))
    bodies={c:literal_body(c)[0] for c in combinations(range(1,14),10)}
    below12={c:v for c,v in bodies.items() if max(c)<=12}
    minimum=min(below12.values())
    bad13=[(c,v) for c,v in bodies.items() if v<Q(6,77)]
    need(len(below12)==66 and minimum==Q(16277,194040)>Q(6,77),"no ten-body violation below13")
    need(len(bodies)==286 and len(bad13)==12,"first body head has12 failures")
    print("BODY_HEAD",len(below12),"THROUGH12_MINIMUM",minimum,"THROUGH13",len(bodies),"VIOLATIONS",len(bad13))
    full=tuple(sorted(tuple(3*c for c in C)+(1,5,11)))
    fullmass,fullcomponents,fullisolated=literal_body(full)
    need(len(full)==len(set(full))==13,"full velocity set has13 distinct positive speeds")
    need(fullmass==Q(25331,756756),"full consumer safe measure")
    need(body_predicate(full,1,14),"fixed denominator14 weak witness")
    print("FULL_CONSUMER",full,"MASS",fullmass,"INTERVALS",len(fullcomponents),"ISOLATED",list(map(str,fullisolated)),"SAFE_AT_1_OVER_14",True)
    table=adjacent_proof()
    manifest={"tails":[[w,str(v)] for w,v in tails.items()],
              "body":[str(mass),[[str(a),str(b)] for a,b in components],list(map(str,isolated))],
              "bodies":[[c,str(v)] for c,v in bodies.items()],"adjacent":table}
    print("SEMANTIC_SHA256",sha256(json.dumps(manifest,separators=(",",":")).encode()).hexdigest())
    print("PASS independent literal predicates, endpoints, finite heads and unbounded family proof")
    print("GATES",GATES)


if __name__=="__main__":
    main()
