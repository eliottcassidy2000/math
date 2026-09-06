#!/usr/bin/env python3
"""Two independent exact safe-body constructions and endpoint synchronization.

Universe: every ten-element subset of1..13; a separate stated small-clock
control inside1..14; named first-exit and singleton tests. No LRC claim.
"""
from fractions import Fraction as Q
from itertools import combinations
from hashlib import sha256

CHECKS=0
R=Q(1,14)
TARGET=Q(6,77)


def need(test,detail):
    global CHECKS
    CHECKS+=1
    if not test:
        raise RuntimeError(detail)


def safe(C,x):
    return all(min((c*x)%1,(-c*x)%1)>=R for c in C)


def meet(A,B):
    out=[];i=j=0
    while i<len(A) and j<len(B):
        lo,hi=max(A[i][0],B[j][0]),min(A[i][1],B[j][1])
        if lo<=hi:
            out.append((lo,hi))
        ar,br=A[i][1],B[j][1]
        i+=ar<=br
        j+=br<=ar
    return out


def body_intersection(C):
    out=[(Q(0),Q(1))]
    for c in C:
        out=meet(out,[(Q(14*k+1,14*c),Q(14*k+13,14*c)) for k in range(c)])
    return out


def body_sweep(C):
    ends=sorted({Q(0),Q(1)}|{Q(14*k+e,14*c)%1 for c in C for k in range(c) for e in (-1,1)})
    pieces=[(x,x) for x in ends if safe(C,x)]
    pieces += [(lo,hi) for lo,hi in zip(ends,ends[1:]) if safe(C,(lo+hi)/2)]
    out=[]
    for lo,hi in sorted(pieces):
        if out and lo<=out[-1][1]:
            out[-1]=(out[-1][0],max(hi,out[-1][1]))
        else:
            out.append((lo,hi))
    return out


def mass(parts):
    return sum((hi-lo for lo,hi in parts),Q(0))


def owners(T,y):
    answer=[]
    for w in T:
        z=w*y
        n=(z+Q(1,2)).numerator//(z+Q(1,2)).denominator
        answer.append((n,(-n*pow(w,-1,3))%3 if abs(z-n)<3*R else None))
    return answer


def spoiled(T,y):
    return {o for n,o in owners(T,y)}=={0,1,2}


def first_exit(I,T):
    lo,hi=I
    state=owners(T,lo)
    labels={o for n,o in state if o is not None}
    if labels!={0,1,2}:
        return lo,min(set(range(3))-labels),'left endpoint'
    event,index=min(((Q(n)+3*R)/w,i) for i,(w,(n,o)) in enumerate(zip(T,state)))
    if event<=hi:
        return event,state[index][1],'first effective-tooth exit'
    return None


def literal_component_escape(I,T):
    lo,hi=I
    ends={lo,hi}|{Q(14*k+e,14*w)%1 for w in T for k in range(w) for e in (-3,3)}
    ends=sorted(x for x in ends if lo<=x<=hi)
    probes=ends+[(a+b)/2 for a,b in zip(ends,ends[1:])]
    return any(not spoiled(T,x) for x in probes)


def main():
    failures=[];digest=sha256();leader=(Q(1),None)
    for C in combinations(range(1,14),10):
        parts=body_intersection(C)
        need(parts==body_sweep(C),('complete safe components',C))
        value=mass(parts)
        digest.update(repr((C,parts,value)).encode())
        if value<leader[0]:leader=value,C
        if value<TARGET:failures.append((max(C),C,value))
    need(len(failures)==12,'complete286-body floor hostile count')
    first=min(failures)
    need(first==(13,(1,2,3,4,5,7,8,9,11,13),Q(21514,315315)),'first height/lex hostile')
    star=(1,2,3,5,7,8,9,11,12,13)
    need(leader==(Q(14249,252252),star),'finite bank minimum')
    print('COMPLETE BODY HEAD286; below6/77',len(failures),'FIRST',first,'LEADER',leader)
    print('SEMANTIC SHA256',digest.hexdigest())
    parts=body_intersection(star)
    print('HOSTILE COMPONENTS',parts)
    need((Q(3,14),Q(3,14)) in parts,'zero-measure safe component retained')
    T=(1,5,11)
    S=tuple(3*c for c in star)+T
    need(safe(S,Q(1,14)),'Haar-gate failure is not a counterexample')
    need(mass(parts)<TARGET,'even equality-tail-adaptive Haar gate fails')
    # Complete fixed-denominator sieve in this bounded bank. For q<=14,
    # a reduced p/q is safe iff no row speed is divisible by q.
    smallclock=[]
    for other in combinations((1,2,3,4,5,6,7,8,9,11,12),7):
        C=tuple(sorted(other+(10,13,14)))
        row=tuple(3*c for c in C)+T
        if all(any(s%q==0 for s in row) for q in range(2,15)):
            smallclock.append((mass(body_intersection(C)),C))
    need(len(smallclock)==209,'exact no-clock<=14 control bank')
    need(min(smallclock)==(Q(4163,51480),(1,2,3,8,9,10,11,12,13,14)),'filtered bank minimum')
    need(all(v>TARGET for v,C in smallclock),'bounded filtered bank passes Haar gate')
    print('FILTERED CONTROL209 maxbody14,no reduced clock<=14,T1,5,11; minimum',min(smallclock))
    bodies=[(6,),(99,),tuple(range(1,11)),star,(1,2,3,4,5,7,8,9,11,13)]
    tails=[(1,5,11),(1,19,79),(5,37,43),(41,47,49)]
    decisions=0
    for C in bodies:
        for T in tails:
            for I in body_intersection(C):
                result=first_exit(I,T)
                need((result is not None)==literal_component_escape(I,T),('exact endpoint/event decision',C,I,T,result))
                if result:
                    y,j,kind=result
                    need(I[0]<=y<=I[1] and safe(tuple(3*c for c in C)+T,(y+j)/3),
                         ('literal recovered physical witness',C,I,T,result))
                decisions+=1
    inside=(Q(5,28),Q(9,28))
    endpoint=(Q(89,462),Q(31,154))
    need(inside in body_intersection((6,)),'literal interior-event component')
    need(endpoint in body_intersection((99,)),'literal equality-event component')
    for C,I in (((6,),inside),((99,),endpoint)):
        event=first_exit(I,(1,5,11))
        need(event[:2]==(Q(31,154),2),'first exit/free owner exact')
        print('FIRST-EXIT CONTROL',C,I,event,'physical',(event[0]+event[1])/3)
    blocked=(Q(5,28),Q(13,70))
    need(blocked in body_intersection(tuple(range(1,11))) and first_exit(blocked,(1,5,11)) is None,
         'inherited body witness component entirely spoiled')
    need(spoiled((1,5,11),Q(2,11)),'canonical synchronization hostile retained')
    print('COMPONENT DECISIONS',decisions,'CHECKS',CHECKS)
    print('REFUTED universal ten-body6/77 Haar floor; PROVED inherited-owner first-exit decision compiler')
    print('OPEN actual arbitrary-row entry and synchronization; no new family closure claimed from fixed clocks')


if __name__=='__main__':
    main()
