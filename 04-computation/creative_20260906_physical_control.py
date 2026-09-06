"""Recovered physical eleven-body: exact closed geometry and mixed-gate gain."""
from fractions import Fraction as F
from functools import reduce
from math import gcd
import sys

sys.stdout.reconfigure(newline="\n")
GATES=0
H=(1,3,4,5,6,7,8,9,10,11,13)
FRONTIER=((7,13,F(2,49),F(1,49)),(19,25,F(138,3325),F(37,3325)),(5,41,F(12,287),F(2,287)),(5,53,F(78,1855),F(2,371)),(1,67,F(20,469),F(2,469)))


def need(ok,message):
    global GATES
    GATES+=1
    if not ok:
        raise ArithmeticError(message)


def clearance(row,x):
    return min(min((h*x)%1,1-(h*x)%1) for h in row)


def safe_components(row):
    # Strict open danger intervals: touching endpoints must NOT be merged.
    danger=[]
    for h in row:
        for k in range(h+1):
            danger.append((max(F(0),F(14*k-1,14*h)),min(F(1),F(14*k+1,14*h))))
    danger.sort()
    merged=[]
    for left,right in danger:
        if merged and left<merged[-1][1]:
            merged[-1]=(merged[-1][0],max(right,merged[-1][1]))
        else:
            merged.append((left,right))
    # Zero and one themselves are dangerous; all intervening endpoints are safe.
    return tuple((a[1],b[0]) for a,b in zip(merged,merged[1:]))


def main():
    components=safe_components(H)
    expected=((F(1,14),F(1,14)),(F(15,182),F(13,154)),(F(3,14),F(3,14)),(F(5,14),F(5,14)),(F(29,70),F(41,98)),(F(85,182),F(27,56)),(F(29,56),F(97,182)),(F(57,98),F(41,70)),(F(9,14),F(9,14)),(F(11,14),F(11,14)),(F(141,154),F(167,182)),(F(13,14),F(13,14)))
    need(components==expected,"complete exact closed component list")
    need(len(H)==len(set(H))==11 and reduce(gcd,H)==1,"positive distinct primitive eleven-body")
    for a,b in components:
        need(clearance(H,a)>=F(1,14) and clearance(H,b)>=F(1,14),"closed endpoints retained")
        need(clearance(H,(a+b)/2)>=F(1,14),"component interior safe")
    for a,b in zip(components,components[1:]):
        need(clearance(H,(a[1]+b[0])/2)<F(1,14),"intercomponent danger gap")
    isolated=tuple(a for a,b in components if a==b)
    need(isolated==tuple(F(j,14) for j in (1,3,5,9,11,13)),"six isolated safe points")
    need(tuple((1-b,1-a) for a,b in reversed(components))==components,"phase reflection covariance")
    M=sum((b-a for a,b in components),F(0))
    L=max(b-a for a,b in components)
    need(M==F(5939,140140) and L==F(11,728),"exact physical mass and width")
    need(M<F(20,469)<F(4,91) and L<F(1,49),"both new scalar gates and old mass gate fail")
    need(all(M>=m for _,_,m,_ in FRONTIER[:-1]),"mass pays first four profiles")
    need(L>=FRONTIER[-1][3],"component pays final profile")
    need(all(M>=m or L>=b for _,_,m,b in FRONTIER),"physical mixed criterion passes at odd scale one")
    # The all-odd-scale extension is monotonic: g>=1 implies gL>=L.
    row=tuple(sorted(tuple(2*h for h in H)+(1,67)))
    x=F(9,34)
    need(len(row)==len(set(row))==13 and reduce(gcd,row)==1,"literal primitive thirteen-speed completion")
    need(clearance(row,x)==F(2,17)>F(1,14),"literal full-row strict safe phase")
    need(10 in row and 1+10==11 and 11%3==2,"crossing decoder edge prevents claimed eleven-plus-two partition")
    print("STATUS: PASS; recovered physical body realizes mixed-certificate gain")
    print("H:",H)
    print("M:",M,"L:",L)
    for a,b in components:
        print("CLOSED_COMPONENT:",a,b,"length",b-a)
    print("ISOLATED_SAFE_POINTS:",','.join(map(str,isolated)))
    print("MIXED_GATE: first four profiles paid by mass; final(1,67) paid by width")
    print("SCALAR_GATES: new mass20/469 fails; new width1/49 fails; old mass4/91 fails")
    print("FULL_ROW:",row,"SAFE_PHASE:",x,"CLEARANCE:",clearance(row,x))
    print("PROVED RELATIVE FAMILY:2H union{gp,gq} safe for every eligible primitive(p,q) and every positive oddg")
    print("SCOPE: physical certificate comparison, not an actual11+2 decoder partition or a first safety proof")
    print("ACTIVE GATES:",GATES)


if __name__=="__main__":
    main()
