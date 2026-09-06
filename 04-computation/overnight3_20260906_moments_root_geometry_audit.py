#!/usr/bin/env python3
"""Independent exact replay of the complete-channel root/sign certificate.

No producer imports, SymPy, root approximation or root isolation. Literal
balanced-count enumeration rebuilds both coefficient rows. Rational Sturm
chains certify the supplied first-root intervals; independent polynomial
division and interval Horner certify every true Laurent second-row sign.
"""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
from math import factorial,gcd,lcm
from pathlib import Path
import json
import sys

sys.stdout.reconfigure(newline="\n")
GATES=0


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def trim(poly):
    result=list(poly)
    while result and not result[-1]:
        result.pop()
    return result


def remainder(a,b):
    a,b=list(map(Q,a)),list(map(Q,b))
    while len(a)>=len(b):
        shift=len(a)-len(b)
        factor=a[-1]/b[-1]
        for j,v in enumerate(b):
            a[j+shift]-=factor*v
        a=trim(a)
    return a


def primitive(poly):
    poly=list(map(Q,poly))
    if not poly:
        return []
    denominator=lcm(*(v.denominator for v in poly))
    integers=[int(v*denominator) for v in poly]
    content=gcd(*integers)
    return [v//content for v in integers]


def value(poly,x):
    result=Q(0)
    for c in reversed(poly):
        result=result*x+c
    return result


def sturm(poly):
    derivative=[j*poly[j] for j in range(1,len(poly))]
    seq=[primitive(poly),primitive(derivative)]
    while seq[-1]:
        rem=remainder(seq[-2],seq[-1])
        if not rem:
            break
        seq.append(primitive([-v for v in rem]))
    return seq


def variation(seq,x):
    signs=[]
    for poly in seq:
        v=value(poly,x)
        if v:
            signs.append(1 if v>0 else -1)
    return sum(a!=b for a,b in zip(signs,signs[1:]))


def variation_negative_infinity(seq):
    signs=[(1 if poly[-1]>0 else -1)*(-1 if (len(poly)-1)%2 else 1) for poly in seq]
    return sum(a!=b for a,b in zip(signs,signs[1:]))


def enclosure(poly,left,right):
    lower=upper=Q(0)
    for c in reversed(poly):
        endpoints=(lower*left,lower*right,upper*left,upper*right)
        lower,upper=min(endpoints)+c,max(endpoints)+c
    return lower,upper


def universe():
    rows=set()
    for A,B in ((1,2),(2,3),(3,5),(4,7),(5,9),(5,12)):
        for h in (2,3,4,6,8):
            for r,t in {(0,0),(B-1,0),(0,A-1),(B-1,A-1)}:
                a=A*(B*h+r)+B*t
                count=0
                for n in range(1,500):
                    g=n+B*h+r+t
                    if A*n>(B-A)*t and gcd(a,g)==1:
                        rows.add((n,A,B,h,r,t))
                        count+=1
                        if count==2:
                            break
                need(count==2,"complete declared corner bank")
    rows.update(((2,1,2,2,1,0),(16,5,9,3,6,4)))
    return sorted(rows)


def literal_row(a,b,c,m,A,t,k):
    counts={}
    for zcount in range(m+1):
        numerator=a*m-(a+c)*zcount
        if numerator%(a+b):
            continue
        ycount=numerator//(a+b)
        xcount=m-ycount-zcount
        if min(xcount,ycount)<0:
            continue
        need((zcount-k*t)%A==0,"literal charge solution in scalar residue")
        j=(zcount-k*t)//A
        counts[j]=factorial(m)//(factorial(xcount)*factorial(ycount)*factorial(zcount))
    need(bool(counts),"nonempty balanced row")
    low,high=min(counts),max(counts)
    need(len(counts)==high-low+1,"complete scalar support")
    content=gcd(*counts.values())
    return [counts[j]//content for j in range(low,high+1)],low


def main():
    path=Path("04-computation/overnight3_20260906_moments_root_geometry_certificates.json")
    data=json.loads(path.read_text())
    cases=data["cases"]
    need([tuple(c["parameters_n_A_B_h_r_t"]) for c in cases]==universe(),"exact221-support universe")
    carries,degrees,signs=Counter(),Counter(),Counter()
    roots=second_roots=0
    for case in cases:
        n,A,B,h,r,t=case["parameters_n_A_B_h_r_t"]
        a,b,c=case["support"]
        g=case["g"]
        need((a,b,c)==(A*(B*h+r)+B*t,A*g-a,B*g-a),"declared support identity")
        need(gcd(a+b,a+c)==g and gcd(a,g)==1,"primitive first-return scale")
        P,first_low=literal_row(a,b,c,g,A,t,1)
        R,second_low=literal_row(a,b,c,2*g,A,t,2)
        need(P==case["first_coefficients_ascending"] and first_low==0,"first literal multinomial row")
        need(R==case["second_coefficients_ascending"] and second_low==case["second_laurent_low"],"second literal multinomial row and shift")
        seq=sturm(P)
        need(len(seq[-1])==1,"first polynomial squarefree")
        second_seq=sturm(R)
        need(len(second_seq[-1])==1,"second polynomial squarefree")
        need(variation_negative_infinity(second_seq)-variation(second_seq,Q(0))==len(R)-1,
             "every second root real and strictly negative")
        second_roots+=len(R)-1
        rem=remainder(R,P)
        certificates=case["certificates"]
        need(len(certificates)==len(P)-1,"certificate count equals first degree")
        previous=None
        for cert in certificates:
            left,right=map(Q,cert["interval"])
            need(left<=right<0,"negative rational interval")
            need(previous is None or previous<left,"strictly disjoint ordered root intervals")
            previous=right
            if left==right:
                need(value(P,left)==0 and value(seq[1],left)!=0,"exact simple rational root")
            else:
                need(value(P,left)!=0 and value(P,right)!=0,"nonroot isolating endpoints")
                need(variation(seq,left)-variation(seq,right)==1,"independent Sturm count exactly one")
            lower,upper=enclosure(rem,left,right)
            need((lower,upper)==tuple(map(Q,cert["remainder_enclosure"])),"independent division and interval enclosure")
            need(lower>0 or upper<0,"strict nonzero root-value enclosure")
            compressed=1 if lower>0 else -1
            raw=compressed*(-1 if second_low%2 else 1)
            need(compressed==cert["compressed_sign"] and raw==cert["raw_laurent_sign"]==-1,"true second Laurent sign")
            signs[(second_low,compressed,raw)]+=1
            roots+=1
        carries[(2*t//A,2*r//B)]+=1
        degrees[h]+=1
    need(len(cases)==221 and roots==1015 and second_roots==2242,"full finite root universe")
    need(carries==Counter({(0,0):60,(0,1):60,(1,0):50,(1,1):51}),"all four carry signatures")
    need(degrees==Counter({2:44,3:45,4:44,6:44,8:44}),"first-degree bank counts")
    print("PASS independent literal balanced rows, rational Sturm intervals, division, and raw signs")
    print("SUPPORTS",len(cases),"BALANCED_ROWS",2*len(cases),"FIRST_ROOTS",roots,"SECOND_ROOTS",second_roots)
    print("DEGREES",sorted(degrees.items()))
    print("CARRIES",sorted(carries.items()))
    print("LOW_COMPRESSED_RAW_SIGNS",sorted(signs.items()))
    print("PRODUCER_CERTIFICATE_SHA256",sha256(path.read_bytes()).hexdigest())
    print("GATES",GATES)


if __name__=="__main__":
    main()
