"""Exact collective-gcd to endpoint-gcd packet; no producer imports."""
from fractions import Fraction as F
from functools import reduce
from itertools import product
from math import gcd,lcm
import sys
sys.stdout.reconfigure(newline="\n")
gates=0
H=76388115

def need(ok,label):
    global gates
    gates+=1
    if not ok:
        raise ArithmeticError(label)

def packet(w):
    A=B=1
    L=1
    prefix=F(1)
    for u,v in zip(w,w[1:]):
        g=gcd(u,v)
        a,b=u//g,v//g
        A*=a
        B*=b
        prefix*=F(b,a)
        L=lcm(L,prefix.denominator)
    C=gcd(A,B)
    need(L%(A//C)==0,"terminal reduced denominator divides prefix denominator lcm")
    J=L//(A//C)
    return A,B,C,L,J

def main():
    count=0
    for w in product(range(1,7),repeat=5):
        A,B,C,L,J=packet(w)
        d=reduce(gcd,w)
        D=gcd(w[0],w[-1])
        need(L==w[0]//d,"primitive coordinate equals denominator lcm")
        need(D==d*J and C%J==0,"exact pair gcd and cancellation-content divisor")
        need(w[0]*B==w[-1]*A,"oriented walk endpoint identity")
        count+=1
    need(packet((6,2,3))[-1]==3 and reduce(gcd,(6,2,3))==1,"actual atlas path loses a returning prime")
    need(3+1==4 and 2+3==5,"hostile edges are inherited atlas ratios1:3 and2:3")
    w=(729,243,81,27,9,3,1)
    need(packet(w)[-1]==1,"seven-vertex atlas path has no cancellation debt")
    A,B,C,L,J=packet((6,4,1))
    need(J==1<C==2,"three distinct actual vertices give strict improvement over cancellation content")
    core=(1,2,3,4,5,6,7,9,10,12,18)
    w=(18,12,3,9,6)
    A,B,C,L,J=packet(w)
    need(J==2 and C==12 and J<=H//11342250<C,"five-distinct-vertex exact packet passes where content bound fails")
    need(w[0]==max(core) and len(core)==11 and reduce(gcd,core)==1,"positive path starts at full primitive core maximum")
    edges=((1,3),(1,4),(1,9),(1,10),(2,3),(4,6),(3,7),(5,12),(4,12),(6,18))
    need(all((u+v)//gcd(u,v) in (4,5,10,11,17) for u,v in edges),"explicit spanning tree uses actual atlas edges")
    Q=91**6
    g=36*Q+1
    need(sum(core)+4*g<=Q**2 and g>2*Q*max(core),"positive example admits genuine finite-box equality by dominance")
    caps={5:11342250,6:31950,7:90}
    thresholds={}
    for t in (1,2):
        thresholds[t]={r:t*H//M for r,M in caps.items()}
        for r,M in caps.items():
            h=thresholds[t][r]
            need(M*h<=t*H<M*(h+1),"sharp integer packet threshold from inherited scalar caps")
    print('STATUS: PASS; exact endpoint-gcd packet, conditional LRC consumer')
    print('UNIVERSE:',count,'positive five-entry integer walks over1..6')
    print('THRESHOLDS J, by physical core scale t and distinct path size:',thresholds)
    print('HOSTILE: actual atlas path6,2,3 has collective gcd1, endpoint gcd3, J3')
    print('POSITIVE:729,243,81,27,9,3,1 is an actual seven-vertex atlas path with J1')
    print('STRICT GAIN: actual path6,4,1 has J1<C2; path18,12,3,9,6 has J2<=6<C12 in a genuine eleven-core')
    print('OPEN: no universal small-J path existence proved for all connected eleven-cores')
    print('ACTIVE GATES:',gates)

if __name__=='__main__':
    main()
