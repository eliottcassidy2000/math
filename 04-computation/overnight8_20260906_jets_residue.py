#!/usr/bin/env python3
"""Odd four-jet precision, full p7 equilateral Smith, and a residue packet.

Analytic proofs are in the companion report. This exact control program
uses no repository imports. It checks the cubic polynomial identities,
complete small residue fibers, literal full Hasse Smith matrices, and the
general coefficient packet independently against reciprocal products.
Every gate remains active under Python -O. Output is deterministic LF.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import comb
import sys
from sympy import Matrix, Poly, ZZ, symbols, expand, resultant
from sympy.matrices.normalforms import smith_normal_form

sys.stdout.reconfigure(newline='\n')
GATES=0


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def vp(x,p):
    x=Q(x)
    if not x:
        return float('inf')
    a,b=abs(x.numerator),x.denominator
    value=0
    while a%p==0:
        a//=p;value+=1
    while b%p==0:
        b//=p;value-=1
    return value


def pmul(a,b,mod=None):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            c[i+j]+=x*y
    if mod:
        c=[x%mod for x in c]
    return c


def ppow(a,n,mod=None):
    b=[1]
    for _ in range(n):
        b=pmul(b,a,mod)
    return b


def peval(a,x,mod):
    out=0
    for c in reversed(a):
        out=(out*x+c)%mod
    return out


def shift(a,c,mod):
    return [sum(a[j]*comb(j,r)*c**(j-r) for j in range(r,len(a)))%mod
            for r in range(len(a))]


def inverse_jets(nodes,m):
    """Literal rational coefficients of each Q_i(x_i+T)^(-m)."""
    alljets=[]
    for i,x in enumerate(nodes):
        b=[Q(1)]+[Q(0)]*(m-1)
        for j,y in enumerate(nodes):
            if i==j:
                continue
            f=[Q((-1)**r*comb(m+r-1,r),(x-y)**(m+r)) for r in range(m)]
            b=[sum(b[s]*f[r-s] for s in range(r+1)) for r in range(m)]
        alljets.append(b)
    return alljets


def loss(nodes,m,p):
    return max(-vp(c,p) for jet in inverse_jets(nodes,m) for c in jet if c)


def full_smith(nodes,m,p):
    M=len(nodes)*m
    J=Matrix([[comb(j,r)*x**(j-r) if j>=r else 0 for j in range(M)]
              for x in nodes for r in range(m)])
    D=smith_normal_form(J,domain=ZZ)
    answer=tuple(vp(D[i,i],p) for i in range(M))
    need(tuple(sorted(answer))==answer,('Smith ordering',nodes,m,p))
    return answer


def fourjet_prediction(nodes,p):
    need(p>=3 and p%2,('odd-prime predictor domain',p))
    depths=sorted(vp(x-y,p) for x,y in combinations(nodes,2))
    e=depths[0]
    need(depths[1]==e,('ultrametric depth shape',nodes,p,depths))
    d=depths[2]-e
    if d:
        return 11*e+7*d-int(p==5)
    if not e:
        return 0
    scale=p**e
    normalized=[Q(x-nodes[0],scale) for x in nodes]
    # An unordered residue triple is an arithmetic progression iff one
    # of its points is the average of the other two.
    arithmetic=any(vp(2*normalized[i]-normalized[j]-normalized[k],7)>=1
                   for i,j,k in ((0,1,2),(1,0,2),(2,0,1))) if p==7 else False
    beta=int(p in (3,5) or (p==7 and arithmetic))
    return 11*e-beta


def p7_partition(e,a):
    need(e>=0 and a%7 not in (0,1),'p7 equilateral partition domain')
    if not e:
        return (0,)*12
    kappa=int(a%7 in (2,4,6))
    return (0,0,0,0,e,2*e,4*e,5*e+1,7*e,8*e,10*e-1+kappa,11*e-kappa)


def minimal_weight_certificate():
    a=symbols('a')
    rows=tuple((x,r) for x in (1,a) for r in range(4))
    slopes=(1,3,7,12,19,27)
    factors={
        (3,):4,
        (7,):4*a,
        (3,7):40*a*(a-1),
        (2,3,7):200*a*(a-1)*(2*a-1),
        (3,6,7):200*a**4*(a-1)*(a-2),
        (2,3,6,7):2100*a**4*(a-1)**4,
        (1,2,3,6,7):560*a**4*(a-1)**4*(7*a*a-7*a+2),
        (2,3,5,6,7):560*a**9*(a-1)**4*(2*a*a-7*a+7),
        (1,2,3,5,6,7):1680*a**9*(a-1)**9,
    }
    seen=set()
    counts=[]
    for r,slope in enumerate(slopes,1):
        count=0
        for I in combinations(range(8),r):
            for J in combinations(range(4,12),r):
                W=sum(J)-sum(rows[i][1] for i in I)
                need(W>=slope,('universal residual minor weight',I,J,W))
                count+=1
                if W==slope:
                    need(J==tuple(range(4,4+r)) and I in factors,('complete equality weight class',I,J))
                    P=expand(Matrix([[comb(j,rows[i][1])*rows[i][0]**(j-rows[i][1]) for j in J] for i in I]).det())
                    need(expand(P-factors[I])==0,('exact minimal-weight factorization',I,J))
                    if r>=4:
                        need(all(int(c)%7==0 for c in Poly(P,a).all_coeffs()),('one-digit coefficient floor',I))
                    seen.add(I)
        counts.append(count)
    need(counts==[64,784,3136,4900,3136,784] and sum(counts)==12804,'whole rank1..6 weight universe')
    need(seen==set(factors),'all9 equality-weight minors and no missing witness')
    for a0 in range(2,49):
        if a0%7 in (0,1):
            continue
        for r in range(1,7):
            values=[vp(int(expand(P).subs(a,a0)),7) if hasattr(P,'subs') else vp(P,7)
                    for I,P in factors.items() if len(I)==r]
            need(min(values)==int(r>=4),('minimal-weight joint attainment control',r,a0))
    count=0
    for a0 in range(2,49):
        if a0%7 in (0,1):
            continue
        for e in (0,1,2,4,7):
            need(full_smith((0,7**e,a0*7**e),4,7)==p7_partition(e,a0),('full all-depth formula control',e,a0))
            count+=1
    print('P7_MINIMAL_WEIGHT ranks1..6',counts,'total12804; only9 equality-weight minors expanded')
    print('P7_FULL_PARTITION controls',count,'residue bit AP; e0 mask separately retained')


def cubic_identities():
    a,u,v=symbols('a u v')
    C0=(1+a)*(1+a+a*a)
    C1=(2-a)*(a*a-3*a+3)
    Ca=(2*a-1)*(3*a*a-3*a+1)
    need(expand(C0+C1-7*(a*a-a+1))==0,'cubic companion sum')
    need(expand(Ca-6*C0+7*(3*a*a+a+1))==0,'third cubic companion')
    need(expand(C0-((a+3)*(a*a-a+1)+2*(2*a-1)))==0,'quadratic reduction')
    need(resultant(C0,C1,a)==4116,'cubic resultant 2^2*3*7^3')
    direct=sum(comb(3+j,j)*comb(6-j,3-j)*u**j*v**(3-j) for j in range(4))
    need(expand(direct-20*(u+v)*(u*u+u*v+v*v))==0,'top reciprocal factorization')
    # At p=3 the two remaining quadratic factors are units on a=2.
    need((2*2+2+1)%3==1 and (2*2-3*2+3)%3==1,'ternary linear cap factors')
    roots=[]
    for r in range(2,7):
        cs=[int(f.subs(a,r))%7 for f in (C0,C1,Ca)]
        if cs==[0,0,0]:
            roots.append(r)
            need((r*r-r+1)%7!=0,('septenary one-digit cap',r))
    need(roots==[2,4,6],'septenary exact common roots')
    for p in (3,5,7,11,13,17,19):
        for r in range(2,p*p):
            if r%p in (0,1):
                continue
            actual=min(vp(20*int(f.subs(a,r)),p) for f in (C0,C1,Ca))
            beta=int(p in (3,5) or (p==7 and r%7 in (2,4,6)))
            need(actual==beta,('complete small lift cubic cap',p,r,actual,beta))
    print('CUBIC_IDENTITIES resultant4116; exact joint cap beta=1 atp3,p5,or p7 AP, otherwise0')


def fourjet_controls():
    direct=0
    literal=0
    for p in (3,5,7,11,13):
        for a in range(2,p*p):
            if a%p in (0,1):
                continue
            for e in (0,1,3):
                nodes=(0,p**e,p**e*a)
                predicted=fourjet_prediction(nodes,p)
                need(loss(nodes,4,p)==predicted,('equilateral inverse loss',p,e,a))
                direct+=1
            if 2<=a<=min(p-1,4):
                nodes=(0,p,p*a)
                need(full_smith(nodes,4,p)[-1]==fourjet_prediction(nodes,p),('equilateral literal Smith',p,a))
                literal+=1
        for e in (0,1,3,5):
            for d in (1,2,4):
                for unit in (1,1+p,-1):
                    nodes=(0,p**e,p**(e+d)*unit)
                    predicted=fourjet_prediction(nodes,p)
                    need(loss(nodes,4,p)==predicted,('close inverse loss',p,e,d,unit))
                    direct+=1
                    if unit==1 and e in (0,1,3) and d in (1,2):
                        need(full_smith(nodes,4,p)[-1]==predicted,('close literal Smith',p,e,d))
                        literal+=1
    twins=((0,7,14),(0,7,21))
    lists=((0,0,0,0,1,2,4,6,7,8,10,10),(0,0,0,0,1,2,4,6,7,8,9,11))
    for nodes,expected in zip(twins,lists):
        actual=full_smith(nodes,4,7)
        literal+=1
        need(actual==expected,('minimal displayed twins',nodes))
        need(sum(actual)==48,'same determinant valuation')
        for unit,translation in ((-1,4),(2,-11)):
            changed=tuple(unit*x+translation for x in reversed(nodes))
            need(full_smith(changed,4,7)==expected,('signed affine/permutation control',changed))
            literal+=1
    need(sum(min(10,s) for s in lists[0])==48 and sum(min(10,s) for s in lists[1])==47,
         'finite-precision kernel exponents48versus47')
    for a,last in ((2,(30,32)),(3,(29,33))):
        nodes=(0,343,343*a)
        full=full_smith(nodes,4,7)
        literal+=1
        need(full==(0,0,0,0,3,6,12,16,21,24)+last,('higher-depth full Smith control',a))
    try:
        fourjet_prediction((0,2,4),2)
    except RuntimeError:
        pass
    else:
        raise RuntimeError('dyadic misuse must be rejected')
    print('FOURJET_CONTROLS rational',direct,'literal12x12Smith',literal)
    print('P7_TWINS (0,7,14)',lists[0],'; (0,7,21)',lists[1])
    print('P7_TWINS_KERNEL modulo7^10 exponent48versus47; determinant48both')


def mod_jet(nodes,m,p):
    answer=[]
    for i,x in enumerate(nodes):
        b=[1]+[0]*(m-1)
        for j,y in enumerate(nodes):
            if i==j:
                continue
            v=pow((x-y)%p,-1,p)
            f=[(-1)**r*comb(m+r-1,r)*pow(v,m+r,p)%p for r in range(m)]
            b=[sum(b[s]*f[r-s] for s in range(r+1))%p for r in range(m)]
        answer.append(b)
    return answer


def coefficient_packet(nodes,p):
    k=(p-1)//2
    F=[1]
    for x in nodes:
        F=pmul(F,[-x,1],p)
    A=ppow(F,k,p)
    packet=[A[j] for j in range(p-1,len(A),p)]
    return F,A,packet


def general_packet_controls():
    count=0
    saturated=0
    cancelled=0
    for p in (3,5,7,11,13):
        k=(p-1)//2;m=k+1
        # Complete subsets through p7; deterministic prefix/affine controls
        # thereafter. This is a control universe, not the theorem premise.
        if p<=7:
            universe=[S for n in range(2,p+1) for S in combinations(range(p),n)]
        else:
            universe=sorted(set(tuple(sorted({(u*j+c)%p for j in range(n)}))
                                for n in range(2,p+1) for u,c in ((1,0),(2,1),(3,2))))
        for nodes in universe:
            F,A,packet=coefficient_packet(nodes,p)
            jets=mod_jet(nodes,m,p)
            q=len(packet)
            need(q<len(nodes),('packet evaluation injectivity degree',p,nodes,q))
            for i,x in enumerate(nodes):
                Qi=1
                for y in nodes:
                    if y!=x:
                        Qi=Qi*(x-y)%p
                direct=shift(A,x,p)[p-1]
                predicted=peval(packet,x,p)
                need(direct==predicted,('translation sparse coefficient identity',p,nodes,x))
                need(jets[i][k]==pow(Qi,-p,p)*direct%p,('reciprocal versus coefficient packet',p,nodes,x))
            live=any(packet)
            need(live==any(j[k] for j in jets),('simultaneous packet detection',p,nodes))
            saturated+=int(live);cancelled+=int(not live);count+=1
            if len(nodes) in (2,3,4) and nodes==tuple(range(len(nodes))):
                actual=loss(tuple(p*x for x in nodes),m,p)
                ceiling=(len(nodes)*m-1)
                need((actual==ceiling)==live,('actual precision ceiling control',p,nodes,actual))
        _,_,packet=coefficient_packet(tuple(range(p)),p)
        need(not any(packet),('complete residue class cancellation hostile',p))
    print('GENERAL_PACKET controls',count,'saturated',saturated,'cancelled',cancelled)
    print('GENERAL_PACKET complete-residue F=X^p-X has zero packet; arbitrary exact loss not inferred')


def deuring_controls():
    for p in (3,5,7,11,13,19,23,31):
        k=(p-1)//2;m=k+1
        H=[comb(k,j)**2%p for j in range(k+1)]
        need(H==list(reversed(H)),('reciprocal S3 identity',p))
        reflected=[sum(H[j]*comb(j,r)*(-1)**r for j in range(r,len(H)))%p for r in range(len(H))]
        need(reflected==[((-1)**k*c)%p for c in H],('reflection S3 identity',p))
        roots=[]
        ordinary=[]
        for a in range(2,p):
            _,_,packet=coefficient_packet((0,1,a),p)
            ha=peval(H,a,p)
            need(packet==[((-1)**k*ha)%p],('Deuring scalar coefficient',p,a))
            jets=mod_jet((0,1,a),m,p)
            need((not any(row[k] for row in jets))==(ha==0),('three-node simultaneous top jet',p,a))
            # Independent finite elliptic point-count check; interpretation
            # is cited separately and never used by the algebraic proof.
            points=1+sum(1+int((x*(x-1)*(x-a))%p!=0)*(1 if pow((x*(x-1)*(x-a))%p,k,p)==1 else -1)
                         for x in range(p))
            need((points==p+1)==(ha==0),('finite Legendre point-count control',p,a,points))
            (roots if not ha else ordinary).append(a)
        if p%4==3:
            need(p-1 in roots,'odd-degree arithmetic-progression Deuring root')
            if p>=7:
                need(bool(ordinary),'ordinary residue beyond degree bound')
        print('DEURING',p,'m',m,'roots',roots,'ordinary_count',len(ordinary))


if __name__=='__main__':
    cubic_identities()
    fourjet_controls()
    minimal_weight_certificate()
    general_packet_controls()
    deuring_controls()
    print('PASS_OPTIMIZATION_LIVE_GATES',GATES)
