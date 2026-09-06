#!/usr/bin/env python3
"""Exact tests for normalized reciprocal jets, precision pruning and dilation envelopes.

No repository or earlier producer mathematical imports. The written proof
supplies unbounded scope; these exact tests retain hostile observer changes.
"""
from collections import defaultdict
from fractions import Fraction as Q
from math import comb, prod
import sys
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

sys.stdout.reconfigure(newline="\n")
GATES=0


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def val(x,p):
    x=Q(x)
    if not x:
        return float("inf")
    a,b=abs(x.numerator),x.denominator
    out=0
    while a%p==0:
        a//=p
        out+=1
    while b%p==0:
        b//=p
        out-=1
    return out


def mul(a,b,K):
    out=[Q(0)]*(K+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<=K:
                out[i+j]+=x*y
    return out


def factor(alpha,m,K):
    return [(-1)**r*comb(m+r-1,r)*alpha**r for r in range(K+1)]


def states(nodes,mult,p):
    result=[]
    for i,x in enumerate(nodes):
        K=mult[i]-1
        depths={j:val(x-y,p) for j,y in enumerate(nodes) if j!=i}
        f=max(depths.values(),default=0)
        D=sum(mult[j]*h for j,h in depths.items())
        groups=defaultdict(list)
        for j,h in depths.items():
            groups[h].append(j)
        b=[Q(1)]+[Q(0)]*K
        shells=[]
        for h,js in sorted(groups.items()):
            shell=[Q(1)]+[Q(0)]*K
            for j in js:
                shell=mul(shell,factor(Q(p**h,x-nodes[j]),mult[j],K),K)
            need(all(val(c,p)>=0 for c in shell),("integral shell jet",i,h))
            attenuated=[c*p**((f-h)*r) for r,c in enumerate(shell)]
            b=mul(b,attenuated,K)
            shells.append((h,shell))
        result.append((D,f,b,shells))
    return result


def direct_coefficients(nodes,mult,i):
    # Independently expand Q_i(x_i+T), then invert it coefficient by coefficient.
    K=mult[i]-1
    q=[Q(1)]+[Q(0)]*K
    for j,y in enumerate(nodes):
        if j!=i:
            m=mult[j]
            q=mul(q,[Q(comb(m,r)*(nodes[i]-y)**(m-r)) for r in range(min(m,K)+1)],K)
    a=[1/q[0]]
    for r in range(1,K+1):
        a.append(-sum(q[s]*a[r-s] for s in range(1,r+1))/q[0])
    return q[0],a


def residue(c,p,K):
    modulus=p**K
    need(c.denominator%p!=0,("integral coefficient residue",c,p))
    return c.numerator*pow(c.denominator,-1,modulus)%modulus


def modular_shell_product(state,own_order,p,K):
    D,f,b,shells=state
    if K==0:
        return [0]*(own_order+1)
    modulus=p**K
    out=[1]+[0]*own_order
    for h,shell in shells:
        coefficients=[]
        for r,c in enumerate(shell):
            attenuation=r*(f-h)
            coefficients.append(0 if attenuation>=K else
                                residue(c,p,K-attenuation)*p**attenuation)
        new=[0]*(own_order+1)
        for i,x in enumerate(out):
            for j,y in enumerate(coefficients):
                if i+j<=own_order:
                    new[i+j]=(new[i+j]+x*y)%modulus
        out=new
    return out


def state_loss(state,p):
    return max(D+r*f-val(c,p) for D,f,b,shells in state for r,c in enumerate(b))


def profile(state,mult,p):
    M=sum(mult)
    out={}
    for m,(D,f,b,shells) in zip(mult,state):
        for r,c in enumerate(b):
            if c:
                d=M-m+r
                intercept=D+r*f-val(c,p)
                out[d]=max(out.get(d,-float("inf")),intercept)
    need(len(out)<=max(mult),"at most max multiplicity distinct dilation slopes")
    return dict(sorted(out.items()))


def exact_case(name,nodes,mult,p,matrix=True):
    st=states(nodes,mult,p)
    for i,(D,f,b,shells) in enumerate(st):
        q0,a=direct_coefficients(nodes,mult,i)
        need(val(q0,p)==D,("metric constant denominator",name,i))
        need(b==[c*q0*p**(r*f) for r,c in enumerate(a)],("annular jet versus direct inverse",name,i))
        moments={s:sum(mult[j]*Q(p**f,nodes[i]-y)**s for j,y in enumerate(nodes) if j!=i)
                 for s in range(1,mult[i])}
        for r in range(1,mult[i]):
            need(r*b[r]==sum((-1)**s*moments[s]*b[r-s] for s in range(1,r+1)),
                 ("exact power-moment recurrence with division precision visible",name,i,r))
    L=state_loss(st,p)
    baseline=max(s[0] for s in st)
    recovered=baseline
    for i,(D,f,b,shells) in enumerate(st):
        K=max(0,D+(mult[i]-1)*f-baseline)
        modular=modular_shell_product(st[i],mult[i]-1,p,K)
        if K:
            need(modular==[residue(c,p,K) for c in b],("attenuated modular shell product",name,i))
        for r in range(1,mult[i]):
            threshold=max(0,D+r*f-baseline)
            if threshold:
                z=modular[r]%(p**threshold)
                if z:
                    recovered=max(recovered,D+r*f-val(z,p))
    need(recovered==L,("metric-baseline adaptive precision",name))
    envelope=profile(st,mult,p)
    for e in (0,1,3,7):
        moved=tuple(p**e*x for x in nodes)
        direct=max(-val(c,p) for i in range(len(nodes))
                   for c in direct_coefficients(moved,mult,i)[1])
        predicted=max(intercept+d*e for d,intercept in envelope.items())
        need(direct==predicted,("all-order dilation envelope control",name,e))
    if matrix:
        M=sum(mult)
        H=Matrix([[comb(j,r)*x**(j-r) if j>=r else 0 for j in range(M)]
                  for x,m in zip(nodes,mult) for r in range(m)])
        S=smith_normal_form(H,domain=ZZ)
        largest=abs(int(S[M-1,M-1]))
        need(val(largest,p)==L,("literal integer Hasse Smith",name))
        determinant=prod(abs(nodes[i]-nodes[j])**(mult[i]*mult[j])
                         for i in range(len(nodes)) for j in range(i))
        need(prod(abs(int(S[i,i])) for i in range(M))==determinant,("Hasse determinant",name))
    if len(nodes)>=3:
        old=states(nodes[:-1],mult[:-1],p)[0]
        new=st[0]
        K=mult[0]-1
        rescaled=[c*p**((new[1]-old[1])*r) for r,c in enumerate(old[2])]
        expected=mul(rescaled,factor(Q(p**new[1],nodes[0]-nodes[-1]),mult[-1],K),K)
        need(expected==new[2],("node-adjoining update including new nearest depth",name))
    print('CASE',name,'p',p,'loss',L,'profile',envelope)
    return st,envelope


def terminal_b2(nodes,x):
    reciprocal=[Q(2,x-y) for y in nodes if y!=x]
    return Q(3,2)*(3*sum(reciprocal)**2+sum(c*c for c in reciprocal))


def terminal_pair_controls():
    nodes=(0,2,-7,-5,1,7,9)
    coeffs=[terminal_b2(nodes,x) for x in (0,2)]
    need(coeffs==[Q(1452032,33075)]*2,"seven-node symmetric finite hostile")
    need([val(c,2) for c in coeffs]==[11,11],"three-node pair budget five fails with outside nodes")
    print('TERMINAL_PAIR_HOSTILE nodes',nodes,'b2',coeffs[0],'both valuations11')
    for h in range(2,13):
        r=2**h-1
        K=h+2
        outside=tuple(1+2**K*j for j in range(r))
        nodes=(0,2)+outside
        precision=h+4
        modulus=2**(precision+1)
        expected=6*(r+1)*(3*r+1)
        for x in (0,2):
            rs=[-1 if x==0 else 1]+[2*pow((x-z)%modulus,-1,modulus)%modulus for z in outside]
            numerator=3*(3*sum(rs)**2+sum(c*c for c in rs))
            need(numerator%2==0,"quadratic normalization division is exact")
            b2=(numerator//2)%(2**precision)
            need(val(b2,2)==h+2,("unbounded simultaneous budget family",h,x))
            need((b2-expected)%(2**(h+3))==0,("proved simultaneous congruence",h,x))
            need((-3*sum(rs))%2==1 and max(3,4,3-h)==4,
                 ("actual pair-local loss is exactly four at the displayed depth",h,x))
            if h<=5:
                need(val(terminal_b2(nodes,x),2)==h+2,("literal rational family control",h,x))
        need(3*(r-1)*K>5,"target firewall: outsider order-zero metric baseline dominates the pair")
    print('UNBOUNDED_PAIR_BUDGET h>=2 r=2^h-1 outsiders1+2^(h+2)j; v2(b2)=h+2; controls h2..12')
    print('PAIR_FIREWALL displayed-depth local loss4, uncorrected cost<=5; outside baseline >=3(r-1)(h+2)>5')


def formal_order_control():
    A=[Q(1,11),Q(1,231)]
    B=[Q(1,15),Q(1,35)]
    ja=mul(factor(A[0],3,2),factor(A[1],3,2),2)
    jb=mul(factor(B[0],3,2),factor(B[1],3,2),2)
    need(ja[:2]==jb[:2] and ja[2]!=jb[2],"first-order branch state does not determine the second-order branch")
    need(ja==[1,Q(-2,7),Q(947,17787)] and jb==[1,Q(-2,7),Q(179,3675)],"exact geometric branch witness")
    print('FORMAL_ORDER_HOSTILE same b0,b1',ja[:2],'different b2',ja[2],jb[2])
    packets=[]
    for outsider in (1,3):
        rs=[Q(-1),Q(-2,outsider)]
        P1=3*sum(rs)
        P2=3*sum(c*c for c in rs)
        b2=(P1*P1+P2)/2
        packets.append((val(P1,2),val(P2,2),val(b2,2)))
    need(packets==[(0,0,4),(0,0,2)],"moment valuations alone destroy the quadratic cancellation")
    print('MOMENT_VALUATION_HOSTILE (vP1,vP2,vb2)',packets)


def main():
    print('SCOPE arbitrary positive multiplicities; fixed-anchor tree-shell jets; exact adaptive precision and all-depth envelope')
    cases=[('uniform-A',(0,1,2),(3,3,3),2),('uniform-B',(1,0,3),(3,3,3),2),
           ('mixed-A',(0,2,1),(2,2,1),2),('mixed-B',(1,3,0),(2,2,1),2),
           ('deep-pair',(0,1,8),(3,3,3),2),
           ('ternary-twojet-A',(0,9,27,81),(2,2,2,2),3),
           ('ternary-twojet-B',(0,9,54,81),(2,2,2,2),3),
           ('unequal-prime5',(-3,2,7,12),(1,3,2,1),5),
           ('seven-node-pair',(0,2,-7,-5,1,7,9),(3,3,3,3,3,3,3),2),
           ('singleton',(-7,),(5,),2)]
    envelopes={}
    for name,nodes,mult,p in cases:
        st,envelope=exact_case(name,nodes,mult,p,matrix=sum(mult)<=12)
        envelopes[name]=envelope
    need(envelopes['uniform-A']=={6:3,7:4,8:1} and envelopes['uniform-B']=={6:3,7:4,8:3},
         "dyadic uniform hostile Gamma recovered as order-two intercept")
    need(envelopes['mixed-A']=={3:2,4:1} and envelopes['mixed-B']=={3:2,4:0},
         "heterogeneous warning retained in slope merge")
    a=envelopes['uniform-A']
    need(max(c+4*d for d,c in a.items())==33 and max(c+4*d for d,c in a.items() if d<8)==32,
         "current-depth pruning cannot be made permanent under dilation")
    print('PRUNING_FIREWALL uniform A order2 is inactive now but gives33 at e4 versus lower-order32')
    terminal_pair_controls()
    formal_order_control()
    print('PASS_OPTIMIZATION_LIVE_GATES',GATES)


if __name__=='__main__':
    main()
