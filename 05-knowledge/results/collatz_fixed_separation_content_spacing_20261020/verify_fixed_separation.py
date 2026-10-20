#!/usr/bin/env python3
"""Independent exact physical/gcd controls for h>=3. No repository imports.
Exhaustive: 6<=p<=40, every primitive r and 3<=h<=floor(p/2).
Additional: p=43,61,89, every primitive r, h in {3,4,5,7,11,floor(p/2)}.
Resultants: complete exhaustive universe through p=18.
"""
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from itertools import product
from math import comb, gcd
from pathlib import Path
import json

checks = 0

def check(ok, why):
    global checks
    checks += 1
    if not ok:
        raise ArithmeticError(why)


def bezout(a,b):
    old_r,r,old_s,s,old_t,t=a,b,1,0,0,1
    while r:
        k=old_r//r
        old_r,r=r,old_r-k*r
        old_s,s=s,old_s-k*s
        old_t,t=t,old_t-k*t
    check(old_r==1,'primitive Bezout clock')
    return old_s,old_t


@lru_cache(None)
def gate(h,t,d1,d2):
    T=2**(h+t);W=3**h
    alpha=lambda d:Q(1) if d==1 else Q(-1,2)
    K=-1/(alpha(d1)+alpha(d2)*Q(W,T))
    four={(1,1):Q(-T,T+W),(-1,-1):Q(2*T,T+W),
          (1,-1):Q(-2*T,2*T-W),(-1,1):Q(-T,W-T//2)}
    check(K==four[d1,d2],'four gate formulas')
    eps=1 if K>0 else -1
    delta=abs(K.numerator).bit_length()-1;J=K.denominator
    check(abs(K.numerator)==2**delta,'exact dyadic numerator')
    if h>=3:
        check(J>1 and gcd(J,6)==1,'external denominator prime')
        check(h<=delta<=2*h+1,'dyadic exponent range')
    vertices=[]
    for a,b,c in product((-2,2),(-1,1),(0,1)):
        e=a-b+delta*c
        vertices.append(2**(max(-e,0)+c)*3**max(-b,0)*J**c
                        +2**max(e,0)*3**max(b,0))
    C=6*J+2**(delta+3)
    if delta>=3:
        check(max(vertices)==C,'exact coefficient-box maximum')
    if h>=3:
        check(C<=34*4**h,'uniform coefficient constant')
    return eps,delta,J,C


def physical(source,positions,moves):
    """Move named initial source ones, checking the current target bit."""
    bits=list(source);pos=list(positions);q=len(bits);edges=[]
    for label,direction in moves:
        old=pos[label];new=(old+direction)%q
        if bits[old]!=1 or bits[new]!=0:return None
        bits[old],bits[new]=0,1;pos[label]=new
        edges.append(old if direction==1 else new)
    return tuple(bits),tuple(edges)


def carries(word):
    numerator=0
    for i,bit in enumerate(word):
        numerator=(3 if bit else 1)*numerator+bit*2**i
    count=0;direct=0
    for i in range(len(word)-1,-1,-1):
        if word[i]:direct+=2**i*3**count;count+=1
    check(numerator==direct,'affine versus positional carry')
    return numerator


@lru_cache(None)
def short_pair(p,r,u):
    first=None
    for H in range(1,p+1):
        for b,c,m in product(range(-H,H+1),repeat=3):
            v=(b,c,m)
            if max(map(abs,v))!=H or (m-b*r-c*u)%p:continue
            if next(x for x in v if x)<0:continue
            if first is None:first=v;continue
            if any(first[i]*v[j]!=first[j]*v[i] for i,j in ((0,1),(0,2),(1,2))):
                H1=max(map(abs,first))
                check((H1*H)**3<=p*p,'successive-minimum product')
                return first,v
    raise ArithmeticError('missing independent pair')


def polynomial(p,r,u,eps,delta,J,C,v):
    b,c,m=v
    if c<0:b,c,m=-b,-c,-m
    a,remainder=divmod(m-b*r-c*u,p)
    check(remainder==0,'relation integrality')
    H=max(map(abs,v));e=a-b+delta*c
    A=2**max(-e,0)*3**max(-b,0)*J**c
    B=eps**c*2**max(e,0)*3**max(b,0)
    out=[0]*(max(max(m,0)+c,max(-m,0))+1)
    for j in range(c+1):out[max(m,0)+j]+=A*comb(c,j)*(-1)**(c-j)
    out[max(-m,0)]-=B
    while len(out)>1 and out[-1]==0:out.pop()
    check(abs(a)<=2*H,'coefficient a bound')
    check(out[0]!=0 and sum(out)!=0,'boundary roots absent')
    check((len(out)==1)==(m==0 and c==0),'constant characterization')
    check(len(out)-1<=2*H and sum(map(abs,out))<=C**H,'degree and l1 bounds')
    return out


def remainder(f,g):
    out=list(map(Q,f));g=list(map(Q,g))
    while len(out)>=len(g) and any(out):
        s=len(out)-len(g);k=out[-1]/g[-1]
        for j,x in enumerate(g):out[s+j]-=k*x
        while len(out)>1 and out[-1]==0:out.pop()
    return out


def resultant(f,g):
    """Exact rational Euclidean recurrence, unlike inherited Sylvester code."""
    m,n=len(f)-1,len(g)-1
    if n==0:return Q(g[0])**m
    if m<n:return (-1)**(m*n)*resultant(g,f)
    rem=remainder(f,g)
    if not any(rem):return Q(0)
    k=len(rem)-1
    return (-1)**(m*n)*Q(g[-1])**(m-k)*resultant(g,rem)


def evaluate(f,x,M):
    out=0
    for a in reversed(f):out=(out*x+a)%M
    return out


def run():
    counts=dict(shapes=0,separation_shapes=0,legal_presentations=0,
                illegal_presentations=0,ordered_physical_paths=0,half_turn_presentations=0,
                half_turn_ordered_paths=0,cut_crossing_paths=0,positive_clock_presentations=0,
                negative_clock_presentations=0,resultant_checks=0,integral_presentations=0)
    by_h={};max_content=None;prime_power=None;illegal_controls=[];half_turn_control=None
    digest=sha256();resultant_digest=sha256();largest_bound_ratio=Q(0)
    for h in range(3,65):
        for t in range(h+1):
            for d1,d2 in product((1,-1),repeat=2):
                eps,delta,J,C=gate(h,t,d1,d2)
                largest_bound_ratio=max(largest_bound_ratio,Q(C,34*4**h))
    # h=2 is an explicit hostile to a no-smooth-gate claim without h>=3.
    smooth_controls=[]
    for t in range(3):
        for d1,d2 in product((1,-1),repeat=2):
            eps,delta,J,C=gate(2,t,d1,d2)
            if J==1:smooth_controls.append(dict(h=2,t=t,directions=[d1,d2],K=eps*2**delta))
    check(smooth_controls==[{'h':2,'t':0,'directions':[1,-1],'K':8},
                           {'h':2,'t':2,'directions':[-1,1],'K':-16}], 'smooth boundary controls')
    for p in list(range(6,41))+[43,61,89]:
        hs=range(3,p//2+1) if p<=40 else sorted({3,4,5,7,11,p//2})
        for r in range(1,p):
            if gcd(p,r)!=1:continue
            counts['shapes']+=1;q=p+r;D=2**q-3**p;M=abs(D)
            source=tuple(((i+1)*p)//q-(i*p)//q for i in range(q))
            # Independently enumerate physical ones in order, then attach labels 1,...,p.
            positions=tuple(i for i,bit in enumerate(source) if bit)
            residues=tuple((-label*r)%p for label in range(1,p+1))
            check(len(positions)==p,'source weight')
            A,B=bezout(p,r)
            z=(pow(2,A,M)*pow(3*pow(2,-1,M)%M,B,M))%M
            check(pow(z,p,M)==2%M and (2*pow(z,r,M)-3)%M==0,'whole clock')
            check(gcd(M,z*(z-1)*6)==1,'clock units')
            for h in hs:
                counts['separation_shapes']+=1;physical_paths=set();presentations=0
                for j in range(p):
                    k=(j-h)%p;u,v=residues[j],residues[k];t=(u+h*r)//p
                    check(0<=t<=h and v==u+h*r-t*p,'residue transport')
                    check((pow(z,v,M)*2**(h+t)-pow(z,u,M)*3**h)%M==0,'whole clock ratio')
                    for d1,d2 in product((1,-1),repeat=2):
                        forward=physical(source,positions,[(j,d1),(k,d2)])
                        reverse=physical(source,positions,[(k,d2),(j,d1)])
                        legal=(u<r if d1==1 else u>=p-r) and (v<r if d2==1 else v>=p-r)
                        check((forward is not None)==legal and (reverse is not None)==legal,'legal iff both physical orders')
                        if not legal:
                            counts['illegal_presentations']+=1
                            if len(illegal_controls)<4:illegal_controls.append(dict(p=p,r=r,h=h,u=u,v=v,directions=[d1,d2]))
                            continue
                        output,edges=forward;out2,edges2=reverse
                        check(output==out2 and edges2==edges[::-1],'commutation and edge order')
                        e1,e2=edges
                        check(not ({e1,(e1+1)%q}&{e2,(e2+1)%q}),'actual edges disjoint')
                        physical_paths.update((edges,edges2));presentations+=1
                        counts['legal_presentations']+=1
                        counts['positive_clock_presentations' if D>0 else 'negative_clock_presentations']+=1
                        by_h[h]=by_h.get(h,0)+1
                        carry=carries(output);g=gcd(M,carry)
                        eps,delta,J,C=gate(h,t,d1,d2)
                        raw=(2+(z-1)*((2 if d1==1 else -1)*pow(z,u,M)+(2 if d2==1 else -1)*pow(z,v,M)))%M
                        reduced=(J*pow(z,u,M)*(z-1)-eps*2**delta)%M
                        check(gcd(M,raw)==gcd(M,reduced)==g,'whole content iff identities')
                        check(gcd(g,6*J)==1,'cleared denominators remain units')
                        row=dict(p=p,r=r,h=h,u=u,v=v,t=t,directions=[d1,d2],g=g,edges=list(edges),J=J,delta=delta,epsilon=eps)
                        if max_content is None or g>max_content['g']:max_content=row
                        if prime_power is None:
                            for ell in [5,7,11,13,17,19,23,29,31]:
                                if g%(ell**3)==0:
                                    prime_power=row|{'prime_cube_divisor':ell**3};break
                        counts['integral_presentations']+=(g==M)
                        digest.update(json.dumps([p,r,h,j,d1,d2,edges,carry,g,J,delta,eps],separators=(',',':')).encode())
                        if p<=18:
                            pair=short_pair(p,r,u)
                            f,F=[polynomial(p,r,u,eps,delta,J,C,w) for w in pair]
                            check(evaluate(f,z,g)==evaluate(F,z,g)==0,'polynomial roots modulo full gcd')
                            R=resultant(f,F)
                            check(R.denominator==1 and R!=0 and R.numerator%g==0,'nonzero full-content resultant')
                            heights=[max(map(abs,w)) for w in pair]
                            check(abs(R)<=C**(4*heights[0]*heights[1]),'resultant l1 bound')
                            counts['resultant_checks']+=1
                            resultant_digest.update(json.dumps([p,r,h,u,d1,d2,pair,f,F,R.numerator,g],separators=(',',':')).encode())
                expected=presentations if 2*h==p else 2*presentations
                check(len(physical_paths)==expected,'half-turn duplication exactly two, otherwise one')
                counts['ordered_physical_paths']+=len(physical_paths)
                counts['cut_crossing_paths']+=sum(q-1 in e for e in physical_paths)
                if 2*h==p:
                    counts['half_turn_presentations']+=presentations
                    counts['half_turn_ordered_paths']+=len(physical_paths)
                    if half_turn_control is None and presentations:
                        half_turn_control=dict(p=p,r=r,h=h,oriented_presentations=presentations,ordered_paths=len(physical_paths))
    check(prime_power is not None,'must include a prime-cube control')
    check(counts['integral_presentations']==0,'finite universe unexpectedly integral')
    # Separate full-prime-cube, constant-polynomial, and dependent-projection controls.
    p,r,u,h,t,d1,d2=19,14,2,3,2,1,-1
    eps,delta,J,C=gate(h,t,d1,d2)
    M=abs(2**(p+r)-3**p);g=125
    A,B=bezout(p,r)
    z=pow(2,A,M)*pow(3*pow(2,-1,M)%M,B,M)%M
    controlled_pairs=[short_pair(p,r,u),((0,1,u),(p,1,u)),((p,0,0),(0,1,u))]
    special_resultants=[]
    for pair in controlled_pairs:
        f,F=[polynomial(p,r,u,eps,delta,J,C,w) for w in pair]
        R=resultant(f,F)
        check(evaluate(f,z,g)==evaluate(F,z,g)==0,'special full-prime-cube roots')
        check(R.denominator==1 and R!=0 and R.numerator%g==0,'special whole-gcd resultant')
        special_resultants.append(dict(vectors=pair,polynomials=[f,F],
                                       resultant=R.numerator,gcd_with_clock=gcd(M,R.numerator),
                                       actual_g=g))
    check(len(special_resultants[-1]['polynomials'][0])==1,'nonzero constant polynomial retained')
    check(special_resultants[-1]['gcd_with_clock']==M>g,'resultant gate is not an iff')
    return dict(status='FINITE-EXACT independent controls, not an all-length finite-residual closure',
                universe={'exhaustive_p':[6,40],'exhaustive_h':'3<=h<=floor(p/2)',
                          'extra_p':[43,61,89],'extra_h':'3,4,5,7,11,floor(p/2)',
                          'r':'1<=r<p, gcd(p,r)=1','all_four_directions':True,
                          'resultant_exhaustive_p':[6,18],'algebraic_gate_checks_h':[3,64]},
                counts=counts,legal_presentations_by_h=by_h,maximum_content=max_content,
                prime_cube_control=prime_power,illegal_physical_controls=illegal_controls,
                half_turn_control=half_turn_control,h_two_smooth_controls=smooth_controls,
                special_resultant_controls=special_resultants,
                maximum_tested_coefficient_ratio=str(largest_bound_ratio),
                physical_digest_sha256=digest.hexdigest(),resultant_digest_sha256=resultant_digest.hexdigest(),
                checks=checks,script_sha256=sha256(Path(__file__).read_bytes()).hexdigest())

if __name__=='__main__':
    print(json.dumps(run(),sort_keys=True,indent=2))
