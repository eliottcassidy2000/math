#!/usr/bin/env python3
"""Exact replay: the dyadic fruit lift and sorted physical row envelopes."""
from fractions import Fraction as Q
from math import gcd,isqrt,lcm
from random import Random

def need(ok, message):
    if not ok:
        raise ArithmeticError(message)

def F(q,t):
    return t**3-t**2+q*t-(6*q-1)/7

def D(q):
    return (q+1)*(1+29*q-196*q*q)/49

def additions(a,b):
    if a is None: return b
    if b is None: return a
    x,y=a; xx,yy=b
    if x==xx and y==-yy: return None
    s=(yy-y)/(xx-x) if x!=xx else (3*x*x+842*x+43904)/(2*y)
    X=s*s-421-x-xx
    return X,s*(x-X)-y

def convergents(a,b):
    p0,p1,q0,q1=0,1,1,0
    while b:
        t,d=divmod(a,b)
        p0,p1=p1,t*p1+p0
        q0,q1=q1,t*q1+q0
        yield p1,q1
        a,b=b,d

def lift(coeff,r,bits):
    den=lcm(*(c.denominator for c in coeff))
    need(den%2==1,"nonintegral normalized polynomial")
    coeff=tuple(int(c*den) for c in coeff)
    def H(t,p):
        out=0
        for c in coeff:
            out=(out*t+c)%p
        return out
    need(H(r,2)==0,"missing residue root")
    for k in range(1,bits):
        if H(r,1<<(k+1)):
            r+=1<<k
            need(H(r,1<<(k+1))==0,"nonsimple Hensel digit")
    return r

def invert(q):
    a,d=q.numerator,q.denominator
    if d%2 and (a+d)%16:
        return None
    if D(q)==0: return None
    M,B=14*(d+abs(a)),7*d
    bits=(2*M*B).bit_length()
    P=1<<bits
    if d%2==0:
        s=lift((1/q,-1/q,Q(1),-(6-1/q)/7),0,bits)
    else:
        e=q+1
        s=(1+2*lift((Q(1),Q(1),e/4,e/56),1,bits))%P
    for j,b in convergents(s,P):
        p=s*b-P*j
        if abs(p)>M or b>B or b%2==0:
            continue
        t=Q(p,b)
        if F(q,t): continue
        c=t-1
        h=q+c*t
        v=c*c-4*h
        root=Q(isqrt(v.numerator),isqrt(v.denominator))
        need(root*root==v,"rational quadratic")
        return tuple(sorted((t,(-c+root)/2,(-c-root)/2)))
    return None

def fruit():
    H=(Q(-196),Q(196))
    P=None
    for n in range(1,37):
        P=additions(P,H)
        for e in (0,1):
            pt=additions(P,(Q(0),Q(0))) if e else P
            q=-pt[0]/196-1
            aa,dd=q.numerator,q.denominator
            need(D(q)==(pt[1]/(196*7))**2,"square discriminant")
            dyad=(dd%2==0 or (aa+dd)%16==0)
            triad=(dd%3==0 or (aa+dd)%3==0)
            need(dyad==triad==(n%3==0),"dyadic/ternary classification")
            result=invert(q)
            need((result is not None)==dyad,"absolute root inverse")
            if result:
                need(sum(result)==1 and all(F(q,t)==0 for t in result),"cubic replay")
                need(sum(t/(1-t) for t in result)==4,"fruit equation")
                if n==27 and not e:
                    need(min(result)>0,"positive height control")
    print("fruit: 72 exact quotient points and both label translates")

def p(n):
    z=0
    while n:
        v=n>>1
        if n%2==0:
            z+=v+v%2
        n=v
    return z

def G(n,c,k,phi,b):
    return c*n+k*p(n)+phi*(n%2)+b

def words():
    random=Random(559)
    for L in range(2,12):
        N=1<<L
        for _ in range(6):
            c,k,phi,b=(random.randint(1,300) for _ in range(4))
            for e in (0,1):
                brute=max(G(n,c,k,phi,b) for n in range(e,N,2))
                sorted_max=max(G(N-(1<<j)+e,c,k,phi,b) for j in range(1,L+1))
                need(brute==sorted_max,"wrong sorted maximum")
        if L>=8:
            M=1<<(L//2)
            want0=200*N-(94 if L%2==0 else 138)*M+33
            want1=(389*N)//2-(741 if L%2==0 else 1093)*M//8+233
            need(max(G(n,189,11,11,11) for n in range(0,N,2))==want0,"physical even formula")
            need(max(G(n,189,11,11,11) for n in range(1,N,2))==want1,"physical odd formula")
    print("row: exact parity maxima across 60 varied physical-style arithmetic controls and the inherited mask-zero row")
if __name__=="__main__":
    fruit()
    words()
