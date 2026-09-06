#!/usr/bin/env python3
"""Exact bounded controls for the improved global duplication stability constant.

All dimensions are handled by the companion analytic envelope proof.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from math import isqrt, comb
from hashlib import sha256

@dataclass(frozen=True)
class QR:
    d: F
    a: F = F(0)
    b: F = F(0)

    def __post_init__(self):
        d,a,b=F(self.d),F(self.a),F(self.b)
        sn,sd=isqrt(d.numerator),isqrt(d.denominator)
        if sn*sn==d.numerator and sd*sd==d.denominator:
            a,b=a+b*F(sn,sd),F(0)
        for key,value in [('d',d),('a',a),('b',b)]:
            object.__setattr__(self,key,value)

    def cast(self,x):
        if not isinstance(x,QR): return QR(self.d,x)
        if x.d==self.d: return x
        if x.b==0: return QR(self.d,x.a)
        raise ValueError('mixed quadratic fields')

    def __add__(self,x):
        x=self.cast(x);return QR(self.d,self.a+x.a,self.b+x.b)
    __radd__=__add__
    def __neg__(self):return QR(self.d,-self.a,-self.b)
    def __sub__(self,x):return self+-self.cast(x)
    def __rsub__(self,x):return self.cast(x)+-self
    def __mul__(self,x):
        x=self.cast(x)
        return QR(self.d,self.a*x.a+self.d*self.b*x.b,self.a*x.b+self.b*x.a)
    __rmul__=__mul__
    def __truediv__(self,x):
        x=self.cast(x);den=x.a*x.a-self.d*x.b*x.b
        if den==0:raise ZeroDivisionError
        return self*QR(self.d,x.a/den,-x.b/den)
    def __rtruediv__(self,x):return self.cast(x)/self
    def __pow__(self,n):
        if n<0:return (1/self)**(-n)
        out=QR(self.d,1)
        for _ in range(n):out*=self
        return out
    def sign(self):
        sg=lambda x:(x>0)-(x<0)
        if self.b==0:return sg(self.a)
        if self.a==0 or sg(self.a)==sg(self.b):return sg(self.b)
        value=self.a*self.a-self.d*self.b*self.b
        return sg(value) if self.a>0 else -sg(value)
    def __str__(self):
        return str(self.a) if self.b==0 else f'({self.a})+({self.b})sqrt({self.d})'


checks=0
trace=[]


def gate(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(label)


sqrt2=QR(2,0,1)
one=QR(2,1)
zero=QR(2,0)
B=sqrt2-1
m=B*B
h=sqrt2/2
A=2-sqrt2
K0=(12*sqrt2-16)/3
cstar=(13-8*sqrt2)/3


# Formal bivariate polynomials over Q(sqrt2), in a and b.
def clean(p):return {ij:c for ij,c in p.items() if c!=zero}
def const(c):return clean({(0,0):zero+c})
def add(p,q):
    out=p.copy()
    for ij,c in q.items():out[ij]=out.get(ij,zero)+c
    return clean(out)
def neg(p):return {ij:-c for ij,c in p.items()}
def sub(p,q):return add(p,neg(q))
def mul(p,q):
    out={}
    for (i,j),c in p.items():
        for (k,l),d in q.items():
            ij=(i+k,j+l);out[ij]=out.get(ij,zero)+c*d
    return clean(out)
def scale(c,p):return mul(const(c),p)
def power(p,n):
    out=const(1)
    for _ in range(n):out=mul(out,p)
    return out
def derivative(p):
    return clean({(i,j-1):j*c for (i,j),c in p.items() if j})
def diag(p):
    out={}
    for (i,j),c in p.items():out[(i+j,0)]=out.get((i+j,0),zero)+c
    return clean(out)


a={(1,0):one};b={(0,1):one};u=const(1)
t=add(a,b)
C=sub(const(2*B),scale(m,t))
delta=sub(const(2),scale(sqrt2,t))
p3=add(power(a,3),mul(sub(u,power(a,2)),b))
p4=add(power(a,4),mul(sub(u,power(a,2)),power(b,2)))
Fenv=add(sub(scale(m,add(u,t)),p3),mul(C,p4))
target=add(sub(sub(const(B),scale(3*K0/8,delta)),p3),
           mul(add(const(A),scale(3*K0/8,delta)),p4))
gate(Fenv==target,'full target objective identity')
secant=sub(sub(u,scale(2,mul(C,b))),const(m))
secant_rhs=add(scale(2*m,mul(b,sub(a,b))),
               scale(4*m,mul(sub(const(h),b),sub(const(1+h),b))))
gate(secant==secant_rhs,'entire signed-tail secant factorization')
T=sub(add(scale(m,sub(sub(add(u,power(a,2)),scale(2,mul(a,b))),scale(3,power(b,2)))),
          scale(4*B,b)),u)
U=sub(add(scale(m,sub(const(2),scale(6,power(b,2)))),scale(4*B,b)),u)
gate(derivative(Fenv)==mul(sub(u,power(a,2)),T),'full derivative identity')
gate(sub(U,T)==scale(m,add(sub(sub(u,power(a,2)),power(b,2)),scale(2,mul(b,sub(a,b))))),
     'derivative domain remainder')
gate(U==scale(-6*m,mul(sub(b,const(h)),sub(b,const((4+sqrt2)/6)))),
     'derivative upper quadratic factorization')
gate(diag(Fenv)==scale(2*m,mul(sub(u,a),power(sub(const(h),a),2))),
     'diagonal boundary factorization')
gate(((4+sqrt2)/6-h).sign()>0 and K0.sign()>0,'strict factor signs')
gate(m*(1+sqrt2)**2==one and 3*K0/4==sqrt2*m,'two-atom final constant')
gate(K0==2*(6*sqrt2-8)/3,'exact doubling of prior constant')


def conv(a,b,cap=4):
    out=[a[0]*0]*(min(cap,len(a)+len(b)-2)+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<len(out):out[i+j]+=x*y
    return out
def literal(roots,cap=4):
    out=[roots[0]*0+1]
    for r in roots:out=conv(out,[r*0+1,r],cap)
    out += [roots[0]*0]*(cap+1-len(out))
    return out
def factor_power(r,n,cap=4):
    return [comb(n,j)*r**j for j in range(min(cap,n)+1)]
def nested_sign(x,y):
    sx,sy=x.sign(),y.sign()
    if sy==0:return sx
    if sx==0 or sx==sy:return sy
    z=(x*x-2*y*y).sign()
    return z if sx>0 else -z


# The same bounded signed prefix universe as the inherited source, with
# the stronger constant checked by literal coefficients and complete energy.
grid=[F(-2),F(-1),F(-1,2),F(0),F(1,2),F(1),F(2)]
rows=0;exceptions=[]
for x,y in product(grid,repeat=2):
    prefix=[F(1),x,y]
    S=sum(prefix)
    if S==0:
        exceptions.append((str(x),str(y),'zero prefix sum'));continue
    ep=(S*S-sum(r*r for r in prefix))/2
    raw=prefix+[-ep/S]
    r=[v/sum(raw) for v in raw]
    E=literal([v*v for v in r],2)[2]
    if E==0:
        exceptions.append((str(x),str(y),'zero energy'));continue
    coeff=literal(r)
    D=-conv(coeff,coeff)[4]
    aa,bb=(sorted([v for v in r if v>0],reverse=True)+[F(0),F(0)])[:2]
    dist2=2-sqrt2*(aa+bb)
    p3v,p4v=(sum(v**k for v in r) for k in [3,4])
    Cv=2*B-m*(aa+bb)
    envelope=aa**3-Cv*aa**4+(1-aa*aa)*(bb-Cv*bb*bb)
    gate(sum(r)==1 and sum(v*v for v in r)==1 and coeff[2]==0,'signed normalization and cancellation')
    gate(D==(5-8*p3v+3*p4v)/6 and E==(1-p4v)/2,'literal / Newton agreement')
    gate((envelope-(p3v-Cv*p4v)).sign()>=0,'signed tail envelope control')
    gate((D/E-cstar-K0*dist2).sign()>0,'stronger global stability control')
    gate(literal(r+[F(0)]*3)==coeff,'zero padding')
    trace.append(('signed',str(x),str(y),str(D/E),str(dist2)))
    rows+=1


families=[]
for N in range(2,9):
    for L in [2,8,100]:
        d=F(L+N-1,N*L)
        unit=QR(d,1)
        q=(N-1)/QR(d,L,L)
        S=N-L*q
        gate(q.sign()>0 and S.sign()>0,'uniform-atom family signs')
        gate(N*(N-1)-2*N*L*q+L*(L-1)*q*q==QR(d,0),'uniform-atom exact tuning')
        coeff=conv(factor_power(unit,N),factor_power(-q,L))
        squares=conv(factor_power(unit,N,2),factor_power(q*q,L,2),2)
        E=squares[2]/S**4
        D=-conv(coeff,coeff)[4]/S**4
        J=D/E
        tp=2/S
        gate(coeff[2]==QR(d,0) and E.sign()>0,'uniform-atom actual cancellation')
        # J-cstar-K0*(2-sqrt2*tp)
        # =J+19/3+8tp - (16/3+16tp/3)*sqrt2.
        gate(nested_sign(J+F(19,3)+8*tp,-F(16,3)*(1+tp))>0,
             'uniform-atom exact global estimate')
        if N+L<=16:
            gate(literal([unit]*N+[-q]*L)==coeff,'independent uniform-atom literal multiplication')
        trace.append(('uniform',N,L,str(J),str(tp)))
        families.append((N,L,str(J)))

# Explicit unequal-three-atom limiting shapes, normalized in square norm.
# They are legitimate closure data; the note gives the negative-dust lift.
atom_rows=[]
for raw in [[F(1),F(1),F(1)], [F(1),F(1),F(1,2)],
            [F(1),F(3,4),F(1,2)], [F(1),F(1,2),F(1,2)],
            [F(1),F(1),-F(1,2)], [F(1),F(1,2),-F(1,4)]]:
    norm=sum(x*x for x in raw)
    r=[QR(norm,0,x/norm) for x in raw]
    p3v,p4v=(sum(v**k for v in r) for k in [3,4])
    J=(5-8*p3v+3*p4v)/(3*(1-p4v))
    top=sorted([x for x in raw if x>0],reverse=True)+[F(0),F(0)]
    tp=QR(norm,0,(top[0]+top[1])/norm)
    gate(nested_sign(J+F(19,3)+8*tp,-F(16,3)*(1+tp))>0,
         'unequal signed macro-shape stronger estimate')
    atom_rows.append(tuple(str(x) for x in raw))
    trace.append(('macro',tuple(str(x) for x in raw),str(J),str(tp)))

# One-atom E->0 control keeps the divisor: J is exactly one.
for tval in [F(1,2),F(1,10),F(1,1000)]:
    raw=[F(1),tval,-tval/(1+tval)]
    r=[x/sum(raw) for x in raw]
    coeff=literal(r)
    E=literal([x*x for x in r],2)[2]
    D=-conv(coeff,coeff)[4]
    gate(E>0 and D==E and coeff[2]==0,'one-atom boundary exact ratio one')
    gate((1-cstar-K0*(2-sqrt2*(r[0]+r[1]))).sign()>0,'one-atom boundary stronger bound')
    trace.append(('one_atom',str(tval),str(E)))

# The relaxed equal-root diagonal has a limiting quotient exactly K0.
gate(4/(3*sqrt2*(1+h)**2)==K0,'sharp envelope obstruction limit')
for x in [F(1,2),F(3,5),F(7,10)]:
    J=F(5-3*x,3*(1+x))
    quotient=(J-cstar)/(2-2*sqrt2*x)
    gate(quotient==4/(3*sqrt2*(1+x)*(1+h)),'relaxed diagonal quotient identity')
    gate((quotient-K0).sign()>0,'relaxed diagonal strict finite gap')
    trace.append(('relaxed',str(x),str(quotient)))

print('PASS improved dimension-independent global duplication stability')
print('K0=(12sqrt2-16)/3; exactly twice the inherited bound; best global constant remains OPEN')
print('Proof certificates: energy-adjusted objective, signed secant, derivative/domain remainder, diagonal and two-atom boundaries')
print('Declared signed universe: 49 prefixes;',rows,'eligible;',len(exceptions),'exceptions')
for row in exceptions:print('EXCEPTION',*row)
print('Actual uniform-atom families:',len(families),'with N=2,...,8 and L=2,8,100')
print('Other controls:',len(atom_rows),'signed macro shapes; 3 one-atom boundary rows; 3 relaxed diagonal rows')
print('Relaxation obstruction: fractional number of equal roots between two and three; K0 sharp for this envelope only')
print('Exact gates:',checks)
print('Trace SHA256:',sha256(repr(trace).encode()).hexdigest())
