#!/usr/bin/env python3
"""Exact algebra controls for the sharp regional three-root packing theorem."""
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


@dataclass(frozen=True)
class K4:
    # Basis 1,sqrt2,sqrt3,sqrt6.
    c: tuple
    def __post_init__(self):object.__setattr__(self,'c',tuple(map(F,self.c)))
    def cast(self,x):return x if isinstance(x,K4) else K4((x,0,0,0))
    def __add__(self,x):
        x=self.cast(x);return K4(tuple(a+b for a,b in zip(self.c,x.c)))
    __radd__=__add__
    def __neg__(self):return K4(tuple(-x for x in self.c))
    def __sub__(self,x):return self+-self.cast(x)
    def __rsub__(self,x):return self.cast(x)+-self
    def __mul__(self,x):
        x=self.cast(x);out=[F(0)]*4
        for i,a in enumerate(self.c):
            for j,b in enumerate(x.c):
                repeat=i&j
                out[i^j]+=a*b*(2 if repeat&1 else 1)*(3 if repeat&2 else 1)
        return K4(tuple(out))
    __rmul__=__mul__
    def __truediv__(self,x):return K4(tuple(a/F(x) for a in self.c))
    def __pow__(self,n):
        out=self.cast(1)
        for _ in range(n):out*=self
        return out
    def sign(self):
        a,b,c,d=self.c
        x,y=QR(3,a,c),QR(3,b,d)
        sx,sy=x.sign(),y.sign()
        if sy==0:return sx
        if sx==0 or sx==sy:return sy
        sg=(x*x-2*y*y).sign()
        return sg if sx>0 else -sg
    def __str__(self):return str(self.c)


checks=0
trace=[]
def gate(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(label)

zero=K4((0,0,0,0));one=K4((1,0,0,0))
u=K4((0,1,0,0));v=K4((0,0,1,0))
h=u/2;z=v/3;A=2-u;B=u-1
gamma=A*(3-v)/4
C0=(1+u+v-u*v)/2
C3=(3-v)/2
K3=4*u*gamma/3
cstar=(13-8*u)/3

# Formal bivariate polynomials over Q(sqrt2,sqrt3).
def clean(p):return {ij:c for ij,c in p.items() if c!=zero}
def cn(c):return clean({(0,0):zero+c})
def add(p,q):
    out=p.copy()
    for ij,c in q.items():out[ij]=out.get(ij,zero)+c
    return clean(out)
def neg(p):return {ij:-c for ij,c in p.items()}
def sub(p,q):return add(p,neg(q))
def mul(p,q):
    out={}
    for (i,j),a in p.items():
        for (k,l),b in q.items():
            ij=(i+k,j+l);out[ij]=out.get(ij,zero)+a*b
    return clean(out)
def sc(c,p):return mul(cn(c),p)
def pw(p,n):
    out=cn(1)
    for _ in range(n):out=mul(out,p)
    return out
def der(p):return clean({(i-1,j):i*c for (i,j),c in p.items() if i})

t={(1,0):one};c={(0,1):one};unit=cn(1)
C=sub(cn(C0),sc(gamma,t))
p3=add(sc(F(1,2),sub(sub(sc(3,t),pw(t,3)),sc(3,mul(t,pw(c,2))))),pw(c,3))
p4=sc(F(1,2),add(add(sub(sub(sub(sc(3,pw(c,4)),sc(2,mul(pw(c,2),pw(t,2)))),sc(2,pw(c,2))),pw(t,4)),sc(2,pw(t,2))),unit))
obj=add(sub(sub(unit,C),p3),mul(C,p4))
ftt=add(sub(sc(3,t),sc(2*C0,sub(add(sc(3,pw(t,2)),pw(c,2)),unit))),
        sc(2*gamma,mul(t,sub(add(sc(5,pw(t,2)),sc(3,pw(c,2))),cn(3)))))
gate(der(der(obj))==ftt,'full three-root second derivative')
gate(C0==A+u*gamma and C0-2*gamma*z==C3,'energy coefficient and regional endpoint')
gate((3*v/8-C3).sign()>0 and (7*v/4-3).sign()>0,'strict tail convexity constants')
gate((F(3,16)-gamma).sign()>0 and gamma.sign()>0,'concavity gamma bound')
gate((C0-F(27,32)).sign()>0,'concavity C0 bound')
gate((2*z-F(8,7)).sign()>0 and (F(10,7)-u).sign()>0,'rational t interval')
gate(F(15)*F(8,7)**2-27*F(8,7)+5==F(-307,49),'left convex quadratic endpoint')
gate(F(15)*F(10,7)**2-27*F(10,7)+5==F(-145,49),'right convex quadratic endpoint')
gate(140**2<81**2*3,'strict Q at initial t')
gate(F(35,187)<F(3,16) and F(61,72)>F(27,32),'rational parameter refinements')
for n,d,val,orientation in [(24,17,2,-1),(17,12,2,1),(19,11,3,-1),(26,15,3,1)]:
    gate(((n*n>val*d*d)-(n*n<val*d*d))==orientation,'radical rational bracket')

def envelope(ap,bp):
    cp=sub(cn(C0),sc(gamma,add(ap,bp)))
    ep3=add(pw(ap,3),mul(sub(unit,pw(ap,2)),bp))
    ep4=add(pw(ap,4),mul(sub(unit,pw(ap,2)),pw(bp,2)))
    return add(sub(sub(unit,cp),ep3),mul(cp,ep4))

# b=z boundary: full polynomial identity, then exact positive endpoint.
polyP=add(add(sub(pw(t,3),sc(1+u,pw(t,2))),sc(F(2,3),t)),cn(F(2,3)))
edge=mul(sc(gamma,mul(sub(unit,t),sub(t,cn(z)))),polyP)
gate(envelope(t,cn(z))==edge,'complete b=z boundary factorization')
aa=u*v/3
Pend=aa**3-(1+u)*aa**2+2*aa/3+F(2,3)
gate(Pend==2*u*(2*z-1)/3 and Pend.sign()>0,'edge cubic positive endpoint')
gate((F(8,3)-2*(1+u)*z).sign()<0,'edge cubic monotonicity')

# a=b boundary, with the square relation c^2=1-2x^2 retained.
diag=sc(2*gamma,mul(mul(sub(unit,t),sub(t,cn(z))),sub(t,cn(h))))
gate(envelope(t,t)==diag,'complete equal-top envelope factorization')
def reduce_square(p):
    out={}
    relation=sub(unit,sc(2,pw(t,2)))
    for (i,j),cc in p.items():
        term=mul({(i,j%2):cc},pw(relation,j//2))
        out=add(out,term)
    return clean(out)
Cdiag=sub(cn(C0),sc(2*gamma,t))
actual_diag=add(sub(sub(unit,Cdiag),add(sc(2,pw(t,3)),pw(c,3))),
                mul(Cdiag,add(sc(2,pw(t,4)),pw(c,4))))
correction=mul(mul(pw(c,2),sub(t,c)),sub(unit,mul(Cdiag,add(t,c))))
gate(reduce_square(sub(actual_diag,add(diag,correction)))=={},'exact diagonal tail correction')
den=mul(add(t,c),add(t,cn(h)))
bracket_num=sub(sc(3,mul(mul(add(t,cn(z)),add(t,cn(h))),sub(unit,mul(Cdiag,add(t,c))))),
                sc(gamma,mul(sub(unit,t),add(t,c))))
right=mul(mul(pw(c,2),sub(t,cn(z))),bracket_num)
gate(reduce_square(sub(mul(actual_diag,den),right))=={},'full cleared positive diagonal factorization')
gate((2-v-F(1,4)).sign()>0 and (F(1,4)-gamma).sign()>0,'positive diagonal brace constants')
gate((2*z+h-1).sign()>0,'diagonal quotient at most one')

def conv(a,b,cap=4):
    out=[a[0]*0]*(min(cap,len(a)+len(b)-2)+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<len(out):out[i+j]+=x*y
    return out
def literal(r,cap=4):
    out=[r[0]*0+1]
    for x in r:out=conv(out,[x*0+1,x],cap)
    out += [r[0]*0]*(cap+1-len(out))
    return out
def factor_power(r,n,cap=4):return [comb(n,j)*r**j for j in range(min(n,cap)+1)]

# Bounded actual signed cancellation rows, with the regional filter explicit.
eligible=0;rejected=[]
for x,y in product([F(1,2),F(3,4),F(1)],
                   [F(-1),F(-3,4),F(-1,2),F(-1,4),F(0),F(1,4),F(1,2),F(3,4)]):
    pre=[F(1),x,y];S=sum(pre)
    ep=(S*S-sum(a*a for a in pre))/2
    raw=pre+[-ep/S];rawS=sum(raw)
    r=[a/rawS for a in raw]
    co=literal(r);E=literal([a*a for a in r],2)[2]
    a,b=(sorted([a for a in r if a>0],reverse=True)+[F(0),F(0)])[:2]
    if E==0 or b*b<F(1,3):
        rejected.append((str(x),str(y),'outside region or zero energy'));continue
    D=-conv(co,co)[4]
    gate(sum(r)==1 and sum(a*a for a in r)==1 and co[2]==0,'actual regional normalization')
    gate((D/E-cstar-K3*(2-u*(a+b))).sign()>0,'strict regional target via literal coefficients')
    gate(sum(a for a in r if a<0)<0,'negative tail retained')
    gate(literal(r+[F(0)]*3)==co,'actual zero padding')
    trace.append(('actual',str(x),str(y),str(D/E),str(a+b)))
    eligible+=1

# Explicit finite family lies inside the region, unlike equal-three dust.
family=[]
for n in [3,4,8,20]:
    L=n*n
    Q=F(3)+F(6,n*n)
    dd=Q+(9-Q)/L
    un=QR(dd,1)
    q=(9-Q)/QR(dd,3*L,L)
    S=3-L*q
    x,y=F(n+1,n),F(n-2,n)
    coeff=conv(conv(factor_power(QR(dd,x),2),factor_power(QR(dd,y),1)),
               factor_power(-q,L))
    E=conv(conv(factor_power(QR(dd,x*x),2,2),factor_power(QR(dd,y*y),1,2),2),
           factor_power(q*q,L,2),2)[2]/S**4
    gate(9-Q-6*L*q+L*(L-1)*q*q==QR(dd,0),'sharpness exact tuning')
    gate(q.sign()>0 and (F(2,L)-q).sign()>0 and S.sign()>0,'sharpness sign and dust bound')
    gate(S*S==Q+L*q*q and (3*x*x-S*S).sign()>0,'finite sharpness family actually inside region')
    gate(coeff[2]==QR(dd,0) and E.sign()>0,'sharpness literal cancellation and energy')
    if n<=4:
        gate(literal([QR(dd,x)]*2+[QR(dd,y)]+[-q]*L)==coeff,'independent literal sharpness factors')
    family.append((n,L,str(S),str(E)))
    trace.append(('sharpness',n,L,str(S),str(E)))

# The two surrogate equality types, and a strict rational three-root point.
for name,r,sg in [('three_equal',[z,z,z],0),('two_equal',[h,h,zero],0),
                  ('strict',[one*F(2,3),one*F(2,3),one*F(1,3)],1)]:
    p3v,p4v=(sum(a**k for a in r) for k in [3,4])
    cc=C0-gamma*(r[0]+r[1])
    val=1-p3v-cc*(1-p4v)
    gate(sum(a*a for a in r)==one and val.sign()==sg,'surrogate equality and strictness control')
    trace.append(('surrogate',name,str(val)))
gate(K3*(2-2*u*z)==3-4*v/3-cstar,'exact sharp limiting ratio')

print('PASS sharp regional K3 duplication stability')
print('Region: largest two positive normalized roots >=1/sqrt3; all finite rows have strict inequality')
print('Tail packing preserves p2 and the selected top roots, not p1; used only as an objective comparison')
print('Declared rational prefix universe24: eligible',eligible,'rejected',len(rejected))
for row in rejected:print('REJECTED',*row)
print('Actual sharpness family controls:',family)
print('Proof controls: full three-root curvature and both boundary factorizations; exact radical sign certificates')
print('Exact gates:',checks)
print('Trace SHA256:',sha256(repr(trace).encode()).hexdigest())
