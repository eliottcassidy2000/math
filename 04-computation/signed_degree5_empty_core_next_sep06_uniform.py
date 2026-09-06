#!/usr/bin/env python3
"""Small exact controls for all-degree k=2 duplication constants.

All two-value configurations at n=4,...,12; exact quadratic arithmetic.
The unbounded conclusion is proved analytically in the companion note.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import isqrt, factorial
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

def mul(a,b):
    out=[a[0]*0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out

def elementary(roots):
    one=roots[0]*0+1
    out=[one]
    for r in roots:out=mul(out,[one,r])
    return out

def rho(t):return (5-3*t)/(3*(1+t))

orbits=0
table=[]
for n in range(4,13):
    d2=F(n-1,2*(n-2))
    tn=QR(d2,F(2,n),F(n-4,n))
    cn=rho(tn)
    # Exact squared certificate tn<1/sqrt(2), keeping both sides positive.
    bound=F(n*n,2)+4-(n-4)**2*d2
    gate(bound>0 and bound*bound>8*n*n,'finite uniform root-sum bound')
    gate(tn.sign()>0 and (1-tn).sign()>0,'maximal root sum in (0,1)')
    for m in range(1,n//2+1):
        d=F(n-1,m*(n-m))
        for branch in [-1,1]:
            x=QR(d,F(1,n),F(branch*(n-m),n))
            y=QR(d,F(1,n),-F(branch*m,n))
            roots=[x]*m+[y]*(n-m)
            a=elementary(roots);q=mul(a,a)[4]
            e=sum(u*u*v*v for u,v in combinations(roots,2))
            zero=QR(d,0);one=QR(d,1)
            gate(sum(roots)==one and sum(u*u for u in roots)==one,'two moments')
            gate(a[2]==zero,'two-value cancellation')
            t=x+y;xy=x*y
            gate(xy==(t-1)/n,'two-value product')
            p3,p4=sum(u**3 for u in roots),sum(u**4 for u in roots)
            gate(p3==((n-1)*t+1)/n,'third moment elimination')
            gate(p4==((n-1)*t*t+1)/n,'fourth moment elimination')
            if e==zero:
                gate(m==1 and branch==1 and q==zero,'only excluded zero-energy orbit')
            else:
                gate(-q/e==rho(t),'fractional-linear margin')
                if branch==1:
                    gate(m>=2 and (n-2*m)**2*d <= (n-4)**2*d2,'integer multiplicity maximizer')
                else:
                    gate((QR(d,F(2,n))-t).sign()>=0,'minus branch below center')
            trace.append((n,m,branch,str(t),str(e),str(q)))
            orbits+=1
    # Direct unnormalized extremal factors, separate from normalized orbits.
    length=n-2
    qn=1/QR(d2,length,length)
    eq_roots=[QR(d2,1)]*2+[-qn]*length
    a=elementary(eq_roots)
    e=sum(u*u*v*v for u,v in combinations(eq_roots,2))
    gate(a[2]==QR(d2,0),'literal extremal cancellation')
    gate(-mul(a,a)[4]/e==cn,'literal extremal exact constant')
    gate(qn.sign()>0 and (1-qn).sign()>0,'extremal sign partition')
    if n==4:gate(cn==QR(d2,F(7,9)),'quartic recovery')
    if n==5:gate(cn==QR(d2,F(27,29),-F(8,29)),'quintic recovery')
    table.append((n,str(tn),str(cn)))

# The limiting polynomial times exponential keeps the linear tail drift.
d=F(2);one=QR(d,1);zero=QR(d,0)
beta=QR(d,2,-1)
exp_coeff=[(-beta)**j/factorial(j) for j in range(5)]
limit=mul([one,2*one,one],exp_coeff)[:5]
cstar=QR(d,F(13,3),-F(8,3))
gate(limit[2]==zero,'linear-drift limit cancellation')
gate(-mul(limit,limit)[4]==cstar,'linear-drift limiting margin with energy one')
gate((cstar-F(1,3)).sign()>0,'strict improvement over terminal square')
gate(rho(QR(d,0,F(1,2)))==cstar,'fractional-linear limiting constant')

# n=3 is outside the Hessian lemma and is checked through the literal map.
r=[F(1),F(1),-F(1,2)]
a=elementary(r);e=sum(x*x*y*y for x,y in combinations(r,2))
gate(a[2]==0 and -mul(a,a)[4]==e,'degree-three exact ratio one')

print('PASS all-degree k=2 signed duplication constants')
print('t_n=[2+(n-4)sqrt((n-1)/(2(n-2)))]/n; c_n=(5-3t_n)/(3(1+t_n))')
print('Uniform strict bound: (13-8sqrt(2))/3; sharp infimum, no finite positive-energy equality')
print('Equality: two reciprocal roots 1; n-2 reciprocal roots -q_n; q_n=1/[n-2+sqrt((n-2)(n-1)/2)]')
print('Declared exact universe: n=4,...,12;',orbits,'oriented two-value configurations; 9 literal extremal rows')
for n,t,c in table:print('ROW',n,'t=',t,'c=',c)
print('Limit: (1+s)^2 exp[-(2-sqrt(2))s]; first coefficient at index2 zero; energy1')
print('Exact gates:',checks)
print('Trace SHA256:',sha256(repr(trace).encode()).hexdigest())
