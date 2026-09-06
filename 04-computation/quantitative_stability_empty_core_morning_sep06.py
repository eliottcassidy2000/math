#!/usr/bin/env python3
"""Bounded exact controls for an explicit dimension-independent stability bound.

The analytic proof, not the finite controls, establishes all dimensions.
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


checks = 0
trace = []


def gate(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)


def mul(a, b, cap=4):
    out = [a[0]*0]*(min(cap, len(a)+len(b)-2)+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j < len(out):
                out[i+j] += x*y
    return out


def literal(roots, cap=4):
    out = [roots[0]*0+1]
    for r in roots:
        out = mul(out, [r*0+1,r], cap)
    out += [roots[0]*0]*(cap+1-len(out))
    return out


def factor_power(r, n, cap=4):
    return [comb(n,j)*r**j for j in range(min(n,cap)+1)]


def nested_sign(x, y):
    """Sign of x+y*sqrt(2) with x,y in the same exact quadratic field."""
    sx,sy = x.sign(),y.sign()
    if sy == 0:
        return sx
    if sx == 0 or sx == sy:
        return sy
    z = (x*x-2*y*y).sign()
    return z if sx > 0 else -z


sqrt2=QR(2,0,1)
one=QR(2,1)
A=2-sqrt2
B=sqrt2-1
m=3-2*sqrt2
K=(6*sqrt2-8)/3
cstar=(13-8*sqrt2)/3
gate(K.sign()>0 and m.sign()>0, 'positive explicit global constant')
gate(m*(1+sqrt2)**2 == one, 'reciprocal square constant')
gate(K == 4*m/(3*sqrt2), 'two-case constants coincide')
gate(1-A*sqrt2 == m, 'uniform secant lower bound')

# Entire two-atom identity, as formal polynomials in t over Q(sqrt2).
lhs=[B+A/2,-3*one/2,A,one/2,-A/2]
rhs=mul(mul([one,-2*one,one],[sqrt2,-one]),[one,A])
gate(lhs == [x/2 for x in rhs], 'exact two-atom gap factorization')

# Rational signed cancellation controls. For every one of the 49 prefixes
# (1,x,y), x,y in this fixed seven-point set, tune the last entry to cancel
# e2. Retain and declare exceptional zero-sum/zero-energy prefixes.
grid=[F(-2),F(-1),F(-1,2),F(0),F(1,2),F(1),F(2)]
rational_rows=0
skipped=[]
for x,y in product(grid,repeat=2):
    prefix=[F(1),x,y]
    S=sum(prefix)
    if S == 0:
        skipped.append((str(x),str(y),'zero prefix sum'))
        continue
    ep=(S*S-sum(r*r for r in prefix))/2
    raw=prefix+[-ep/S]
    total=sum(raw)
    gate(total != 0, 'nonzero normalized cancellation sum')
    roots=[r/total for r in raw]
    coeff=literal(roots)
    E=literal([r*r for r in roots],2)[2]
    if E == 0:
        skipped.append((str(x),str(y),'zero energy'))
        continue
    D=-mul(coeff,coeff)[4]
    positives=sorted([r for r in roots if r>0],reverse=True)+[F(0),F(0)]
    a,b=positives[:2]
    p3,p4=(sum(r**k for r in roots) for k in (3,4))
    H=p3-A*p4
    f=lambda u:u-A*u*u
    dist2=2-sqrt2*(a+b)
    gap=D/E-cstar-K*dist2
    gate(sum(roots)==1 and sum(r*r for r in roots)==1, 'rational normalized moments')
    gate(coeff[2]==0 and E>0, 'literal rational cancellation and energy')
    gate(E==(1-p4)/2 and D==(5-8*p3+3*p4)/6, 'rational Newton/literal agreement')
    gate((a*a*f(a)+(1-a*a)*f(b)-H).sign()>=0, 'signed tail envelope')
    gate(dist2.sign()>0 and gap.sign()>0, 'explicit global stability on signed control')
    gate(literal(roots+[F(0)]*3)==coeff, 'zero padding preserves actual coefficients')
    trace.append(('rational',str(x),str(y),str(D/E),str(dist2),str(gap)))
    rational_rows+=1


def compressed_control(name, groups, expected_a, expected_b):
    """Check actual coefficients, all root energy, and global bound."""
    d=groups[0][0].d
    u=QR(d,1)
    S=sum(n*r for r,n in groups)
    gate(S.sign()>0, 'positive family normalization')
    coeff=[u]
    squares=[u]
    for r,n in groups:
        coeff=mul(coeff,factor_power(r,n))
        squares=mul(squares,factor_power(r*r,n,2),2)
    p={k:sum(n*r**k for r,n in groups)/S**k for k in range(1,5)}
    E=squares[2]/S**4
    D=-mul(coeff,coeff)[4]/S**4
    J=D/E
    t=(expected_a+expected_b)/S
    gate(coeff[2]==QR(d,0) and p[1]==u and p[2]==u, 'family literal cancellation')
    gate(E.sign()>0 and E==(1-p[4])/2, 'family complete positive energy')
    gate(D==(5-8*p[3]+3*p[4])/6, 'family literal square / moment agreement')
    # J-cstar-K*(2-sqrt2*t) = J+1+4t-(4/3)(1+2t)sqrt2.
    gate(nested_sign(J+1+4*t,-F(4,3)*(1+2*t))>0,
         'global quadratic bound in exact biquadratic comparison')
    if sum(n for _,n in groups)<=33:
        expanded=[r for r,n in groups for _ in range(n)]
        gate(literal(expanded)==coeff, 'independent factor-by-factor family product')
    trace.append((name,str(S),str(J),str(t)))
    return str(J)


family_rows=[]
for L in [2,3,8,30,100]:
    d=F(L+1,2*L)
    u=QR(d,1)
    q=1/QR(d,L,L)
    gate(q.sign()>0 and (1-q).sign()>0, 'finite-degree extremizer root signs')
    J=compressed_control('two_equal_'+str(L),[(u,2),(-q,L)],u,u)
    family_rows.append(('two_equal',L,J))

    d=F(L+2,3*L)
    u=QR(d,1)
    q=2/QR(d,L,L)
    gate(6-6*L*q+L*(L-1)*q*q==QR(d,0), 'three-positive tuning quadratic')
    gate(q.sign()>0 and (1-q).sign()>0, 'three-positive family root signs')
    J=compressed_control('three_equal_'+str(L),[(u,3),(-q,L)],u,u)
    family_rows.append(('three_equal',L,J))

for v in [F(1,2),F(2,3),F(9,10)]:
    for L in [10,100]:
        d=1+v*v+2*v/L
        u=QR(d,1)
        q=2*v/QR(d,L*(1+v),L)
        gate(q.sign()>0, 'unequal two-atom negative dust')
        J=compressed_control('two_unequal_'+str(v)+'_'+str(L),
                             [(u,1),(QR(d,v),1),(-q,L)],u,QR(d,v))
        family_rows.append(('two_unequal_'+str(v),L,J))

# Four exact positive-energy controls approaching the one-atom boundary.
for t in [F(1,2),F(1,10),F(1,100),F(1,1000)]:
    raw=[F(1),t,-t/(1+t)]
    roots=[r/sum(raw) for r in raw]
    coeff=literal(roots)
    D=-mul(coeff,coeff)[4]
    E=literal([r*r for r in roots],2)[2]
    dist2=2-sqrt2*(roots[0]+roots[1])
    gate(coeff[2]==0 and E>0 and D==E, 'one-atom boundary exact ratio one')
    gate((1-cstar-K*dist2).sign()>0, 'one-atom boundary global estimate')
    trace.append(('one_atom',str(t),str(E),str(dist2)))

# Exact first-order constants and the limiting global obstruction.
C_dust=(28*sqrt2-32)/3
gate(C_dust == 4*(sqrt2-F(1,2))-2*cstar, 'dust derivative of the exact quotient')
K_two=(64-44*sqrt2)/3
gate(K_two == 8*(A*sqrt2+1)/(3*sqrt2*(sqrt2+1)**2), 'local unequal-atom constant')
gate((K_two-F(59,100)).sign()>0, 'local constant exceeds 0.59')

# K_three numerator x+y sqrt2, denominator u+v sqrt2, with coefficients
# in Q(sqrt3). Bound its exact value between 0.350 and 0.351.
x,y=QR(3,-F(4,3),-F(4,3)),QR(3,F(8,3))
u,v=QR(3,2),QR(3,0,-F(2,3))
gate(nested_sign(u,v)>0, 'three-atom distance positive')
for bound,expected in [(F(7,20),1),(F(351,1000),-1)]:
    gate(nested_sign(x-bound*u,y-bound*v)==expected, 'three-atom sharp-constant obstruction')
gate(F(351,1000)<F(59,100), 'three-atom obstruction defeats local two-atom guess')

print('PASS explicit dimension-independent quadratic signed-duplication stability')
print('Global K=(6sqrt2-8)/3; distance squared=2-sqrt2*(largest two positive roots sum)')
print('Declared universe: 49 rational prefixes; eligible rows=',rational_rows,'exceptions=',len(skipped))
for row in skipped:print('EXCEPTION',*row)
print('Families: 5 actual finite-degree minimizers; 5 three-positive controls; 6 unequal-two controls')
for row in family_rows:print('FAMILY',*row)
print('Boundary: 4 exact ratio-one rows; zeros retained; both root signs retained')
print('Power2 optimal; best global constant OPEN, between (6sqrt2-8)/3 and K_three in (0.350,0.351)')
print('Exact gates:',checks)
print('Trace SHA256:',sha256(repr(trace).encode()).hexdigest())
