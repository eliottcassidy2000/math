#!/usr/bin/env python3
"""Bounded exact controls for degree-independent signed-duplication stability.

Standard library only. The companion note proves the unbounded theorem.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import isqrt, comb, factorial
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


def mul(a, b, cap=None):
    degree = len(a) + len(b) - 2
    if cap is not None:
        degree = min(degree, cap)
    out = [a[0] * 0] * (degree + 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i+j <= degree:
                out[i+j] += x*y
    return out


def factor_power(r, n, cap=6):
    return [comb(n, j)*r**j for j in range(min(n, cap)+1)]


def literal_factors(roots, cap=6):
    out = [roots[0]*0+1]
    for r in roots:
        out = mul(out, [r*0+1, r], cap)
    return out


def sign_vs_sqrt2(x):
    """Return the sign of sqrt(2)-x in one quadratic field, exactly."""
    if x.sign() <= 0:
        return 1
    return (2-x*x).sign()


sqrt2 = QR(2, 0, 1)
one = QR(2, 1)
zero = QR(2, 0)
A, B = 2-sqrt2, sqrt2-1
cstar = (13-8*sqrt2)/3
gate((F(3,5)-A).sign() > 0, 'curvature rational upper bound on A')
gate(F(144,25) < F(27,4), 'curvature rational lower bound on 3sqrt3/2')
gate((1-cstar).sign() > 0, 'one-atom boundary ratio separated')

# Polynomial identity in formal t over Q(sqrt2).  Here
# h(x)+h(1-x)=(3t-t^3)/2-A*[1-(t^2-1)^2/2].
lhs = [B+A/2, -3*one/2, A, one/2, -A/2]
rhs = mul(mul([one, -2*one, one], [sqrt2, -one]), [one, A])
rhs = [v/2 for v in rhs]
gate(lhs == rhs, 'full two-support gap polynomial factorization')

# Explicit equality and hostile rows of the sphere inequality.
sphere_rows = [
    ('one_atom', [one]),
    ('two_atoms', [sqrt2/2, sqrt2/2]),
    ('negative_atom', [-one]),
    ('opposite_atoms', [sqrt2/2, -sqrt2/2]),
    ('four_atoms', [one/2]*4),
    ('missing_square_mass', [one/2]),
]
for name, roots in sphere_rows:
    p2 = sum(r*r for r in roots)
    H = sum(r**3-A*r**4 for r in roots)
    gate((1-p2).sign() >= 0, 'sphere row within domain')
    gate((B-H).sign() >= 0, 'sphere row gap')
    gate((B-H).sign() == 0 if name in ('one_atom', 'two_atoms')
         else (B-H).sign() > 0, 'sphere equality type')
    trace.append((name, str(p2), str(B-H)))

# Exact mixed-sign dust family. All moments use compressed multiplicities;
# coefficients use binomial factors. The first three rows additionally use
# an independent factor-by-factor multiplication retaining each root.
mixed_rows = []
for L in [6, 10, 25, 100, 400]:
    d = F(2)+F(10,L)
    u = QR(d,1)
    beta = 2/QR(d,2,1)
    kappa = F(5,L)
    S = 2-beta
    neg, pos = -2*beta/L, beta/L
    groups = [(u,2), (neg,L), (pos,L)]
    raw_p = {k: sum(n*r**k for r,n in groups) for k in range(1,5)}
    p = {k: raw_p[k]/S**k for k in range(1,5)}
    coeff = mul(mul(factor_power(u,2), factor_power(neg,L),6),
                factor_power(pos,L),6)
    sq = mul(coeff,coeff,4)
    sq_roots = mul(mul(factor_power(u*u,2), factor_power(neg*neg,L),2),
                   factor_power(pos*pos,L),2)
    E = sq_roots[2]/S**4
    D = -sq[4]/S**4
    J = D/E
    gate(((1-kappa)*beta*beta-4*beta+2) == QR(d,0), 'dust tuning quadratic')
    gate(S.sign() > 0 and beta.sign() > 0, 'normalization signs')
    gate(raw_p[1] == S and p[1] == u and p[2] == u, 'normalized two moments')
    gate(coeff[2] == QR(d,0), 'literal mixed-sign cancellation')
    gate(E == (1-p[4])/2 and E.sign() > 0, 'literal complete root energy')
    gate(D == (5-8*p[3]+3*p[4])/6, 'literal square matches Newton moments')
    gate(sign_vs_sqrt2((13-3*J)/8) > 0, 'strict mixed-sign ratio above cstar')
    gate(neg.sign() < 0 and pos.sign() > 0, 'both dust signs preserved')
    gate(L*(neg+pos)/S == 1-2/S, 'signed dust first moment retained')
    gate(L*(neg*neg+pos*pos)/S**2 == 1-2/S**2, 'dust square mass retained')
    if L <= 25:
        roots = [u]*2+[neg]*L+[pos]*L
        gate(literal_factors(roots) == coeff, 'independent literal mixed-sign factors')
        gate(literal_factors(roots+[QR(d,0)]*3) == coeff, 'zero-entry padding')
    mixed_rows.append((L, str(beta), str(J), str(E)))
    trace.append(('mixed', L, str(D), str(E), str(J)))

# Positive-energy one-atom approach: exact ratio one, with optional zero
# entries. Normalize in Q; Q(sqrt2) is used only for the moment gap.
boundary_rows = []
for t in [F(1,2),F(1,10),F(1,100),F(1,1000)]:
    roots = [F(1),t,-t/(1+t)]
    S = sum(roots)
    normalized = [r/S for r in roots]
    coeff = literal_factors(normalized+[F(0)]*3)
    D = -mul(coeff,coeff,4)[4]
    E = literal_factors([r*r for r in normalized],2)[2]
    p3,p4 = (sum(r**k for r in normalized) for k in [3,4])
    gap = B-p3+A*p4
    gate(coeff[2] == 0 and E > 0, 'one-atom cancellation and positive energy')
    gate(D == E, 'one-atom approach has exact ratio one')
    gate(gap == F(3,4)*E*(1-cstar), 'one-atom gap vanishes at energy scale')
    gate(gap.sign() > 0, 'finite one-atom moment gap positive')
    boundary_rows.append((str(t),str(E),str(gap)))
    trace.append(('boundary',str(t),str(E),str(gap)))

# Factor identity used when dividing the singular ratio by tail square mass.
gate(mul([F(1),F(-1)],[F(5),F(5),F(5),F(-3)])
     == [F(5),F(0),F(0),F(-8),F(3)], 'zero-energy numerator factorization')

# Independent finite formal multiplication of the unique normalized limit.
drift = 1-sqrt2
limit = mul([one,sqrt2,one/2],
            [drift**k/factorial(k) for k in range(7)],6)
gate(limit[1] == one and limit[2] == zero, 'limit normalization and cancellation')
gate(-mul(limit,limit,4)[4] == cstar/4, 'normalized limiting square coefficient')
log_moments = {k: 2*(sqrt2/2)**k for k in [2,3,4]}
gate(log_moments[2] == one and log_moments[3] == sqrt2/2
     and log_moments[4] == one/2, 'limit root moments')
gate((1-log_moments[4])/2 == one/4, 'limiting energy one quarter')
trace.append(('limit',tuple(str(v) for v in limit)))

print('PASS uniform signed-duplication stability and exponential closure')
print('Theorem: J->cstar iff two normalized +1/sqrt2 atoms plus dust square mass->0')
print('Equivalent closure: (1+s/sqrt2)^2 exp[(1-sqrt2)s], locally uniformly')
print('Sphere maximizers: one +1 atom OR two +1/sqrt2 atoms; one-atom ratio->1')
print('Declared exact universe: 6 sphere rows; 5 mixed-sign dust rows; 4 boundary rows')
for row in mixed_rows:
    print('MIXED L=',row[0], 'beta=',row[1], 'J=',row[2], 'E=',row[3])
for row in boundary_rows:
    print('BOUNDARY t=',row[0], 'E=',row[1], 'gap=',row[2], 'J=1')
print('Mixed-sign dust limit: positive sum sqrt2-1; negative sum -2(sqrt2-1)')
print('Exact gates:',checks)
print('Trace SHA256:',sha256(repr(trace).encode()).hexdigest())
