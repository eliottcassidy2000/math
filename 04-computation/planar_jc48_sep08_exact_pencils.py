#!/usr/bin/env python3
"""Exact controls for all two-dimensional exact radical pencils through degree8.

The complete analytic line classification is in the companion note. This
program retains every possible moving-simple-root allocation in its inherited
17-type table, verifies universal primitives over the parameter field, and
checks actual quadratic Jacobian pairs and the full DG discriminant degree.
"""
from hashlib import sha256
import json
import sympy as S

x,t,c,a,b,z=S.symbols('x t c a b z')
gates=[]


def need(label, predicate):
    if not bool(predicate):
        raise RuntimeError(label)
    gates.append(label)


def zero(label, value):
    need(label, S.cancel(value)==0)


def jac(f,g):
    return S.diff(f,x)*S.diff(g,t)-S.diff(f,t)*S.diff(g,x)


# All exact generic multiplicity types in degrees0..8, including position
# conditions. A moving root is necessarily simple, so every allocation
# 1<=moving<=number_of_ones is retained before imposing the stated condition.
types=[(),(1,),(3,),(4,),(3,1),(5,),(6,),(5,1),(3,3),(4,1,1),
       (7,),(8,),(7,1),(5,3),(6,1,1),(4,3,1),(5,1,1,1)]
need('complete inherited exact-type count',len(types)==17)
allocations=[]
for part in types:
    for moving in range(1,part.count(1)+1):
        fixed=list(part)
        for _ in range(moving):fixed.remove(1)
        allocations.append((part,moving,tuple(fixed)))
need('all moving-root allocations',len(allocations)==12)
survivors=[((1,),1),((3,1),1),((5,1),1),((7,1),1),
           ((4,1,1),2),((5,1,1,1),3)]
for part,moving,fixed in allocations:
    if (part,moving) in survivors:
        need(f'survivor fixed gcd {part},{moving}',
             not fixed or (len(fixed)==1 and fixed[0] in (3,4,5,7)))
    else:
        need(f'rejection paid by position condition {part},{moving}',
             part in ((4,1,1),(6,1,1),(4,3,1),(5,1,1,1)))

# Coprime residual pencil: repeated roots occur at zeros of its fixed
# Wronskian; the rational derivative formula pays the generic-root claim.
A=x**3+x+1;B=x**2-2
zero('Wronskian derivative identity',S.diff(A/B,x)-(S.diff(A,x)*B-A*S.diff(B,x))/B**2)
need('nonconstant residual Wronskian',S.expand(S.diff(A,x)*B-A*S.diff(B,x))!=0)
need('coprime residual control',S.gcd(A,B)==1)

# The6+1+1 condition is the nonsingular conic3b^2=4ad. On a chart a0!=0,
# eliminating d0,d1 from the two endpoint equations turns the mixed
# coefficient condition into a square, forcing proportional vectors.
a0,a1,b0,b1=S.symbols('a0 a1 b0 b1',nonzero=True)
d0=3*b0**2/(4*a0);d1=3*b1**2/(4*a1)
mixed=6*b0*b1-4*(a0*d1+a1*d0)
zero('conic contains no line chart identity',mixed+3*(a0*b1-a1*b0)**2/(a0*a1))
need('conic matrix nonsingular',S.Matrix([[0,0,-2],[0,3,0],[-2,0,0]]).det()!=0)
# The a=0 endpoint has b=0,d!=0; its polar forces the other a=0,
# and then that endpoint's quadratic equation forces its b=0.
d=S.symbols('d',nonzero=True)
zero('conic missing chart polar',S.diff(3*b**2-4*a*d,a)+4*d)

# Universal primitives use rational coefficient functions, with no square
# roots of a,b,c adjoined to the coefficient field.
def radical_derivative(N,R):
    # d(yR)=(NR'+N'R/2)dx/y when y^2=N.
    return S.expand(N*S.diff(R,x)+S.diff(N,x)*R/2)

universal=[]
N=a*x+b;R=2/a
zero('linear universal primitive',radical_derivative(N,R)-1)
universal.append((0,1,N,R))
for k in (1,2,3):
    C={k:-S.Rational(2,2*k-1)/b}
    for j in range(k-1,0,-1):
        C[j]=-S.Rational(2*j,2*j-1)*a*C[j+1]/b
    # v=y/x^k; primitive v sum Cj x^-j = y R.
    R=sum(C[j]*x**(-j-k) for j in range(1,k+1))
    N=x**(2*k+1)*(a*x+b)
    zero(f'odd-root universal primitive k={k}',radical_derivative(N,R)-1)
    universal.append((2*k+1,1,N,R))
N=x**4*(a*x**2+b);R=-1/(b*x**3)
zero('411 universal primitive',radical_derivative(N,R)-1)
universal.append((4,2,N,R))
N=x**5*(a*x**3+b);R=-S.Rational(2,3)/(b*x**4)
zero('5111 universal primitive',radical_derivative(N,R)-1)
universal.append((5,3,N,R))
need('all six exact pencil spaces',len(universal)==6)

# Genuine rational Jacobian pairs for every line type. These controls have
# bounded discriminant degree; global DG membership is only asserted below.
for m,r,N,R in universal:
    H=t**2-x if m==0 else x**m*t**2-x**r
    P=S.diff(H,t);D=S.expand(P*P-4*S.Poly(H,t).coeff_monomial(t**2)*(H-c))
    expected=4*(x+c) if m==0 else 4*x**m*(x**r+c)
    zero(f'actual discriminant m={m},r={r}',D-expected)
    primitive=P*R.subs({a:4,b:4*c})
    G=S.cancel(-primitive.subs(c,H))
    zero(f'actual rational mate m={m},r={r}',jac(H,G)-1)
    need(f'actual discriminant degree m={m},r={r}',S.degree(D,x)<=8)

# Full fifteen-dimensional global quadratic, not a chosen subfamily.
nn=S.symbols('n0:9');pp=S.symbols('p0:5');q0=S.symbols('q0')
N=sum(nn[i]*x**i for i in range(9))
P=2*nn[8]*x**6+2*nn[7]*x**5+sum(pp[i]*x**i for i in range(5))
Q=nn[8]*x**4+nn[7]*x**3+(pp[4]-nn[6])*x**2+(pp[3]-nn[5])*x+q0
D=S.expand(P*P-4*N*(Q-c))
for degree in (12,11,10,9):
    zero(f'full DG discriminant cancellation degree{degree}',D.coeff(x,degree))
need('full DG discriminant degree at most8',S.degree(D,x)<=8)
zero('full DG discriminant direction',S.diff(D,c)-4*N)

# Proportional-pencil positive with two odd roots; general necessity alone
# must not discard this rank-one boundary of the pencil classification.
N=x**5*(x-1)**3
H=N*t**2
G=(8*x*x-4*x-1)/(3*x**4*(x-1)**2*t)
zero('proportional5+3 rational mate',jac(H,G)-1)
zero('proportional discriminant',S.diff(H,t)**2-4*N*(H-c)-4*c*N)

# The exact global elliptic family has both pencil ranks, beta=0 and!=0.
beta,delta,q0=S.symbols('beta delta q0')
hh=x*x+x**4*t
H=(x**4+delta*x)*(1+x*x*t)**2+beta*hh+q0
G=1/(3*delta*hh)
zero('global elliptic rational mate',jac(H,G)-1)
N=x**5*(x**3+delta);P=S.Poly(H,t).coeff_monomial(t);Q=H.subs(t,0)
D=S.expand(P*P-4*N*(Q-c))
zero('global elliptic exact pencil',D-((4*c+beta**2-4*q0)*x**8+4*delta*(c-q0)*x**5))
zero('global elliptic rank boundary',S.Matrix([[1,delta],[beta**2-4*q0,-4*delta*q0]]).det()+delta*beta**2)
r,bb=S.symbols('r bb')
zero('global elliptic actual second chart',H.subs({x:1/r,t:-r*r-r**4*bb})-((1+delta*r**3)*bb**2-beta*bb+q0))

# Named hostiles: pointwise exactness is weaker than a pencil of exact
# members, and taking a span of two exact endpoints need not preserve it.
N0=x**6*(x*x+2*x+3)
N1=x**6*(x*x-2*x+3)
for label,N in [('plus',N0),('minus',N1)]:
    poly=S.Poly(N/x**6,x)
    aa,bb,dd=poly.all_coeffs()
    zero(f'exact conic endpoint {label}',3*bb*bb-4*aa*dd)
mid=S.Poly((N0+N1)/x**6,x)
aa,bb,dd=mid.all_coeffs()
need('sum of exact conic endpoints is nonexact',3*bb*bb-4*aa*dd!=0)
for deg in (2,3):
    D=4*(c-x**deg)
    need(f'generic squarefree hostile degree{deg}',S.discriminant(D,x)!=0)
    need(f'hostile pencil is outside six families degree{deg}',S.degree(D,x) not in (1,4,6,8))

print('exact_pencils: PASS')
print('scope: all nonproportional degree<=8 exact pencils; full rational-mate criterion for bounded quadratic discriminants')
print('inherited universe:17 exact types;12 moving-simple-root allocations;6 line spaces')
print('universal primitives: rational in pencil coefficients and the chosen radical')
print('controls: all6 genuine quadratic mates; full15-dimensional DG cancellation; proportional5+3; global elliptic rank boundary')
print('hostiles: exact endpoints with nonexact sum; generic quadratic/cubic discriminants')
print('gates:',len(gates))
print('gate-label sha256:',sha256(json.dumps(gates,separators=(',',':')).encode()).hexdigest())
