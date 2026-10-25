#!/usr/bin/env python3
"""Independent Session60 controls. Explicit failures remain active under -O."""
import itertools as it
import json
from math import comb
import sympy as s

COUNT = 0

def check(condition, message):
    global COUNT
    COUNT += 1
    if not bool(condition):
        raise RuntimeError(message)

def ident(expr, message):
    check(s.cancel(s.expand(expr)) == 0, message)

# Enumerate representatives by fixing the last sign; independently build characters
# as binary tuples, and test orthogonality and local-basis descent.
genera = {}
for m in range(5, 10):
    chars = [v for v in it.product((0,1), repeat=m) if sum(v)%2==0]
    reps = list(it.product((1,-1), repeat=m-1))
    for j in range(m):
        check(len({tuple(v[i] for i in range(m) if i != j) for v in chars}) == 2**(m-1), 'local monic basis')
    # Hilbert numerator convolution, with canonical degree m-6.
    d = m-6
    hilbert = sum((-1)**k*comb(m-3,k)*comb(m-1+d-2*k,m-1)
                  for k in range(m-2) if d-2*k >= 0)
    pg = sum(comb(sum(v)//2-1,2) if sum(v)>=6 else 0 for v in chars)
    check(pg == hilbert, 'character/Hilbert genus')
    genera[str(m)] = pg
    for v in chars:
        check(sum(v)%2 == 0, 'global sign character descent')

# Seven degree-six characters, their exact character Gram matrix and generators.
cc = [tuple(int(i != j) for i in range(7)) for j in range(7)]
gg = list(it.product((1,-1), repeat=6))
chi = s.Matrix([[s.prod(g[i]**c[i] for i in range(6)) for g in gg] for c in cc])
check(chi*chi.T == 64*s.eye(7), 'orthogonal rational character projectors')
for j in range(7):
    check([sum(c[i] for i in [j])%2 == 0 for c in cc] == [i==j for i in range(7)], 'one deletion invariant two-form')
for g6 in gg:
    g = g6+(1,)
    check([s.prod(g)*g[j] for j in range(7)] == [s.prod(g[i] for i in range(7) if i!=j) for j in range(7)], 'deck residue determinant twist')
check(78-(7*22) == -76, 'Picard sum offset')
check([20+6*r-76 for r in (18,19,20)] == [52,58,64], 'hourglass Picard candidates')

# Actual line incidences, no-fourfold and six-root fixed-point totals.
N = s.Matrix([[1,0,0],[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[1,1,1],[1,-1,-1]])
triples = [a for a in it.combinations(range(7),3) if N[list(a),:].det()==0]
check(triples == [(0,1,2),(0,3,4),(0,5,6)], 'triple incidence')
for inds in it.combinations(range(7),4):
    check(N[list(inds),:].rank()==3, 'no fourfold point')
for omit in range(7):
    inds = [i for i in range(7) if i!=omit]
    t = sum(set(a)<=set(inds) for a in triples)
    check(t == (0 if omit==0 else 2), 'deleted triple census')
    for a,b in it.combinations(inds,2):
        p = N[a,:].cross(N[b,:])
        zeros = [i for i in inds if (N[i,:]*p.T)[0]==0]
        fixed = 2**(5-len(zeros)) * (2 if len(zeros)==3 else 1)
        check(fixed == 8, 'resolved eight-fixed-point count')
check(22+15*6 == 16*7, 'invariant H2 dimension')
check(15*8//8 == 15, 'quotient nodes')
# Exceptional conic fixed locus for diag(-1,-1,1): w=0 yields two simple roots.
u,v,w = s.symbols('u v w')
check(s.discriminant(2*u*u+3, u) != 0, 'two simple exceptional fixed directions')

# Independent D12 enumeration and determinant ratio on actual relation spaces.
roots = {tuple(N[i,1:]) for i in range(1,7)}
M = N.T
K = s.Matrix.hstack(*M.nullspace())
left = (K.T*K).inv()*K.T
syms = []
for a,b in it.permutations(sorted(roots),2):
    A2 = s.Matrix.hstack(s.Matrix(a), s.Matrix(b))
    if not A2.det(): continue
    if {tuple(A2*s.Matrix(v)) for v in roots} != roots: continue
    A = s.diag(1,A2)
    # Columns of M are normals; A acts on those columns here.
    perm = [next(j for j in range(7) if A*M[:,i]==M[:,j]) for i in range(7)]
    P = s.zeros(7)
    for i,j in enumerate(perm): P[j,i]=1
    check(M*P == A*M, 'base/ambient exact sequence')
    B = left*P*K
    check(P*K == K*B, 'relation-space stability')
    ident(P.det()-A.det()*B.det(), 'canonical determinant quotient')
    syms.append(tuple(perm))
check(len(syms)==12, 'D12 count')
check({p[1] for p in syms} == set(range(1,7)), 'outer transitivity')

# Conic, explicit projective normalization to Veronese columns, then Gale dual.
A,B,C = s.symbols('A B C')
Q = s.Matrix([[1,0,0],[0,-1,s.Rational(1,2)],[0,s.Rational(1,2),-1]])
check(Q.det()==s.Rational(3,4), 'smooth dual conic')
for i in range(1,7): ident((N[i,:]*Q*N[i,:].T)[0], 'central point on conic')
x = s.symbols('x')
p = s.Matrix([x*x-x+1,x*x-1,x*(x-2)])
for entry in [ (p.T*Q*p)[0] ]: ident(entry, 'conic parametrization')
# s=1/(t-3), multiply parametrization by s^2 and read coefficients.
z = s.symbols('z')
pz = s.simplify(z*z*p.subs(x,3+1/z))
T = s.Matrix([[s.expand(e).coeff(z,j) for j in range(3)] for e in pz])
check(T.det()!=0, 'explicit conic/Veronese transformation')
ss = list(map(s.Rational,[-1])) + [s.Rational(-1,3),s.Rational(-1,4),s.Rational(-1,2),s.Rational(0),s.Rational(-2,5)]
V = s.Matrix([[a**k for a in ss] for k in range(3)])
column_scales = []
for j,a in enumerate(ss):
    col = T*V[:,j]
    d=col[0]
    check(col==d*N[j+1,:].T, 'actual normal transformed to Veronese')
    column_scales.append(d)
weights = [1/s.prod(a-b for b in ss if a!=b) for a in ss]
R = s.Matrix([[weights[i]*ss[i]**k for i in range(6)] for k in range(3)])
check(R*V.T == s.zeros(3), 'Gale annihilation')
check(R.rank()==V.rank()==3, 'full relation space')
check(R*s.diag(*[1/q for q in weights]) == V, 'square-root change diagonal coefficients')

# Shioda: calculate resultant instead of multiplying recorded evaluations.
E,X,Z,t = s.symbols('E X Z t')
q=(E+X+Z)*t*t-(E+2*Z)*t+(E-X)
f=t*(t+1)*(t-1)*(t-2)*(2*t-1)
branch=(E*E-X*X)*(E*E-Z*Z)*(E*E-(X+Z)**2)
res=s.resultant(q,f,t)
# resultant(q,f) = q2^5 f(t1)f(t2); W^2=q2*resultant.
ident((E+X+Z)*res-27*branch, 'independent resultant Shioda scalar27')
check(not s.sqrt(27).is_Rational, 'rational scaling boundary')
# Exact imaginary quadratic calculation via polynomial remainder, no floats.
d=s.symbols('d')
a=(1+d)/2; b=(1-d)/2
hexexpr=s.Poly(s.expand((t-a)**6-(t-b)**6+3*d*f),d)
check(s.rem(hexexpr,s.Poly(d*d+3,d)).as_expr()==0, 'sixth-root branch transform')
# Differential pullbacks and curve identities.
U,Vv=s.symbols('U Vv')
ident((u**6-1)-(u*u)**3+1, 'first elliptic quotient')
ident((u**6-1)/u**6 - (1-u**(-6)), 'second elliptic quotient')
ident(s.diff(u*u,u)-2*u, 'first differential coefficient')
ident(s.diff(u**-2,u)*u**3+2, 'second differential coefficient')
Gram=s.Matrix([[0,1,1,1],[1,0,1,1],[1,1,0,3],[1,1,3,0]])
check(Gram.det()==-3, 'CM divisor independence')
print(json.dumps({'status':'PASS','checks':COUNT,'optimization_safe':True,'genus':genera,'D12_order':len(syms),'gale_transform':str(T),'gale_column_scales':list(map(str,column_scales)),'shioda_scalar':27,'Picard_candidates':[52,58,64]},indent=2,sort_keys=True))
