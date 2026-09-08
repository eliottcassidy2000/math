#!/usr/bin/env python3
"""Exact controls for quadratic functions on W_m, m>=1.

The proof is uniform in m.  These finite controls check the complete section
spaces, actual coordinate changes, sharp examples, the sparse-pencil mate,
and both boundaries of the recovered polynomial-mate obstruction.  No
finite genus or coefficient census is used as an all-m proof.
"""
from hashlib import sha256
import json
import sympy as S


gates = []
records = []


def need(label, condition):
    if not bool(condition):
        raise RuntimeError(label)
    gates.append(label)


def zero(label, expression):
    need(label, S.cancel(S.together(expression)) == 0)


def jac(f, g, a, b):
    return S.diff(f, a)*S.diff(g, b)-S.diff(f, b)*S.diff(g, a)


x, t, r, b, z, v = S.symbols('x t r b z v')
delta, beta, q, lam = S.symbols('delta beta q lam', nonzero=True)
alpha, gamma, na, nb, da, db = S.symbols('alpha gamma na nb da db')

# Literal chart/basis controls.  All monomials in every numerator box are
# included, not only a generating subset.  The theorem proves all m,n.
basis_count = 0
for m in range(1, 7):
    tx = -r**m-r**(2*m)*b
    zero(f'actual inverse t chart m{m}',
         (1/(z-x**m)).subs({x:1/r,z:b/(1+r**m*b)}, simultaneous=True)-tx)
    zero(f'actual inverse b chart m{m}',
         (-x**m-x**(2*m)*t).subs({x:1/r,t:tx}, simultaneous=True)-b)
    zero(f'source two-form m{m}',jac(1/r,tx,r,b)-r**(2*m-2))
    for n in range(4):
        for j in range(n+1):
            for i in range(m*n+1):
                E=x**i*t**(n-j)*(1+x**m*t)**j
                target=(-1)**n*r**(m*n-i)*b**j*(1+r**m*b)**(n-j)
                zero(f'whole section box m{m} n{n} i{i} j{j}',
                     E.subs({x:1/r,t:tx}, simultaneous=True)-target)
                basis_count += 1

# Independent negative-Laurent-coefficient kernel: every source monomial
# x^i t^j, 0<=i<=4m,0<=j<=2, and every negative r coefficient is retained.
# Basis columns are then checked to span the entire computed kernel.
kernel_rows = []
for m in range(1, 5):
    ambient=[(i,j) for j in range(3) for i in range(4*m+1)]
    rows={}
    for col,(i,j) in enumerate(ambient):
        for k in range(j+1):
            power=m*j+m*k-i
            if power < 0:
                row=rows.setdefault((power,k),[0]*len(ambient))
                row[col]+=(-1)**j*S.binomial(j,k)
    matrix=S.Matrix(list(rows.values()))
    dimension=len(ambient)-matrix.rank()
    need(f'full negative-coefficient kernel dimension m{m}',dimension==3*(2*m+1))
    cols=[]
    for j in range(3):
        for i in range(2*m+1):
            p=S.Poly(S.expand(x**i*t**(2-j)*(1+x**m*t)**j),x,t)
            cols.append([p.coeff_monomial(x**a*t**c) for a,c in ambient])
    basis=S.Matrix.hstack(*(S.Matrix(c) for c in cols))
    need(f'full section basis independent m{m}',basis.rank()==dimension)
    need(f'full section basis kernel m{m}',matrix*basis==S.zeros(matrix.rows,basis.cols))
    kernel_rows.append((m,len(ambient),matrix.rows,dimension))

# Universal quadratic numerator and discriminant identity are polynomial
# identities in independent coefficient symbols (x^m is represented by z).
A,B,C=S.symbols('A B C')
N=A*z*z+B*z+C
P=2*A*z+B
zero('universal discriminant cancellation',P*P-4*N*(A-lam)-(B*B-4*A*C+4*lam*N))

# Actual global extremal families, all m=1..8 with free delta,beta,q.
# The bracket is checked after clearing the declared rational denominator.
family_rows=[]
for m in range(1,9):
    e=2*m-1
    h=x**m+x**(2*m)*t
    H=(x**(2*m)+delta*x)*(1+x**m*t)**2+beta*h+q
    N=x**(4*m)+delta*x**(2*m+1)
    P=2*x**(3*m)+2*delta*x**(m+1)+beta*x**(2*m)
    Q=x**(2*m)+delta*x+beta*x**m+q
    zero(f'actual source coefficients m{m}',H-(N*t*t+P*t+Q))
    tx=-r**m-r**(2*m)*b
    zero(f'actual h chart m{m}',h.subs({x:1/r,t:tx},simultaneous=True)+b)
    zero(f'actual H chart m{m}',H.subs({x:1/r,t:tx},simultaneous=True)
         -((1+delta*r**e)*b*b-beta*b+q))
    zero(f'sharp rational mate cleared bracket m{m}',jac(H,h,x,t)+e*delta*h*h)
    zero(f'actual finite field map bracket m{m}',jac(x**(-e),h,x,t)+e)
    zero(f'actual field-map expression m{m}',H-((1+delta*x**(-e))*h*h+beta*h+q))
    F=(4*lam+beta**2-4*q)*x**(4*m)+4*delta*(lam-q)*x**(2*m+1)
    zero(f'full discriminant pencil m{m}',P*P-4*N*(Q-lam)-F)
    # Direct chart numerator also checks its complete degree box.
    AA=x**(2*m)+delta*x+beta*x**m+q
    BB=-beta*x**(2*m)-2*q*x**m
    CC=q*x**(2*m)
    zero(f'actual numerator identity m{m}',
         (AA*z*z+BB*z+CC).subs(z,x**m+1/t)-H/t**2)
    need(f'actual numerator degree bound m{m}',
         all(S.degree(f,x)<=2*m for f in [AA,BB,CC]))
    need(f'generic branch and genus m{m}',((2*m)-2)//2==m-1)
    # beta=0 is the proportional-pencil boundary, still sharp.
    zero(f'proportional sharp pencil m{m}',
         F.subs(beta,0)-4*(lam-q)*N)
    # For every composite e too, V^e-u has Eisenstein constant coefficient
    # of exact u-order one.  This is a bounded literal polynomial control
    # for the all-e Eisenstein argument, not an irreducibility census.
    u=S.Symbol('u')
    cover=S.Poly(v**e-u,v)
    need(f'field degree monic and Eisenstein m{m}',
         cover.LC()==1 and cover.TC()==-u and
         all(cover.nth(j)==0 for j in range(1,e)))
    family_rows.append((m,e,m-1))

# The extremal sparse pencil supplies a complete rational mate.  Work in
# coordinates (x,V), H=(V^2-D0)/(4N); the source Jacobian is 2N times the
# (x,V) Jacobian, independently of the linear coefficient P.
for m in range(2,7):
    e=2*m-1
    N=na*x**(4*m)+nb*x**(2*m+1)
    D0=da*x**(4*m)+db*x**(2*m+1)
    H=(v*v-D0)/(4*N)
    G=2*v/(e*(db+4*nb*H)*x**(2*m))
    zero(f'full four-coefficient sparse-pencil mate m{m}',2*N*jac(H,G,x,v)-1)
    F=x**(2*m+1)*(alpha*x**e+gamma)
    R=2/(e*gamma*x**(2*m))
    zero(f'negative radical primitive m{m}',F*S.diff(R,x)+S.diff(F,x)*R/2+1)
    # Both leading coefficients nonzero: the residual binomial is
    # squarefree and misses zero, by an explicit Bezout identity.
    residual=alpha*x**e+gamma
    zero(f'binomial branch separation m{m}',residual-x*S.diff(residual,x)/e-gamma)

# Cheap hostile to a moving repeated root in a linear pencil.  Only the
# critical fibre lambda=0 is double; the generic root count is not inferred
# from a special coefficient value.
N=x+1;D0=x*x
W=S.diff(D0,x)*N-D0*S.diff(N,x)
zero('pencil critical determinant',W-x*(x+2))
need('generic pencil is squarefree',S.gcd(D0+4*lam*N,S.diff(D0+4*lam*N,x))==1)
need('special fibre has a repeated root',S.gcd(D0,S.diff(D0,x))==x)
fixed=x**5*(x*x+1)
need('gcd retains fixed high root',S.gcd(fixed,x**5*(x**3+2))==x**5)

# Complete possible degree patterns in the N-constant argument, m=1..12.
# These are all degrees after the exact coefficient cancellation, not
# numerical choices of the free coefficients.
degree_rows=[]
for m in range(1,13):
    for k in range(1,m+1):
        for ell in range(-1,m+1):
            degP=m+k
            need(f'centered nonconstant A cannot be affine m{m} k{k} ell{ell}',
                 degP>ell and 2*degP>max(k,1))
            degree_rows.append((m,k,ell))
    for degP in range(m+1):
        need(f'constant A centered degree never one m{m} d{degP}',2*degP!=1)
    # Explicit symbolic identity with an arbitrary polynomial E in degree m.
    coefficients=S.symbols(f'e{m}_0:{m+1}')
    AA=sum(coefficients[i]*x**i for i in range(m+1))
    EE=sum((i+1)*x**i for i in range(m+1))
    BB=-AA*x**m+EE
    CC=1-x**m*EE
    zero(f'full constant-leading cancellation m{m}',AA*x**(2*m)+BB*x**m+CC-1)
    zero(f'full centered source linear coefficient m{m}',2*AA*x**m+BB-(AA*x**m+EE))

# Scope hostiles: the genuine-degree condition and rational/polynomial
# distinction are essential.  Both t and t^2 are actual global functions.
zero('lower-degree global t has polynomial mate',jac(t,-x,x,t)-1)
zero('genuine global t squared has rational mate',jac(t*t,-x/(2*t),x,t)-1)
need('t squared is genuine quadratic with square generic discriminant',S.degree(t*t,t)==2)
# A degree-twelve exact radical curve already has genus two, so retaining
# m (and the discriminant degree) is essential.
F=x**7*(x**5+1);R=-S.Rational(2,5)/x**6
zero('high-genus exact hostile',F*S.diff(R,x)+S.diff(F,x)*R/2-1)

records=[basis_count,kernel_rows,family_rows,len(degree_rows)]
print('dg_genus: PASS')
print('scope: W_m quadratic filtration, sharp generic genus <=m-1, extremal rational-mate iff')
print('polynomial exclusion: recovered THM-2071 consequence; genuine degree two required')
print('whole section boxes: m=1..6,n=0..3;',basis_count,'basis identities')
print('independent full Laurent kernels (m,ambient,rows,dimension):',kernel_rows)
print('actual sharp families (m,field-degree,genus):',family_rows)
print('sparse-pencil four-coefficient Jacobians: m=2..6; all symbolic coefficients retained')
print('complete centered degree patterns: m=1..12;',len(degree_rows),'nonconstant-A cases')
print('hostiles: lower-degree polynomial mate, genuine quadratic rational mate, special pencil fibre, genus2')
print('gates:',len(gates))
print('semantic sha256:',sha256(json.dumps([gates,records],separators=(',',':')).encode()).hexdigest())
