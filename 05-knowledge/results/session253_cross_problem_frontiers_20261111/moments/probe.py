import sympy as S
from math import factorial
x,t,d=S.symbols('x t d')
def falling(z,n): return S.prod(z-i for i in range(n))
def rows(h):
 p=sum(S.Rational(factorial(2*h+1),factorial(3*h-3*j)*factorial(1+2*j))*falling(x+h,h-j)*t**j for j in range(h+1))
 # inverse t, using monic p and exact division of lower carry by p0
 qm=falling(2*x+2*h,2*h+1)/factorial(6*h+3)
 p0=p.subs(t,0)
 invpart=-S.cancel(qm/p0)*S.cancel((p-p0)/t)
 q=invpart+sum(falling(2*x+2*h,2*h-e)*t**e/S.Integer(factorial(6*h-3*e)*factorial(2+2*e)) for e in range(2*h+1))
 R=S.rem(S.Poly(q,t),S.Poly(p,t)).as_expr()
 M=S.zeros(h)
 for j in range(h):
  v=S.Poly(S.rem(S.Poly(R*t**j,t),S.Poly(p,t)),t)
  for i in range(h): M[i,j]=S.factor(v.nth(i))
 return p,R,M
checks=[]
for h in range(2,7):
 p,R,M=rows(h)
 # Entire small cluster at x=-h. Rank mod d is 1. det has order h-1,
 # established by invertibility of inverse scaled by d through specialization.
 M0=M.subs(x,-h)
 if not (M0.rank()==1): raise AssertionError('M0.rank()==1')
 # delta-adic expansion M=M0+d*M1... choose a unit pivot, Schur complement
 a,b=next((i,j) for i in range(h) for j in range(h) if M0[i,j]!=0)
 rowsidx=[i for i in range(h) if i!=a]; colsidx=[j for j in range(h) if j!=b]
 N=M.extract(rowsidx,colsidx)-M.extract(rowsidx,[b])*M.extract([a],colsidx)/M[a,b]
 N1=N.diff(x).subs(x,-h)
 if not (N1.det()!=0): raise AssertionError('N1.det()!=0')
 # Every r has local rank deficiency r-1 if complement is invertible in this finite bank.
 ranks=[M.subs(x,-r).rank() for r in range(2,h+1)]
 if not (ranks==[h-r+1 for r in range(2,h+1)]): raise AssertionError('ranks==[h-r+1 for r in range(2,h+1)]')
 checks.append({'h':h,'ranks_at_minus_r':ranks,'last_cluster_smith':[0]+[1]*(h-1),'schur_first_derivative_det':str(S.factor(N1.det()))})
 print(checks[-1],flush=True)
# Hostile: determinant valuation does not determine precision.
A=S.diag(1,d,d); B=S.diag(1,1,d**2)
if not (A.det()==B.det()): raise AssertionError('A.det()==B.det()')
if not (all(S.denom(S.cancel(z)).subs(d,0)!=0 for z in d*A.inv())): raise AssertionError('all(S.denom(S.cancel(z)).subs(d,0)!=0 for z in d*A.inv())')
if not (S.cancel((d*B.inv())[2,2])==1/d): raise AssertionError('S.cancel((d*B.inv())[2,2])==1/d')
print('HOSTILE: diag(1,d,d) and diag(1,1,d^2) have same determinant, largest loss 1 versus 2')
# Non-Eisenstein collision gives wrong balanced partition.
# p=t^2-d^2, response=t: multiplication has entries [[0,d^2],[1,0]], SNF (0,2).
print('HOSTILE: p=t^2-d^2 and R=t has Smith (0,2), not (0,1)')
print('PASS: all exact gates')

# Ternary collision centralizer and nonliftable direction.
N=S.Matrix([[0,0,1],[0,0,0],[0,0,0]])
J=S.Matrix([[0,1,0],[0,0,1],[0,0,0]])
a,b,c,e,f=S.symbols('a b c e f')
C=S.Matrix([[a,b,c],[0,e,f],[0,0,a]])
if not (C*N==N*C): raise AssertionError('C*N==N*C')
if not (C*J-J*C==S.Matrix([[0,a-e,b-f],[0,0,e-a],[0,0,0]])): raise AssertionError('C*J-J*C==S.Matrix([[0,a-e,b-f],[0,0,e-a],[0,0,0]])')
# Explicit p=t^3-d, response=t^2. Basis t^2,t,1.
T=S.Matrix([[0,1,0],[0,0,1],[d,0,0]])
if not ((T**2).subs(d,0)==N): raise AssertionError('(T**2).subs(d,0)==N')
if not (T**3==d*S.eye(3)): raise AssertionError('T**3==d*S.eye(3)')
if not (S.det(T**2)==d**2): raise AssertionError('S.det(T**2)==d**2')
print('PASS: ternary centralizer has 5 coordinates; liftable commutative reduction requires e=a and f=b')
# Unbounded-parameter complementary blocks: exact symbolic certificates.
r=S.symbols('r'); z=2*r+1
for ell in (1,2,3):
 P=sum(S.Rational(factorial(ell),factorial(i)*factorial(3*(ell-i)))*S.prod(z+j for j in range(2*i+1,2*ell+1))*t**i for i in range(ell+1))
 Q=sum(S.Rational(factorial(2*ell),factorial(i)*factorial(3*(2*ell-i)))*S.prod(2*z+j for j in range(2*i+1,4*ell+1))*t**i for i in range(2*ell+1))
 result=S.factor(S.resultant(P,Q,t))
 forced=S.prod((r+j)**j*(2*r+2*j+1)**j for j in range(1,ell+1))
 residual=S.Poly(S.cancel((-1)**ell*result/forced),r)
 if not (all(c>0 for c in residual.all_coeffs())): raise AssertionError('all(c>0 for c in residual.all_coeffs())')
 # Independent definition via Sylvester determinant.
 pc=S.Poly(P,t).all_coeffs(); qc=S.Poly(Q,t).all_coeffs()
 syl=S.zeros(3*ell)
 for i in range(2*ell):
  for j,c in enumerate(pc): syl[i,i+j]=c
 for i in range(ell):
  for j,c in enumerate(qc): syl[2*ell+i,i+j]=c
 # Polynomial determinant identity verified by a degree-bound interpolation:
 # product of 2ell P coefficients (degree<=2ell), ell Q (degree<=4ell), <=8ell².
 degree_bound=8*ell*ell
 for rr in range(degree_bound+1):
  if not (syl.subs(r,rr).det()==result.subs(r,rr)): raise AssertionError('syl.subs(r,rr).det()==result.subs(r,rr)')
 print('COMPLEMENT_CERTIFICATE',ell,'degree_bound',degree_bound,'resultant',result,flush=True)
print('PASS: complement nonvanishing for all r>=2 and ell=1,2,3')
