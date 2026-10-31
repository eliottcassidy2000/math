import itertools,json
import sympy as S
from check import MS,D
P=S.Matrix([[1,0,0],[-1,1,1],[0,-1,0]]) # e0-e1,e1-e2,e1
Q=S.Matrix([[0,0,1],[1,1,0],[-1,0,0]]) # e1-e2,e1,e0
R=S.Matrix([[0,1],[1,0]]) # e1,e0
B18=S.kronecker_product(P,Q,R);B=S.diag(B18,B18);BI=B.inv()
expected=None
for c,M in enumerate(MS):
 T=BI*M*B
 if any(T[i,j]!=0 for i in range(D) for j in range(i)):raise RuntimeError(('not triangular',c))
 dg=[int(T[i,i]) for i in range(D)]
 if expected is None:expected=dg
 elif dg!=expected:raise RuntimeError('letter-dependent diagonal')
counts={str(k):expected.count(k) for k in set(expected)}
if counts!={'0':8,'1':16,'2':10,'3':2}:raise RuntimeError('multiplicity')
checks=0
for word in itertools.product(range(4),repeat=2):
 M=MS[word[1]]*MS[word[0]];T=BI*M*B
 if any(T[i,j] for i in range(D) for j in range(i)):raise RuntimeError('word lower triangle')
 if any(T[i,i]!=expected[i]**2 for i in range(D)):raise RuntimeError('word diagonal')
 checks+=1
print(json.dumps(dict(common_triangular_letters=4,diagonal_multiplicities=counts,two_letter_products=checks)))
from check import INIT,OUT,state
I=S.eye(D);M=MS[1]*MS[0];v=M*INIT
if (M-4*I)*(M-I)**2*v!=S.zeros(D,1):raise RuntimeError('01 recurrence')
for k in range(1,4):
 if (OUT*M**k*INIT)[0]!=S.Rational(4**k+6*k-1,9):raise RuntimeError('01 initial terms')
A=MS[0]
if A*(A-I)*(A-2*I)*(A-3*I)!=S.zeros(D):raise RuntimeError('0 minimal polynomial')
v=state((2,3,3))
for k in range(1,4):
 if (OUT*A**k*v)[0]!=14*3**k-17:raise RuntimeError('cancelled2mode initial terms')
print(json.dumps(dict(polynomial_mode_hostile='H((01)^k)=(4^k+6k-1)/9 for k>=1',cancelled2mode_hostile='H(2330^k)=14*3^k-17 for k>=1')))
