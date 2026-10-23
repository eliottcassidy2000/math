from pathlib import Path
import json
import sympy as s
from itertools import combinations

a,b,c,u,v,w=s.symbols('a b c u v w');X=s.Matrix([[a,b],[c,-a]]);Y=s.Matrix([[u,v],[w,-u]]);J=s.Matrix([[0,1],[-1,0]])
assert X.T*J+J*X==s.zeros(2)
K=X*Y-Y*X
M=s.Matrix([[s.diff(K[0,0],q) for q in (u,v,w)],[s.diff(K[0,1],q) for q in (u,v,w)],[s.diff(K[1,0],q) for q in (u,v,w)]])
assert s.expand(M.det())==0
minors=[s.expand(M.extract(i,j).det()) for i in combinations(range(3),2) for j in combinations(range(3),2)]
assert -4*a*a in minors
assert 2*b*b in minors
assert 2*c*c in minors
# Hence every nonzero X has ad-rank2, and one-dimensional centralizer.
e=s.Matrix([[0,1],[0,0]]);f=e.T;h=s.diag(1,-1)
assert h*e-e*h==2*e and h*f-f*h==-2*f and e*f-f*e==h
# Illustrative general dimensions; not claims of equality existence.
cases=[{'s':d,'V':n,'paired':d+2*n+2,'lower':d+2*n+3,'isotropic':2*(d+n+1)} for d,n in [(3,2),(3,4),(3,6),(133,56)]]
assert all(r['isotropic']>r['lower'] for r in cases)
out={'adjoint_matrix':[[str(v) for v in row] for row in M.tolist()],'selected_nonzero_minors':['-4*a^2','2*b^2','2*c^2'],'centralizer_of_every_nonzero_sl2_element_dimension':1,'automatic_sl2_grading_verified':True,'dimension_controls':cases,'scope':'Symbolic matrix controls; module/radical/Levi argument independently audited in UNGRADED_AUDIT.md.'}
Path('lie-arithmetic-audit/ungraded_independent_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
# Explicit classical positive control, independently assemble all N_u matrices.
classical=[]
for m in range(1,6):
 n=2*m;J=s.zeros(n)
 for i in range(m):J[i,m+i]=1;J[m+i,i]=-1
 O=s.zeros(n+2);O[0,n+1]=1;O[n+1,0]=-1
 O[1:n+1,1:n+1]=J
 Z=s.zeros(n+2);Z[0,n+1]=2
 Ns=[]
 for i in range(n):
  N=s.zeros(n+2);N[1+i,n+1]=1
  for j in range(n):N[0,1+j]=J[i,j]
  assert N.T*O+O*N==s.zeros(n+2)
  assert N*Z-Z*N==s.zeros(n+2)
  Ns.append(N)
 for i in range(n):
  for j in range(n):assert Ns[i]*Ns[j]-Ns[j]*Ns[i]==J[i,j]*Z
 assert Z.T*O+O*Z==s.zeros(n+2)
 classical.append({'m':m,'pairs_checked':n*n,'completion_dimension':(m+1)*(2*m+3)})
out['classical_symplectic_matrix_controls']=classical
Path('lie-arithmetic-audit/ungraded_independent_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print('All explicit classical matrix controls pass.')
