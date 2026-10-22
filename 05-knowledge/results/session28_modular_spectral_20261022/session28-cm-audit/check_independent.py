"""Independent consumer audit of CM certificate and harmonic/Hessian algebra.
Does not execute producer scripts. Requires sympy.
"""
import json, math, hashlib
from pathlib import Path
import sympy as s
from sympy.matrices.normalforms import hermite_normal_form
HERE=Path(__file__).resolve().parent
SRC=HERE.parent/'session28-rank88'/'check_cm163_glue.json'
data=json.loads(SRC.read_text()); p=163
Y=s.Matrix(data['overlattice_basis_numerator']); B=Y/p
G=s.diag(*[s.Matrix([[2,1],[1,82]]) for _ in range(4)])
H=s.Matrix(data['overlattice_gram']); M=s.Matrix(data['alpha_integral_matrix'])
assert H==B.T*G*B and H.det()==1
assert all(t.q==1 for t in H) and all(H[i,i]%2==0 for i in range(8))
assert all(H[:i,:i].det()>0 for i in range(1,9))
assert all(t.q==1 for t in B.inv()) and B.det()==s.Rational(1,p*p)
for v in data['isotropic_generators']:
 g=s.Matrix([t for a in v for t in [82*a,-a]])/p
 assert all(t.q==1 for t in B.inv()*g)
assert M*M-M+41*s.eye(8)==s.zeros(8)
assert M.T*H+H*M==H and M.T*H*M==41*H
J=2*M-s.eye(8)
assert J*J==-p*s.eye(8) and J.T*H==-H*J
# Direct norm-2 enumeration, exact LDL branch bounds.
L,D=H.LDLdecomposition(hermitian=False); x=[0]*8; roots=[]; visited=[0]
def rec(i,used):
 visited[0]+=1
 if i<0:
  if used==2:roots.append(tuple(x))
  return
 c=sum(L[j,i]*x[j] for j in range(i+1,8)); R=(2-used)/D[i,i]
 if R<0:return
 bound=math.isqrt(int(s.ceiling(R)))+1
 lo=int(s.floor(-c))-bound; hi=int(s.ceiling(-c))+bound
 for k in range(lo,hi+1):
  value=used+D[i,i]*(k+c)**2
  if value<=2:
   x[i]=k;rec(i-1,value)
rec(7,s.Integer(0))
assert len(roots)==240 and s.Matrix(roots).rank()==8
# Additional composite discriminants, outside producer's prime-only census.
composites=[]
for d in [15,27,35,39,51,55,63,75,87,91,99]:
 a,b=next((a,b) for a in range(d) for b in range(d) if (1+a*a+b*b)%d==0)
 m=(d+1)//4; A=s.Matrix([[2,1],[1,2*m]]); g0=s.diag(A,A,A,A)
 V=[(1,0,a,b),(0,1,-b,a)]
 gs=[s.Matrix([t for z in v for t in [2*m*z,-z]]) for v in V]
 B0=hermite_normal_form(s.Matrix.hstack(d*s.eye(8),*gs))/d
 h=B0.T*g0*B0; alpha=B0.inv()*s.diag(*[s.Matrix([[0,-m],[1,1]]) for _ in range(4)])*B0
 assert h.det()==1 and all(t.q==1 for t in h) and all(h[i,i]%2==0 for i in range(8))
 assert all(t.q==1 for t in alpha) and alpha.T*h+h*alpha==h
 assert alpha**2-alpha+m*s.eye(8)==s.zeros(8)
 composites.append({'d':d,'a':a,'b':b,'determinant':1,'norm_alpha':m})
# Direct polynomial Laplacian audit for non-diagonal trace-free symmetric H.
harmonic=[]
for n in [2,3,4,5,6]:
 xs=s.Matrix(s.symbols('x:'+str(n))); r2=(xs.T*xs)[0]
 A=s.Matrix(n,n,lambda i,j: ((i+1)*(j+1)+i+j)%5-2); A=(A+A.T)/2
 A[n-1,n-1]-=s.trace(A)
 tr=s.trace(A*A); h=(xs.T*A*xs)[0]; h2=(xs.T*A*A*xs)[0]
 P=s.expand(h*h-s.Rational(4,n+4)*r2*(h2-tr*r2/n)-s.Rational(2,n*(n+2))*r2*r2*tr)
 assert s.expand(sum(s.diff(P,z,2) for z in xs))==0
 # Compare expanded standard harmonic formula.
 Palt=s.expand(h*h-s.Rational(4,n+4)*r2*h2+s.Rational(2,(n+4)*(n+2))*r2*r2*tr)
 assert P==Palt
 harmonic.append(n)
assert s.Rational(168498000*2*8*8,88*90)==2723200
out={'status':'PASS','producer_input_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),'cm163':{'certificate_identities':'direct matrix arithmetic, no producer execution','norm2_roots':len(roots),'root_span_rank':8,'enumeration_nodes':visited[0],'positive_principal_minors':True},'composite_discriminant_controls':composites,'harmonic_laplacian_dimensions':harmonic,'minimal_shell_isotropic_constant':2723200,'limits':['Does not enumerate rank88 short vectors or certify its minimum.','Harmonic theta modularity and rank8 classification are classical cited inputs.','No claim that all rank88 Gaussian Hessians are positive.']}
text=json.dumps(out,indent=2)+'\n';print(text);(HERE/'check_independent.json').write_text(text)
