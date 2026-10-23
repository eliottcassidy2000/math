from pathlib import Path
import json
import sympy as s
z=s.symbols('z');u=(13-z)/12
P=lambda m:s.expand(13-12*s.chebyshevt(m,u))
D=lambda v:(v-1)*(v-25)
for m in range(1,13):
 assert s.expand(D(P(m))-D(z)*s.chebyshevu(m-1,u)**2)==0
for m in range(1,6):
 for n in range(1,6):assert s.expand(P(m).subs(z,P(n))-P(m*n))==0
A=s.Matrix([[2,1],[1,1]]);B=s.Matrix([[3,-1],[1,0]]);U=s.Matrix([[1,-1],[0,1]])
assert A*U==U*B
rows=[]
for m in range(1,21):
 L=s.zeros(m)
 for i in range(m-1):
  L[i,i]+=1;L[i+1,i+1]+=1;L[i,i+1]-=1;L[i+1,i]-=1
 forests=int((s.eye(m)+L).det());tr=int(s.trace(A**m));c=13-6*tr
 assert forests==s.fibonacci(2*m)
 assert tr==s.lucas(2*m)
 assert int(D(c))==180*forests**2
 assert P(m).subs(z,-5)==c
 rows.append({'vertices':m,'direct_forest_determinant':forests,'trace':tr,'c':c})
assert P(2).subs(z,13)==25
out={'symbolic_discriminant_identities':12,'symbolic_compositions':25,'direct_path_laplacian_determinants':rows,'scope':'Symbolic finite polynomial checks and direct determinants independent of producer rational-grid/Fibonacci controls; all-m proof separately audited.'}
Path('lie-arithmetic-audit/powering_independent_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
