"""Non-DRT controls by disjoint cyclic-triangle reversals; no fitted parameters."""
import json, sys, math
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'session30-decorated-engine'))
from predict import predict

def paley_deleted(q):
 r={i*i%q for i in range(1,q)}
 return np.array([[0 if i==j else 1 if (j-i)%q in r else -1 for j in range(q-1)] for i in range(q-1)],dtype=np.int64)

def reverse_disjoint_triangles(S):
 T=S.copy(); unused=set(range(len(S))); triangles=[]
 while len(unused)>=3:
  found=None
  for a in sorted(unused):
   for b in sorted(unused-{a}):
    if T[a,b]!=1:continue
    for c in sorted(unused-{a,b}):
     if T[b,c]==1 and T[c,a]==1:found=(a,b,c);break
    if found:break
   if found:break
  if not found:break
  a,b,c=found
  for i,j in ((a,b),(b,c),(c,a)):T[i,j]*=-1;T[j,i]*=-1
  triangles.append(found);unused.difference_update(found)
 return T,triangles

def scalar_A1(S,z=1.):
 n=len(S);Q=z*S/n;I=np.eye(n);Q2=Q@Q
 G=-Q2@np.linalg.inv(I+Q2);H=Q2@np.linalg.inv(I-Q2);O=Q@np.linalg.inv(I-Q2)
 D=np.diag(H-G);d=S.sum(axis=1);delta=d/n
 x=np.sum(2*np.diag(G)*np.diag(H)-1.5*np.diag(G)**2)+np.sum(G*(H*H-O*O-(2/3)*G*G))-D@G@D
 path=-z*z*delta@np.linalg.solve(I+Q2,delta)
 lol=-2*D@G@np.ones(n)
 return float(x+path+lol)

rows=[]
for q in (7,11,19,23,31,43,67,103,139):
 S=paley_deleted(q);T,tris=reverse_disjoint_triangles(S);n=len(T);d=T.sum(axis=1);E=T-S;K=T@T
 assert np.array_equal(d,S.sum(axis=1));assert np.all(np.abs(d)==1)
 assert np.all(E.sum(axis=1)==0);assert np.max(np.abs(E).sum(axis=1))<=4
 B=K+(n-1)*np.eye(n,dtype=np.int64)
 defect=K-(np.ones((n,n),dtype=np.int64)-(n+1)*np.eye(n,dtype=np.int64)+np.outer(d,d))
 a1=scalar_A1(T)
 row={'n':n,'triangles_reversed':tris,'arcs_flipped':3*len(tris),'max_row_l1_E':int(np.max(np.abs(E).sum(axis=1))),'max_offdiag_S2':int(np.max(np.abs(B))),'max_abs_Sd':int(np.max(np.abs(T@d))),'deleted_DRT_identity_defect_max':int(np.max(np.abs(defect))),'A1':a1,'n2_A1_plus_4_5_over_n':n*n*(a1+4.5/n),'A1_predicted_subleading':28/3,'n3_A1_remainder':n**3*(a1+4.5/n-(28/3)/n**2)}
 if n<=10:
  p=predict(T);assert abs(p['A1']-a1)<1e-10
  row['exact_kernel_predictor']=p
 rows.append(row)
 print(json.dumps(row),flush=True)
Path(__file__).with_name('finite-controls.json').write_text(json.dumps({'scope':'Finite numerical checks, not a proof of asymptotic bounds or exact Hamiltonian asymptotics.','controls':rows},indent=2)+'\n')
