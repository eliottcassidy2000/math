from pathlib import Path
import sympy as s,json
p=str(Path(__file__).with_name('transitive_hessian_probe.json'));rows=json.load(open(p))
def inertia(A):
 A=s.Matrix(A);pos=neg=zero=0;steps=[]
 while A.rows:
  n=A.rows;k=next((i for i in range(n) if A[i,i]),None)
  if k is not None:
   ix=[k]+[i for i in range(n) if i!=k];A=A.extract(ix,ix);v=A[0,0]
   pos+=int(bool(v>0));neg+=int(bool(v<0));steps.append({'block':1,'pivot':str(v)})
   A=A[1:,1:]-A[1:,0]*A[0,1:]/v
  else:
   pair=next(((i,j) for i in range(n) for j in range(i) if A[i,j]),None)
   if pair is None:zero+=n;steps.append({'zero_block':n});break
   i,j=pair;ix=[i,j]+[k for k in range(n) if k not in pair];A=A.extract(ix,ix);P=A[:2,:2];steps.append({'block':2,'offdiag':str(P[0,1])});pos+=1;neg+=1
   A=A[2:,2:]-A[2:,:2]*P.inv()*A[:2,2:]
 return [pos,zero,neg],steps
for r in rows:
 r['signature_exact'],r['rational_congruence_steps']=inertia(r['matrix'])
 assert r['signature_exact']==r['signature_float']
print(json.dumps(rows,indent=2))
