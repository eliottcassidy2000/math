from pathlib import Path
from itertools import permutations
from collections import defaultdict
import json,sympy as s
HERE=Path(__file__).resolve().parent;BASE=HERE.parent/'session29-profiles'
D=json.loads((BASE/'seamed-results.json').read_text());C=json.loads((BASE/'compressed-matrices.json').read_text());W=json.loads((HERE/'wiring-results.json').read_text())
pairs=[(i,j) for i in range(5) for j in range(i+1,5)]
def require(b,msg):
 if not b:raise RuntimeError(msg)
def adj(m):
 a=[[0]*5 for _ in range(5)]
 for k,(i,j) in enumerate(pairs):a[i][j]=(m>>k)&1;a[j][i]=1-a[i][j]
 return a
def op(a):return list(map(list,zip(*a)))
def canon(a,b,c):
 best=None;auto=0
 for mid in permutations([i for i in range(5) if i not in (b,c)]):
  p=(b,)+mid+(c,);m=sum(a[p[i]][p[j]]<<k for k,(i,j) in enumerate(pairs))
  if best is None or m<best:best=m;auto=1
  elif m==best:auto+=1
 return best,auto
flip=1023^sum(1<<k for k,(i,j) in enumerate(pairs) if j==i+1)
bi=defaultdict(lambda:{'A':[],'B':[],'mirror_fixed_tiles':{'A':[],'B':[]}});tilebi={}
for seed,tiles in D['tile_alphabet'].items():
 typ='A' if seed=='40' else 'B'
 for t in tiles:
  x=adj(t^flip);b,auto=canon(x,0,4);mb,_=canon(op(x),4,0);bi[b][typ].append(t);bi[b]['mirror_class']=mb;bi[b]['automorphism_order']=auto;tilebi[t]=b
  tm=sum((1-adj(t)[4-i][4-j])<<k for k,(i,j) in enumerate(pairs))
  if tm==t:bi[b]['mirror_fixed_tiles'][typ].append(t)
for b,d in bi.items():
 ts=d['A']+d['B'];d['K']=C[str(ts[0])]
 require(all(C[str(t)]==d['K'] for t in ts),'bi-root transfer descends')
 require(bi[d['mirror_class']]['mirror_class']==b,'mirror involution')
outer={}
for q,d in W['root_classes'].items():
 ts=d['A']+d['B'];v=[row[0] for row in C[str(ts[0])]]
 require(all([row[0] for row in C[str(t)]]==v for t in ts),'exit-root column descends')
 outer[q]={'A':len(d['A']),'B':len(d['B']),'v':v}
B=sum((s.Matrix(d['K']) for d in bi.values()),s.zeros(3));v=sum((s.Matrix(d['v']) for d in outer.values()),s.zeros(3,1));J=s.diag(1,-1,1);z,w=s.symbols('z w')
require(J*B==B.T*J,'summed distinct-class adjoint')
ordinary=s.cancel(z*z*(v.T*J*(s.eye(3)-z*B).inv()*v)[0]);r2=int((v.T*J*v)[0]);table=W['incidence_table'];mass2=sum(x['H'] for x in table);fix2=sum(x['H'] for x in table if x['self_converse'])
require(r2==2*mass2-fix2,'unweighted r2 H mass')
# Six raw second moments in coordinate order H²,S²,R²,HS,HR,SR.
x=s.symbols('h t r');mons=[x[0]**2,x[1]**2,x[2]**2,x[0]*x[1],x[0]*x[2],x[1]*x[2]]
L=s.zeros(6)
for d in bi.values():
 y=s.Matrix(d['K'])*s.Matrix(x)
 for i,p in enumerate([y[0]**2,y[1]**2,y[2]**2,y[0]*y[1],y[0]*y[2],y[1]*y[2]]):
  poly=s.Poly(p,*x)
  for j,m in enumerate(mons):L[i,j]+=poly.coeff_monomial(m)
m0=s.Matrix([sum(m.subs(dict(zip(x,d['v']))) for d in outer.values()) for m in mons]);even=s.Matrix([[1,-1,1,0,0,0]])
F=sum((s.Matrix(d['K']) for b,d in bi.items() if d['mirror_class']==b),s.zeros(3));poly=s.Poly((s.Matrix(x).T*J*F*s.Matrix(x))[0],*x);odd=s.Matrix([[poly.coeff_monomial(m) for m in mons]])
E=s.cancel((even*(s.eye(6)-w*L).inv()*m0)[0]);O=s.cancel((odd*(s.eye(6)-w*L).inv()*m0)[0]);fixed=z*z*E.subs(w,z*z)+z**3*O.subs(w,z*z);merged=s.cancel((ordinary+fixed)/2)
a=len(outer);b=len(bi);f=sum(d['mirror_class']==k for k,d in bi.items());fblue=sum(any(d['mirror_fixed_tiles'].values()) for k,d in bi.items() if d['mirror_class']==k)
counts=[]
for r in range(2,9):
 total=a*a*b**(r-2);fixedcount=a*b**((r-2)//2)*(f if r%2 else 1);counts.append({'r':r,'ordinary':total,'fixed':fixedcount,'merged':(total+fixedcount)//2})
result={'outer_classes':outer,'birooted_classes':{str(k):v for k,v in sorted(bi.items())},'a':a,'b':b,'f':f,'fixed_classes_with_fixed_tile':fblue,'B':list(map(list,B.tolist())),'v':list(v),'ordinary_H_GF':str(ordinary),'ordinary_H_r2':r2,'merged_H_r2':mass2,'fixed_H_r2':fix2,'moment_L':L.tolist(),'moment_initial':list(m0),'fixed_even_H_GF_in_w':str(E),'fixed_odd_H_GF_in_w':str(O),'merged_H_GF':str(merged),'target_class_counts':counts}
(HERE/'all-length-results.json').write_text(json.dumps(result,indent=2,default=int)+'\n')
print(json.dumps({k:result[k] for k in ['a','b','f','fixed_classes_with_fixed_tile','B','v','ordinary_H_GF','ordinary_H_r2','merged_H_r2','fixed_H_r2','target_class_counts']},indent=2,default=int))
