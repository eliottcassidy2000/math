from pathlib import Path
import json, sympy as S, hashlib, subprocess
root=Path(__file__).resolve().parents[2]/'session220-ribbon'
if not root.exists(): root=Path(__file__).resolve().parents[2]/'ribbon'

def load(w):
 data=list(map(int,(root/f'span{w}_intrep.txt').read_text().split()));r,L=data[:2]
 init=S.Matrix(data[2:2+r]);out=S.Matrix([data[2+r:2+2*r]]);pos=2+2*r; mats={}
 for _ in range(L):
  letter=data[pos];pos+=1;mats[letter]=S.Matrix(r,r,data[pos:pos+r*r]);pos+=r*r
 assert pos==len(data)
 return mats,init,out

def literal(word):
 n=len(word);pre=[0]*n
 # A predecessor u of vertex v points to v, following the literal back-arc grammar.
 for u in range(n):
  for v in range(u+1,n):
   back=(v-u<=2 and (word[v]>>(v-u-1)&1))
   if back:pre[u]|=1<<v
   else:pre[v]|=1<<u
 dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if not mask>>v&1:continue
   old=mask^(1<<v);bits=old&pre[v]
   while bits:
    bit=bits&-bits;bits-=bit;dp[mask][v]+=dp[old][bit.bit_length()-1]
 return sum(dp[-1])

A,v,h=load(2);families=[]
for word in ((0,2,0,3,2,3,3,3),(0,2,1,3,3,2,3,2)):
 Q=S.eye(14)
 for c in word:Q=A[c]*Q
 x=S.Symbol('x')
 assert Q.rank()==2
 assert Q.charpoly(x).as_expr()==x**14-267*x**13+530*x**12
 assert Q**3-267*Q**2+530*Q==S.zeros(14)
 assert int((h*v)[0])==0
 counts=[literal(word),literal(word*2)]
 assert counts==[65,16825]
 for k in range(1,13):
  numerator=63*265**k+200*2**k
  assert numerator%263==0
  assert int((h*Q**k*v)[0])==numerator//263
 families.append({'word':''.join(map(str,word)),'literal_counts_k1_k2':counts,'algebraic_k_checked':12})
 # Cyclic rotations preserve characteristic polynomial but need not preserve physical initial counts.
 for i in range(8):
  R=S.eye(14)
  for c in word[i:]+word[:i]:R=A[c]*R
  assert (R-2*S.eye(14)).det()==0 and (R-S.eye(14)).det()!=0
bounds={}
for w in (2,3):
 mats,_,_=load(w);r=mats[0].rows;maxentry=max(abs(int(e)) for M in mats.values() for e in M)
 bounds[w]={'dimension':r,'max_letter_entry':maxentry,'product_abs_bound':r*maxentry*1012,'elimination_abs_bound':1012**2+1012}
 assert max(bounds[w]['product_abs_bound'],bounds[w]['elimination_abs_bound'])<2**63
report={'status':'PASS','families':families,'modular_int64_bounds':bounds,'full_census_rerun':False,'scope':'Independent exact matrix algebra; literal subset DP at 8 and 16 vertices; sixteen rotations tested for eigenvalues 2 and 1; all recorded modular candidates rechecked; full census code reviewed but not rerun.'}
# Recheck all recorded first-prime candidates with independent modular elimination.
def singular(M,p):
 rows=[[int(t)%p for t in M.row(i)] for i in range(M.rows)]
 n=len(rows)
 for j in range(n):
  pivot=next((i for i in range(j,n) if rows[i][j]),None)
  if pivot is None:return True
  rows[j],rows[pivot]=rows[pivot],rows[j]
  inv=pow(rows[j][j],-1,p)
  for i in range(j+1,n):
   coeff=rows[i][j]*inv%p
   if coeff:
    for k in range(j,n):rows[i][k]=(rows[i][k]-coeff*rows[j][k])%p
 return False
candidate_counts=[]
for w,d in ((2,9),(3,5)):
 mats,_,_=load(w);r=mats[0].rows
 candidates=json.loads((root/f'screen{w}_{d}.json').read_text()); survivors=[]
 for word in candidates:
  Q=S.eye(r)
  for c in word:Q=mats[c]*Q
  assert singular(Q-2*S.eye(r),1009)
  if singular(Q-2*S.eye(r),1013):survivors.append(word)
 assert len(candidates)==(50 if w==2 else 30)
 assert len(survivors)==(16 if w==2 else 0)
 candidate_counts.append({'span':w,'first_prime_candidates':len(candidates),'second_prime_survivors':len(survivors)})
report['candidate_recheck']=candidate_counts
print(json.dumps(report,indent=2))
