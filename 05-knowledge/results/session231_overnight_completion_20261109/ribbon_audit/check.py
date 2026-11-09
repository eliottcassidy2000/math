def require(condition, message):
    if not condition:
        raise RuntimeError(message)

import os,sys,json,hashlib,itertools
from pathlib import Path
import numpy as np
import sympy as S
ROOT=Path(os.environ.get('MATH_REPO',str(Path(__file__).resolve().parents[4])))
SRC=ROOT/'04-computation/experiments/ribbon_spectrum_20261104'
sys.path.insert(0,str(SRC))
P=1000003
require(S.isprime(P), 'check failed: S.isprime(P)')
def emit(**x):print(json.dumps(x),flush=True)
def read(w):
 path=SRC/f'span{w}_intrep.txt'; data=list(map(int,path.read_text().split())); it=iter(data)
 n,q=next(it),next(it); y=[next(it) for _ in range(n)]; h=[next(it) for _ in range(n)]
 A={}
 for _ in range(q):
  mask=next(it); require(mask not in A, 'check failed: mask not in A')
  A[mask]=np.array([[next(it) for _ in range(n)] for _ in range(n)],dtype=np.int64)
 require(next(it,None) is None and set(A)==set(range(2**w)), 'check failed: next(it,None) is None and set(A)==set(range(2**w))')
 emit(input=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),dimension=n,mask_order=list(A))
 return A,y,h

def quotient(w,A,y,h):
 from ribbon import transfer_system
 order,states,mats,final=transfer_system(w)
 fulls={}
 for letter,M in mats.items():
  full=S.zeros(len(order))
  for (j,i),v in M.items():full[j,i]=v
  fulls[sum(1<<(s-1) for s in letter)]=full
 phi=S.zeros(1,len(order))
 for j in final:phi[j]=1
 # Recover an intertwiner in the ACTUAL stored basis, without assuming a
 # platform-dependent LLL basis is reproduced byte-for-byte.
 n=len(y); pairs=[]; echelon={}
 def offer(c,d):
  v=[int(z)%P for z in c]
  for j,row in sorted(echelon.items()):
   if v[j]:
    factor=v[j];v=[(x-factor*z)%P for x,z in zip(v,row)]
  j=next((j for j,x in enumerate(v) if x),None)
  if j is None:return
  inv=pow(v[j],-1,P);echelon[j]=[(x*inv)%P for x in v]
  pairs.append((c,d))
 offer(np.array(h,dtype=object),np.array(list(phi),dtype=object))
 fullnp={m:np.array(F.tolist(),dtype=object) for m,F in fulls.items()}
 cursor=0
 while len(pairs)<n and cursor<len(pairs):
  c,d=pairs[cursor];cursor+=1
  for m in A:offer(c@A[m].astype(object),d@fullnp[m])
 require(len(pairs)==n, 'check failed: len(pairs)==n')
 C=S.Matrix([list(c) for c,d in pairs]);D=S.Matrix([list(d) for c,d in pairs])
 B=C.inv()*D
 require(list(B[:,0])==y, 'check failed: list(B[:,0])==y')
 require(S.Matrix([h])*B==phi, 'check failed: S.Matrix([h])*B==phi')
 for mask,full in fulls.items():require(S.Matrix(A[mask])*B==B*full, 'check failed: S.Matrix(A[mask])*B==B*full')
 emit(quotient=w,profile_dimension=len(order),quotient_dimension=B.rows,identities='PASS')

def nonsingular_batch(arr):
 # Batched modular row reduction, no producer elimination routine imported.
 R=arr.copy()%P; size,n,_=R.shape; ok=np.ones(size,dtype=bool)
 for c in range(n):
  eligible=np.flatnonzero(ok)
  if not len(eligible):break
  sub=R[eligible,c:,c]; has=np.any(sub!=0,axis=1)
  ok[eligible[~has]]=False; ids=eligible[has]
  if not len(ids):continue
  piv=c+np.argmax(sub[has]!=0,axis=1)
  tmp=R[ids,c,:].copy(); R[ids,c,:]=R[ids,piv,:];R[ids,piv,:]=tmp
  inv=np.array([pow(int(v),-1,P) for v in R[ids,c,c]],dtype=np.int64)
  factors=(R[ids,c+1:,c]*inv[:,None])%P
  R[ids,c+1:,c+1:]=(R[ids,c+1:,c+1:]-factors[:,:,None]*R[ids,c,c+1:][:,None,:])%P
 return ok

def scan(w,A,alphabet,maxlen):
 n=len(next(iter(A.values()))); front=np.eye(n,dtype=np.int64)[None,:,:]; total=0
 maxA=max(int(np.max(np.abs(A[m]))) for m in alphabet)
 for length in range(1,maxlen+1):
  maxprev=int(np.max(np.abs(front)))
  multiplication_bound=n*maxA*maxprev
  require(multiplication_bound+8<2**63, 'check failed: multiplication_bound+8<2**63')
  front=np.concatenate([A[m]@front for m in alphabet],axis=0)
  require(len(front)==len(alphabet)**length, 'check failed: len(front)==len(alphabet)**length')
  counts={a:0 for a in (2,4,8)}
  for start in range(0,len(front),256):
   batch=front[start:start+256]
   for a in counts:
    ok=nonsingular_batch(batch-a*np.eye(n,dtype=np.int64))
    counts[a]+=int(np.sum(~ok))
  require(not any(counts.values()), counts)
  total+=len(front)
  emit(span=w,alphabet=alphabet,length=length,level_words=len(front),cumulative=total,modular_singular=counts,max_entry=int(np.max(np.abs(front))),certified_dot_product_bound=multiplication_bound)
 return total

A2,y2,h2=read(2); A3,y3,h3=read(3)
quotient(2,A2,y2,h2);quotient(3,A3,y3,h3)
require(scan(2,A2,[0,1,2,3],6)==5460, 'check failed: scan(2,A2,[0,1,2,3],6)==5460')
require(scan(3,A3,list(range(8)),3)==584, 'check failed: scan(3,A3,list(range(8)),3)==584')
require(scan(2,A2,[1,2,3],9)==29523, 'check failed: scan(2,A2,[1,2,3],9)==29523')
word=[3,1,2,3,3,3,1,2,1]; M=S.eye(14)
for a in word:M=S.Matrix(A2[a])*M
x=S.Symbol('x'); cp=M.charpoly(x).as_expr()
require(S.expand(cp-x**11*(x**3-971*x**2-1302*x+96))==0, 'check failed: S.expand(cp-x**11*(x**3-971*x**2-1302*x+96))==0')
require((M**14).rank()==3, 'check failed: (M**14).rank()==3')
# rank over F2 using Python integer bit rows
rows=[sum((int(v)%2)<<j for j,v in enumerate(row)) for row in (M**14).tolist()]
r=0
for c in range(14):
 pivot=next((j for j in range(r,14) if rows[j]>>c&1),None)
 if pivot is None:continue
 rows[r],rows[pivot]=rows[pivot],rows[r]
 for j in range(r+1,14):
  if rows[j]>>c&1:rows[j]^=rows[r]
 r+=1
require(r==1, 'check failed: r==1')
emit(hostile_word=word,charpoly=str(cp),stable_rank_Q=3,stable_rank_F2=r,status='PASS')
# Positive control: known odd-radix tripling block.
B=S.Matrix(A2[2])*S.Matrix(A2[1]); require(S.expand(B.charpoly(x).as_expr()-x**13*(x-3))==0, 'check failed: S.expand(B.charpoly(x).as_expr()-x**13*(x-3))==0')
require(not nonsingular_batch(np.array(B.tolist(),dtype=np.int64)[None,:,:]-3*np.eye(14,dtype=np.int64)).all(), 'check failed: not nonsingular_batch(np.array(B.tolist(),dtype=np.int64)[None,:,:]-3*np.eye(14,dtype=np.int64)).all()')
emit(positive_control='tripling eigenvalue 3 detected',status='PASS')
