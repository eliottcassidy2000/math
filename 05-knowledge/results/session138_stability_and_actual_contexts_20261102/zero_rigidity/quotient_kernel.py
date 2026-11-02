"""Exact run-word quotient compiler, Python arbitrary-precision integers.

No inherited path counter is used here. The module's optional audit imports
one only to compare all coefficients after independently compiling them.
"""
from functools import lru_cache

# Local ordered path-cover/run weights; ratios f(r)/f(r-1), f(0)=1.
WEIGHTS={'C3':(3,2,1),'TT2':(1,2),'K1':(1,)}

def partition(q,kinds,gradient=False):
 """Return Z and, optionally, all ordered-edge derivatives at q.

 Quotient edges are independently weighted indeterminates; diagonal edges
 are forbidden because consecutive identical block runs must be merged.
 Counts include all nonempty-run counts up to each block's capacity.
 """
 k=len(q); steps=[WEIGHTS[t] for t in kinds]; bases=[len(t)+1 for t in steps]; strides=[1]
 for b in bases[:-1]:strides.append(strides[-1]*b)
 N=strides[-1]*bases[-1] if k else 1
 digits=[tuple((s//strides[i])%bases[i] for i in range(k)) for s in range(N)]
 @lru_cache(None)
 def suffix(s,i):
  r=digits[s]; z=int(all(r))
  for j in range(k):
   if j!=i and r[j]<len(steps[j]) and q[i]>>j&1:
    z+=steps[j][r[j]]*suffix(s+strides[j],j)
  return z
 if not k:return (1,[]) if gradient else 1
 Z=sum(steps[i][0]*suffix(strides[i],i) for i in range(k))
 if not gradient:return Z
 F=[[0]*k for _ in range(N)]
 for i in range(k):F[strides[i]][i]=steps[i][0]
 G=[[0]*k for _ in range(k)]
 for s in range(1,N):
  r=digits[s]
  for i,f in enumerate(F[s]):
   if not f:continue
   for j in range(k):
    if j==i or r[j]==len(steps[j]):continue
    sn=s+strides[j]; c=steps[j][r[j]]
    G[i][j]+=f*c*suffix(sn,j)
    if q[i]>>j&1:F[sn][j]+=f*c
 if sum(sum(row) for s,row in enumerate(F) if all(digits[s]))!=Z:raise RuntimeError('forward/backward mismatch')
 return Z,G

def kernel(q,root):
 """Deleted-triangle parent H and full symmetric capacities, block labels.

 Parent labels: each block in quotient order; root's two surviving vertices
 are u,w with u->w; intact blocks have the standard cyclic labels0,1,2.
 """
 k=len(q); kinds=['C3']*k;kinds[root]='TT2';h,G=partition(q,kinds,True)
 Cblock=[[G[i][j]+G[j][i] for j in range(k)] for i in range(k)]
 others=[i for i in range(k) if i!=root]
 D=tuple(sum(1<<j for j,v in enumerate(others) if q[u]>>v&1) for u in others)
 _,DG=partition(D,['C3']*(k-1),True)
 DC=[[0]*(k-1) for _ in range(k-1)]
 for i in range(k-1):
  for j in range(i):
   z=DG[i][j]+DG[j][i]
   if z%9:raise RuntimeError('intact cross divisor')
   DC[i][j]=DC[j][i]=z//9
 # Within-block capacity in D for its local ear-gradient.
 DI=[]
 for i in range(k-1):
  types=['C3']*(k-1);types[i]='TT2';DI.append(2*partition(D,types))
 d={}
 for i,b in enumerate(others):
  # Core ear gradient at an all/none block cut R. For fixed vertex a in
  # block i, its two blockmates contribute DI_i*(2 r_i-1).
  ri=int(q[root]>>b&1)
  v=(2*ri-1)*DI[i]
  for j,c in enumerate(others):
   if j!=i:v+=3*DC[i][j]*(int(q[root]>>c&1)-int(q[b]>>c&1))
  d[b]=v
 sizes=[3]*k;sizes[root]=2;offset=[];z=0
 for m in sizes:offset.append(z);z+=m
 C=[[0]*z for _ in range(z)]
 for i in range(k):
  types=kinds[:];types[i]='K1' if i==root else 'TT2';inside=2*partition(q,types)
  for a in range(sizes[i]):
   for b in range(a):C[offset[i]+a][offset[i]+b]=C[offset[i]+b][offset[i]+a]=inside
  for j in range(i):
   if i==root or j==root:
    b=j if i==root else i
    v=Cblock[root][b]
    if v%3:raise RuntimeError('root sum divisor')
    ss=v//3
    if (ss+d[b])%2:raise RuntimeError('root split parity')
    for t in range(3):
     for a,cc in enumerate(((ss+d[b])//2,(ss-d[b])//2)):
      C[offset[root]+a][offset[b]+t]=C[offset[b]+t][offset[root]+a]=cc
   else:
    v=Cblock[i][j]
    if v%9:raise RuntimeError('cross divisor')
    for a in range(3):
     for b in range(3):C[offset[i]+a][offset[j]+b]=C[offset[j]+b][offset[i]+a]=v//9
 return h,C

if __name__=='__main__':
 import sys,itertools,json,hashlib
 from reference_paths import tournament,c3_blowup,H,endpoint_Q_dp
 rows=[]
 for k in range(1,5):
  for code in range(1<<(k*(k-1)//2)):
   q=tournament(k,code);T=c3_blowup(q);z=partition(q,['C3']*k)
   if z!=H(T):raise RuntimeError(('full count',k,code,z,H(T)))
   for root in range(k):
    vs=[i for i in range(3*k) if i!=3*root];B=tuple(sum(1<<j for j,v in enumerate(vs) if T[u]>>v&1) for u in vs)
    h,C=kernel(q,root);s,e,K=endpoint_Q_dp(B)
    CC=[[K[i][j]+K[j][i] for j in range(len(B))] for i in range(len(B))]
    if h!=sum(s) or C!=CC:raise RuntimeError(('kernel mismatch',k,code,root,h,sum(s)))
    rows.append((k,code,root,h,sum(map(sum,C))))
 print(json.dumps({'quotient_orders':[1,2,3,4],'labelled_quotients':75,'deletion_rows':len(rows),'all_coefficients_match':True,'sha256':hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest()},sort_keys=True))
