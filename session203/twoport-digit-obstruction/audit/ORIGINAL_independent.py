import sys,pathlib,itertools,json
from collections import Counter
from fractions import Fraction as Q
PTH=pathlib.Path(__file__).resolve().parents[1]/'session203-twoport-digit-obstruction';sys.path.insert(0,str(PTH));import verify as author
P=author.P;N=len(P);ix={p:i for i,p in enumerate(P)}

def trans(p,mask,drop,swap):
 out=Counter();allowed=set(i for i in range(2)if mask>>i&1)
 for pre in [None]+list(range(len(p))):
  if pre is not None and p[pre][1] in allowed:continue
  for suc in [None]+list(range(len(p))):
   if suc is not None and p[suc][0]not in allowed:continue
   if pre is not None and pre==suc:continue
   parts=[x for k,x in enumerate(p)if k!=pre and k!=suc]+[(2 if pre is None else p[pre][0],2 if suc is None else p[suc][1])]
   keep=[i for i in range(3)if i!=drop]
   if swap:keep=keep[::-1]
   rel={-1:-1,drop:-1,keep[0]:0,keep[1]:1}
   q=tuple(sorted((rel[a],rel[b])for a,b in parts))
   if sum(a==-1 for a,b in q)<=2:out[q]+=1
 return out
M=[]
for op in author.OPS:
 mat=[[0]*N for _ in range(N)]
 for j,p in enumerate(P):
  for q,c in trans(p,*op).items():mat[ix[q]][j]+=c
 M.append(mat)
for k,mat in enumerate(M):
 for j in range(N):assert [(i,mat[i][j])for i in range(N)if mat[i][j]]==sorted(author.M[k][j])

def row(v,mat):return [sum(v[i]*mat[i][j]for i in range(N)if v[i] and mat[i][j])for j in range(N)]
def determinant(mat):
 a=[list(map(Q,r))for r in mat];d=Q(1)
 for k in range(N):
  hit=next((i for i in range(k,N)if a[i][k]),None)
  if hit is None:return 0
  if hit!=k:a[k],a[hit]=a[hit],a[k];d=-d
  v=a[k][k];d*=v
  for i in range(k+1,N):
   if not a[i][k]:continue
   z=a[i][k]/v
   for j in range(k+1,N):a[i][j]-=z*a[k][j]
   a[i][k]=0
 return d
D=[]
for mat in M:D.append(determinant([[v-2*(i==j)for j,v in enumerate(r)]for i,r in enumerate(mat)]))
valid=[i for i,d in enumerate(D)if not d];av,ad=author.eigen2();assert valid==av and D==ad
cert=json.loads((PTH/'CERTIFICATE.json').read_text());seen=set();maxden=1
for case in cert['cases']:
 ci,a,b=case['ci'],case['a'],case['b'];key=(ci,a,b);assert key not in seen;seen.add(key)
 ell=[3*x+2*y for x,y in zip(*author.C[ci])];u=row(ell,M[a]);v=row(ell,M[b]);r0=[x-2*y for x,y in zip(u,ell)];delta=[y-x for x,y in zip(u,v)];H=[]
 for root in [r0,delta]:
  for depth in range(1,6):
   for word in itertools.product((a,b),repeat=depth):
    rr=root
    for o in word:rr=row(rr,M[o])
    H.append([x-y for x,y in zip(rr,root)])
 H += [[int(j==i)for j in range(N)]for i in author.BAD]
 assert len(H)==128
 G=[6*int(j==author.BULK)-delta[j]for j in range(N)];dual={int(i):Q(v)for i,v in case['dual']};z=dual.pop(128,Q(0));assert z>=0 and all(0<=i<128 for i in dual)
 for q in dual.values():maxden=max(maxden,q.denominator)
 slack=[sum(y*H[i][j]for i,y in dual.items())-z*G[j]-delta[j]for j in range(N)];assert min(slack)>=0
assert seen==set(itertools.product(range(6),valid,valid))
print(json.dumps(dict(status='PASS',method='Independent direct endpoint rewrite, dense integer rows, Fraction Gaussian determinant, explicit word enumeration and rational slack verification',dimension=N,matrices=len(M),eigen2_moves=valid,exact_cases=len(seen),max_denominator=maxden),indent=2))
