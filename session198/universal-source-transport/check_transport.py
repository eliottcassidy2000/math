"""Finite exact certificate for an all-word universal-transport obstruction."""
import json,itertools,hashlib

def need(b,s):
 if not b:raise RuntimeError(s)
GEN={
'O':((1,0,1,0),(1,2,1,2),(0,0,0,0),(0,0,1,1)),
'Z':((2,1,1,1),(0,1,0,0),(0,0,1,1),(0,0,0,0)),
'J':((0,1,0,1),(0,0,0,0),(1,0,1,1),(0,0,1,0))}
P={
'N':((6,4,10,9),(4,3,7,6),(10,7,15,14),(9,6,14,13)),
'M':((11,8,20,18),(8,6,14,13),(20,14,34,31),(18,13,31,29)),
'A':((14,10,21,20),(10,7,16,15),(21,16,33,31),(20,15,31,29)),
'B':((14,10,26,23),(11,8,20,18),(20,15,35,32),(19,14,34,31)),
'C':((23,18,44,40),(18,14,34,31),(44,34,80,73),(40,31,73,67)),
'T':((14,11,20,19),(10,8,15,14),(26,20,35,34),(23,18,32,31)),
'J0':((12,12,24,24),(12,12,21,21),(15,12,30,27),(15,12,30,27)),
'J1':((12,12,21,21),(12,12,18,18),(18,15,32,30),(18,15,31,29)),
'J2':((12,12,18,18),(12,12,15,15),(21,18,32,31),(21,18,30,29)),
'J3':((12,12,15,15),(12,12,12,12),(24,21,30,30),(24,21,27,27))}
K={'O':(1,-1,-1,1),'Z':(0,0,1,-1),'J':(1,1,0,-1)}
EX={'N':[(1,0,2,2),(1,1,1,1),(1,1,3,2)],'M':[(1,1,3,3),(2,1,3,2),(1,1,3,2)],'A':[(3,2,3,3),(1,1,2,2),(4,2,6,6)],'B':[(1,1,2,2),(3,2,3,3),(1,1,3,2)],'C':[(1,1,3,3),(4,3,7,6),(1,1,5,4)],'T':[(2,1,5,4),(1,1,1,1),(6,4,12,10)],'J1':[(0,0,1,1),(0,0,2,2),(3,6,3,4)],'J2':[(0,0,2,2),(0,0,1,1),(6,9,8,10)]}
def tr(A):return tuple(zip(*A))
def mv(A,v):return tuple(sum(a*b for a,b in zip(row,v))for row in A)
def mm(A,B):return tuple(tuple(sum(a*b for a,b in zip(row,col))for col in tr(B))for row in A)
def profile(T,r):
 n=len(T);F=(1<<n)-1
 H=[0]*(F+1);H[0]=1;st=[[0]*n for _ in range(F+1)];en=[[0]*n for _ in range(F+1)]
 for S in range(1,F+1):
  for v in range(n):
   if not S>>v&1:continue
   rest=S^(1<<v)
   if not rest:st[S][v]=en[S][v]=1
   else:
    st[S][v]=sum(st[rest][x]for x in range(n)if rest>>x&1 and T[v]>>x&1)
    en[S][v]=sum(en[rest][x]for x in range(n)if rest>>x&1 and T[x]>>v&1)
  H[S]=sum(en[S])
 D=H[F];E=en[F][r];d=H[F^(1<<r)];U=sum(H[S]*st[F^S][r]for S in range(F+1))
 need(st[F][r]==1,'literal start-one source')
 return (D-E-d,U-D+E-d,E,d-E),(D,E,U,d)
def build(word):
 T=[2,0];r=0
 for x in word:
  n=len(T);T.append(0)
  for v in range(n):
   if x in 'JZ' and v==r:T[n]|=1<<v
   else:T[v]|=1<<n
  if x=='J':r=n
 return T,r
vectors=[];physical=[]
for word,expected in [('',(0,0,0,1)),('Z',(1,0,1,0)),('O',(0,2,0,1)),('ZJ',(0,0,2,1)),('ZJZ',(3,0,3,0))]:
 T,r=build(word);y,old=profile(T,r)
 need(y==expected,('physical span',word,y,expected));physical.append(dict(word=word,rows=T,root=r,positive=y,profile=old))
 if word!='ZJZ':vectors.append(y)
# Leibniz determinant uses columns given by the actual source vectors.
A=tr(vectors);det=0
for p in itertools.permutations(range(4)):
 term=(-1)**sum(p[i]>p[j]for i in range(4)for j in range(i+1,4))
 for i in range(4):term*=A[i][p[i]]
 det+=term
need(abs(det)==4,'physical linear span')
for letter,G in GEN.items():need(mv(G,K[letter])==(0,0,0,0),'generator kernel')
for family,expected in EX.items():
 for idx,letter in enumerate(K):
  v=mv(P[family],K[letter]);need(v==expected[idx],'kernel table')
  need(all(v[i]>=(0,0,1,1)[i]for i in range(4)),'physical strict positivity')
 need(tr(P[family])in [P[x]for x in EX],'transpose closure')
need(mv(P['J0'],K['O'])==(0,0,0,0),'J0 zero kernel')
need(mv(P['J3'],K['Z'])==(0,0,0,0),'J3 zero kernel')
for a,F,G,b in [('Z','J0','J0','O'),('J','J0','J3','J'),('O','J3','J3','Z')]:
 need(mm(tr(GEN[a]),P[F])==mm(P[G],GEN[b]),'surviving transport')
out=dict(status='PASS',kernel_checks=3,strict_table_vectors=24,physical_spanning_determinant=det,physical_witnesses=physical,surviving_identities=3,scope='All-word universal-identity obstruction only; target-dependent source fibers and exact coverage remain open.')
print(json.dumps(out,sort_keys=True,indent=2))
