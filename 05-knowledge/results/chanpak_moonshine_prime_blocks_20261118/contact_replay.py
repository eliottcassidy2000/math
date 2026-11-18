from itertools import combinations
import json

def require(x,msg='check failed'):
 if not x:raise AssertionError(msg)
def dp(out):
 n=len(out);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for m in range(1,1<<n):
  for v in range(n):
   if m>>v&1:
    for u in range(n):
     if m>>u&1 and out[u]>>v&1:D[m][v]+=D[m^(1<<v)][u]
 return D

def reverse(out):return [sum(((out[j]>>i)&1)<<j for j in range(len(out))) for i in range(len(out))]
def profile(out):
 D=dp(out);E=dp(reverse(out));a=[sum(row) for row in D];a[0]=1;F=len(a)-1
 return {'H':a[-1],'P2':sum(a[s]*a[F^s] for s in range(F+1)),'start':E[-1],'end':D[-1]}
def build(A,C):
 r,s=len(A),len(C);n=r+s+2;b=r;x=n-1;out=[0]*n
 for i in range(r):out[i]=A[i] | sum(1<<j for j in range(r,r+s+1))
 out[b]=sum(1<<j for j in range(r+1,n))
 for i in range(s):out[r+1+i]=C[i]<<(r+1)
 out[x]=sum(1<<j for j in list(range(r))+list(range(r+1,r+s+1)))
 return out

def arc_cover(out,u,v):
 # Prefix ending u and complementary suffix beginning v; no seam arc required.
 D=dp(out);E=dp(reverse(out));F=(1<<len(out))-1
 return sum(D[s][u]*E[F^s][v] for s in range(F+1) if s>>u&1 and not s>>v&1)

def strong(out):
 for v in range(len(out)):
  seen=1<<v
  while True:
   z=seen
   for u in range(len(out)):
    if seen>>u&1:z|=out[u]
   if z==seen:break
   seen=z
  if seen!=(1<<len(out))-1:return False
 return True

def inflate(C,v,k):
 sizes=[k if i==v else 1 for i in range(len(C))];offset=[];n=0
 for s in sizes:offset.append(n);n+=s
 out=[0]*n
 for i in range(len(C)):
  for a in range(sizes[i]):
   out[offset[i]+a]|=sum(1<<(offset[i]+b) for b in range(a))
   for j in range(len(C)):
    if C[i]>>j&1:out[offset[i]+a]|=((1<<sizes[j])-1)<<offset[j]
 return out,offset
A=[0,1];C=[16,1,3,7,14];T=build(A,C);records=[]
require(profile(T)['H']==45)
for i,j in [(1,3),(0,0)]:
 u=i;v=3+j;old=arc_cover(T,u,v);new=arc_cover(T,v,u);U=T.copy();U[u]^=1<<v;U[v]^=1<<u
 p=profile(U);require(p['H']==45+new-old);require(strong(U))
 records.append({'contact':[u,v],'removed_paths':old,'added_paths':new,'delta':new-old,'rows':U,'profile':p})
require([r['profile']['H'] for r in records]==[131,133])
twins=[]
for v in range(5):
 values=[]
 for k in range(1,5):
  D,o=inflate(C,v,k);hs=[]
  for i,j in [(1,3),(0,0)]:
   U=build(A,D);u=i;w=3+o[j];U[u]^=1<<w;U[w]^=1<<u;hs.append(profile(U)['H'])
  values.append(hs)
 twins.append({'vertex':v,'counts_for_lengths_1_to_4':values})
print(json.dumps({'baseline':profile(T),'contacts':records,'twins':twins,'status':'PASS'},indent=2))
