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
from itertools import permutations

def odd_cycles(T,k):
 count=0
 for V in combinations(range(len(T)),k):
  for rest in permutations(V[1:]):
   w=(V[0],)+rest
   count+=all((T[w[i]]>>w[(i+1)%k])&1 for i in range(k))
 return count
A=[0,1]
Cs=[[16,1,3,7,14],[6,16,26,19,1]]
require(Cs[1]==[sum(1<<j for j in range(5) if j!=i and ((i<j)!=(i+j+2 in {5,6}))) for i in range(5)])
contacts=[[(1,3),(0,0)],[(1,3),(1,0)]]
records=[]
for C,ports in zip(Cs,contacts):
 require(odd_cycles(C,3)==3 and odd_cycles(C,5)==1)
 T=build(A,C); require(profile(T)['H']==45)
 row={'C':C,'C_profile':profile(C),'cycles3':3,'cycles5':1,'host_rows':T,'host_profile':profile(T),'returns':[]}
 for i,j in ports:
  u=i;v=3+j;old=arc_cover(T,u,v);new=arc_cover(T,v,u);U=T.copy();U[u]^=1<<v;U[v]^=1<<u
  P=profile(U);require(P['H']==45+new-old);require(strong(U))
  row['returns'].append({'arc':[u,v],'removed':old,'added':new,'rows':U,'profile':P})
 records.append(row)
require([[x['profile']['H'] for x in r['returns']] for r in records]==[[131,133],[141,143]])
require(sorted(records[0]['C_profile']['start'])!=sorted(records[1]['C_profile']['start']))
require(sorted(records[0]['C_profile']['end'])!=sorted(records[1]['C_profile']['start']))
print(json.dumps({'records':records,'common_C_slot_polynomial_coefficients':[0,2,0,6,0,1],'common_host_slot_polynomial_coefficients':[0,0,0,8,0,26,0,10,0,1],'status':'PASS'},indent=2))
