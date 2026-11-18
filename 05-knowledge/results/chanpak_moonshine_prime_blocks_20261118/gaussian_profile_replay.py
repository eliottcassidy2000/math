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
from math import gcd,isqrt
from itertools import permutations

def ordinal(A,B):
 n=len(A);m=len(B)
 return [a+(((1<<m)-1)<<n) for a in A]+[b<<n for b in B]
def shell(T):
 n=len(T)
 return [a|(1<<(n+1)) for a in T]+[(1<<n)-1,1<<n]
def vp(n,p):
 z=0
 while n%p==0:z+=1;n//=p
 return z
def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
seeds={'one':[0],'TT2':[0,1],'C3':[2,4,1],'H5':[8,1,3,6],'H11':[24,1,3,6,14],'rigidH9':[16,1,3,7,14]}
profiles={k:profile(t) for k,t in seeds.items()}
records=[]
for ka,kb in [('TT2','TT2'),('H5','H5'),('one','H11'),('TT2','H5'),('rigidH9','rigidH9')]:
 A,B=seeds[ka],seeds[kb];ha,pa=profiles[ka]['H'],profiles[ka]['P2'];hb,pb=profiles[kb]['H'],profiles[kb]['P2'];N=ha*hb+pa*pb;D=ha*pb-hb*pa;T=shell(ordinal(A,B));actual=profile(T)['H']
 require(actual==N);require(strong(T));require(N*N+D*D==(ha*ha+pa*pa)*(hb*hb+pb*pb));require(2*abs(D)<N)
 if min(len(A),len(B))>=2:require(D%4==0)
 records.append({'sources':[ka,kb],'N':N,'D':D,'prime':prime(N),'rows':T})
for k,t in seeds.items():
 h,p=profiles[k]['H'],profiles[k]['P2'];require(p>=2*h)
 if len(t)>=2:require(p%4==0 and (h*h+p*p)%8==1)
 for ell in [3,7,11,19,23,31,43,47]:require(vp(h*h+p*p,ell)==2*min(vp(h,ell),vp(p,ell)))
C=seeds['rigidH9'];n=len(C);aut=sum(all(((C[i]>>j)&1)==((C[p[i]]>>p[j])&1) for i in range(n) for j in range(n)) for p in permutations(range(n)))
require(aut==1 and gcd(profiles['rigidH9']['H'],profiles['rigidH9']['P2'])==3)
family=[]
for m in range(1,7):
 A=[(1<<i)-1 for i in range(m)];T=shell(ordinal(A,seeds['H11']));H=profile(T)['H'];require(H==11+100*2**m);family.append({'m':m,'H':H,'prime':prime(H)})
print(json.dumps({'profiles':profiles,'examples':records,'H11_dyadic_family':family,'rigidH9_automorphisms':aut,'status':'PASS'},indent=2))
