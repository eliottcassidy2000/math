"""Independent unit-interface verification; assertions are not used."""
import json,math,hashlib
from pathlib import Path

def need(v,msg):
 if not v:raise RuntimeError(msg)

def graph(L):
 s=(L-1).bit_length();rows=[0]*(s+2);a=0;b=s+1
 for i in range(s):
  x=i+1
  for y in range(x+1,s+1):rows[x]|=1<<y
  rows[x if (L-1)>>i&1 else b]|=1<<(b if (L-1)>>i&1 else x)
 rows[b]|=1<<a
 for x in range(1,s+1):rows[a]|=1<<x
 F=(1<<(s+2))-1
 return [F^(1<<i)^r for i,r in enumerate(rows)],b

def count(T,b):
 n=len(T);F=(1<<n)-1
 end=[[0]*n for _ in range(F+1)];start=[[0]*n for _ in range(F+1)];hs=[1]+[0]*F
 for S in range(1,F+1):
  for v in range(n):
   if not S>>v&1:continue
   R=S^(1<<v)
   if R==0:end[S][v]=start[S][v]=1
   else:
    end[S][v]=sum(end[R][u]for u in range(n)if R>>u&1 and T[u]>>v&1)
    start[S][v]=sum(start[R][u]for u in range(n)if R>>u&1 and T[v]>>u&1)
  hs[S]=sum(end[S])
 U=sum(hs[S]*start[F^S][b]for S in range(F+1))
 L=sum(end[S][b]*hs[F^S]for S in range(F+1))
 # Removing b cuts cycles at b. Count paths of the deletion whose first
 # and last incidences meet b, independently of source-cover formulas.
 rootless=F^(1<<b);cyc=0
 for a in range(n):
  if a==b or not T[b]>>a&1:continue
  dp=[[0]*n for _ in range(F+1)];dp[1<<a][a]=1
  for S in range(1,F+1):
   if S&~rootless:continue
   for v in range(n):
    if dp[S][v]:
     for w in range(n):
      if rootless>>w&1 and not S>>w&1 and T[v]>>w&1:dp[S|1<<w][w]+=dp[S][v]
  cyc+=sum(dp[rootless][v]for v in range(n)if T[v]>>b&1)
 return hs[F],end[F][b],U,hs[rootless],start[F][b],L,cyc

def g(L):
 if L==1:return 1
 if L%2==0:return g(L//2)
 t=(L+1)//2;return g(t)+t-t%2

def eval_bits(w,mod=None):
 n=0;h=1
 for bit in w:
  if bit=='0':h+=n+n%2;n*=2
  else:n=2*n+1
  if mod:n%=mod;h%=mod
 return n,h

def compile_word(M,u,z):
 q=M;k=0
 while q%2==0:q//=2;k+=1
 a=2
 while pow(2,a,q)!=1%q:a+=1
 delta=((1<<a)-1)*((1<<k)-1)+1
 need(math.gcd(delta,M)==2,'loop content')
 r=(u-1)%M;m=M//2
 j=0 if m==1 else (((z-g(r+1))//2)*pow(delta//2,-1,m))%m
 j+=m
 return ('1'*a+'0'*k)*j+format(r,f'0{(M-1).bit_length()}b')

physical=[]
for L in range(2,129):
 T,b=graph(L);p=count(T,b)
 expected=(g(L)+L-L%2,1,2*L-1,1,1)
 need(p[:5]==expected and p[6]==1,('physical interface',L,p,expected))
 # Independently: the whole root-deleted graph is a transitive chain.
 deletion=[i for i in range(len(T))if i!=b]
 need(all(T[deletion[j]]>>deletion[i]&1 for i in range(len(deletion))for j in range(i+1,len(deletion))),'converse transitive deletion')
 physical.append((L,*p))
T,b=graph(1);p=count(T,b)
need(p[4]==0 and p[6]==0,'excluded singleton boundary')
checks=0;exact=0;max_bits=0;digest=hashlib.sha256()
for M in range(2,42,2):
 for H in range(1,M,2):
  for U in range(1,M,2):
   u=(U+1)//2;z=(H-u+u%2)%M
   w=compile_word(M,u,z);n,h=eval_bits(w,M)
   L=(n+1)%M
   need(L==u and h==z,'joint compiler')
   need((h+L-L%2)%M==H and (2*L-1)%M==U,'interface targets')
   if M<=10:
    nn,hh=eval_bits(w);actual=nn+1
    need(hh==g(actual),'two independent integer recurrences')
    need((hh+actual-actual%2)%M==H and (2*actual-1)%(2*M)==U,'exact strengthened cover lift');exact+=1
   digest.update(f'{M}:{H}:{U}:{w}\n'.encode());max_bits=max(max_bits,len(w));checks+=1
out=dict(status='PASS',literal_L=[2,128],physical_graphs=len(physical),physical_profile_hash=hashlib.sha256(json.dumps(physical).encode()).hexdigest(),joint_interface_tests=checks,even_moduli=[2,40],exact_big_integer_replays=exact,max_word_bits=max_bits,word_digest=digest.hexdigest(),singleton_boundary='correctly excluded',scope='All-order unit-interface congruence theorem only, not exact-value coverage.')
print(json.dumps(out,sort_keys=True,indent=2))
