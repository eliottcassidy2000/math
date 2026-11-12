"""Exact all-order inverse at prescribed receiver V; incoming mask means vertex beats root."""
from functools import lru_cache

def require(p,msg):
 if not p:raise ValueError(msg)
def endpoints(T):
 n=len(T);N=1<<n;dp=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for v in range(n):dp[1<<v][v]=1
 for s in range(1,N):
  for v in range(n):
   c=dp[s][v]
   if not c:continue
   H[s]+=c;t=T[v]&~s
   while t:
    b=t&-t;t-=b;dp[s|b][b.bit_length()-1]+=c
 return H,dp
def profile(T):
 n=len(T);F=(1<<n)-1;R=F^1;H,dp=endpoints(T)
 rev=tuple(sum(1<<j for j in range(n) if T[j]>>i&1) for i in range(n));Hr,dr=endpoints(rev)
 V=e=U=d=0
 for s in range(1,F+1,2):
  sg=1 if s.bit_count()%2 else -1;v=dp[s][0]*H[F^s];u=dr[s][0]*Hr[F^s]
  V+=v;e+=sg*v;U+=u;d+=sg*u
 Q=k=0;s=R
 while True:
  t=H[s]*H[R^s];Q+=t;k+=t*((-1)**s.bit_count())
  if not s:break
  s=(s-1)&R
 return dict(H=H[F],h=H[R],Q=Q,U=U,V=V,d=d,e=e,k=k)
def fuse(a,b):
 return dict(H=a['h']*b['H']+(a['H']-a['h'])*b['h']+(a['U']-a['h'])*(b['V']-b['h'])+(a['d']-a['h'])*(b['e']-b['h']),h=a['h']*b['h'],Q=a['Q']*b['Q'],U=a['h']*b['U']+(a['U']-a['h'])*b['Q'],V=a['Q']*b['V']+(a['V']-a['Q'])*b['h'],d=a['h']*b['d']+(a['d']-a['h'])*b['k'],e=a['k']*b['e']+(a['e']-a['k'])*b['h'],k=a['k']*b['k'])
I=dict(H=1,h=1,Q=1,U=1,V=1,d=1,e=1,k=1)
S=dict(H=1,h=1,Q=2,U=2,V=1,d=0,e=1,k=0)
E=dict(H=1,h=1,Q=2,U=1,V=2,d=1,e=0,k=0)
def unit(p,n=None):
 if n is None:n=p.bit_length()
 require(n>=0 and 0<=p<1<<n,'unit bit range');z=I.copy()
 for i in range(n):z=fuse(z,E if p>>i&1 else S)
 return z
def rooted(D,x):
 n=len(D);require(0<=x<1<<n,'cut range')
 return tuple([(((1<<n)-1)^x)<<1]+[(D[i]<<1)|(x>>i&1) for i in range(n)])
@lru_cache(None)
def cores(D):return tuple(profile(rooted(D,x)) for x in range(1<<len(D)))
def graph_fuse(A,B):
 n=len(A);m=len(B);T=[0]*(n+m-1)
 for Z,ids in ((A,list(range(n))),(B,[0]+list(range(n,n+m-1)))):
  for i in range(len(Z)):
   for j in range(len(Z)):
    if Z[i]>>j&1:T[ids[i]]|=1<<ids[j]
 for i in range(1,n):
  for j in range(n,n+m-1):T[i]|=1<<j
 return tuple(T)
def unit_graph(p,n):return rooted(tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n)),p)
def witness(D,w):return graph_fuse(graph_fuse(unit_graph(w['p'],w['a']),rooted(D,w['cut'])),unit_graph(w['q'],w['b']))
def fiber(D,V):
 C=cores(tuple(D));h=C[0]['h'];Q=C[0]['Q'];require(h>0 and h%2==1,'odd positive h')
 if V<h:return []
 z=V-h;a=p=0;ans=[]
 while z>=0:
  pre=unit(p,a)
  for x,c in enumerate(C):
   t=z-(c['V']-h)
   if t<0 or t%Q:continue
   q=t//Q;b=q.bit_length();out=fuse(fuse(pre,c),unit(q,b));require(out['h']==h and out['V']==V,'fiber invariant')
   ans.append(dict(a=a,p=p,cut=x,q=q,b=b,profile=out))
  if z==0:break
  bit=z&1;p+=bit<<a;z=(z-h*bit)//2;a+=1
 return ans
def triangle_deltas(r):
 require(r>=0,'receiver code');out={0}
 for a in range(r.bit_length()+1):
  k=(r>>a)&3;t=r>>(a+2);p=((r>>(a-1))&1) if a else 1;f=3*t+2*p-(t&1)
  if k==1:out.add(f)
  elif k==2:out.add(-f-2)
 return out
def ordinal(*parts):
 out=[];off=0;n=sum(map(len,parts))
 for D in parts:
  tail=((1<<n)-1)^((1<<(off+len(D)))-1)
  out.extend((r<<off)|tail for r in D);off+=len(D)
 return tuple(out)
TRI=(2,4,1)
S4=(2,4,9,3)
A5=(6,12,24,17,3)
E6=(44,1,2,6,15,30)
BC=ordinal(TRI,S4)
# Fixed literal representatives discovered by the inverse; no random search needed.
EIGHT=[(54365,A5,8,57,7,4),(54367,A5,8,57,11,4),(54369,BC,10,313,96,0),(54371,A5,2,1,3,291),(54373,A5,2,1,5,291),(54375,E6,4,9,57,40),(54377,A5,0,0,7,1165),(54379,A5,0,0,11,1165)]
def eight_witnesses():
 out=[]
 for H,D,a,p,x,q in EIGHT:
  w=dict(a=a,p=p,cut=x,q=q,b=q.bit_length());pr=fuse(fuse(unit(p,a),cores(D)[x]),unit(q));require((pr['H'],pr['h'],pr['V'])==(H,15,139878),'eight exact packet');w['profile']=pr
  out.append((D,w))
 return out
