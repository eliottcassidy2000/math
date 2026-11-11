from fractions import Fraction
from itertools import combinations
from explore import H,fans,tt,flip

def packet(a,r):
 n=len(a); ids=[u for u in range(n) if u!=r]; k=len(ids)
 b=[sum(((a[u]>>v)&1)<<j for j,v in enumerate(ids)) for u in ids]
 A=sum(((a[r]>>u)&1)<<j for j,u in enumerate(ids));allmask=(1<<k)-1
 h=[0]*(1<<k);h[0]=1
 st=[[0]*k for _ in h];en=[[0]*k for _ in h]
 for mask in range(1,1<<k):
  for v in range(k):
   if mask>>v&1:
    rest=mask^(1<<v)
    st[mask][v]=1 if not rest else sum(st[rest][u] for u in range(k) if b[v]>>u&1)
    en[mask][v]=1 if not rest else sum(en[rest][u] for u in range(k) if b[u]>>v&1)
  h[mask]=sum(st[mask])
 q=sum(h[s]*h[allmask^s] for s in range(1<<k))
 R=sum(sum(st[s][u] for u in range(k) if A>>u&1)*h[allmask^s] for s in range(1<<k))
 L=sum(sum(en[s][u] for u in range(k) if A>>u&1)*h[allmask^s] for s in range(1<<k))
 s=sum(st[allmask][u] for u in range(k) if A>>u&1)
 e=sum(en[allmask][u] for u in range(k) if A>>u&1)
 F=H(a);hh=h[allmask]
 beta=-L-hh;gamma=R-q;delta=F+L+2*hh-s-R+e
 return q,beta,gamma,delta

def predict(p,m,n):
 q,b,c,d=p;return q*2**(m+n)+b*2**m+c*2**n+d

def curvature(p):
 q,b,c,d=p;return Fraction(q*d-b*c,q)
if __name__=='__main__':
 for n in range(1,6):
  edges=list(combinations(range(n),2));vs=[];menus={}; seen={};collision=None
  for code in range(1<<len(edges)):
   a=tt(n)
   for k,(i,j) in enumerate(edges):
    if code>>k&1:a=flip(a,i,j)
   for r in range(n):
    p=packet(a,r);v=curvature(p);vs.append(v)
    # single-fan H at k=0,1,2 completely determines scalar continuation
    sig=(H(a),H(fans(a,r,0,1)),H(fans(a,r,0,2)))
    if sig in seen and seen[sig][2]!=p and collision is None:
     aa,rr,pp=seen[sig]
     if predict(pp,1,1)!=predict(p,1,1):collision=(sig,aa,rr,pp,a,r,p)
    seen[sig]=(a,r,p)
    if n<=3:
     for m in range(1,4):
      for nn in range(1,4):assert H(fans(a,r,m,nn))==predict(p,m,nn)
  print(n,'curvature',min(vs),max(vs),'distinct',len(set(vs)),'collision',collision)
