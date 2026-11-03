"""Controls for unconditional first-step identity; no bare assertions."""
from itertools import combinations

def need(p,m):
 if not p:raise AssertionError(m)
def H(a,start=None):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):
  if start is None or start==i:d[1<<i][i]=1
 for m in range(1,1<<n):
  for i in range(n):
   c=d[m][i]
   if c:
    z=a[i]&~m
    while z:
     v=z&-z;j=v.bit_length()-1;d[m|v][j]+=c;z-=v
 return sum(d[-1])
def strong(a):
 n=len(a)
 for root in range(n):
  seen=1<<root;old=0
  while seen!=old:
   old=seen
   for i in range(n):
    if seen>>i&1:seen|=a[i]
  if seen!=(1<<n)-1:return False
 return True
def query(a,p):
 n=len(a);return tuple(a[i]|((1<<n) if i!=p else 0) for i in range(n))+(1<<p,)
def sub(a,q):
 ids=[i for i in range(len(a)) if i!=q];return tuple(sum(((a[i]>>j)&1)<<k for k,j in enumerate(ids)) for i in ids),ids
def g(a,q,r):
 b,ids=sub(a,q);return H(b,ids.index(r))
def ear(a,q):
 n=len(a);return tuple(a[i]|((1<<n) if i==q else 0) for i in range(n))+((((1<<n)-1)^(1<<q))|(1<<(n+1)),(1<<n)-1)
def obs(a,p):
 h=H(a);a=query(a,p);ha=H(a);a=query(a,len(a)-1);haa=H(a)
 return h,(ha-h)//2,(haa-ha)//2
def gs(a,q,r,p):
 out=[]
 for i in range(3):
  out.append(g(a,q,r))
  if i<2:
   a=query(a,p);p=len(a)-1
 return tuple(out)
def check(a,q,r,p):
 G=gs(a,q,r,p);v=obs(a,p);b=ear(a,q);w=obs(b,p)
 target=(3*v[0]+2*G[0],3*v[1]+G[1]-G[0],3*v[2]+G[2]-G[1])
 need(w==target,'observer correction')
 need(G[0]<=G[1]<=G[2],'query injection')
 need(gs(b,len(a),q,p)==(0,0,0),'automatic zero guard for all queries')
 if strong(a):need(strong(b),'strongness')
 return v,G,b,w
for a,q,r,p,expect in [((8,1,3,6),2,3,3,((5,3,4),(1,1,1),(17,9,12))),((20,1,2,7,14),3,4,4,((9,10,11),(2,2,2),(31,30,33)))]:
 v,G,b,w=check(a,q,r,p);need((v,G,w)==expect,'named example');need(strong(a) and strong(b),'named examples strong')
 guards=[]
 for q2 in range(len(b)):
  pred=[i for i in range(len(b)) if b[i]>>q2&1]
  if len(pred)==1:guards.append((q2,pred[0],g(b,q2,pred[0])))
 need(guards==[(len(a),q,0),(len(a)+1,len(a),3 if G[0]==1 else 6)],'next guard census')
 need(guards[1][2]==H(a)-H(a,q),'alternative guard formula')
 print('example',v,'G=',G,'ear',b,'observer',w,'next guards',guards)
# Exhaust all labeled tournaments of orders 3 and 4, all valid guards and all p!=q.
tests=0
for n in (3,4):
 e=list(combinations(range(n),2))
 for bits in range(1<<len(e)):
  a=[0]*n
  for k,(i,j) in enumerate(e):
   if bits>>k&1:a[i]|=1<<j
   else:a[j]|=1<<i
  for q in range(n):
   r=[i for i in range(n) if a[i]>>q&1]
   if len(r)!=1:continue
   for p in range(n):
    if p!=q:check(tuple(a),q,r[0],p);tests+=1
print('all labeled tournaments n=3,4: marked guard cases',tests)
# Sharp size boundary of monotonicity: X−q singleton, p=r.
a=(0,1);need(gs(a,0,1,1)==(1,0,1),'order-two edge case')
print('order-two edge case G=(1,0,1); PASS')
