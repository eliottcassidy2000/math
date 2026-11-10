#!/usr/bin/env python3
"""Exact scoped verification, no optimization-removable gates."""
import itertools,json,math,sys

def need(v,msg):
 if not v: raise RuntimeError(msg)

def table(rows):
 n=len(rows);top=(1<<n)-1
 need(all(not(r>>i&1) and r<=top for i,r in enumerate(rows)), 'row shape')
 need(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(n)for j in range(i+1,n)), 'tournament')
 dp=[[0]*n for _ in range(top+1)];h=[0]*(top+1);h[0]=1
 for m in range(1,top+1):
  for v in range(n):
   if m>>v&1:
    t=m^(1<<v);dp[m][v]=1 if not t else sum(dp[t][w] for w in range(n) if rows[w]>>v&1)
  h[m]=sum(dp[m])
 return h,dp

def state(rows,u):
 h,dp=table(rows);n=len(rows);top=(1<<n)-1
 sp=[[0]*n for _ in range(top+1)];sp[1<<u][u]=1
 for m in range(1,top+1):
  if not m>>u&1:continue
  for v in range(n):
   if m>>v&1 and v!=u:
    t=m^(1<<v);sp[m][v]=sum(sp[t][w] for w in range(n) if rows[w]>>v&1)
 return (h[-1],sum(sum(sp[m])*h[top^m] for m in range(top+1)),dp[-1][u],h[top^(1<<u)])

def brute(rows,u):
 n=len(rows);top=(1<<n)-1;h=[0]*(top+1);h[0]=1;s=[0]*(top+1);e=0
 for k in range(1,n+1):
  for p in itertools.permutations(range(n),k):
   if all(rows[a]>>b&1 for a,b in zip(p,p[1:])):
    m=sum(1<<v for v in p);h[m]+=1
    if p[0]==u:s[m]+=1
    if k==n and p[-1]==u:e+=1
 return h,(h[-1],sum(s[m]*h[top^m] for m in range(top+1)),e,h[top^(1<<u)])

def strong(rows):
 top=(1<<len(rows))-1
 for u in range(len(rows)):
  seen=1<<u
  while True:
   nxt=seen
   for v,r in enumerate(rows):
    if seen>>v&1:nxt|=r
   if nxt==seen:break
   seen=nxt
  if seen!=top:return False
 return True

def graft(a,u,b):
 n=len(a);m=len(b);rows=list(a)+[r<<n for r in b]
 for i in range(n):
  for j in range(m):
   if i==u:rows[n+j]|=1<<i
   else:rows[i]|=1<<(n+j)
 return rows

def sink(a):
 n=len(a);return [r|(1<<n) for r in a]+[0]

def HQ(a):
 h,_=table(a);t=len(h)-1;return h[-1],sum(x*h[t^m] for m,x in enumerate(h))

G=[22,60,8,81,4,29,55]
Z=[126,180,120,2,168,72,154,45]
C=[2,4,1]
need(HQ(G)==(49,812),'gadget')
hb,zb=brute(G,0);need(hb==table(G)[0] and zb==state(G,0),'independent gadget path enumeration')
hb,zb=brute(Z,3);need(hb==table(Z)[0] and zb==state(Z,3),'independent root path enumeration')
need(zb==(259,266,70,77),'root state')
need(strong(G) and strong(Z),'strong seed checks')
need(HQ(Z)==(259,4456),'rooted versus gadget content distinction')
W=graft(C,0,G);need(HQ(W)[0]==1673 and strong(W),'strong content-seven graft')
V=graft(sink(Z),3,[0]);need(HQ(V)[0]==791 and strong(V),'strong SB')
need(all(239%d for d in range(2,math.isqrt(239)+1)) and all(113%d for d in range(2,math.isqrt(113)+1)),'prime checks')
# Exhaust all labelled tournaments through five vertices; compare the exact
# packing-cover sum against the cut injection and complement-parity bound.
census=[]
for n in range(1,6):
 es=list(itertools.combinations(range(n),2));count=0
 for bits in range(1<<len(es)):
  a=[0]*n
  for k,(i,j) in enumerate(es):
   x,y=(i,j) if bits>>k&1 else(j,i);a[x]|=1<<y
  h,q=HQ(a);need(q>=(n+1)*h,'cut lower bound')
  if n>=2:need(q%4==0,'four divisibility')
  count+=1
 census.append([n,count])
# Independent integer obstruction: EVERY relaxed physical source with a
# possible normalized quotient seven is considered, not only listed graphs.
# q>=2h, 4|q, h odd, h+q<=7 forces h=1,q=4.
pairs=[(h,q)for h in range(1,8,2)for q in range(4,8,4)if q>=2*h and h+q<=7]
need(pairs==[(1,4)],'forced pair')
# d+4(D-E)+3f=7 implies D<=3. Include all feasible d,E and f<=2.
hostile=[]
for D in [1,3]:
 for d in range(1,D+1,2):
  for E in range(min(d,(D-1)//2)+1):
   for f in range(3):
    # U=D+d-E+f and U=E mod2 force f even.
    if f%2:continue
    if d+4*(D-E)+3*f==7:hostile.append([D,E,d,f])
need(not hostile,'quotient seven exclusion')
if '--fail' in sys.argv:need(False,'deliberate negative control')
print(json.dumps({'status':'PASS','labelled_HQ_census':census,'gadget':{'rows':G,'H':49,'Q':812,'gcd':7,'normalized_pair':[7,116]},'source_C3_state_h_R_e_d':state(C,0),'graft':{'vertices':len(W),'H':1673,'prime_cofactor':239,'strong':True},'root_content_seven':{'rows':Z,'root':3,'state':zb,'content':7},'SB':{'vertices':len(V),'H':791,'prime_cofactor':113,'strong':True},'quotient_seven_forced_pairs':pairs,'relaxed_source_solutions':hostile},indent=2))
