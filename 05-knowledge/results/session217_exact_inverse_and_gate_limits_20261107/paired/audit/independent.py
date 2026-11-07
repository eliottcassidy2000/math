#!/usr/bin/env python3
import itertools,json

def tourn(n,code):
 r=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   x,y=(i,j) if code>>k&1 else(j,i);r[x]|=1<<y;k+=1
 return tuple(r)
def add(r,c):
 n=len(r);out=list(r)+[0]
 for i in range(n):
  newbeats=(i!=0) if c=='Q' else(i==0)
  if newbeats:out[n]|=1<<i
  else:out[i]|=1<<n
 return tuple(out)
def rename(r,order):return tuple(sum(1<<j for j,v in enumerate(order) if r[u]>>v&1) for u in order)
def normal(r,q,z):
 n=len(r);labels=[('r',0)]+[('q',i) for i in range(q)]+[('d',i) for i in range(1,n)]+[('z',i) for i in range(z)]
 out=[0]*len(labels);ranks={'q':0,'d':1,'z':2}
 for a,(ta,ia) in enumerate(labels):
  for b,(tb,ib) in enumerate(labels):
   if a==b:continue
   if ta=='r':win=tb=='q' or(tb=='d' and bool(r[0]>>ib&1))
   elif tb=='r':win=ta=='z' or(ta=='d' and bool(r[ia]&1))
   elif ta==tb:win=(bool(r[ia]>>ib&1) if ta=='d' else ia<ib)
   else:win=ranks[ta]<ranks[tb]
   if win:out[a]|=1<<b
 return tuple(out)
def profile(r):
 n=len(r);M=(1<<n)-1;H=[0]*(M+1);S=[0]*(M+1);E=[0]*(M+1);H[0]=1
 for mask in range(1,M+1):
  ids=[i for i in range(n) if mask>>i&1]
  for p in itertools.permutations(ids):
   if all(r[a]>>b&1 for a,b in zip(p,p[1:])):
    H[mask]+=1;S[mask]+=p[0]==0;E[mask]+=p[-1]==0
 L=sum(E[m]*H[M^m] for m in range(M+1));R=sum(S[m]*H[M^m] for m in range(M+1));D=M^1
 return H[M],S[M],E[M],L,R,H[D],sum(H[m]*H[D^m] for m in range(M+1) if not m&1)
def f(p,c):
 h,s,e,L,R,d,b=p
 return (h+L-s,d,e,2*L-d,R+b,d,2*b) if c=='Q' else(h+R-e,s,d,L+b,2*R-d,d,2*b)
graphchecks=literal=stablechecks=axes=0;degenerate=0
for n in range(1,5):
 for code in range(1<<(n*(n-1)//2)):
  r=tourn(n,code)
  for length in range(7):
   for word in itertools.product('QZ',repeat=length):
    g=r
    for c in word:g=add(g,c)
    qids=[n+i for i,c in enumerate(word) if c=='Q'][::-1];zids=[n+i for i,c in enumerate(word) if c=='Z']
    if rename(g,[0]+qids+list(range(1,n))+zids)!=normal(r,len(qids),len(zids)):raise RuntimeError('graph normal form')
    graphchecks+=1
  if n>3:continue
  p=profile(r)
  if n==1:
   if profile(add(r,'Q'))==f(p,'Q'):raise RuntimeError('singleton hostile did not fail')
   continue
  for c in 'QZ':
   if profile(add(r,c))!=f(p,c):raise RuntimeError('transfer')
   literal+=1
  sr=add(add(r,'Q'),'Z');ps=profile(sr);h,s,e,L,R,d,b=ps;a=L-d;rho=R-d
  if s!=e or e!=d or a<d or rho<d or h<d:raise RuntimeError('stable bounds')
  for q in range(5):
   for z in range(5):
    u=2**q;v=2**z
    pred=(h+(u-1)*a+(v-1)*rho+(u-1)*(v-1)*b,d,d,u*(a+(v-1)*b)+d,v*(rho+(u-1)*b)+d,d,u*v*b)
    cur=ps
    for c in 'Q'*q+'Z'*z:cur=f(cur,c)
    if pred!=cur or cur[0]<d*max(u,v):raise RuntimeError('stable formula or bound')
    for j in range(9):
     x=cur[0]+cur[3]-d+2*j
     xp=2*u*v*b+2*u*(a-b)+v*(rho-b)+h-a-rho+b+2*j
     if x!=xp:raise RuntimeError('compiler base formula')
    if q+z<=2:
     g=sr
     for c in 'Z'*z+'Q'*q:g=add(g,c)
     if profile(g)!=pred:raise RuntimeError('literal closed form')
     literal+=1
    stablechecks+=1
  h,s,e,L,R,d,b=p
  for k in range(1,7):
   for c in 'QZ':
    cur=p
    for _ in range(k):cur=f(cur,c)
    if c=='Q':hh=h+L-s+(2**k-2)*(L-d);ll=d+2**k*(L-d);ss=d
    else:hh=h+R-e+(2**k-2)*(R-d);ll=L+(2**k-1)*b;ss=s
    if cur[0]!=hh or cur[3]!=ll or cur[1]!=ss:raise RuntimeError('axis')
    axes+=1
  degenerate+=L==d
print(json.dumps(dict(status='PASS',canonical_graph_checks=graphchecks,literal_permutation_profile_checks=literal,stable_exponent_pairs=stablechecks,axis_checks=axes,degenerate_Q_sources=degenerate,singleton_transfer_hostile='fails as expected; graph commutation still holds'),indent=2))
