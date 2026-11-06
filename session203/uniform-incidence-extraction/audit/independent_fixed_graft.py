#!/usr/bin/env python3
from random import Random
from functools import lru_cache
import json
rng=Random(203)
def need(c,m):
 if not c:raise RuntimeError(m)
def data(a,root=None):
 n=len(a);full=(1<<n)-1;inc=[sum(((a[j]>>i)&1)<<j for j in range(n)) for i in range(n)]
 @lru_cache(None)
 def end(mask,v):
  if mask==(1<<v):return 1
  old=mask^(1<<v);q=old&inc[v];s=0
  while q:
   b=q&-q;s+=end(old,b.bit_length()-1);q-=b
  return s
 @lru_cache(None)
 def start(mask,v):
  if mask==(1<<v):return 1
  old=mask^(1<<v);q=old&a[v];s=0
  while q:
   b=q&-q;s+=start(old,b.bit_length()-1);q-=b
  return s
 def H(mask):
  return sum(end(mask,v) for v in range(n) if mask>>v&1) if mask else 1
 out={'h':H(full)}
 if root is not None:
  out.update(d=H(full^(1<<root)),e=end(full,root),s=start(full,root),R=sum(start(s,root)*H(full^s) for s in range(1,full+1) if s>>root&1),L=sum(end(s,root)*H(full^s) for s in range(1,full+1) if s>>root&1))
 return out
def induced(a,ids):return [sum(((a[v]>>w)&1)<<k for k,w in enumerate(ids)) for v in ids]
def make(u,B,S,bc,sc):
 bits=list(map(int,bin(u-1)[2:]));m=len(bits);nb=len(B);ns=len(S);n=1+nb+m+ns;a=[0]*n
 # c=0, then B, private chain, S (r is first S vertex)
 for i in range(nb):
  for j in range(nb):
   if B[i]>>j&1:a[1+i]|=1<<(1+j)
  for j in range(1+nb,n):a[1+i]|=1<<j
 for i in range(m):
  for j in range(i+1,m):a[1+nb+i]|=1<<(1+nb+j)
  for j in range(ns):a[1+nb+i]|=1<<(1+nb+m+j)
 for i in range(ns):
  for j in range(ns):
   if S[i]>>j&1:a[1+nb+m+i]|=1<<(1+nb+m+j)
 cwin=list(bc)+bits+list(sc)
 for v,win in enumerate(cwin,1):
  if win:a[0]|=1<<v
  else:a[v]|=1
 return a,bits
def randgraph(n):
 a=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   if rng.randrange(2):a[i]|=1<<j
   else:a[j]|=1<<i
 return a
rows=[]
for trial in range(45):
 nb=trial%4;ns=1+(trial//4)%4;B=randgraph(nb);S=randgraph(ns)
 bc=[rng.randrange(2) for _ in B];sc=[0]+[rng.randrange(2) for _ in S[1:]] # old r->c
 k=2+trial%2;vals=[]
 for u in range(2**k+1,2**k+5):
  a,bits=make(u,B,S,bc,sc);m=len(bits);nb=len(B);ns=len(S)
  L=data(induced(a,[0]+list(range(1,nb+1))),0)
  R=data(induced(a,[0]+list(range(1+nb+m,1+nb+m+ns))),0)
  W=data(induced(a,[0]+list(range(1+nb,1+nb+m))),0)
  x=L['R']-L['d'];y=R['L']-R['d'];ell=1+sum((1-b)*2**i for i,b in enumerate(bits))
  need(W['s']==1 and W['e']==u%2 and W['R']==u and W['L']==ell,'unit data')
  K=R['d']*L['h']-R['d']*x-L['d']*R['d']-L['d']*y+L['d']*(R['h']-R['s'])
  expected=L['d']*R['d']*W['h']+L['d']*y*u+R['d']*x*ell+x*y*2**m+L['d']*(R['s']-R['d'])*(u%2)+K
  actual=data(a)['h'];need(expected==actual,'general physical graft formula');vals.append(actual)
 need(vals[3]-vals[2]-vals[1]+vals[0]==-2*data(B)['h']*data(S)['h'],'mixed difference')
 rows.append({'trial':trial,'B_order':nb,'S_order':ns,'k':k,'values':vals,'mixed_difference':vals[3]-vals[2]-vals[1]+vals[0]})
print(json.dumps({'status':'PASS','context_count':len(rows),'graphs_checked':4*len(rows),'largest_graph_order':12,'scope':'independent direct subset recursion; finite controls for separately proved all-order formula','rows':rows},indent=2,sort_keys=True))
