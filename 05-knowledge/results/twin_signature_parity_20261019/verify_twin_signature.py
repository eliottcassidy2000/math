#!/usr/bin/env python3
"""Exact controls; run normally and with python3 -O."""
from functools import lru_cache
from math import comb
from collections import Counter

def adj_of(edges,n):
 a=[[] for _ in range(n)]
 for u,v in edges:a[u].append(v);a[v].append(u)
 return a

def centers(a):
 n=len(a)
 if n<3:return list(range(n))
 deg=[len(x) for x in a];leaves=[i for i,x in enumerate(deg) if x==1];remain=n
 while remain>2:
  remain-=len(leaves);new=[]
  for v in leaves:
   for w in a[v]:
    deg[w]-=1
    if deg[w]==1:new.append(w)
  leaves=new
 return leaves

def root(a,v,p=-1):return tuple(sorted(root(a,w,v) for w in a[v] if w!=p))
def canon(a):
 c=centers(a)
 return ('V',root(a,c[0])) if len(c)==1 else ('E',tuple(sorted((root(a,c[0],c[1]),root(a,c[1],c[0])))))

@lru_cache(None)
def rooted_aut(t):
 z=1
 for s,m in Counter(t).items():
  from math import factorial
  z*=rooted_aut(s)**m*factorial(m)
 return z

def aut(a):
 c=centers(a)
 if len(c)==1:return rooted_aut(root(a,c[0]))
 x,y=root(a,c[0],c[1]),root(a,c[1],c[0])
 return rooted_aut(x)*rooted_aut(y)*(2 if x==y else 1)

@lru_cache(None)
def ar(t):
 z=1
 for s,m in Counter(t).items():z*=comb(2*ar(s),m) if m<=2*ar(s) else 0
 return z

def afree(a):
 c=centers(a)
 if len(c)==1:return ar(root(a,c[0]))
 x,y=root(a,c[0],c[1]),root(a,c[1],c[0])
 return ar(x)*ar(y)*(1 if x==y else 2)

def hereditary(a):
 n=len(a)
 for mask in range(1,1<<n):
  verts=[i for i in range(n) if mask>>i&1]
  if len(verts)<2:continue
  edges=[(i,j) for i in verts for j in a[i] if i<j and mask>>j&1]
  if len(edges)!=len(verts)-1:continue
  d={v:i for i,v in enumerate(verts)}
  b=adj_of([(d[i],d[j]) for i,j in edges],len(verts))
  if aut(b)==1:return False
 return True

def trees(maxn):
 out={1:[()]}
 for n in range(2,maxn+1):
  cur={}
  for old in out[n-1]:
   for v in range(n-1):
    e=old+((v,n-1),);a=adj_of(e,n);cur.setdefault(canon(a),e)
  out[n]=list(cur.values())
 return out

def orient_data(edges,n,mask):
 a=adj_of(edges,n);out=[0]*n
 for k,(u,v) in enumerate(edges):
  if mask>>k&1:out[v]|=1<<u
  else:out[u]|=1<<v
 return a,out

def oroot(a,out,v,p=-1):return tuple(sorted((int(bool(out[v]>>w&1)),oroot(a,out,w,v)) for w in a[v] if w!=p))
@lru_cache(None)
def oroot_aut(t):
 from math import factorial
 z=1
 for (sgn,s),m in Counter(t).items():z*=oroot_aut(s)**m*factorial(m)
 return z

def ocanon_aut(a,out):
 c=centers(a)
 if len(c)==1:
  t=oroot(a,out,c[0]);return ('V',t),oroot_aut(t)
 u,v=c
 if not (out[u]>>v&1):u,v=v,u
 x,y=oroot(a,out,u,v),oroot(a,out,v,u)
 return ('E',x,y),oroot_aut(x)*oroot_aut(y)

def embedding_count(edges,n,mask,host):
 a,ori=orient_data(edges,n,mask)
 rootv=max(range(n),key=lambda v:len(a[v]))
 order=[rootv];parent={rootv:-1}
 for v in order:
  for w in a[v]:
   if w!=parent[v]:parent[w]=v;order.append(w)
 allbits=(1<<n)-1;pos=[-1]*n
 def go(i,unused):
  if i==n:return 1
  v=order[i];p=parent[v]
  if p<0:cands=unused
  elif ori[p]>>v&1:cands=unused&host[pos[p]]
  else:cands=unused&(~host[pos[p]])&~(1<<pos[p])
  if i==n-1:return cands.bit_count()
  total=0
  while cands:
   bit=cands&-cands;cands-=bit;pos[v]=bit.bit_length()-1
   total+=go(i+1,unused^bit)
  return total
 return go(0,allbits)

def transitive(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]
def random_host(n,rng):
 host=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   if rng.randrange(2):host[i]|=1<<j
   else:host[j]|=1<<i
 return host

def require(condition,message):
 if not condition:raise RuntimeError(message)

def brute_emb(edges,n,mask,host):
 from itertools import permutations
 arcs=[(v,u) if mask>>k&1 else (u,v) for k,(u,v) in enumerate(edges)]
 return sum(all(host[p[u]]>>p[v]&1 for u,v in arcs) for p in permutations(range(n)))

def v2fact(n):
 z=0
 while n:n//=2;z+=n
 return z

def v_budget(m,k):
 total=0;power=2
 while power<=2*(m+k):
  total+=max(0,(m+k+power-1)//power-k);power*=2
 return total

def main():
 import json,pathlib
 ts=trees(11)
 expected=[(1,1,0,0),(1,1,0,0),(1,1,0,0),(2,2,1,0),(3,3,1,0),(6,6,2,0),(11,10,4,0),(23,16,9,0),(47,24,19,3),(106,35,46,18),(235,49,106,64)]
 rows=[];orientation_checks=0
 for n,es in ts.items():
  hs=zs=bad=0
  for e in es:
   a=adj_of(e,n);h=hereditary(a);z=afree(a)==0
   hs+=h;zs+=z;bad+=z and not h
   if n<=9:
    observed=set()
    for mask in range(1<<max(0,n-1)):
     code,group=ocanon_aut(*orient_data(e,n,mask));orientation_checks+=1
     if group==1:observed.add(code)
    require(len(observed)==afree(a),f'capacity at n={n}: {e}')
  require((len(es),hs,zs,bad)==expected[n-1],f'census mismatch n={n}')
  rows.append({'n':n,'unlabelled_trees':len(es),'hereditarily_symmetric':hs,'automorphism_forced_even':zs,'nonhereditary_automorphism_forced_even':bad})
 # The actual nine-vertex witness and the retained asymmetric subtree.
 e=((0,1),(0,2),(0,3),(0,4),(0,5),(1,6),(2,7),(6,8));a=adj_of(e,9)
 require(aut(a)==6 and afree(a)==0 and not hereditary(a),'nine-vertex witness')
 for mask in range(256):
  _,group=ocanon_aut(*orient_data(e,9,mask))
  require(group%2==0,'oriented witness unexpectedly asymmetric')
 # Certificates cover every nonhereditary tree on at most eight vertices.
 data=json.loads(pathlib.Path(__file__).with_name('minimality_witnesses.json').read_text())
 checked=0
 for row in data['rows']:
  n=row['n'];need={canon(adj_of(e,n)) for e in ts[n] if not hereditary(adj_of(e,n))};have=set()
  require(not row['unresolved'],'unresolved minimality entry')
  for w in row['witnesses']:
   e=tuple(map(tuple,w['edges']));mask=w['orientation_mask'];h=w['host_out_bits']
   require(all(not(h[i]>>i&1) for i in range(n)),'loop in host')
   require(all(((h[i]>>j&1)+(h[j]>>i&1))==1 for i in range(n) for j in range(i+1,n)),'host not tournament')
   b=embedding_count(e,n,mask,transitive(n));c=embedding_count(e,n,mask,h)
   require(b==w['transitive_count'] and c==w['host_count'],'certificate count mismatch')
   require(b==brute_emb(e,n,mask,transitive(n)) and c==brute_emb(e,n,mask,h),'independent permutation count mismatch')
   require((b-c)%2==1,'certificate not nonconstant');have.add(canon(adj_of(e,n)));checked+=1
  require(need==have,f'minimality does not cover n={n}')
 require(checked==8,'wrong number of nonconstancy certificates')
 # P5 with only its first edge reversed: asymmetric but universally even on five-tournaments.
 p5=tuple((i,i+1) for i in range(4));mask=1
 a,out=orient_data(p5,5,mask)
 require(ocanon_aut(a,out)[1]==1,'P5 control symmetric')
 for bits in range(1<<10):
  host=[0]*5;j=0
  for u in range(5):
   for v in range(u+1,5):
    if bits>>j&1:host[u]|=1<<v
    else:host[v]|=1<<u
    j+=1
  require(embedding_count(p5,5,mask,host)%2==0,'P5 tournament parity failure')
 require(brute_emb(p5,5,mask,out)==1,'all-host P5 converse control')
 # Dynamic optimization independent of closed factorial budget.
 factorial_budget_cases=0
 cost=[v2fact(i) for i in range(65)]
 prev=[0]+[10**9]*64
 for k in range(1,9):
  cur=[min(prev[m-i]+cost[i] for i in range(m+1)) for m in range(65)]
  for m in range(65):
   require(cur[m]==v_budget(m,k),f'factorial budget k={k},m={m}');factorial_budget_cases+=1
  prev=cur
 result={'status':'PASS','scope':'Exact capacity/filter census, not a classification of good trees beyond eight vertices.','orientation_masks_checked':orientation_checks,'minimality_nonconstancy_pairs':checked,'P5_tournaments_checked':1024,'factorial_budget_cases':factorial_budget_cases,'census':rows}
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
