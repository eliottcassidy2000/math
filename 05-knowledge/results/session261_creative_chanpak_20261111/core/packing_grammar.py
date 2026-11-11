#!/usr/bin/env python3
"""Exact control for the core-resource theorem, not a Chan-Pak census."""
import itertools,json,random,pathlib

def tournament(core,word):
 k=len(core); out=list(core)+[0]*len(word)
 for i,a in enumerate(word):
  v=k+i
  for c in range(k):
   if a>>c&1:out[c]|=1<<v
   else:out[v]|=1<<c
  for j in range(i+1,len(word)):out[v]|=1<<(k+j)
 return out

def triangles(t):
 out=[]
 for a,b,c in itertools.combinations(range(len(t)),3):
  if ((t[a]>>b&1) and (t[b]>>c&1) and(t[c]>>a&1)) or ((t[b]>>a&1) and(t[c]>>b&1) and(t[a]>>c&1)):
   out.append((1<<a)|(1<<b)|(1<<c))
 return out

def matching(edges):
 dp={0}
 for e in set(edges):
  dp|={u|e for u in tuple(dp) if not(u&e)}
 # Hyperedges have varying sizes. Need matched edge count, not union bitcount.
 best={0:0}
 for e in set(edges):
  for u,v in tuple(best.items()):
   if not(u&e):best[u|e]=max(best.get(u|e,0),v+1)
 return max(best.values(),default=0)

def resources(core,A):
 k=len(core); ed={m:[] for m in triangles(core)}
 for c in range(k):
  lo=[a for a in A if a>>c&1];hi=[a for a in A if not(a>>c&1)]
  if lo and hi:ed[1<<c]=[lo[0],hi[0]]
 for c,d in itertools.combinations(range(k),2):
  for a in A:
   t=tournament(core,[a])
   if ((1<<c)|(1<<d)|(1<<k)) in triangles(t):
    ed.setdefault((1<<c)|(1<<d),[a]);break
 return ed

def check(core,A):
 ed=resources(core,A); allsets=[(0,[])]
 for m,w in ed.items():
  allsets += [(used|m,ws+[w]) for used,ws in tuple(allsets) if not(used&m)]
 used,parts=max(allsets,key=lambda p:len(p[1])); w=sum(parts,[])
 got=matching(triangles(tournament(core,w)))
 if got!=len(parts):raise RuntimeError(('attainment',core,A,ed,w,got,len(parts)))
 # A fixed word's packing cannot exceed the core-resource matching bound.
 for n in range(4):
  for w in itertools.product(A,repeat=n):
   actual=matching(triangles(tournament(core,w)))
   if actual>len(parts):raise RuntimeError(('upper',core,A,w))
 return len(parts)

def all_tournaments(k):
 pairs=list(itertools.combinations(range(k),2))
 for mask in range(1<<len(pairs)):
  t=[0]*k
  for i,(a,b) in enumerate(pairs):
   u,v=(a,b) if mask>>i&1 else(b,a);t[u]|=1<<v
  yield t

def main():
 checked=0; distribution={}
 for k in range(1,4):
  for core in all_tournaments(k):
   for mask in range(1,1<<(1<<k)):
    A=[a for a in range(1<<k) if mask>>a&1]
    r=check(core,A);checked+=1;distribution[str((k,r))]=distribution.get(str((k,r)),0)+1
 rng=random.Random(261112)
 for _ in range(50):
  k=4;core=list(all_tournaments(4))[rng.randrange(64)];A=rng.sample(range(16),rng.randrange(1,6));check(core,A);checked+=1
 edges=[(0,1),(0,2),(0,4),(1,5),(2,3),(3,7),(4,6),(5,7),(6,7)]
 core=[58,60,57,16,32,8]; safe=set(); full=[]
 for a,b in itertools.product(edges,repeat=2):
  A=frozenset(x+(y<<3) for x in a for y in b)
  for n in range(1,5):
   for B in itertools.combinations(A,n):
    if matching(resources(core,B))<=2:safe.add(frozenset(B))
 maximal=[A for A in safe if not any(A<B for B in safe)]
 squares=[A for A in maximal if len(A)==4]; small=sorted([sorted(A) for A in maximal if len(A)<4])
 if len(squares)!=45 or small!=[[8,40],[16,24],[32,48],[57,61],[58,59],[60,62]]:raise RuntimeError(('tensor',len(squares),small))
 print(json.dumps({'tensor':{'safe_subalphabets':len(safe),'maximal_squares':len(squares),'maximal_edges':small},'status':'PASS' ,'core_alphabet_pairs':checked,'complete':'all core orientations k<=3 and nonempty alphabets; every word length<=3; explicit maximizing word for every pair','random':'50 k4 pairs, seed261112','distribution':distribution},indent=2))
if __name__=='__main__':main()
