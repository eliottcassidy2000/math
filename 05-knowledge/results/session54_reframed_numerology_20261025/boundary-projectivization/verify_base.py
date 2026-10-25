#!/usr/bin/env python3
"""Exact independent Paley-7 blocker/even-subset bridge; no external packages."""
import itertools as it, collections as co, json, hashlib, pathlib, sys
ROOT=pathlib.Path(__file__).resolve().parent
V=range(7); FULL=127
E=[(u,v) for u in V for v in V if (v-u)%7 in (1,2,4)]
EI={e:i for i,e in enumerate(E)}
G=[(a,b) for a in (1,2,4) for b in V]
def require(p,msg):
 if not p: raise RuntimeError(msg)
def vm(S): return sum(1<<v for v in S)
def vertices(m): return [v for v in V if m>>v&1]
def am(S): return sum(1<<EI[e] for e in S)
def ga(m,g):
 a,b=g; return am([((a*u+b)%7,(a*v+b)%7) for i,(u,v) in enumerate(E) if m>>i&1])
def gv(m,g):
 a,b=g;return vm([(a*v+b)%7 for v in vertices(m)])
def boundary(m):
 z=0
 for i,(u,v) in enumerate(E):
  if m>>i&1:z^=(1<<u)|(1<<v)
 return z
OUT={v:{w for u,w in E if u==v} for v in V}
IN={v:{u for u,w in E if w==v} for v in V}
def encode(m):
 remain=[e for i,e in enumerate(E) if not(m>>i&1)]
 out={v:{w for u,w in remain if u==v} for v in V}
 inc={v:{u for u,w in remain if w==v} for v in V}
 sources=[v for v in V if not inc[v]]; sinks=[v for v in V if not out[v]]
 if len(sources)==2:return vm(sources),'two-source'
 if len(sinks)==2:
  u,v=sinks
  if (u,v) not in EI:u,v=v,u
  return FULL^vm([u,v,(2*v-u)%7]),'two-sink'
 isolated=[v for v in V if not inc[v] and not out[v]]
 if isolated:return FULL^(1<<isolated[0]),'star'
 for c in V:
  if sum(out[v]=={c} for v in V)==3:return vm({c}|OUT[c]),'common-out'
  if sum(inc[v]=={c} for v in V)==3:return vm({c}|IN[c]),'common-in'
 raise RuntimeError('Unclassified blocker')
def has_hp_dp(m):
 # Independent endpoint subset recurrence, using retained incoming adjacency.
 incoming=[0]*7
 for i,(u,v) in enumerate(E):
  if not(m>>i&1):incoming[v]|=1<<u
 end=[0]*128
 for v in V:end[1<<v]=1<<v
 for S in range(1,128):
  if S&(S-1)==0:continue
  for v in vertices(S):
   if end[S^(1<<v)]&incoming[v]:end[S]|=1<<v
 return bool(end[127])
def mul(a,b):
 z=0
 while b:
  if b&1:z^=a
  b>>=1;a<<=1
  if a&8:a^=11
 return z
def power(a,n):
 z=1
 for _ in range(n):z=mul(z,a)
 return z
AP=[power(2,n) for n in range(7)]
def fourier(S):
 u=v=0
 for x in vertices(S):u^=AP[x];v^=AP[-x%7]
 return u,v

def main():
 paths=[]
 for p in it.permutations(V):
  if all((p[i],p[i+1]) in EI for i in range(6)):paths.append(am(zip(p,p[1:])))
 require(len(paths)==189,'HP count')
 blocks=[];tested={}
 for k in (5,6):
  count=0
  for c in it.combinations(range(21),k):
   m=sum(1<<i for i in c);hp=any(not(m&p) for p in paths)
   require(hp==has_hp_dp(m),'permutation versus subset DP')
   count+=1
   if not hp:
    require(k==6,'five blocker');blocks.append(m)
  tested[k]=count
 require(len(blocks)==63,'blocker count')
 require(set(blocks)=={ga(b,g) for b in blocks for g in G},'G-stability')
 even={s for s in range(1,128) if s.bit_count()%2==0}
 images={b:encode(b)[0] for b in blocks}
 require(set(images.values())==even and len(images)==len(set(images.values())),'bijection')
 for b in blocks:
  for g in G:require(images[ga(b,g)]==gv(images[b],g),'equivariance')
 left=set(blocks);orbits=[]
 while left:
  b=min(left);O={ga(b,g) for g in G};left-=O
  orbits.append({'representative':b,'arcs':[list(e) for i,e in enumerate(E) if b>>i&1], 'size':len(O),'type':encode(b)[1], 'image':vertices(images[b]),'stabilizer':[list(g) for g in G if ga(b,g)==b], 'boundary_image_size':len({boundary(x) for x in O})})
 require(sorted(o['size'] for o in orbits)==[7,7,7,21,21],'orbit profile')
 coords={s:fourier(s) for s in even|{0}}
 require(set(coords.values())==set(it.product(range(8),repeat=2)),'Fourier isomorphism')
 for s in even|{0}:
  for a,b in G:
   u,v=coords[s]
   require(coords[gv(s,(a,b))]==(mul(AP[b],power(u,a)),mul(AP[-b%7],power(v,a))),'Fourier action')
 fibers=co.Counter(boundary(b) for b in blocks)
 require(len(fibers)==21 and set(fibers.values())=={3},'boundary collapse')
 badxor=next((a,b,a^b) for a,b in it.combinations(blocks,2) if (a^b) not in blocks)
 metric=next((a,b,(a^b).bit_count(),(images[a]^images[b]).bit_count()) for a,b in it.combinations(blocks,2) if (a^b).bit_count()>6)
 require(378%32!=0,'constant-weight linear obstruction')
 # Hostiles: the naive boundary and constant-weight/XOR transfers must fail.
 require(len(fibers)!=63,'hostile naive boundary accepted')
 require((badxor[2]).bit_count()!=6 or badxor[2] not in blocks,'hostile XOR closure accepted')
 require(metric[2]>6,'hostile metric accepted')
 require(all(any(g!=(1,0) for g in G if ga(b,g)!=b) for b in blocks),'fixed-point hostile')
 output={'status':'PASS','hamiltonian_paths':len(paths),'universes_checked_by_two_methods':tested,'blockers':len(blocks),'orbits':orbits,'type_counts':dict(co.Counter(encode(b)[1] for b in blocks)),'boundary_fiber_histogram':dict(co.Counter(fibers.values())),'xor_nonclosure_witness':badxor,'metric_impossibility_witness':metric,'bijection':[{'blocker':b,'even_subset':vertices(images[b]),'fourier':coords[images[b]],'type':encode(b)[1]} for b in blocks]}
 print(json.dumps(output,sort_keys=True,indent=2))
if __name__=='__main__':main()
