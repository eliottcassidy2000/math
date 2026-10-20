#!/usr/bin/env python3
"""Independent complete n=5 incidence audit; imports no producer code."""
from itertools import combinations, permutations
from collections import Counter, defaultdict
from pathlib import Path
import hashlib,json
checks=0
def check(v,msg):
 global checks
 checks+=1
 if not v: raise RuntimeError(msg)
n=5
E=tuple(combinations(range(n),2)); bit={e:1<<i for i,e in enumerate(E)}
B=lambda a,b:bit[tuple(sorted((a,b)))]
full=(1<<10)-1
perms=tuple(permutations(range(n)))
cut=lambda s:sum(w for (u,v),w in bit.items() if bool(s>>u&1)!=bool(s>>v&1))
cuts={cut(s) for s in range(32)}
cls={x:min(x^c for c in cuts) for x in range(1024)}
classes=sorted(set(cls.values()))
check(len(classes)==64 and len(cuts)==16,'quotient dimensions')
cycles={}
for tail in permutations(range(1,n)):
 word=(0,)+tail
 if word[1]>word[-1]: continue
 mask=sum(B(u,v) for u,v in zip(word,word[1:]+word[:1]))
 cycles[mask]=word
check(len(cycles)==12,'Hamilton cycle census')
cyclepairs=sorted({tuple(sorted((C,full^C))) for C in cycles})
check(len(cyclepairs)==6,'complementary cycle-pair census')
check(all(C in cycles and D in cycles for C,D in cyclepairs),'cycle complements')
# Use a canonical ordering discovered from all cycles, not the producer row masks.
encode=lambda x:sum(((x&C).bit_count()%2)<<i for i,(C,D) in enumerate(cyclepairs))
evenclasses=[x for x in classes if x.bit_count()%2==0]
check(len(evenclasses)==32,'even component size')
for x in range(1024):
 for C,D in cyclepairs:
  if x.bit_count()%2==0: check((x&C).bit_count()%2==(x&D).bit_count()%2,'pair representative independence')
 for c in cuts:
  check(encode(x)==encode(x^c),'cut representative independence')
images={encode(x) for x in evenclasses}
check(images=={y for y in range(64) if y.bit_count()%2==0},'bijective even cycle-coordinate map')
# Literal path incidence across all 1024 tournaments, not generator aggregation.
rows=defaultdict(Counter)
incidences=set()
for x in range(1024):
 for p in perms:
  if all(bool(u<v)^bool(x&B(u,v)) for u,v in zip(p,p[1:])):
   pm=sum(B(u,v) for u,v in zip(p,p[1:]))
   y=x^full^pm
   check((x,y) not in incidences,'path witness collision')
   incidences.add((x,y)); rows[cls[x]][cls[y]]+=1
check(len(incidences)==7680,'total path incidence')
for q in classes:
 check(sum(rows[q].values())==120,'aggregate row sum')
 for r in classes:
  expected=8 if (q^r).bit_count()%2==0 and encode(q^r).bit_count()==4 else 0
  check(rows[q][r]==expected,'all-pair weighted adjacency')
# Recover literal matrix and compare spectra via a matrix polynomial, not Fourier.
V=evenclasses
A=[[int(rows[x][y]>0) for y in V] for x in V]
Y=sorted(images)
H=[[int((x^y).bit_count()==2) for y in Y] for x in Y]
def multiply(M,N):
 return [[sum(a*b for a,b in zip(row,col)) for col in zip(*N)] for row in M]
def exact_spectrum(M,expected):
 size=len(M)
 powers=[[[int(i==j) for j in range(size)] for i in range(size)]]
 for k in range(1,5): powers.append(multiply(powers[-1],M))
 for k in range(5): check(sum(powers[k][i][i] for i in range(size))==sum(mult*eigen**k for eigen,mult in expected.items()),'spectral trace moment')
 value=powers[0]
 for root in expected:
  factor=[[M[i][j]-root*int(i==j) for j in range(size)] for i in range(size)]
  value=multiply(value,factor)
 check(all(x==0 for row in value for x in row),'spectral polynomial annihilator')
 return [sum(powers[k][i][i] for i in range(size)) for k in range(5)]
actual_spectrum={15:1,3:10,-1:15,-5:6}
half_spectrum={15:1,5:6,-1:15,-3:10}
traceA=exact_spectrum(A,actual_spectrum);traceH=exact_spectrum(H,half_spectrum)
# Explicit Sylow subgroups, conjugation, and the six-coordinate action.
def compose(p,q): return tuple(p[q[i]] for i in range(n))
def inverse(p): return tuple(p.index(i) for i in range(n))
def rotations(word):
 g=tuple(word[(word.index(i)+1)%n] for i in range(n))
 power=tuple(range(n)); group=set()
 for k in range(5): group.add(power);power=compose(g,power)
 return frozenset(group)
subgroups=[rotations(cycles[C]) for C,D in cyclepairs]
check(len(set(subgroups))==6,'distinct Sylow5 groups')
for i,(C,D) in enumerate(cyclepairs): check(rotations(cycles[D])==subgroups[i],'complementary rotation groups')
def permute_edges(x,p): return sum(B(p[u],p[v]) for (u,v),w in bit.items() if x&w)
actions=set();normalizer_sizes=[]
for i,G in enumerate(subgroups):
 count=0
 for p in perms:
  pinv=inverse(p)
  conjugate=frozenset(compose(compose(p,g),pinv) for g in G)
  if conjugate==G: count+=1
 normalizer_sizes.append(count)
check(normalizer_sizes==[20]*6,'normalizer orders')
for p in perms:
 action=tuple(cyclepairs.index(tuple(sorted((permute_edges(C,p),permute_edges(D,p))))) for C,D in cyclepairs)
 actions.add(action)
 pinv=inverse(p)
 for i,G in enumerate(subgroups):
  check(frozenset(compose(compose(p,g),pinv) for g in G)==subgroups[action[i]],'cycle/Sylow equivariance')
 for q in evenclasses:
  expected=sum(((encode(q)>>i)&1)<<action[i] for i in range(6))
  check(encode(permute_edges(q,p))==expected,'coordinate equivariance')
check(len(actions)==120,'faithful action')
check({(a[0],a[1]) for a in actions}==set(permutations(range(6),2)),'double transitivity')
# Literal formula (2) and formula (3), independently tested against edge cycles.
words=[(0,4,1,2,3),(0,4,1,3,2),(0,3,1,2,4),(0,3,1,4,2),(0,1,3,2,4),(0,1,4,2,3)]
Cs=[sum(B(u,v) for u,v in zip(w,w[1:]+w[:1])) for w in words]
free=list(combinations(range(1,5),2))
for q in range(64):
 if q.bit_count()%2: continue
 a,b,c,d,e,f=[q>>i&1 for i in range(6)]
 y=[a^c^d,b^c^d,a^b^e,b^c^e,b^d^e,c^d^e]
 z=[y[0]^y[3]^y[4],y[1]^y[3]^y[4],y[0]^y[2]^y[4],y[0]^y[2]^y[3],y[0]^y[1]^y[2],y[2]^y[3]^y[4]]
 check(z==[a,b,c,d,e,f],'explicit inverse')
 field=sum(B(u,v) for k,(u,v) in enumerate(free) if q>>k&1)
 check(y==[(field&C).bit_count()%2 for C in Cs],'explicit cycles match formula')
# Check local adjacency relation and exact clique sizes by local backtracking.
def max_clique(adj):
 best=0
 def search(chosen,candidates):
  nonlocal best
  best=max(best,chosen)
  while candidates:
   if chosen+len(candidates)<=best: return
   v=candidates.pop();search(chosen+1,candidates&adj[v])
 search(0,set(adj));return best
local={g:{h for h in range(64) if h.bit_count()==4 and (g^h).bit_count()==4} for g in range(64) if g.bit_count()==4}
localhalf={g:{h for h in range(64) if h.bit_count()==2 and (g^h).bit_count()==2} for g in range(64) if g.bit_count()==2}
for g,neighbors in local.items():
 check(neighbors=={h for h in local if ((63^g)&(63^h))==0},'local Kneser relation')
for g,neighbors in localhalf.items(): check(neighbors=={h for h in localhalf if (g&h).bit_count()==1},'local triangular relation')
maxA=max_clique(local)+1;maxH=max_clique(localhalf)+1
check((maxA,maxH)==(4,6),'clique maxima')
result={'status':'PASS','checks':checks,'labelled_tournaments':1024,'ordered_incidence_pairs':7680,'switching_classes':64,'component_vertices':32,'weighted_adjacency_comparisons':4096,'canonical_cycle_pairs':cyclepairs,'six_point_action_order':len(actions),'normalizer_orders':normalizer_sizes,'actual_spectrum':actual_spectrum,'halved_cube_spectrum':half_spectrum,'actual_trace_moments_0_through_4':traceA,'halved_cube_trace_moments_0_through_4':traceH,'actual_max_clique':maxA,'halved_cube_max_clique':maxH,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
print(json.dumps(result,indent=2,sort_keys=True))
