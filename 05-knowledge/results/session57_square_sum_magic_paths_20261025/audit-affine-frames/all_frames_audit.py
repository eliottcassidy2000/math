"""Brute unordered triples, with only necessary cross-sum membership pruning."""
import itertools as it, json, pathlib

def require(c,m):
    if not c:raise RuntimeError(m)
def add(a,b):return tuple((x+y)%3 for x,y in zip(a,b))
def neg(a):return tuple(-x%3 for x in a)
def sub(a,b):return add(a,neg(b))
def key(B,D):return tuple(sorted((tuple(sorted(B)),tuple(sorted(D)))))
results=[]
for n in (2,3,4):
 V=list(it.product(range(3),repeat=n));zero=(0,)*n;U={v for v in V if not any(v[2:])};C=(0,0)+(1,)*(n-2);A={add(C,u) for u in U}
 triples=list(it.combinations(U,3)); nonlinear=0
 for B in triples:
  delta={sub(b,c) for b in B for c in B if b!=c};ap=add(add(B[0],B[1]),B[2])==zero
  require(len(delta)==(2 if ap else 6),'three-set difference cardinalities')
  if not ap:
   nonlinear+=1; remainder=(U-{zero})-delta
   require(len(remainder)==2 and {neg(x) for x in remainder}==remainder,'nonlinear complement is punctured line')
 frames=set()
 for B in it.combinations(V,3):
  # B+d contained in A implies all pair differences of B belong to U.
  if any(sub(b,B[0]) not in U for b in B):continue
  candidates={sub(a,B[0]) for a in A}
  for D in it.combinations(sorted(candidates),3):
   if set(B)&set(D):continue
   if {add(b,d) for b in B for d in D}==A:frames.add(key(B,D))
 apap=mixed=neither=0
 for B,D in frames:
  apB=add(add(B[0],B[1]),B[2])==zero;apD=add(add(D[0],D[1]),D[2])==zero
  apap+=apB and apD;mixed+=apB!=apD;neither+=not(apB or apD)
 q=3**n
 require((len(frames),apap,mixed,neither)==(30*(q-9),6*(q-9),24*(q-9),0),'all frame classification/counts')
 paths=set()
 for B,D in frames:
  for b in it.permutations(B):
   for d in it.permutations(D):
    p=tuple(x for pair in zip(b,d) for x in pair);p=min(p,p[::-1]);require(p not in paths,'path injection');paths.add(p)
 require(len(paths)==1080*(q-9),'all-frame path count')
 results.append(dict(q=q,C=C,all_frames=len(frames),AP_AP=apap,mixed=mixed,neither=neither,paths=len(paths),nonlinear_plane_triples=nonlinear))
p=pathlib.Path(__file__).parent/'all-frames-results.json';p.write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results))
