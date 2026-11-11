"""Independent finite controls for even-component domination expansion and endpoint lemma."""
import itertools,json,random,math
from pathlib import Path

def components(n,edges,mask):
 unseen={v for v in range(n) if mask>>v&1};cs=[]
 while unseen:
  q=[min(unseen)];unseen.remove(q[0]);c=set(q)
  while q:
   v=q.pop()
   for w in list(unseen):
    if (min(v,w),max(v,w)) in edges:unseen.remove(w);q.append(w);c.add(w)
  cs.append(c)
 return cs

def d_direct(n,es):
 return sum(all(s>>v&1 or any(s>>w&1 and (min(v,w),max(v,w)) in es for w in range(n) if v!=w) for v in range(n)) for s in range(1<<n))

def proof_controls():
 graphs=0
 for n in range(0,6):
  pairs=list(itertools.combinations(range(n),2))
  for b in range(1<<len(pairs)):
   es={e for i,e in enumerate(pairs) if b>>i&1}
   expansion=0
   for s in range(1<<n):
    cs=components(n,es,s)
    if all(len(c)%2==0 for c in cs):expansion+=2**len(cs)
   if expansion!=d_direct(n,es):raise RuntimeError(('domination identity',n,b))
   graphs+=1
 # Independently enumerate actual directed permutations on selected vertex sets.
 rng=random.Random(2611711);cases=0;fibres=0;strict=0
 for n in range(1,7):
  # All bases up to 4, then small seeded targeted controls, not a spectrum search.
  bits=range(1<<(n*(n-1)//2)) if n<=4 else [rng.getrandbits(n*(n-1)//2) for _ in range(35)]
  for b in bits:
   out=[set() for _ in range(n)]
   for i,(u,v) in enumerate(itertools.combinations(range(n),2)):
    if b>>i&1:out[u].add(v)
    else:out[v].add(u)
   for am in range(1<<n):
    A={v for v in range(n) if am>>v&1};B=set(range(n))-A
    es={(min(a,bb),max(a,bb)) for a in A for bb in B if bb in out[a]}
    P=0; factorial_sum=0;dommass=0
    for mask in range(1,1<<n):
     if mask.bit_count()%2:continue
     vs=[v for v in range(n) if mask>>v&1]
     pc=sum(p[0] in A and p[-1] in B and all(y in out[x] for x,y in zip(p,p[1:])) for p in itertools.permutations(vs))
     cs=components(n,es,mask);admissible=all(len(c)%2==0 for c in cs)
     if admissible:
      floor=math.factorial(len(cs));factorial_sum+=floor;dommass+=2**len(cs)
      if pc<floor:raise RuntimeError(('path fibre lower bound',n,b,am,mask,pc,cs))
      strict+=2*floor>2**len(cs)
     P+=pc;fibres+=1
    d=d_direct(n,es)
    if d!=1+dommass or 2*P<dommass or P<factorial_sum:raise RuntimeError('main inequality')
    cases+=1
 return {'general_graph_identity_scope':'All labelled simple graphs orders 0..5','general_graph_cases':graphs,'endpoint_scope':'All tournaments bases 1..4 and every A/B colouring; 35 seeded bases each of orders5,6 with every colouring; direct permutations','seed':2611711,'coloured_bases':cases,'even_subset_fibres':fibres,'strict_factorial_fibres':strict,'result':'PASS; finite controls support the all-order proofs'}
if __name__=='__main__':
 r=proof_controls();Path(__file__).with_name('extension4_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
