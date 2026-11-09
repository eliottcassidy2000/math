import itertools as it,json,math,hashlib,collections,time,gzip
from pathlib import Path
BASE=Path(__file__).resolve().parent.parent/'fibers'
OUT=Path(__file__).parent

def check(b,msg):
 if not b:raise RuntimeError(msg)
def validate(adj):
 n=len(adj)
 check(all(not(a>>i&1) and 0<=a<1<<n for i,a in enumerate(adj)), 'invalid diagonal')
 check(all(((adj[i]>>j)^(adj[j]>>i))&1 for i in range(n) for j in range(i)), 'not tournament')
def literal(adj):
 validate(adj);n=len(adj);h=0;q=0;aw=[0]*n;gw=[0]*n;sw=[0]*n;ew=[0]*n;kk=[[0]*n for _ in adj]
 # Enumerate each ordered pair of disjoint covering paths by a permutation and split.
 for p in it.permutations(range(n)):
  breaks=[i+1 for i in range(n-1) if not(adj[p[i]]>>p[i+1]&1)]
  if not breaks:
   h+=1;sw[p[0]]+=1;ew[p[-1]]+=1;splits=range(n+1)
  elif len(breaks)==1:splits=breaks
  else:continue
  for k in splits:
   q+=1
   if k<n:aw[p[k]]+=1
   if k:gw[p[k-1]]+=1
   if k and k<n:kk[p[k-1]][p[k]]+=1
 return [(h,q,sum(aw[v] for v in range(n) if m>>v&1),sum(gw[v] for v in range(n) if m>>v&1),sum(kk[u][v] for u in range(n) for v in range(n) if m>>u&1 and m>>v&1),sum(sw[v] for v in range(n) if m>>v&1),sum(ew[v] for v in range(n) if m>>v&1)) for m in range(1<<n)]
def decode(code,n):
 a=[0]*n;b=0
 for i in range(n):
  for j in range(i+1,n):
   if code>>b&1:a[i]|=1<<j
   else:a[j]|=1<<i
   b+=1
 return a

def canonical(a):
 # Independent exhaustive degree-compatible labellings, no refinement heuristic.
 groups=collections.defaultdict(list)
 for i,v in enumerate(a):groups[v.bit_count()].append(i)
 best=None;aut=0
 for blocks in it.product(*(it.permutations(groups[d]) for d in sorted(groups))):
  p=sum(blocks,());code=0
  for i in range(len(p)):
   for j in range(i+1,len(p)):code=2*code+((a[p[i]]>>p[j])&1)
  if best is None or code<best:best=code;aut=1
  elif code==best:aut+=1
 return best,aut

def graph(states,name):
 h,q=states[0][:2];physical={s[2:4] for s in states};tot=0;ret=0;sources=0;short=None;trivial=[];dig=hashlib.sha256()
 cert=gzip.open(OUT/(name+".graph.jsonl.gz"),"wt")
 for src in sorted(physical):
  distance={src:0};todo=collections.deque([src]);hits=0;edges=[]
  while todo:
   a,g=todo.popleft()
   # Use interface equations directly; test both x values, never assume parity rule.
   for x,y in it.product((0,1),repeat=2):
    A=2*a+h*y-q*x;twog=g+q*y-h*x
    if twog%2:continue
    G=twog//2
    if not(-h<A<q and -h<G<q):continue
    dest=A,G;edges.append((a,g,x,y,A,G))
    if dest in physical:
     if src==dest and src in [(0,0),(q-h,q-h)]:
      trivial.append((src,(a,g),x,y,distance[a,g]+1))
      check((a,g)==src and x==y and x==(src!=(0,0)), 'NONTRIVIAL ENDPOINT RETURN')
     else:
      hits+=1;short=min(short or 10**9,distance[a,g]+1)
    if dest not in distance:distance[dest]=distance[a,g]+1;todo.append(dest)
  cert.write(json.dumps(dict(source=src,states=sorted(distance),edges=sorted(edges)),separators=(",",":"))+"\n")
  tot+=len(distance);ret+=hits;sources+=bool(hits)
  dig.update(repr((src,sorted(distance),sorted(edges))).encode())
 cert.close()
 return dict(physical=len(physical),source_counted_reachable=tot,return_sources=sources,return_ingresses=ret,shortest=short,endpoint_ingresses=trivial,all_edges_sha256=dig.hexdigest())

def census(n):
 lines=[s for s in (BASE/f'classes9_{n}.txt').read_text().splitlines() if s and not s.startswith('#')]
 seen=set();mass=0;fibers=collections.defaultdict(set);target=set();witnesses={}
 for line in lines:
  row=line.split();code=int(row[0]);a=decode(code,n);c,aut=canonical(a)
  check(c not in seen,'duplicate isomorphism class');seen.add(c)
  check(aut==int(row[5]),'automorphism count mismatch');mass+=math.factorial(n)//aut
  ss=literal(a)
  check(ss[0][0]==int(row[1]),'Hamilton census mismatch')
  for m,(h,q,A,G,T,S,E) in enumerate(ss):
   for b in (0,1):
    key=h,q,A,G,b;v=T+S+b*E;fibers[key].add(v);witnesses.setdefault((key,v),(code,m))
   if (h,q,A,G)==(45,432,129,129):target.add(T+S)
 check(mass==2**(n*(n-1)//2),'orbit mass incomplete')
 best=ubest=0;uw=None;bw=None
 for key,vals in fibers.items():
  for v in vals:
   if v-2 in vals:continue
   r=1
   while v+2*r in vals:r+=1
   if r>best:best=r;bw=(key,[witnesses[key,v+2*i] for i in range(r)])
   if math.gcd(*key[:2])==1 and r>ubest:ubest=r;uw=(key,[witnesses[key,v+2*i] for i in range(r)])
 return dict(classes=len(lines),independently_verified_orbit_mass=mass,max_run=best,unsaturated_max_run=ubest,number_fibers=len(fibers),target=sorted(target),best_witness=bw,unsaturated_witness=uw)

if __name__ == "__main__":
 out={}
 for name,a in {'B29':[26,44,57,16,34,9],'B95':[64,121,35,37,13,81,28],'B389':[20,85,72,195,204,159,33,71]}.items():
  ss=literal(a);out[name]=graph(ss,name)
  f=collections.defaultdict(set)
  for s in ss:f[s[2:4]].add(s[4:])
  out[name]['max_AG_fiber']=max(map(len,f.values()))
  out[name]['AG_ambiguous']=[list(v) for v in f.values() if len(v)>1]
  print(name,out[name]['source_counted_reachable'],flush=True)
 for n in (6,7):
  out['census'+str(n)]=census(n);print('census',n,out['census'+str(n)],flush=True)
 (OUT/'AUDIT.json').write_text(json.dumps(out,indent=2))
