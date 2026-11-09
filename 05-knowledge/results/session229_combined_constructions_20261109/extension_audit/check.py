import json, math, hashlib
from pathlib import Path
from collections import defaultdict,Counter
SRC=Path(__file__).resolve().parent.parent/'rich'

def require(ok,msg):
 if not ok: raise RuntimeError(msg)

def coords(adj):
 # Independently enumerate directed path sequences, with no producer imports or DP.
 n=len(adj); N=1<<n; F=N-1
 start=[[0]*n for _ in range(N)]; end=[[0]*n for _ in range(N)]
 def walk(first,last,mask):
  start[mask][first]+=1; end[mask][last]+=1
  avail=adj[last]&~mask
  while avail:
   bit=avail&-avail; avail-=bit; walk(first,bit.bit_length()-1,mask|bit)
 for v in range(n): walk(v,v,1<<v)
 h=[sum(row) for row in start];h[0]=1
 require(all(sum(end[s])==h[s] for s in range(1,N)),'endpoint census')
 H=h[F];Q=sum(h[s]*h[F^s] for s in range(N))
 A=[sum(start[s][v]*h[F^s] for s in range(1,N)) for v in range(n)]
 G=[sum(end[s][v]*h[F^s] for s in range(1,N)) for v in range(n)]
 K=[[sum(end[s][v]*start[F^s][w] for s in range(1,F)) for w in range(n)] for v in range(n)]
 result=[]
 for m in range(N):
  vs=[v for v in range(n) if m&(1<<v)]
  sig=(H,Q,sum(A[v] for v in vs),sum(G[v] for v in vs),sum(start[F][v] for v in vs),sum(end[F][v] for v in vs))
  T=sum(K[v][w] for v in vs for w in vs)
  result.append((sig,T))
 return result

CORES=[[[48,49,3,7,12,28],[40,37,17,6,11,28]],[[64,121,35,37,13,81,28],[20,85,72,195,204,159,33,71]]]
def census(cores,filename):
 fibers=defaultdict(lambda:defaultdict(list)); eligible=Counter();candidates=0;cuts=0; lookup={}
 for ci,core in enumerate(cores):
  n=len(core)
  for orientation in range(1<<n):
   candidates+=1
   # Construct arcs explicitly from their meaning, rather than copying producer expression.
   a=list(core)+[0]
   for v in range(n):
    if orientation&(1<<v): a[n]|=1<<v
    else:a[v]|=1<<n
   require(all(not (a[i]>>i&1) for i in range(n+1)),'loop')
   require(all(((a[i]>>j&1)+(a[j]>>i&1))==1 for i in range(n+1) for j in range(i)),'not tournament')
   values=coords(a)
   if math.gcd(*values[0][0][:2])!=1:continue
   eligible[ci]+=1; cuts+=len(values)
   for mask,(sig,T) in enumerate(values):
    witness=(ci,orientation,mask)
    fibers[sig][T].append(witness);lookup[witness]=(sig,T,a)
 runs=[]
 for sig,ts in fibers.items():
  for t in sorted(ts):
   if t-2 in ts:continue
   run=[];u=t
   while u in ts:run.append(u);u+=2
   runs.append((sig,run))
 best=max(len(r) for s,r in runs); maximal=[(s,r) for s,r in runs if len(r)==best]
 receipt=json.loads((SRC/filename).read_text())
 for field,val in [('eligible_extensions',sum(eligible.values())),('interfaces',len(fibers)),('max_run',best),('maximal_run_count',len(maximal))]:require(receipt[field]==val,field)
 if 'max_fiber_size' in receipt: require(receipt['max_fiber_size']==max(map(len,fibers.values())),'fiber max')
 pair_runs=[(s,r) for s,r in runs if len(r)>=2]
 all_pair_witnesses=[w for s,r in pair_runs for t in r for w in fibers[s][t]]
 orientations=sorted(set((c,o) for c,o,m in all_pair_witnesses))
 require(all(o in (0,(1<<len(cores[c]))-1) for c,o,m in all_pair_witnesses),'nonordinal witness to pair')
 if 'maximal_run_seed_orientations' in receipt:require(receipt['maximal_run_seed_orientations']==[list(v) for v in orientations],'orientations')
 seen=set()
 for rec in receipt['witnesses']:
  sig=tuple(rec['signature']);run=[]
  for w in rec['run']:
   key=(w['core'],w['orientation'],w['mask']);s,t,a=lookup[key]
   require((s,t,a)==(sig,w['T'],w['adj']),'physical receipt mismatch');run.append(t)
  require((sig,run) in maximal,'run mismatch');seen.add((sig,tuple(run)))
 if filename=='EXTEND_UNSATURATED.json':require(len(seen)==len(maximal),'missing full runs')
 canonical=[[list(sig),[[t,[list(w) for w in ws]] for t,ws in sorted(ts.items())]] for sig,ts in sorted(fibers.items())]
 digest=hashlib.sha256(json.dumps(canonical,separators=(',',':')).encode()).hexdigest()
 nonordinal_multis=[]
 for sig,ts in fibers.items():
  if len(ts)>1 and any(o not in (0,(1<<len(cores[c]))-1) for ws in ts.values() for c,o,m in ws):
   nonordinal_multis.append({'signature':sig,'values':sorted(ts),'witnesses':{str(t):ws for t,ws in ts.items()}})
 return dict(candidates=candidates,eligible_by_core=dict(eligible),eligible_extensions=sum(eligible.values()),cut_records=cuts,interfaces=len(fibers),fiber_size_histogram=dict(sorted(Counter(map(len,fibers.values())).items())),max_run=best,maximal_run_count=len(maximal),all_pair_physical_witness_count=len(all_pair_witnesses),all_pair_orientations=orientations,nonordinal_multivalue_interface_count=len(nonordinal_multis),nonordinal_multivalue_examples=nonordinal_multis[:3],complete_census_sha256=digest)

out={name:census(cores,name) for cores,name in zip(CORES,['EXTEND_TRIPLES.json','EXTEND_UNSATURATED.json'])}
out['inputs_sha256']={p:hashlib.sha256((SRC/p).read_bytes()).hexdigest() for p in ['extend_triples.py','extend_unsaturated.py','EXTEND_TRIPLES.json','EXTEND_UNSATURATED.json','EXTENSION_PROOF.md','PROOF.md','search.py']}
out['status']='PASS';out['method']='Exhaustive directed path sequence enumeration; source coordinates of ordered path covers; complete full-interface grouping retaining every physical witness, no producer imports; explicit checks (no assert).'
print(json.dumps(out,indent=2,sort_keys=True))
