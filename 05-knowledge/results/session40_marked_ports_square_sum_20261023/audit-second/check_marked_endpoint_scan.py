#!/usr/bin/env python3
"""Independent Python scan, endpoint existence from backwards subset recursion."""
import collections,functools,json,csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def locate(*paths):
 for p in paths:
  if p.exists():return p
 raise RuntimeError('Missing required input: '+str(paths))

def require(b,msg):
 if not b:raise RuntimeError(msg)
def adjacency(code):
 A=[set() for _ in range(8)];k=27
 for i in range(8):
  for j in range(i+1,8):
   if code>>k&1:A[i].add(j)
   else:A[j].add(i)
   k-=1
 return A

def signature(A,s,t):
 colors=[0 if v==s else 1 if v==t else 2 for v in range(8)];history=[]
 while True:
  nc=len(set(colors));keys=[]
  for v in range(8):
   counts=collections.Counter(colors[w] for w in A[v]);keys.append((colors[v],tuple(counts[c] for c in range(nc))))
  history.append(tuple(sorted(keys)));uniq=sorted(set(keys));new=[uniq.index(k) for k in keys]
  if len(uniq)==nc:return tuple(history),tuple(sorted(collections.Counter(new).values())),[[v for v in range(8) if new[v]==c] for c in range(nc)]
  colors=new

def run():
 codes=sorted(int(r['canon']) for r in csv.DictReader(open(locate(ROOT/'session38-wiring/classes_input_n8.tsv',ROOT.parent/'session38_incidence_ports_electrical_20261023/wiring/classes_input_n8.tsv')),delimiter='\t'))
 groups=collections.defaultdict(list);augmented=collections.defaultdict(list);shapes=collections.Counter();total=0;chirality_controls=0
 for code in codes:
  A=adjacency(code);out=[sum(1<<v for v in row) for row in A]
  @functools.lru_cache(None)
  def ends(mask,start):
   if mask==1<<start:return 1<<start
   rest=mask^(1<<start);options=out[start]&rest;answer=0
   while options:
    b=options&-options;options-=b;answer|=ends(rest,b.bit_length()-1)
   return answer
  for s in range(8):
   bits=ends(255,s)
   for t in range(8):
    if bits>>t&1:
     sig,shape,cells=signature(A,s,t);groups[sig].append((code,s,t));shapes[shape]+=1;total+=1
     triples=[cell for cell in cells if len(cell)==3];J=0
     if len(triples)==2:
      X,Y=triples;q=len(A[X[0]]&set(Y))
      if q in (1,2):
       owners={u:[v for v in Y if (v in A[u])==(q==1)] for u in X}
       require(all(len(z)==1 for z in owners.values()),'unique matching owner')
       J=sum(1 if owners[w][0] in A[owners[u][0]] else -1 for u in X for w in X if w in A[u])
       require(abs(J)==3,'chirality magnitude');chirality_controls+=1
     augmented[(sig,J%8)].append((code,s,t))
 collisions=sorted(sorted(g) for g in groups.values() if len({v[0] for v in g})>1)
 expect=[]
 for line in locate(ROOT/'session40-marked-ports/endpoint_scan.txt',ROOT/'marked-ports/endpoint_scan.txt').read_text().splitlines():
  if line.startswith('collision '):expect.append(sorted(tuple(map(int,x.split(':'))) for x in line.split()[1:]))
 require(collisions==sorted(expect),'collision list mismatch')
 require(total==296534,'endpoint total');require(len(groups)==291810,'signature total');require(len(augmented)==291828,'augmented signatures');require(all(len({z[0] for z in g})==1 for g in augmented.values()),'augmented class collision')
 require(dict(shapes)=={(1,)*8:293686,(1,1,1,1,1,3):2776,(1,1,1,5):16,(1,1,3,3):56},'shape mismatch')
 out={'status':'PASS','classes':len(codes),'actual_endpoint_pairs':total,'signatures':len(groups),'augmented_signatures':len(augmented),'chirality_controls':chirality_controls,'augmented_cross_class_collisions':0,'shapes':{','.join(map(str,k)):v for k,v in sorted(shapes.items())},'cross_class_collision_groups':collisions,'method':'backwards subset recursion and exact Python tuple refinement histories'}
 Path(__file__).with_name('marked-endpoint-audit-receipt.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k!='cross_class_collision_groups'})
if __name__=='__main__':run()
