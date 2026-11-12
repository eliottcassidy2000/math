"""Full target-capped closures; compare exact continuation languages, not just scalar values."""
import json,collections
from capped import build,minimize,count

def language_equivalence(N,targets):
 vs,ts,ws,hs=build(N)
 rows=[]
 for target in targets:
  mt,mh,part=minimize(ts,[int(h==target) for h in hs],0)
  # Count all distinct target-continuation types among each observed scalar h.
  types=collections.defaultdict(set)
  for h,p in zip(hs,part):types[h].add(p)
  rows.append({'target':target,'minimal_exact_DFA_states':len(mt),'same_count_types':{h:len(x) for h,x in sorted(types.items())}})
 return rows

def run():
 table=[]
 for N in (7,13,21,23,25,27,31,33,41,43,63,65,81,83,85,87,125,127):
  vs,ts,ws,hs=build(N)
  ct,co,cp=minimize(ts,hs)
  et,eo,ep=minimize(ts,[int(h==N) for h in hs],0)
  tt,to,tp=minimize(ts,[1 for h in hs],0)
  # Every live SCC in these two minimized automata must be trivial (proved normal form).
  for tr,out in ((et,eo),(tt,to)):
   n=len(tr);adj=[set(x for x in r if x!=i) for i,r in enumerate(tr)]
   indeg=[0]*n
   for row in adj:
    for j in row:indeg[j]+=1
   q=collections.deque(i for i,x in enumerate(indeg) if not x);num=0
   while q:
    i=q.popleft();num+=1
    for j in adj[i]:
     indeg[j]-=1
     if not indeg[j]:q.append(j)
   if num!=n:raise RuntimeError(('nontrivial SCC',N))
  by=collections.defaultdict(set)
  for h,p in zip(hs,ep):by[h].add(p)
  targettypes={h:len(by[h]) for h in ((N-1)//2,(N+1)//2) if h in by}
  table.append({'N':N,'capped_live_states':len(vs),'minimal_count_Moore_states':len(ct),'minimal_threshold_DFA_states':len(tt),'minimal_exact_DFA_states':len(et),'half_target_count_types':targettypes,'all_SCCs_singleton':True})
  print(N,len(vs),len(et),flush=True)
 # concrete target-sensitive ambiguity
 N=41;vs,ts,ws,hs=build(N);et,eo,ep=minimize(ts,[int(h==N) for h in hs],0)
 def act(w):
  q=ep[0]
  for c in w:q=et[q][c]
  return q
 a,b=(1,2,0),(2,1,2);q=collections.deque([(act(a),act(b),())]);seen=set();found=None
 while q:
  x,y,w=q.popleft()
  if (x,y) in seen:continue
  seen.add((x,y))
  if eo[x]!=eo[y]:found=w;break
  for c in range(4):q.append((et[x][c],et[y][c],w+(c,)))
 if found is None:raise RuntimeError('expected distinguishing suffix')
 witness={'target':41,'prefixes':['120','212'],'prefix_counts':[count(a),count(b)],'shortest_distinguishing_suffix':''.join(map(str,found)),'final_counts':[count(a+found),count(b+found)]}
 return {'table':table,'continuation_witness':witness}
if __name__=='__main__':
 out=run()
 with open('05-knowledge/results/session276_corridors_20261112/wildcard/COMPARISON.json','w') as f:json.dump(out,f,indent=2)
