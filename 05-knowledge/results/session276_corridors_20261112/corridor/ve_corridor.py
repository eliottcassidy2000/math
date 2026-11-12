"""Search a response-preserving local corridor, not a count-coverage census."""
import sys,json,itertools
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'session275_exact_selection_20261111'/'switches'))
from carrier import profile
A=(52,13,24,49,34,6);n=len(A);ed=list(itertools.combinations(range(n),2));moves=[(e,) for e in ed]+list(itertools.combinations(ed,2))
def flip(t,es):
 t=list(t)
 for i,j in es:t[i]^=1<<j;t[j]^=1<<i
 return tuple(t)
def common(a,b):
 # n6 exact permutations; retain actual shared frames.
 return next((p for p in itertools.permutations(range(n)) if all(a[u]>>v&1 and b[u]>>v&1 for u,v in zip(p,p[1:]))),None)
def run():
 seen={A:None};pdata={};q=[A];edges={};cache={}
 def pr(t):
  if t not in cache:cache[t]=profile(t)
  return cache[t]
 while q:
  a=q.pop();p=pr(a)
  for es in moves:
   b=flip(a,es)
   if b in seen:continue
   v=pr(b)
   if v['V']!=57 or abs(v['H']-p['H'])>2:continue
   frame=common(a,b)
   if frame is None:continue
   seen[b]=a;edges[b]=(es,frame);q.append(b)
 vals=sorted({pr(a)['H'] for a in seen});best=max(seen,key=lambda a:pr(a)['H']);path=[]
 while best is not None:path.append(dict(rows=best,profile=pr(best),entry=edges.get(best)));best=seen[best]
 print(json.dumps(dict(states=len(seen),values=vals,path=path[::-1],tested_states=len(cache)),indent=2))
if __name__=='__main__':run()
