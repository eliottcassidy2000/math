from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'session275_exact_selection_20261111'/'switches'))
from carrier import profile,fuse
from escape_check import context,T41,T43,T45

def reroot(t,r):
 ids=[r]+[i for i in range(len(t)) if i!=r]
 return tuple(sum(((t[u]>>v)&1)<<j for j,v in enumerate(ids)) for u in ids)
out=[]
for n in range(2,8):
 ts=[fuse(tuple(context(n,1)),tuple(t)) for t in [T41,T43,T45]]
 menus=[]
 for t in ts:
  menu={}
  for r in range(len(t)):
   p=profile(reroot(t,r));menu.setdefault(p['V'],[]).append([r,p['h'],p['U'],p['e'],p['d']])
  menus.append(menu)
 common=sorted(set(menus[0])&set(menus[1])&set(menus[2]))
 def odd(v):return v//(v&-v)
 normalized=[{odd(v):v for v in m} for m in menus]
 commonodd=sorted(set(normalized[0])&set(normalized[1])&set(normalized[2]))
 out.append(dict(n=n,H=[profile(t)['H'] for t in ts],common_V=common,common_odd_parts=commonodd,normalized_matches=[{v:m[v] for v in commonodd} for m in normalized],common_roots=[{v:m[v] for v in common} for m in menus]))
 print(json.dumps(out[-1]),flush=True)
