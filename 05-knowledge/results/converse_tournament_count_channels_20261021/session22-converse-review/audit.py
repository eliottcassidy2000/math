from pathlib import Path
import hashlib,itertools,json,math,sys
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'session22-converse-channels'))
from channels import channel_dp,counts
checks=0
def require(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
pins=json.loads((ROOT/'session22-converse-channels/source-pins.json').read_text())
REPO=next(p for p in ROOT.parents if (p/'AGENTS.md').is_file() and (p/'01-canon').is_dir())
def provenance_path(path):
 if path=='math/01-canon/theorems/THM-283-sc-burnside-formula.md':
  return ROOT/'provenance/THM-283-original-historical-error.md'
 if not path.startswith('math/'):
  raise RuntimeError('Unexpected provenance path: '+path)
 return REPO/path.removeprefix('math/')
for p,h in pins.items():require(hashlib.sha256(provenance_path(p).read_bytes()).hexdigest()==h,'pin '+p)
def literal(n,p):
 edges=list(itertools.combinations(range(n),2));idx={e:i for i,e in enumerate(edges)}
 mapped=[(idx[tuple(sorted((p[u],p[v])))],p[u]>p[v]) for u,v in edges]
 full=(1<<len(edges))-1
 fixed=anti=0
 for x in range(full+1):
  y=sum(((((x>>i)&1)^swap)<<j) for i,(j,swap) in enumerate(mapped))
  fixed+=(x==y);anti+=(x==(y^full))
 unseen=set(range(n));cycles=0
 while unseen:
  v=next(iter(unseen));cycles+=1
  while v in unseen:unseen.remove(v);v=p[v]
 return fixed,anti,cycles
vectors=counts(14); rows=[]
for n in range(6):
 entries=[literal(n,p) for p in itertools.permutations(range(n))]
 require(sum(t[0] for t in entries)==math.factorial(n)*vectors['ordinary'][n],'ordinary literal masks')
 require(sum(t[1] for t in entries)==math.factorial(n)*vectors['self_converse'][n],'anti literal masks')
 for a,b in [(1,0),(1,2),(2,0),(2,1),(3,2)]:
  value=channel_dp(n,a,b)[0][n]
  require(sum(f**a*2**(b*k) for f,anti,k in entries)==math.factorial(n)*value,'multi-channel labelled fixed tuples')
 rows.append({'n':n,'ordinary':vectors['ordinary'][n],'self_converse':vectors['self_converse'][n]})
require(literal(2,(1,0))[1]==2,'transposition allowed antipodal length one')
require(literal(4,(1,2,3,0))[1]==0,'four-cycle forbidden antipodal length two')
for n in range(15):
 require(vectors['ordinary'][n]==(int(n==0)+sum(vectors['strong'][k]*vectors['ordinary'][n-k] for k in range(1,n+1))),'A(1-S)=1')
 require(vectors['self_converse'][n]==sum(vectors['ordinary'][j]*(int(n==2*j)+vectors['strong_self_converse'][n-2*j]) for j in range(n//2+1)),'C=A(z2)(1+U)')
for bad in [True,False,1.,'1',None,[],{},-1]:
 for name in ['bound','a','b']:
  kw={'bound':3,'a':1,'b':0};kw[name]=bad
  try:channel_dp(**kw)
  except ValueError:pass
  else:raise RuntimeError('accepted bad '+name)
  checks+=1
for mode in ['core','no_fixed']:
 for bad in [0,1,None,'False',[],{}]:
  try:channel_dp(3,**{mode:bad})
  except ValueError:pass
  else:raise RuntimeError('accepted bad '+mode)
  checks+=1
receipt={'checks':checks,'literal_mask_burnside':rows,'pins_verified':pins,'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_sha256':{p:hashlib.sha256((ROOT/'session22-converse-channels'/p).read_bytes()).hexdigest() for p in ['README.md','channels.py','check.py','incoming_local_scale.py']}}
print(json.dumps(receipt,indent=2))
