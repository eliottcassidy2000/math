import json,itertools
from pathlib import Path
p=Path(__file__).parent

def require(a,b):
 if a!=b:raise RuntimeError((a,b))

def direct(m,w):
 f=g=j=0
 for first in range(1<<m):
  second=((1<<m)-1)^first
  a=bool(first and (w>>(first.bit_length()-1)&1));b=bool(second and (w>>((second&-second).bit_length()-1)&1))
  f+=a;g+=b;j+=a*b
 return (1,1<<m,f,g,j,(w>>(m-1))&1,w&1)
checks=0
for m in range(1,12):
 for w in range(1<<m):
  f=g=j=e=0
  for k in range(m):
   b=(w>>k)&1
   f,g,j,e=f+b*(1<<k),2*g+b,j+b*(g+e),b
  require((1,1<<m,f,g,j,e,w&1),direct(m,w));checks+=1

def edges(names,suffix,ports,prefix,root):
 d={v:set() for v in names}
 for a,b in itertools.combinations(prefix,2):d[a].add(b)
 for a,b in itertools.combinations(suffix,2):d[a].add(b)
 for a in prefix:d[a].update(suffix)
 d[root].update(prefix)
 for a in suffix:
  d[a if a in ports else root].add(root if a in ports else a)
 return d

def strong(d):
 for v in d:
  seen={v}
  while True:
   new=seen|set.union(*(d[x] for x in seen))
   if new==seen:break
   seen=new
  if len(seen)!=len(d):return False
 return True
retention=0
for t in range(1,9):
 prefix=list(range(10,10+t));root=100;sufs=[[0,1,2],[0,1,2,4],[0,1,2,3,4],[0,1,2,3,4]];ws=[{0,1,2},{0,1,2,4},{0,1,2,4},{0,2,4}];stages=[]
 for su,w in zip(sufs,ws):
  d=edges(su+prefix+[root],su,w,prefix,root);require(strong(d),True);stages.append(d)
 for i in (0,1):
  old,new=stages[i:i+2]
  for a in old:
   for b in old:require(b in old[a],b in new[a])
 old,new=stages[2:];changes=[]
 for a,b in itertools.combinations(old,2):
  if (b in old[a]) != (b in new[a]):changes.append({a,b})
 require(changes,[{1,root}]);retention+=1
res={'binary_recurrence_direct_subset_checks':checks,'physical_named_strong_retention_schedules':retention,'TT3_signature':direct(3,7),'TT5_signature':direct(5,21),'pass':True}
(p/'independent.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res))
