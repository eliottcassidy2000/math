"""Independent strong n<=7 catalogue from a fixed directed cycle.
Canonical code compares adjacency rows lexically within invariant degree cells.
"""
from kernel import *
from itertools import product
from collections import Counter
import json,hashlib,pathlib

def canonical(a):
 cells={}
 for v,row in enumerate(a):cells.setdefault(row.bit_count(),[]).append(v)
 perms=[list(permutations(cells[d]))for d in sorted(cells)]
 return min(induced(a,tuple(v for block in choice for v in block))for choice in product(*perms))

levels=[]
for n in range(3,8):
 fixed={(i,(i+1)%n)for i in range(n)}
 rest=[(i,j)for i in range(n) for j in range(i+1,n)if(i,j)not in fixed and(j,i)not in fixed]
 base=[sum(1<<v for u,v in fixed if u==i)for i in range(n)]
 classes=set()
 for bits in product((0,1),repeat=len(rest)):
  a=base[:]
  for bit,(i,j)in zip(bits,rest):
   if bit:a[i]|=1<<j
   else:a[j]|=1<<i
  classes.add(canonical(a))
 hist=Counter();hvals=Counter()
 for a in classes:
  h=h_recursive(a);ds=[h_recursive(induced(a,[u for u in range(n) if u!=v]))for v in range(n)]
  require(max(ds)<h,('strict deletion',n,a));hist[max(ds)]+=1;hvals[h]+=1
 levels.append({'n':n,'strong_classes':len(classes),'birth_histogram':dict(sorted(hist.items())),'H_histogram':dict(sorted(hvals.items()))})
 require(len(classes)=={3:1,4:1,5:6,6:35,7:353}[n],('class count',n))
require(levels[-1]['birth_histogram'][37]==104,'birth37')
atoms=[]
for a in ((16,9,3,5,14),(24,17,3,6,12)):
 ds=[h_literal(induced(a,[u for u in range(5)if u!=v]))for v in range(5)]
 require(h_literal(a)==15 and max(ds)==5,'H15 atoms')
 atoms.append({'rows':a,'H':15,'d':5,'deletions':ds,'scores':sorted(r.bit_count()for r in a)})
require(atoms[0]['scores']!=atoms[1]['scores'],'nonisomorphism')
for B in range(5,100001):
 p=1;r=0
 while p*15<=3*B:p*=15;r+=1
 require(r>=1,'r positive')
 while p<=B:p*=3
 require(B<p<=3*B and p//3<=B,'uniform triangle scaling')
print(json.dumps({'status':'PASS','levels':levels,'H15_atoms':atoms,'uniform_scale_B_range':[5,100000],'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))
