#!/usr/bin/env python3
"""Independent literal check of the shared0011 anchors; additive scratch."""
import itertools,json
from verify import hp,two_apex,require

records=[];totalpaths=0
for m in range(8,12):
 for mid in itertools.product(range(2),repeat=m-8):
  word=(0,0,1,1)+mid+(0,0,1,1)
  mask=sum(b<<i for i,b in enumerate(word))
  for arc,pa,pb in itertools.product(range(2),repeat=3):
   a=two_apex(m,arc,pa,pb,mask)
   h=hp(a)
   require(h>=3**(m-8),'anchor lower bound')
   # Make the path injection literally and ensure that every word is distinct.
   prefix=list(range(4));suffix=list(range(m-4,m))
   fM=next(i for i in prefix if word[i]==1^pa);prefix.remove(fM)
   fR=next(i for i in prefix if word[i]==1^pb);prefix.remove(fR)
   lL=next(i for i in suffix if word[i]==pa);suffix.remove(lL)
   lM=next(i for i in suffix if word[i]==pb);suffix.remove(lM)
   seen=set()
   for assignment in itertools.product(range(3),repeat=m-8):
    bins=[prefix+[lL],[fM,lM],[fR]+suffix]
    for j,c in enumerate(assignment):bins[c].append(4+j)
    bins=[sorted(q) for q in bins]
    p=tuple([i+2 for i in bins[0]]+[0]+[i+2 for i in bins[1]]+[1]+[i+2 for i in bins[2]])
    require(len(set(p))==m+2,'path covers vertices')
    require(all(a[x]>>y&1 for x,y in zip(p,p[1:])),'path edges')
    seen.add(p)
   require(len(seen)==3**(m-8),'distinct injection')
   totalpaths+=len(seen)
   records.append({'m':m,'middle':list(mid),'arc':arc,'polarities':[pa,pb],'H':h,'lower':len(seen)})
# Quantitative Ramsey counting bound needs no enumeration of huge grids.
# Test each possible histogram for up to7colors, proving finite controls.
histograms=0
for q in range(1,8):
 def comps(n,k):
  if k==1:yield(n,);return
  for a in range(n+1):
   for rest in comps(n-a,k-1):yield(a,)+rest
 for sizes in comps(2*q,q):
  require(sum(s*(s-1)//2 for s in sizes)>=q,'Ramsey per-column pairs')
  histograms+=1
require(3**20==3486784401,'growth integer')
require(6*2**29+1==3221225473,'target integer')
require(3**20>6*2**29+1,'m28 growth')
print(json.dumps({'status':'PASS','graphs':len(records),'literal_injected_paths':totalpaths,
 'histograms':histograms,'records':records},indent=2,sort_keys=True))
