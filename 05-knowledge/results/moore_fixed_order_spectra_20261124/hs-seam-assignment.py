from itertools import permutations
from pathlib import Path
from functools import cache
import runpy
x=runpy.run_path(str(Path(__file__).with_name('hs-phase-slice.py')))
make,count=x['make'],cache(x['count'])
S=x['S'];alts=x['alts'];positions=(1,19,20,21,23,24,26)
free=[1<<i for i in positions]
phases=[]
for u0 in permutations(range(1,5)):
 for v0 in permutations((2,3,4)):
  phases.append(((0,*u0),(0,1,*v0)))

def choices(u,v):
 certificates={}
 for which in (None,alts):
  X,m=make(u,v,False,which)
  if m&S:continue
  masks=[m^sum(b for i,b in enumerate(free) if z>>i&1) for z in range(128)]
  heights=[count(8,s) for s in masks]
  for z,h in enumerate(heights):
   for i,b in enumerate(free):
    if z>>i&1:continue
    k=z|(1<<i)
    if abs(heights[k]-h)==2:
     cert=(which,masks[z],masks[k])
     certificates[min(h,heights[k])]=cert
 return certificates

options=[choices(*phase) for phase in phases]
found={}
def add(i,vis):
 for h in sorted(options[i]):
  if h in vis:continue
  vis.add(h)
  if h not in found or add(found[h],vis):
   found[h]=i
   return True
 return False
for i in sorted(range(len(options)),key=lambda x:(len(options[x]),x)):
 add(i,set())
print('two-choice assignment',len(found),'distinct seams',len(set().union(*options)))
assert len(found)==144 and len(set().union(*options))==170
for h,i in found.items():
 order,A,B=options[i][h]
 assert [min(count(8,A),count(8,B)),abs(count(8,A)-count(8,B))]==[h,2]
 assert (A^B).bit_count()==1 and not ((A|B)&S)
 X,_=make(*phases[i],False,order)
 assert x['read'](X,A)==x['read'](X,B)==phases[i]
