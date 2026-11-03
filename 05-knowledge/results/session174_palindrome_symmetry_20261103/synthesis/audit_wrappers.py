import sys,itertools,json
sys.path.insert(0,'../spectrum')
from common import *
records=[]
for r in range(1,13):
 h=paths(wrap((0,),r))[-1]
 require(h==2**r+1,'singleton exact formula')
 records.append({'r':r,'H':h})
count=0
for n in range(1,5):
 pairs=list(itertools.combinations(range(n),2))
 for bits in range(1<<len(pairs)):
  a=[0]*n
  for k,(i,j) in enumerate(pairs):
   if bits>>k&1:a[i]|=1<<j
   else:a[j]|=1<<i
  p=cover(a);q=[p[1],2*p[2] if n>=2 else 0,6*p[3] if n>=3 else 0]
  x,y,z=q;pred=5*x*x+6*x*y+7*y*y+4*y*z+2*z*z
  h=paths(wrap(a,2))[-1]
  require(h==pred,'positive r2 skeleton count')
  require(h==5 if n==1 else h>=45,'r2 lower bound')
  count+=1
feedback={}
for n in range(1,5):
 vals=set();pairs=list(itertools.combinations(range(n),2))
 for bits in range(1<<len(pairs)):
  z=[0]*n
  for k,(i,j) in enumerate(pairs):
   if bits>>k&1:z[i]|=1<<j
   else:z[j]|=1<<i
  g=[2,((1<<n)-1)<<2]+[(v<<2)|1 for v in z]
  p=cover(z);h=paths(g)[-1]
  require(h==evalF(p,2)+p[1],'arbitrary feedback exact formula')
  require(h not in (7,21),'arbitrary small feedback exclusion')
  require(h>=2**n+1,'universal subset lower bound')
  vals.add(h)
 feedback[n]=sorted(vals)
require(feedback=={1:[3],2:[5],3:[9,15],4:[17,27,37]},'small feedback spectra')
print(json.dumps({'status':'PASS','singleton_cases':records,'labelled_r2_cases':count,'arbitrary_feedback_values':feedback},sort_keys=True,indent=2))
