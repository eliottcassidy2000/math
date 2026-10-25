from pathlib import Path
import itertools,json,math
from fractions import Fraction as F
N=6

def req(x,m):
 if not x:raise RuntimeError(m)

def components(n,edges):
 sets=[{i} for i in range(n)]
 for a,b in edges:
  x=next(s for s in sets if a in s);y=next(s for s in sets if b in s)
  if x is not y:sets.remove(x);sets.remove(y);sets.append(x|y)
 return sorted([sorted(s) for s in sets])

cache={}
def data(n,edges):
 key=(n,tuple(edges))
 if key in cache:return cache[key]
 es=set(edges);sm=0;ev=True
 for p in itertools.permutations(range(n)):
  s=(-1)**sum(p[a]>p[b] for a,b in edges);sm+=s
  ys={tuple(sorted((p[a],p[b]))) for a,b in edges}
  if ys==es and s<0:ev=False
 cache[key]=(ev,F(sm,math.factorial(n)));return cache[key]

connected={s:[0]*7 for s in ['all','good','zero_even','odd_edges','even_edges']};direct={'T':[1],'C':[1],'R':[1]};rows=[]
for n in range(1,7):
 E=list(itertools.combinations(range(n),2));ix={e:k for k,e in enumerate(E)};seen=set();ts=cs=rs=0
 for mask in range(1<<len(E)):
  if mask in seen:continue
  es=[e for k,e in enumerate(E) if mask>>k&1]
  for p in itertools.permutations(range(n)):
   seen.add(sum(1<<ix[tuple(sorted((p[a],p[b])))] for a,b in es))
  even,alpha=data(n,es);comps=components(n,es);info=[]
  for verts in comps:
   where={v:i for i,v in enumerate(verts)};ee=[(where[a],where[b]) for a,b in es if a in where];info.append(data(len(verts),ee))
  req(even==all(e for e,a in info),'component signed-admissibility')
  z=F(1)
  for e,a in info:z*=a
  req(alpha==z,'independent component integral')
  if len(es)%2:req(alpha==0,'odd edges zero')
  if not even:continue
  ts+=1;cs+=len(es)%2==0;rs+=alpha!=0
  if len(comps)==1:
   connected['all'][n]+=1;connected['odd_edges' if len(es)%2 else 'even_edges'][n]+=1
   if alpha:connected['good'][n]+=1
   elif len(es)%2==0:connected['zero_even'][n]+=1
  rows.append({'n':n,'mask':mask,'alpha':str(alpha),'component_sizes':[len(c) for c in comps]})
 for k,v in [('T',ts),('C',cs),('R',rs)]:direct[k].append(v)
# Independent bivariate coefficient recurrence, keeping edge parity explicitly.
def multiset(counts):
 a=[[0,0] for _ in range(7)];a[0][0]=1
 for size,par,number in counts:
  for _ in range(number):
   # ascending recurrence for unbounded multiplicities of each type
   for total in range(size,7):
    for parity in range(2):a[total][parity]+=a[total-size][parity^par]
 return a
full=multiset([(n,p,connected['odd_edges' if p else 'even_edges'][n]) for n in range(1,7) for p in range(2)])
good=multiset([(n,0,connected['good'][n]) for n in range(1,7)])
req([sum(x) for x in full]==direct['T'],'all-type Euler')
req([x[0] for x in full]==direct['C'],'converse parity Euler')
req([x[0] for x in good]==direct['R'],'rank Euler')
prod=json.loads((Path(__file__).resolve().parent.parent/'profile-component-series/certificate.json').read_text())
req(connected==prod['connected'],'independent connected census')
for k in direct:req(direct[k]==prod[k],'producer series '+k)
req(connected['good']==[0,1,0,1,1,7,21],'good counts');req(connected['zero_even'][6]==2,'two connected zero types')
out={'connected':connected,'direct':direct,'parity_multiset_coefficients':full,'records':rows}
Path(__file__).with_name('independent.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2))
z=[x[0] for x in multiset([(n,0,connected['zero_even'][n]) for n in range(1,7)])]
o=[x[0] for x in multiset([(n,1,connected['odd_edges'][n]) for n in range(1,7)])]
def convolution(a,b):return [sum(a[j]*b[k-j] for j in range(k+1)) for k in range(7)]
req(z==prod['Z'] and o==prod['O_even'],'zero and odd parity factors')
zo=convolution(z,o);req(convolution(direct['R'],zo)==direct['C'],'C=RZO_even')
zo[0]-=1;req(convolution(direct['R'],zo)==prod['K'],'kernel factorization')
