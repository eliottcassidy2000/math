from itertools import product
from pathlib import Path
import json, subprocess, sys
OUT=Path(__file__).resolve().parent;OUT.mkdir(exist_ok=True)
def require(v,msg):
 if not v: raise RuntimeError(msg)
def rank(a,p):
 a=[r[:] for r in a]; k=0
 for j in range(len(a[0])):
  q=next((q for q in range(k,len(a)) if a[q][j]%p),None)
  if q is None:continue
  a[k],a[q]=a[q],a[k]; inv=pow(a[k][j]%p,-1,p);a[k]=[(v*inv)%p for v in a[k]]
  for q in range(len(a)):
   if q!=k:
    c=a[q][j];a[q]=[(u-c*v)%p for u,v in zip(a[q],a[k])]
  k+=1
 return k
H=[(0,1,2),(6,7,8),(0,4,8),(1,4,7),(2,4,6)]
V=[(0,3,6),(2,5,8),(0,4,8),(3,4,5),(2,4,6)]
L=list(dict.fromkeys(H+V)); hc={0,1,2,4,6,7,8}
addr=[(1,0),(-1,-1),(0,1),(-1,1),(0,0),(1,-1),(0,-1),(1,1),(-1,0)]
primes=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59]
rows=lambda ls:[[int(j in l) for j in range(9)]+[-1] for l in ls]
counts=[]
for p in primes:
 require([rank(rows(ls),p) for ls in [H,V,L]]==[5,5,7],f'ranks {p}')
 q={i*i%p for i in range(p)};nh=nf=0;root=[0]*4; distinct_values=[0]*4
 for x,z in product(range(p),repeat=2):
  h={0,x%p,-x%p,z%p,-z%p,(x+z)%p,(-x-z)%p}
  v={0,x%p,-x%p,z%p,-z%p,(x-z)%p,(-x+z)%p}
  require(not(len(h)==7 and len(v)==7) or len(h|v)==9, 'two-copy injectivity')
 # Three actual shared entry values E,A,C are arbitrary field elements.
 for e,a,c in product(range(p),repeat=3):
  core=[e,a,(2*e-a)%p,c,(2*e-c)%p]
  extra=[(a+c-e)%p,(3*e-a-c)%p,(a-c+e)%p,(e-a+c)%p]
  forms=[core,core+extra[:2],core+extra[2:],core+extra]
  r=[(1 if v==0 else 2 if v in q else 0) if p!=2 else 1 for v in core+extra]
  for i,idx in enumerate([range(5),range(7),[0,1,2,3,4,7,8],range(9)]):
   w=1
   for j in idx:w*=r[j]
   root[i]+=w; distinct_values[i]+=all(v in q for v in forms[i])
  if len(set(forms[1]))==7 and all(v in q for v in forms[1]):
   nh+=1
   nf+=len(set(forms[3]))==9 and all(v in q for v in forms[3])
 counts.append(dict(prime=p,distinct_hourglasses=nh,completing_distinct_magic=nf,root_cones=root,entry_parameters=distinct_values))
expected=json.load(open(OUT.parent/'hourglasses/certificate.json'))['finite_field_census']
for got,want in zip(counts,expected):
 for k in ['prime','distinct_hourglasses','completing_distinct_magic']:require(got[k]==want[k],str((k,got,want)))
fiber=json.load(open(OUT.parent/'hourglass-fiber-product/results.json'))['finite_fields']
for want in fiber:
 got=next(c for c in counts if c['prime']==want['p'])
 require(got['root_cones']==want['cone_roots'],'root census')
 require(got['entry_parameters']==want['affine_square_entry_parameters'],'entry census')
for e,x,z in product(range(-4,5),repeat=3):
 a=[e+i*x+j*z for i,j in addr]
 require(all(sum(a[i] for i in line)==3*e for line in L),'labelled map line sums')
 require({addr[i] for i in hc}==set(product([-1,0,1],repeat=2))-{(1,-1),(-1,1)},'exact H label map')
 require({addr[i] for i in {0,2,3,4,5,6,8}}==set(product([-1,0,1],repeat=2))-{(1,1),(-1,-1)},'exact V label map')
# Explicit guard failure must throw even under -O.
try:require(False,'intentional probe')
except RuntimeError as ex:require(str(ex)=='intentional probe','probe exception')
else:raise RuntimeError('failure guard disabled')
report={'status':'PASS','ranks_modular_all_primes':primes,'complete_counts':counts,'explicit_exception_probe':'PASS','exact_label_map':addr}
print(json.dumps(report,sort_keys=True))
