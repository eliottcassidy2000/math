"""New exact integral-span calculation; standard library only."""
from fractions import Fraction
from pathlib import Path
from itertools import combinations, product
from math import gcd
import importlib.util,json
p=Path(__file__).resolve().parents[2]/'session251_creative_completion_20261111/ternary/compiler.py'
spec=importlib.util.spec_from_file_location('compiler',p);c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
cores=[(0,),(2,4,1),(2,4,9,3),(2,4,9,19,7),(2,4,25,3,11)]
def state(r,cut):
 d,k,km,i,im,o,om,j,jm=c.signature(r,cut)
 return [d+j+jm,d,d+k,i,o,-im,-om]
def det(rows):
 a=[[Fraction(z) for z in row] for row in rows];r=Fraction(1)
 for j in range(len(a)):
  t=next((i for i in range(j,len(a)) if a[i][j]),None)
  if t is None:return 0
  if t!=j:a[t],a[j]=a[j],a[t];r=-r
  q=a[j][j];r*=q
  for i in range(j+1,len(a)):
   f=a[i][j]/q;a[i]=[u-f*v for u,v in zip(a[i],a[j])]
 if r.denominator!=1:raise AssertionError(r)
 return int(r)
idx=[(0,0),(0,1),(1,0),(1,1),(1,3),(2,0),(2,1),(2,2)]
rows=[state(cores[i],cut) for i,cut in idx]
mins=[det([row for j,row in enumerate(rows) if i!=j]) for i in range(8)]
expected=[-288,0,176,0,0,48,0,0]
if mins!=expected or gcd(*mins)!=16:raise AssertionError(mins)
diffs=[[a-b for a,b in zip(row,rows[0])] for row in rows[1:]]
if det(diffs)!=-160:raise AssertionError(diffs)
affine=[]
for r in cores[3:]:
 extra=state(r,0);ds=diffs+[[a-b for a,b in zip(extra,rows[0])]]
 mm=[det([ds[j] for j in jj]) for jj in combinations(range(8),7)]
 if gcd(*mm)!=32:raise AssertionError(mm)
 affine.append({'extra':extra,'minors':mm,'index':gcd(*mm)})
pts=[(0,0),(1,0),(3,0),(2,0),(0,1),(0,3),(1,1)]
E=[[1,c.g(a)+c.g(b)-2,a*b,a,b,a%2,b%2] for a,b in pts]
if det(E)!=-8:raise AssertionError(E)
checks=0
for row in rows+[state(cores[3],0)]:
 p=[sum(a*b for a,b in zip(e,row)) for e in E]
 if len({z%2 for z in p})!=1 or (p[3]-p[2]+p[1]-3*p[0])%4:raise AssertionError(p)
 h=p[0];i=(p[2]-p[1])//2;o=(p[5]-p[4])//2
 d=(p[3]-p[0]-p[2]+p[1])//2;s=p[1]-h-i;e=p[4]-h-o;L=p[6]-p[1]-p[4]+p[0]
 if [h,d,L,i,o,s,e]!=row:raise AssertionError(row)
 checks+=1
report={'scope':'exact integer and affine spans, not physical realizability','vectors':rows,'linear_maximal_minors':mins,'linear_index':gcd(*mins),'small_parent_affine_index':abs(det(diffs)),'h9_affine_extensions':affine,'observer_determinant':det(E),'observer_index':16*abs(det(E)),'observer_smith':[1,2,2,2,2,2,4],'roundtrips':checks}
print(json.dumps(report,indent=2,sort_keys=True))
