"""Exact controls for the five-dimensional interface and seven-state reduction."""
from pathlib import Path
from fractions import Fraction
import importlib.util, itertools, json
p=Path(__file__).resolve().parents[2]/'session251_creative_completion_20261111/ternary/compiler.py'
spec=importlib.util.spec_from_file_location('inherited_compiler',p)
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
CORES=[(0,),(2,4,1),(2,4,9,3),(2,4,9,19,7),(2,4,25,3,11)]
def state(sig):
 d,kp,km,ip,im,op,om,jp,jm=sig
 if km!=-d:raise AssertionError(('antipode',sig))
 return (d,d+kp,op,om,ip,im,jp+jm)
def product(a,b):
 d,e,u,um,v,vm,q=a;D,E,U,Um,V,Vm,Q=b
 return (d*D,e*E,d*U+u*E,d*Um,e*V+v*D,vm*D,d*Q+q*D+u*V+um*Vm)
def response(s,a,b):
 d,e,u,um,v,vm,q=s
 return d*(c.g(a)+c.g(b)-1)+e*a*b+v*a-vm*(a%2)+u*b-um*(b%2)+q
def join(a,b):
 n,m=len(a),len(b)
 return tuple(x|(((1<<m)-1)<<n) for x in a)+tuple(x<<n for x in b)
def rank(rows):
 if not rows:return 0
 a=[[Fraction(x) for x in row] for row in rows];r=0
 for col in range(len(a[0])):
  pivot=next((i for i in range(r,len(a)) if a[i][col]),None)
  if pivot is None:continue
  a[r],a[pivot]=a[pivot],a[r];v=a[r][col];a[r]=[x/v for x in a[r]]
  for i in range(r+1,len(a)):
   v=a[i][col];a[i]=[x-v*y for x,y in zip(a[i],a[r])]
  r+=1
 return r
def det(rows):
 a=[[Fraction(x) for x in row] for row in rows];out=Fraction(1)
 for col in range(len(a)):
  pivot=next((i for i in range(col,len(a)) if a[i][col]),None)
  if pivot is None:return 0
  if pivot!=col:a[col],a[pivot]=a[pivot],a[col];out=-out
  v=a[col][col];out*=v
  for i in range(col+1,len(a)):
   q=a[i][col]/v
   a[i]=[x-q*y for x,y in zip(a[i],a[col])]
 if out.denominator!=1:raise AssertionError(out)
 return out.numerator
signed=responses=0
for n in range(1,5):
 for code in range(1<<(n*(n-1)//2)):
  rows=c.tournament(n,code);h=c.hp_subsets(rows);full=(1<<n)-1
  conv=sum((-1)**mask.bit_count()*h[mask]*h[full^mask] for mask in range(full+1))
  if conv:raise AssertionError((n,code,conv))
  for cut in range(1<<n):
   sig=c.signature(rows,cut);st=state(sig);signed+=1
   for a,b in itertools.product(range(4),repeat=2):
    if response(st,a,b)!=c.count(sig,a,b):raise AssertionError((rows,cut,a,b))
    responses+=1
basis=[];basis_records=[];posbasis=[];posrecords=[]
for idx,rows in enumerate(CORES):
 for cut in range(1<<len(rows)):
  sig=c.signature(rows,cut);st=list(state(sig));d,kp,km,ip,im,op,om,jp,jm=sig
  if rank(basis+[st])>len(basis):basis.append(st);basis_records.append({'core':idx,'cut':cut,'state':st})
  positive=[d,d+kp,op,ip,jp]
  if rank(posbasis+[positive])>len(posbasis):posbasis.append(positive);posrecords.append({'core':idx,'cut':cut,'state':positive})
joins=0
for a,b in itertools.product(CORES[:3],repeat=2):
 for ca,cb in itertools.product(range(1<<len(a)),range(1<<len(b))):
  lhs=product(state(c.signature(a,ca)),state(c.signature(b,cb)))
  rhs=state(c.signature(join(a,b),ca|(cb<<len(a))))
  if lhs!=rhs:raise AssertionError((a,b,ca,cb,lhs,rhs))
  joins+=1
report={'scope':'new identities only; no new integer-coverage claim','signed_core_cuts':signed,'response_comparisons':responses,'seven_state_ordinal_compositions':joins,'seven_basis_determinant':det(basis),'positive_five_basis_determinant':det(posbasis),'seven_basis':basis_records,'positive_five_basis':posrecords}
if len(basis)!=7 or len(posbasis)!=5:raise AssertionError('dimension')
print(json.dumps(report,indent=2,sort_keys=True))
