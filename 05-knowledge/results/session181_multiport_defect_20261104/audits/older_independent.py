"""Independent audit of the defect-slot module packet."""
import itertools,math,json
from functools import lru_cache
from pathlib import Path
OUT=Path(__file__).resolve().parent

def need(v,m):
 if not v:raise RuntimeError(m)
def h(a):
 @lru_cache(None)
 def start(v,mask):
  rest=mask^(1<<v)
  return 1 if not rest else sum(start(j,rest) for j in range(len(a)) if rest>>j&1 and a[v]>>j&1)
 return sum(start(v,(1<<len(a))-1) for v in range(len(a)))
def covers(a):
 f=[0]*(len(a)+1)
 for p in itertools.permutations(range(len(a))):
  for cuts in range(1<<(len(a)-1)):
   if all(cuts>>i&1 or a[p[i]]>>p[i+1]&1 for i in range(len(a)-1)):f[1+cuts.bit_count()]+=1
 return f
def inflate(q,bs):
 offsets=[0]
 for b in bs:offsets.append(offsets[-1]+len(b))
 a=[0]*offsets[-1]
 for i,b in enumerate(bs):
  for v,row in enumerate(b):
   a[offsets[i]+v]=row<<offsets[i]
   for j in range(len(bs)):
    if q[i]>>j&1:a[offsets[i]+v]|=((1<<len(bs[j]))-1)<<offsets[j]
 return tuple(a)
def predict(q,r,b):
 f=covers(b);ans=0
 for p in itertools.permutations([v for v in range(len(q)) if v!=r]):
  eligible=[];mandatory=[]
  for t in range(len(p)+1):
   good=(t==0 or q[p[t-1]]>>r&1) and (t==len(p) or q[r]>>p[t]&1)
   if good:eligible.append(t)
   if 0<t<len(p) and not q[p[t-1]]>>p[t]&1:mandatory.append(t)
  if not set(mandatory)<=set(eligible):continue
  g,d=len(eligible),len(mandatory)
  ans+=sum(math.comb(g-d,k-d)*f[k] for k in range(max(1,d),min(g,len(b))+1))
 return ans
def tournament(n,c):
 a=[0]*n
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  if c>>k&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return tuple(a)
def tt(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
checks=0
for n in range(2,5):
 for code in range(1<<(n*(n-1)//2)):
  q=tournament(n,code)
  for r in range(n):
   for b in [(0,),tt(2),(2,4,1),(8,1,3,6)]:
    bs=[(0,)]*n;bs[r]=b
    need(h(inflate(q,bs))==predict(q,r,b),'root substitution')
    checks+=1
# Full arbitrary multi-module triangle, where two or more module interiors change.
modules=[(0,),tt(2),(2,4,1),(8,1,3,6)]
bounds=0
for bs in itertools.product(modules,repeat=3):
 N=h(inflate((2,4,1),bs))
 need(all(N>=2**len(b)+1 for b in bs),'arbitrary-module bound')
 bounds+=1
# Characteristic product is checked against literal counts, not the slot evaluator.
rec=0
for q in [(2,4,1),(8,1,3,6)]:
 n=len(q);coeff=[1]
 for a in range(1,n+1):
  nxt=[0]*(len(coeff)+1)
  for j,c in enumerate(coeff):nxt[j]-=a*c;nxt[j+1]+=c
  coeff=nxt
 for r in range(n):
  vals=[]
  for k in range(1,n+3):
   bs=[(0,)]*n;bs[r]=tt(k);vals.append(h(inflate(q,bs)))
  for k in range(2):need(sum(coeff[j]*vals[k+j] for j in range(n+1))==0,'recurrence');rec+=1
receipt={'status':'PASS','independent_root_module_checks':checks,'arbitrary_multimodule_triangle_bounds':bounds,'literal_count_recurrence_checks':rec,'scope_correction_requested':'Section3 strong quotient must have order at least3'}
(OUT/'older_independent.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
