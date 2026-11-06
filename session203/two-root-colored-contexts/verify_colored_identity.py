"""Exact colored endpoint-cover tensors and independent literal graph controls.
No inference of an all-order identity from sampled body tournaments is made:
the finite core/gap tensor is recomputed exactly.
"""
from itertools import permutations,product,combinations
from collections import defaultdict
from pathlib import Path
import json,random
ZERO=(0,0,0,0)

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def add(*polys):
 z=defaultdict(int)
 for p in polys:
  for m,c in p.items():z[m]+=c
 return {m:c for m,c in z.items()if c}
def scale(p,k):return {m:k*c for m,c in p.items()if k*c}
def mul(p,q):
 z=defaultdict(int)
 for a,c in p.items():
  for b,d in q.items():z[tuple(x+y for x,y in zip(a,b))]+=c*d
 return {m:c for m,c in z.items()if c}
ONE={ZERO:1}
X=[]
for i in range(4):
 m=[0]*4;m[i]=1;X.append({tuple(m):1})
x,y,z,w=X
W=add(y,scale(z,-1),mul(x,y),mul(x,w),mul(y,y),mul(y,w))
H=add(ONE,x,y,z,w)
F=add(y,mul(x,y),mul(y,y))
J=add(y,mul(y,y),mul(y,w))
FACT_A=add(mul(add(scale(ONE,4),x),W),mul(add(x,z),F),mul(y,add(scale(ONE,4),scale(x,3),scale(y,4),scale(z,3),scale(w,4))))
FACT_B=add(mul(add(scale(ONE,21),scale(y,42),scale(w,20)),W),scale(mul(add(z,w),J),20),scale(mul(y,H),22))

FAMILIES=[
 dict(name='A',parent=[2,24,3,37,13,23],child=[2,8,67,101,15,87,19],pm=[15,63],cm=[79,127],factor=FACT_A),
 dict(name='B',parent=[32,5,1,39,47,70,31],child=[0,5,129,39,47,199,159,27],pm=[0,5],cm=[0,5],factor=FACT_B),
]

def tensor(core,masks):
 n=len(core);total=defaultdict(int)
 for order in permutations(range(n)):
  p=ONE
  for gap in range(n+1):
   left=order[gap-1]if gap else None
   right=order[gap]if gap<n else None
   factor={}
   if left is None or right is None or core[left]>>right&1:factor[ZERO]=1
   for first,last in product(range(2),repeat=2):
    if(left is None or not(masks[first]>>left&1))and(right is None or masks[last]>>right&1):
     factor=add(factor,X[2*first+last])
   p=mul(p,factor)
   if not p:break
  for m,c in p.items():total[m]+=c
 return dict(total)

def hp_data(a,col):
 n=len(a);N=1<<n;D=[[[0]*2 for _ in range(n)]for _ in range(N)]
 for v in range(n):D[1<<v][v][col[v]]=1
 for S in range(1,N):
  for v in range(n):
   if not any(D[S][v]):continue
   q=a[v]&~S
   while q:
    b=q&-q;q-=b;t=b.bit_length()-1
    for f in range(2):D[S|b][t][f]+=D[S][v][f]
 funcs=[[0]*N for _ in range(4)]
 for S in range(1,N):
  for v in range(n):
   for f in range(2):funcs[2*f+col[v]][S]+=D[S][v][f]
 return funcs

def evaluate(poly,funcs):
 N=len(funcs[0]);cache={ZERO:[1]+[0]*(N-1)}
 def monomial(m):
  if m in cache:return cache[m]
  i=next(i for i,a in enumerate(m)if a);sm=list(m);sm[i]-=1
  old=monomial(tuple(sm));out=[0]*N
  for S in range(1,N):
   A=S
   while A:
    out[S]+=funcs[i][A]*old[S^A];A=(A-1)&S
  cache[m]=out;return out
 return sum(c*monomial(m)[-1]for m,c in poly.items())

def hp(a):
 if not a:return 1
 return sum(v[-1]for v in hp_data(a,[0]*len(a)))

def context(core,masks,body,col):
 n=len(core);a=list(core)+[0]*len(body)
 for j,color in enumerate(col):
  v=n+j
  for i in range(n):
   if masks[color]>>i&1:a[v]|=1<<i
   else:a[i]|=1<<v
  for k in range(j):
   u=n+k
   if body[k]>>j&1:a[u]|=1<<v
   else:a[v]|=1<<u
 return a

def encode(poly):return [{'powers':list(m),'coefficient':c}for m,c in sorted(poly.items())]

def main():
 tensors=[]
 for f in FAMILIES:
  P=tensor(f['parent'],f['pm']);C=tensor(f['child'],f['cm']);D=add(C,scale(P,-2),scale(H,-1))
  need(D==f['factor'],('exact tensor factorization',f['name']))
  tensors.append(dict(family=f['name'],parent=encode(P),child=encode(C),excess=encode(D)))
 inputs=[]
 for n in range(4):
  for code in range(1<<(n*(n-1)//2)):
   body=[0]*n
   for k,(i,j)in enumerate(combinations(range(n),2)):
    if code>>k&1:body[i]|=1<<j
    else:body[j]|=1<<i
   for col in product(range(2),repeat=n):inputs.append((body,col))
 rng=random.Random(203220)
 for n in (4,5):
  for _ in range(12):
   body=[0]*n
   for i,j in combinations(range(n),2):
    if rng.randrange(2):body[i]|=1<<j
    else:body[j]|=1<<i
   inputs.append((body,tuple(rng.randrange(2)for _ in range(n))))
 checks=0;strict=0
 for body,col in inputs:
  funcs=hp_data(body,col);hb=hp(body)
  need(evaluate(H,funcs)==hb,'path series including empty body')
  observedW=hp(context([0],[0,1],body,col))-hb
  need(observedW==evaluate(W,funcs),'independent cut defect polynomial')
  bad=any(col[i]==0 and col[j]==1 and body[i]>>j&1 for i in range(len(body))for j in range(len(body)))
  need(observedW>=0 and (observedW>0)==bad,'cut defect equality boundary')
  for f in FAMILIES:
   parent=hp(context(f['parent'],f['pm'],body,col));child=hp(context(f['child'],f['cm'],body,col))
   excess=child-2*parent-hb
   need(excess==evaluate(f['factor'],funcs),'literal graph versus factor tensor')
   need(excess>=0 and (excess>0)==bad,'exact equality boundary')
   checks+=1;strict+=excess>0
 out=dict(status='EXACT_TENSOR_AND_LITERAL_CONTROLS_PASS',colored_inputs=len(inputs),graph_pairs=checks,strict_cases=strict,
          interpretation='Powers refer to commutative subset convolution of endpoint-color path counts00,01,10,11; constant1 is supported only on the empty body',
          tensors=tensors)
 print(json.dumps(out,indent=2))
 Path(__file__).with_name('VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
