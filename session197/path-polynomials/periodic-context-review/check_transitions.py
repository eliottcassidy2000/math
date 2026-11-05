import itertools,json,math
from pathlib import Path
def need(x,s):
 if not x: raise RuntimeError(s)
def trim(a):
 while len(a)>1 and not a[-1]:a.pop()
 return a
def add(a,b,c=1):
 d=list(a)+[0]*max(0,len(b)-len(a))
 for i,x in enumerate(b):d[i]+=c*x
 return trim(d)
def mul(a,b):
 d=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):d[i+j]+=x*y
 return trim(d)
def shift(a,s):
 d=[0]*(len(a))
 for i,x in enumerate(a):
  for j in range(i+1):d[j]+=x*math.comb(i,j)*s**(i-j)
 return trim(d)
def X(a):return [0]+a if any(a) else [0]
def D(a):return add(X(X(a)),X(add(shift(a,1),shift(a,-1),-1)))
def basis(k,p):
 a=[0,1] if p else [0,0,1]
 for j in range(1,k+1):a=mul(a,[-j*j,0,1])
 return a
def coords(a,p):
 a=list(a);c={};lo=1 if p else 2
 while len(a)-1>=lo:
  k=(len(a)-1-lo)//2;b=basis(k,p);need(len(a)==len(b),'parity');c[k]=a[-1];a=add(a,b,-a[-1])
 need(not any(a),'constant remainder');return c
bands=0
for p in (0,1):
 for k in range(31):
  q=basis(k,p)
  want={k+1:1,k:k*k+6*k+(3 if p else 5)}
  if k:want[k-1]=2*k*(k if p else k+1)*(2*k+1)
  need(coords(D(q),p)==want,'exact triangle bands');bands+=1
  xx=coords(X(q),1-p)
  need(xx==({k:1} if p else {k+1:1,k:(k+1)**2}),'cross parity padding');bands+=1
transitions=0
for length in range(1,7):
 for word in itertools.product('DX',repeat=length):
  if 'D' not in word:continue
  for p in (0,1):
   for k in range(9):
    a=basis(k,p)
    for c in word:a=D(a) if c=='D' else X(a)
    v=coords(a,(p+word.count('X'))%2)
    need(v.get(k,0)>=(k+1)**2 and v.get(k+1,0)>=1 and (not k or v.get(k-1,0)>=1),'period three paths');transitions+=1
long=0
# Nonphysical positive states are deliberately included: domination is a cone theorem.
for word in ('D','DX','XD','XXD','DXDX','DXXD','XXDDX'):
 for p in (0,1):
  for vacuum in (1,3):
   a=add(basis(0,p),basis(2,p),2)
   a=[vacuum*x for x in a]
   for m in range(1,21):
    for c in word:a=D(a) if c=='D' else X(a)
    h=sum(a)
    for k in range(m//2+1):need(h>=vacuum*(k+1)**(2*m-4*k),'all path bounds')
    need(h>=vacuum*(1+m//4)**m,'claimed uniform bound');long+=1
p=Path(__file__).resolve().parent
source=json.loads((p/'source-normal.json').read_text())
physical=json.loads((p/'independent-normal.json').read_text())
need([x['counts'] for x in source['physical_context_words']]==physical['counts'],'independent literal word counts')
print(json.dumps({'status':'PASS','exact_bands':bands,'period_transitions':transitions,'long_period_states':long,'literal_histories_compared':len(physical['counts'])},sort_keys=True))
