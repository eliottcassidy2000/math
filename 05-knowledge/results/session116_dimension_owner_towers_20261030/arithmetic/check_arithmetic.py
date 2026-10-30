"""Exact independent arithmetic/convolution checks. No external dependencies."""
from itertools import product
from math import gcd
import json
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok: raise RuntimeError(msg)
def translate(mask,v,N):
 out=0
 for x in range(N):
  if mask>>x&1:out|=1<<(x^v)
 return out
def convolution(a,b,N,cache=None):
 out=0
 while a:
  bit=a&-a;x=bit.bit_length()-1;a^=bit
  out^=cache[b,x] if cache is not None else translate(b,x,N)
 return out
def phase(d,v,r,c):
 N=1<<d;epsilon=c^((r&v).bit_count()&1)
 return sum((epsilon^((r&x).bit_count()&1)^(x==v))<<x for x in range(N))
def normdef(f,x,y):
 N=len(f);xy=[0]*N
 for a in range(N):
  for b in range(N):xy[a^b]+=(-1 if f[a][b] else 1)*x[a]*y[b]
 return sum(z*z for z in xy)-sum(z*z for z in x)*sum(z*z for z in y)
def ideal(f):
 N=len(f);g=0;witness=None
 for t in range(1,N):
  for a in range(N):
   if a>(a^t):continue
   for b in range(N):
    if b>(b^t):continue
    coeff=2*((-1)**(f[a][b]+f[a^t][b^t])+(-1)**(f[a][b^t]+f[a^t][b]))
    need(coeff in (-4,0,4),'coefficient divisibility')
    g=gcd(g,coeff)
    if coeff and witness is None:
     x=[0]*N;y=[0]*N;x[a]=x[a^t]=y[b]=y[b^t]=1
     need(normdef(f,x,y)==coeff,'sharp sparse witness');witness=(x,y)
 need(g in (0,4),'exact defect ideal')
 return g
report={}
# ALL normalized 4-dimensional sign arrays, not just Fano line arrays.
counts={0:0,4:0}
for bits in product(range(2),repeat=9):
 f=[[0]*4 for _ in range(4)]
 for (a,b),bit in zip(product(range(1,4),repeat=2),bits):f[a][b]=bit
 counts[ideal(f)]+=1
need(counts=={0:2,4:510},'all normalized dimension-four classification')
report['normalized_dim4_defect_ideals']=counts
# Entire dimension-eight line family, plus arithmetic squareclass bridge.
lines=sorted({tuple(sorted((a,b,a^b))) for a in range(1,8) for b in range(1,8) if a!=b})
from math import prod,isqrt
coherent=0
for bits in product(range(2),repeat=7):
 f=[[0]*8 for _ in range(8)]
 for a in range(1,8):f[a][a]=1
 for (u,v,w),z in zip(lines,bits):
  for a,b in ((u,v),(v,w),(w,u)):f[a][b]=z;f[b][a]=z^1
 g=ideal(f);coherent+=g==0
 sigma=prod(p+1 for p,z in zip((19,2,59,29,5,89,7),bits) if z)
 square=sigma%30==0 and isqrt(sigma//30)**2==sigma//30
 need((g==0)==square,'arithmetic bank defect ideal')
need(coherent==16,'coherent count');report['dim8_line_tables']=128
for d in (3,4):
 N=1<<d;params=list(product(range(N),range(N),range(2)))
 masks={p:phase(d,*p) for p in params}
 need(len(set(masks.values()))==2*N*N,'unique phase parameters')
 cache={(m,x):translate(m,x,N) for m in masks.values() for x in range(N)}
 for (v,r,c),a in masks.items():
  need(a.bit_count()%4==(1 if c==0 else 3),'arithmetic unit character')
  need(convolution(a,a,N,cache)==1,'unit involution')
  spectrum=[sum((-1)**(((a>>x)&1)^((u&x).bit_count()&1)) for x in range(N)) for u in range(N)]
  need(max(map(abs,spectrum))==N-2,'extreme Walsh family')
  for (w,s,e),b in masks.items():
   need(convolution(a,b,N,cache)==masks[v^w,r^s,c^e],'additive group law')
 if d==3:need(set(masks.values())=={a for a in range(1<<N) if a.bit_count()%2},'all units at d3')
 else:need(7 not in masks.values(),'three-point hostile at d4')
 report['phase_group_d'+str(d)]={'size':len(masks),'pair_products':len(masks)**2,'index_in_units':1<<(N-2*d-2)}
# Universal weight norm: all group-algebra pairs at dimension8, including nonunits.
N=8;cache={(m,x):translate(m,x,N) for m in range(1<<N) for x in range(N)}
for a in range(1<<N):
 for b in range(1<<N):need(convolution(a,b,N,cache).bit_count()%4==a.bit_count()*b.bit_count()%4,'all parity convolution pairs')
report['all_dim8_parity_pairs']=65536
# Sharp boundary: affine functions no longer have square-zero ideal at d2.
a=sum(((x&1)!=0)<<x for x in range(4));b=sum(((x&2)!=0)<<x for x in range(4))
need(convolution(a,b,4)!=0,'d2 square-zero failure')
# Universal moduli of integer Hamming weight on the entire unit group.
unit_gcd={}
for d in (0,1,2,3):
 N=1<<d;units=[a for a in range(1<<N) if a.bit_count()%2];g=0
 for a in units:
  for b in units:g=gcd(g,convolution(a,b,N).bit_count()-a.bit_count()*b.bit_count())
 need(g==({0:0,1:0,2:8,3:4}[d]),'unit-weight defect gcd')
 unit_gcd[str(d)]=g
need(convolution(7,19,8).bit_count()==5,'mod8 hostile weight5')
report['unit_weight_defect_gcd']=unit_gcd
report.update(status='PASS',checks=checks)
print(json.dumps(report,indent=2,sort_keys=True))
