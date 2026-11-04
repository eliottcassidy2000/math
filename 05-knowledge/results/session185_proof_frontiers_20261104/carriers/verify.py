"""Exact fixed-front growing-tail audit; standard library, no asserts."""
import json,itertools
from pathlib import Path
from inherited import direct_signature,tables,strong
ROOT=Path(__file__).parent
def need(b,m):
 if not b:raise RuntimeError(m)
def mul(x,y):
 h,p,f,g,j,e,s=x;H,P,F,G,J,E,S=y
 return h*H,p*P,f*H+p*F,g*P+h*G,j*H+h*J+e*S+g*F,h*E,s*H

def tail(u):
 if not u:return (1,1,0,0,0,0,0)
 z=None
 for i in range(u.bit_length()):
  b=(u>>i)&1;q=(1,2,b,b,0,b,b);z=q if z is None else mul(z,q)
 return z

def graph(r,m,u):
 r=list(r);W={i for i in range(len(r)) if m>>i&1}
 for k in range(u.bit_length()):
  n=len(r);r=[x|(1<<n) for x in r]+[0]
  if u>>k&1:W.add(n)
 return r,W

def attach(r,W,t):
 n=len(r);out=[sum(1<<j for j in range(i+1,t+n)) for i in range(t)]+[x<<t for x in r]+[(1<<t)-1]
 for i in range(n):
  if i in W:out[t+i]|=1<<(t+n)
  else:out[-1]|=1<<(t+i)
 return out

def odd(x):
 while x and x%2==0:x//=2
 return x

# Actual five-vertex SCC whose normalized F has arbitrarily high dyadic valuation.
r=(2,12,17,21,3);m=8;z=direct_signature(r,{3})
need(z==(11,100,12,27,0,1,4),'front signature')
family=[]
for k in range(12):
 u=(63*2**(20*k+8)-3)//25
 need(25*u==63*2**(20*k+8)-3,'integrality')
 zz=mul(z,tail(u));need(zz[2]==63*2**(20*k+10),'unbounded valuation')
 c=11-11*tail(u)[4]-(u&1)-27*u
 need(zz[0]-zz[4]==c and c<=11-27*u,'transported defect')
 family.append({'m':k,'tail_mask':u,'tail_vertices':u.bit_length(),'F':zz[2],'c':c,'J':zz[4],'v2F':20*k+10})
# Full physical check at first member (17 vertices including TT1/root).
u=645;rr,W=graph(r,m,u);need(strong(r),'front strong');need(direct_signature(rr,W)==mul(z,tail(u)),'literal 15-vertex suffix')
a=attach(rr,W,1);need(strong(a),'literal terminal guard');need(tables(a)[0][-1]==family[0]['c']+2*family[0]['F'],'literal full count')
# Exhaustive arbitrary-length transitive tails for every marked front in all strong types <=5.
data=json.loads((ROOT/'bank.json').read_text());checks=0;hits=[]
for i,m,z0 in data['bank']:
 if not m:continue
 z0=tuple(z0);h,p,F,G,J,e,s=z0;need(G>0,'marked front G')
 for c in (-35,-29):
  bound=(h-J-c)//G
  for u in range(0,max(0,bound)+1):
   zz=z0 if u==0 else mul(z0,tail(u));checks+=1
   if zz[0]-zz[4]==c and odd(zz[2])==63:hits.append({'front':i,'mask':m,'tail':u,'sig':zz})
need(not hits,'unexpected fixed-front missing-child return')
# Independent named positive witness checks from closure worker, no bank signature trusted.
witnesses=json.loads((ROOT/'escape_witnesses.json').read_text());audits=[]
for w in witnesses:
 rr=w['rows'];W={i for i in range(len(rr)) if w['mask']>>i&1};zz=direct_signature(rr,W)
 need(list(zz)==w['sig'],'literal escape signature');need(strong(rr),'escape SCC')
 counts=[]
 for t in (1,2,3):
  a=attach(rr,W,t);need(strong(a),'escape attachment guard');H=tables(a)[0][-1];need(H==zz[0]-zz[4]+2**t*zz[2],'escape count');counts.append(H)
 audits.append({'address':w['address'],'counts_at_carrier_clocks_1_2_3':counts,'signature':zz})
result={'status':'PASS','fixed_front_tail_candidates_checked':checks,'fixed_front_returns':hits,'unbounded_valuation_family':family,'literal_witness_audits':audits,'scope':'Marked strong front of order at most five followed by arbitrary transitive tail; not all SCC<=5 component words.'}
(ROOT/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
