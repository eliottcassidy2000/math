"""Reproducible arithmetic inverse and physical controls. Writes local receipts."""
from pathlib import Path
import json,hashlib
from itertools import combinations
from core import *
from pair_families import *
from inverse import *
P=Path(__file__).parent

def profile(T,r):
 I=inventories(T);return (I['H'],I['e'][r],I['R'][r],I['allH'][((1<<len(T))-1)^(1<<r)]),I['s'][r]
def hc(T):
 require(tournament(T),'invalid witness')
 n=len(T);N=1<<n;A=[[0]*n for _ in range(N)];C=[[0]*n for _ in range(N)]
 for v in range(n):A[1<<v][v]=1
 C[1][0]=1
 for mask in range(1,N):
  bits=mask
  while bits:
   ubit=bits&-bits;bits-=ubit;u=ubit.bit_length()-1
   outs=T[u]&~mask
   while outs:
    bit=outs&-outs;outs-=bit;v=bit.bit_length()-1
    A[mask|bit][v]+=A[mask][u];C[mask|bit][v]+=C[mask][u]
 return sum(A[-1]),sum(C[-1][u]for u in range(n)if T[u]&1)

seed_checks=physical=literal=0
for n in range(1,5):
 pairs=list(combinations(range(n),2))
 for bits in range(1<<len(pairs)):
  V=[0]*n
  for i,(a,b)in enumerate(pairs):V[a if bits>>i&1 else b]|=1<<(b if bits>>i&1 else a)
  I=inventories(V)
  for r in range(n):
   if I['s'][r]!=1:continue
   s=(I['H'],I['e'][r],I['R'][r],I['allH'][((1<<n)-1)^(1<<r)])
   require(s[2]>=s[0]and(n==1 or s[0]-s[1]>=s[3]),'physical input inequalities')
   S=seeds(V,r)
   for kind,q in seed_profiles(s).items():
    if n==1 and kind!='P':continue
    T,a=S[kind];got,start=profile(T,a)
    require(got==q and start==2,'generalized seed profile');seed_checks+=1
   for kind in ('N','M','P','Q'):
    if n==1 and kind!='P':continue
    for v in range(2 if kind=='N'else 1,5):
     T,a=construct(kind,V,r,v);h=scalar(kind,s,v)
     require(hc(T)==(h,2),'physical family H/C')
     require(h>=7*s[0]and v<=parameter_limit(kind,s,h),'inverse cutoffs')
     if len(T)<=7:require(literal_h(T)==h,'literal Hamilton paths');literal+=1
     physical+=1

cores=source_bank(101);core_controls=0
for s,recipe in cores.items():
 T,r=build_source(recipe);q,start=profile(T,r)
 require(q==s and start==1,'physical source bank')
 if s[0]>=3:require(cycles(T)==1,'strong C1 source')
 core_controls+=1

bound=1861;S,inputs=support(bound);C1=source_bank(bound);c1_values={s[0]for s in C1}
missing=sorted(set(range(1,bound+1,2))-c1_values-set(S))
require(missing==[7,21,67],'bounded image')
maxorder=0;positive={}
for h,record in sorted(S.items()):
 T,r=realize(record);require(hc(T)==(h,2),'complete support witness')
 maxorder=max(maxorder,len(T))
 if h in (533,821,1861):positive[h]=dict(record,rows=T,root=r,order=len(T),H=h,C=2)

inverse_checks=0
for h in range(1,502):
 got=inverse(h);require(bool(got)==(h in S),'exact targeted membership')
 for record in got:require(scalar(record['family'],record['source'],record['parameter'])==h,'inverse presentation')
 inverse_checks+=1

# Oversized rectangle includes sources and parameters the target algorithm prunes.
rectangle=0
for s,recipe in source_bank(61).items():
 if s==(1,0,2,1):continue
 for k in range(5):
  D,E,U,d=s;q=s if k==0 else(D,0,U<<k,d)
  for kind in ('N','M','P','Q'):
   if q==(1,1,1,1)and kind!='P':continue
   for v in range(2 if kind=='N'else 1,17):
    h=scalar(kind,q,v)
    require(h>=7*D and v<=parameter_limit(kind,q,h),'oversized exact bounds')
    if h<=501:require(h in S,'oversized source image')
    rectangle+=1

# The erased weak deletion response must genuinely be kept as two states.
require((9,0,18,1)in inputs and(9,0,18,3)in inputs,'weak response collision retained')
for h in (533,821,1861):require(h in S and h not in c1_values,'positive C1 repairs')
receipt=dict(status='PASS',bound=bound,source_quadruples=len(source_bank(bound//7)),candidate_sources=len(inputs),family_output_values=len(S),c1_quadruples=len(C1),union_missing=missing,seed_profiles=seed_checks,small_physical_outputs=physical,literal_path_checks=literal,source_graph_controls=core_controls,complete_support_physical_witnesses=len(S),max_witness_order=maxorder,targeted_inverse_checks=inverse_checks,oversized_rectangle=rectangle,repairs=positive)
(P/'SUPPORT.json').write_text(json.dumps(dict(bound=bound,family_support=S),sort_keys=True,indent=2)+'\n')
(P/'VALIDATION.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
print(json.dumps(receipt,sort_keys=True,indent=2))
