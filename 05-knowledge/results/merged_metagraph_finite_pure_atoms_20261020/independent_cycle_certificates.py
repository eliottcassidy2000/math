#!/usr/bin/env python3
"""Check an exhaustive circulant witness list without importing its producer."""
import json,sys,itertools as it,math,gzip
from pathlib import Path
OUT=Path(__file__).resolve().parent
source=Path(sys.argv[1]) if len(sys.argv)>1 else OUT/'cycle_certificates.json'
data=json.loads(source.read_text());rows=data['certificates']
by_n={}
for r in rows:by_n.setdefault(r['n'],[]).append(r)
for extended in source.parent.glob('circulant_n*_certificates.json.gz'):
    packet=json.loads(gzip.decompress(extended.read_bytes()))
    n=packet['n'];assert n>=17 and n%2==1 and n not in by_n
    assert packet['presentations']==packet['excluded']==2**((n-1)//2)
    by_n[n]=[dict(n=n,mask=r['mask'],certificate=r['certificate']) for r in packet['rows']]
assert sorted(by_n)==list(range(3,max(by_n)+1,2))
checked=[]
for n,rr in sorted(by_n.items()):
 m=(n-1)//2
 assert sorted(r['mask'] for r in rr)==list(range(2**m))
 for row in rr:
  mask=row['mask'];S={d if mask&(1<<(d-1)) else n-d for d in range(1,m+1)}
  if 'S' in row:assert S==set(row['S'])
  assert len(S)==m and S.isdisjoint({n-s for s in S})
  def arc(i,j):return i!=j and (j-i)%n in S
  cert=row['certificate']
  if n<7:
   assert cert is None
   continue
  C=cert['cycle'];P=cert['bad_cut_path']
  assert sorted(C)==list(range(n))==sorted(P)
  assert all(arc(C[i],C[(i+1)%n]) for i in range(n))
  assert any(P==C[k:]+C[:k] for k in range(n))
  suc={C[i]:C[(i+1)%n] for i in range(n)}
  rev={P[i]:P[n-1-i] for i in range(n)}
  a,b=cert['successor_bad_pair'];x,y=cert['reversal_bad_pair']
  assert arc(a,b)!=arc(suc[a],suc[b])
  assert arc(x,y)!=arc(rev[y],rev[x])
 checked.append({'n':n,'presentations':len(rr),'negative_certificates':len(rr) if n>=7 else 0})
# Independent bottom-up subset DP (producer uses recursive memoized DFS).
carousel=[]
for n in range(3,16,2):
 S=set(range(1,(n+1)//2));size=2**(n-1)
 D=[[0]*n for _ in range(size)];D[0][0]=1
 for mask in range(size):
  for u in range(n):
   if not D[mask][u]:continue
   for v in range(1,n):
    if not mask&(1<<(v-1)) and (v-u)%n in S:
     D[mask|1<<(v-1)][v]+=D[mask][u]
 HC=sum(D[-1][u] for u in range(1,n) if (-u)%n in S)
 carousel.append({'n':n,'HC':HC})
assert carousel==data['carousel_HC']
# Sharpened effective cutoff, all comparisons exact integers.
bases=[]
for n in [133,135,137]:
 k=(n-7)//6;L=2*4**k*math.factorial(k);U=3**((n-1)//2)-1
 assert L>U and 4*(k+1)>27
 bases.append({'n':n,'k':k,'double_cycle_lower':L,'group_minus_one_upper':U})
extra=[n for n in range(7,133,2) if 2*4**((n-7)//6)*math.factorial((n-7)//6)>3**((n-1)//2)-1]
assert extra==[121,127,129]
# Explicit carousel obstruction for a large check range; all-order proof is algebraic.
for m in range(3,201):
 n=2*m+1;C=list(range(m-1))+list(range(m,2*m-1))+[2*m,m-1,2*m-1]
 assert sorted(C)==list(range(n))
 assert set((C[(i+1)%n]-C[i])%n for i in range(n))<={1,2,m}
 suc={C[i]:C[(i+1)%n] for i in range(n)}
 assert 1<=m-1<=m and (suc[m-1]-suc[0])%n==2*m-2>m
result={'status':'PASS','complete_presentations_checked':checked,'independent_carousel_HC':carousel,'sharpened_cutoff_bases':bases,'additional_excluded_orders':extra,'maximum_remaining_order':131}
(OUT/'independent_cycle_certificates.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result),flush=True)
