"""Independent literal-word truth tables: no DP or polynomial path products."""
from itertools import permutations
import json,hashlib
from pathlib import Path

def require(x,s):
 if not x:raise ValueError(s)
def truth(rows):
 n=len(rows);N=2*n;allbits=(1<<N)-1;values=[0]*(1<<N)
 for word in permutations(range(n+2)):
  ones=zeros=0
  for u,v in zip(word,word[1:]):
   if u<n and v<n:
    if not rows[u]>>v&1:break
   elif u>=n and v>=n:
    if u!=n:break
   elif u<n:ones|=1<<((v-n)*n+u)
   else:zeros|=1<<((u-n)*n+v)
  else:
   require(not ones&zeros,'contradiction in path constraints')
   free=allbits^(ones|zeros);s=free
   while True:
    values[ones|s]+=1
    if not s:break
    s=(s-1)&free
 return values
T=(24,17,19,6,8);N=(8,17,19,6,9)
require(T[0]>>4&1 and not N[0]>>4&1,'actual flipped pair')
require(all((T[u]>>0&1)==(T[u]>>4&1) for u in [1,2,3]),'twins')
a=truth(T);b=truth(N);diff=[x-y for x,y in zip(a,b)];mob=diff[:]
for bit in range(10):
 for mask in range(1024):
  if mask>>bit&1:mob[mask]-=mob[mask^(1<<bit)]
require(mob[586]==2,'quartic witness')
expected={'0':0,'1':4,'2':20,'3':56,'4':32}
counts={str(k):sum(bool(v) and i.bit_count()==k for i,v in enumerate(mob)) for k in range(11)}
require(all(counts[str(k)]==expected[str(k)] for k in range(5)),'degree census')
require(not any(counts[str(k)] for k in range(5,11)),'degree bound')
source=json.load(open('/workspace/shared/session186/structure/normal.json'))
for row in source['strong_hostile']['queries']:
 m=row['mask'];require(a[m]==row['F_T'] and b[m]==row['F_N'],'16-query literal values')
Path('/workspace/shared/session186/structure_review/full_coefficients.json').write_text(json.dumps({'T':T,'N':N,'nonzero':[[i,v] for i,v in enumerate(mob) if v],'bit_order':'a0,...,a4,c0,...,c4; omitted entries zero'},indent=2)+'\n')
print(json.dumps({'status':'PASS','literal_permutations_per_tournament':5040,'truth_queries':2048,'coefficient_586':mob[586],'degree_counts':counts,'truth_sha256':hashlib.sha256(json.dumps([a,b],separators=(',',':')).encode()).hexdigest()},indent=2))
# Independently repeat the finite minimum-order census using literal truth tables.
from itertools import combinations
def skew_via_subsets(rows,removed,z,j):
 vertices=[v for v in range(len(rows)) if v!=removed]
 def H_end(ids,last):
  return sum(p[-1]==last and all(rows[u]>>v&1 for u,v in zip(p,p[1:])) for p in permutations(ids))
 def H_start(ids,first):
  return sum(p[0]==first and all(rows[u]>>v&1 for u,v in zip(p,p[1:])) for p in permutations(ids))
 out=0
 for mask in range(1,1<<len(vertices)):
  A=[v for k,v in enumerate(vertices) if mask>>k&1];B=[v for k,v in enumerate(vertices) if not mask>>k&1]
  if not B:continue
  if z in A and j in B:out+=H_end(A,z)*H_start(B,j)
  if j in A and z in B:out-=H_end(A,j)*H_start(B,z)
 return out
hist={};skew_checks=0
for q in range(1,5):
 h={}
 for code in range(1<<(q*(q-1)//2)):
  rows=[0]*q
  for k,(u,v) in enumerate(combinations(range(q),2)):
   win,lose=(u,v) if code>>k&1 else (v,u);rows[win]|=1<<lose
  a=rows+[0]
  for u in range(1,q):
   if rows[0]>>u&1:a[q]|=1<<u
   else:a[u]|=1<<q
  a[0]|=1<<q;b=a[:];b[0]^=1<<q;b[q]|=1
  t=truth(a);s=truth(b);c=[x-y for x,y in zip(t,s)]
  for k in range(2*(q+1)):
   for mask in range(len(c)):
    if mask>>k&1:c[mask]-=c[mask^(1<<k)]
  for i,j in permutations(range(1,q),2):
   mon=(1<<i)|(1<<j)|(1<<(q+1+i))|(1<<(q+1+q))
   require(c[mon]==skew_via_subsets(rows,i,0,j),'independent deletion-cover skew')
   skew_checks+=1
  degree=max(mask.bit_count() for mask,v in enumerate(c) if v)
  h[str(degree)]=h.get(str(degree),0)+1
 hist[str(q)]=h
require(hist=={'1':{'1':1},'2':{'3':2},'3':{'3':8},'4':{'4':64}},'complete finite minimum-order census')
Path('/workspace/shared/session186/structure_review/census.json').write_text(json.dumps({'status':'PASS','independent_truth_table_census':hist,'deletion_skew_checks':skew_checks},indent=2)+'\n')
