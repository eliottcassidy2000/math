"""Independent exact controls for physical append-only binary-child classification.
No imports of author code or data. This finite checker is not an all-order proof.
"""
import json,functools
from pathlib import Path
ROOT=Path(__file__).resolve().parent
CAP=5000
seed=(1,0,2,1)
def raw(s,c):
 D,E,U,d=s
 if c=='O':return D,0,2*U,d
 if c=='Z':return D-E+U,d,2*U-d,d
 return D-E+U,D-E,D+U,D
def F(s):return s[0]-s[1]+s[2]
def positive(s):
 D,E,U,d=s;return D-E-d,U-D+E-d,E,d-E
def frompos(y):
 a,f,e,p=y;return a+2*e+p,e,a+f+2*e+2*p,e+p
def posstep(y,c):
 a,f,e,p=y
 if c=='O':return a+e,a+2*f+e+2*p,0,e+p
 if c=='Z':return 2*a+f+e+p,f,e+p,0
 return f+p,0,a+e+p,e
def h(s):return s[2]-s[0]+s[1]-2*s[3]
def apply(w):
 s=seed
 for c in w:s=raw(s,c)
 return s
def need(b,msg):
 if not b:raise RuntimeError(msg)
# Symbolic linear identities by checking a basis, including directions not physical.
for i in range(4):
 y=tuple(int(i==j) for j in range(4));s=frompos(y)
 need(F(s)==(2,1,3,3)[i],'F coefficients')
 for c in 'OZJ':
  ns=raw(s,c);need(positive(ns)==posstep(y,c),'transfer equality')
  need(2*F(ns)>=3*F(s),'coefficientwise growth')
  need(F(ns)==(2*F(s)+h(s) if c=='Z' else 2*F(s)-(s[0]-2*s[1])),'child identity')
seen={seed:''};buckets=[set() for _ in range(CAP+1)];buckets[2].add(seed)
for U in range(CAP+1):
 for s in sorted(buckets[U]):
  need(min(positive(s))>=0 and s[3]>=1 and s[0]%2==1,'physical invariants')
  need(F(s)>=s[2]+1,'target cutoff')
  for c in 'OZJ':
   v=raw(s,c);need(v[2]>U,'strict cover growth')
   if v[2]>CAP or v in seen:continue
   seen[v]=seen[s]+c;buckets[v[2]].add(v)
minus_profiles={seed}
j=1
while 2**j+1<=CAP:
 minus_profiles.add((2**j+1,1,2**j+1,1));j+=1
plus_profiles=set();bases=['O','ZO'];j=1
while apply('Z'*j+'JO')[2]<=CAP:bases.append('Z'*j+'JO');j+=1
for w in bases:
 s=apply(w)
 while s[2]<=CAP:plus_profiles.add(s);s=raw(s,'Z')
A1_profiles=set();j=0
while apply('O'*j)[2]<=CAP:A1_profiles.add(apply('O'*j));j+=1
A1_profiles.add(apply('Z'));j=1
while apply('Z'*j+'J')[2]<=CAP:A1_profiles.add(apply('Z'*j+'J'));j+=1
need(minus_profiles=={s for s in seen if h(s)==-1},'negative slice exhaustion')
need(plus_profiles=={s for s in seen if h(s)==1},'positive slice exhaustion')
need(A1_profiles=={s for s in seen if s[0]-2*s[1]==1},'O/J slice exhaustion')
N=10**8
plus=set();minus=set();k=0
while 6*2**k-1<=N:
 plus.add(6*2**k-1)
 if 10*2**k-1<=N:plus.add(10*2**k-1)
 j=1
 while 6*(2**j+1)*2**k-1<=N:plus.add(6*(2**j+1)*2**k-1);j+=1
 k+=1
j=1
while 2**j+1<=N:
 minus.add(2**j+1)
 if 3*(2**j+1)<=N:minus.add(3*(2**j+1))
 j+=1
need(plus&minus=={5,9,17},'arithmetic intersection')
for sign,expected in [(1,plus),(-1,minus)]:
 actual={F(s) for s in seen if any(F(raw(s,c))==2*F(s)+sign for c in 'OZJ') and 5<=F(s)<=CAP+1}
 need(actual=={n for n in expected if 5<=n<=CAP+1},'complete finite child count sets')
# Append-suffix stress controls, independently of slice membership.
suffix_checks=0
for s in seen:
 n=F(s)
 if n<5 or n>101:continue
 for c in 'OZJ':
  for d in 'OZJ':
   need(F(raw(raw(s,c),d))>2*n+1,'two-letter hostile');suffix_checks+=1
# Physical graph compiler and literal DP.
def graph(w):
 rows=[2,0];r=0
 for c in w:
  n=len(rows);rows.append(0)
  for j in range(n):
   if c!='O' and j==r:rows[n]|=1<<j
   else:rows[j]|=1<<n
  if c=='J':r=n
 return rows
def count(rows):
 n=len(rows);full=(1<<n)-1
 need(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(n) for j in range(i)),'tournament')
 @functools.cache
 def f(last,left,closed):
  if left==0:return int(not closed or bool(rows[last]&1))
  bits=rows[last]&left;v=0
  while bits:
   bit=bits&-bits;bits-=bit;v+=f(bit.bit_length()-1,left^bit,closed)
  return v
 return sum(f(i,full^(1<<i),False) for i in range(n)),f(0,full^1,True)
controls=[]
for w,target in [('Z',5),('O',5),('ZZ',9),('ZO',9),('ZZZ',17),('ZJO',17),('OJOO',25),('OJJOO',49),('JZJJZ',49),('OOZZ',51)]:
 need(F(apply(w))==target,'source witness')
 hp,hc=count(graph(w+'J'));need((hp,hc)==(target,1),'literal graph count')
 controls.append({'source_word':w,'count':hp,'cycles':hc,'vertices':len(w)+3})
for minusw,plusw,n in [('Z','O',5),('ZZ','ZO',9),('ZZZ','ZJO',17)]:
 need(F(apply(minusw+'Z'))==2*n-1 and F(apply(plusw+'Z'))==2*n+1,'both controls')
parents25=[s for s in seen if F(s)==25]
child25=sorted({F(raw(s,c)) for s in parents25 for c in 'OZJ'})
need(len(parents25)==2 and child25==[45,63],'25 complete hostile')
need({F(raw(seed,c)) for c in 'OZJ'}=={5},'parent3 boundary')
result={'status':'PASS','scope':'finite exact controls; all-order proof audited separately','source_cover_cap':CAP,'all_physical_profiles':len(seen),'negative_slice_profiles':len(minus_profiles),'positive_slice_profiles':len(plus_profiles),'O_J_minus_slice_profiles':len(A1_profiles),'arithmetic_intersection_cap':N,'both_child_parents':sorted(plus&minus),'two_letter_stress_checks':suffix_checks,'literal_controls':controls,'hostile25_profiles':parents25,'hostile25_one_letter_children':child25,'parent3_one_letter_children':[5],'singleton':'Outside nontrivial source orbit; balanced inverse allows singleton only as a predecessor.'}
print(json.dumps(result,indent=2,sort_keys=True))
