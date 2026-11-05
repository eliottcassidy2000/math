from pathlib import Path
from itertools import combinations,product,permutations
from functools import lru_cache
import hashlib,json
P=Path(__file__).resolve().parent
def need(x,s):
 if not x:raise RuntimeError(s)
@lru_cache(None)
def g(u):
 if u==1:return 1
 if u%2==0:return g(u//2)
 l=(u+1)//2;return g(l)+l-l%2
def arrays(a):
 n=len(a);N=1<<n;st=[[0]*n for _ in range(N)];en=[[0]*n for _ in range(N)];h=[1]+[0]*(N-1)
 for S in range(1,N):
  for v in range(n):
   if not S>>v&1:continue
   R=S^(1<<v)
   if not R:st[S][v]=en[S][v]=1;continue
   for w in range(n):
    if R>>w&1:
     if a[v]>>w&1:st[S][v]+=st[R][w]
     else:en[S][v]+=en[R][w]
  h[S]=sum(st[S]);need(h[S]==sum(en[S]),'start/end total')
 return st,en,h
def profile(a,r,arr=None):
 st,en,h=arrays(a)if arr is None else arr;F=(1<<len(a))-1
 R=sum(st[S][r]*h[F^S]for S in range(1,F+1)if S>>r&1)
 return h[F],en[F][r],R,st[F][r],h[F^(1<<r)]
def replace_triangle(a,r):
 n=len(a);b=a[:]+[0,0]
 for v in range(n):
  if v==r:continue
  if a[r]>>v&1:b[n]|=1<<v;b[n+1]|=1<<v
  else:b[v]|=(1<<n)|(1<<(n+1))
 b[r]|=1<<n;b[n]|=1<<(n+1);b[n+1]|=1<<r
 return b
def cycle_count(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)];d[1][0]=1
 for S in range(1,1<<n):
  for u in range(n):
   if d[S][u]:
    q=a[u]&~S
    while q:
     b=q&-q;q-=b;d[S|b][b.bit_length()-1]+=d[S][u]
 return sum(d[-1][v]for v in range(n)if a[v]&1)
def physical(word):
 a=[2,0]
 for mask in word:
  n=len(a);a.append(0)
  for v in range(n):
   if v<2 and mask>>v&1:a[n]|=1<<v
   else:a[v]|=1<<n
 return a
def word(u,v):return [1 if x=='1'else 3 for x in bin(u-1)[2:]]+[2,3]+[(v-1)>>i&1 for i in range((v-1).bit_length())]
def E(u,v):
 D=g(u);p=u%2;t=v%2
 return (D+3*u-p-2)*g(v)+(2*D+8*u-p-4)*v-u*(1+t)-p
# Every ordered assignment of triangle vertices to three path slots.
a=[2,4,1];partition_counts=[0]*4
for w in permutations(range(3)):
 for i in range(4):
  for j in range(i,4):
   A,B,C=w[:i],w[i:j],w[j:]
   if not all(all(a[x]>>y&1 for x,y in zip(s,s[1:]))for s in (A,B,C)):continue
   if B and B[0]==0:
    partition_counts[0]+=1
    if not C:partition_counts[1]+=1
   if A and C and C[0]==0:
    partition_counts[2]+=1
    if not B:partition_counts[3]+=1
need(partition_counts==[7,3,4,2],'literal module partition constants')
general=0;start_checks=0
for n in range(2,6):
 edges=list(combinations(range(n),2))
 for bits in range(1<<len(edges)):
  a=[0]*n
  for k,(u,v)in enumerate(edges):a[u if bits>>k&1 else v]|=1<<(v if bits>>k&1 else u)
  arr=arrays(a)
  for r in range(n):
   if a[r].bit_count()!=1:continue
   b=a[r].bit_length()-1;h,ea,Ra,sa,d=profile(a,r,arr);_,eb,Rb,_,_=profile(a,b,arr)
   out=replace_triangle(a,r);hp,ep,Rp,sp,dp=profile(out,r)
   need((hp,ep,Rp,dp)==(3*h+3*Rb-3*eb,ea+Rb-eb,2*h+Ra+4*Rb-2*eb,h+Rb-eb),'general triangle law')
   if n>=3:need(sp==3*sa,'general start law');start_checks+=1
   general+=1
need(general==1390,'complete labelled incidence universe')
family=0
for u in range(2,18):
 D=g(u);p=u%2
 seed=physical(word(u,1));h,e,R,s,d=profile(seed,0)
 need((h,e,R,d)==(3*D+9*u-3*p-6,D+2*u-p-2,3*D+11*u-2*p-6,D+3*u-p-2),'physical gate seed')
 need(s==cycle_count(seed)==3,'seed start/cycles')
 for v in range(1,18):
  a=physical(word(u,v));h,e,R,s,d=profile(a,0)
  need(h==E(u,v) and s==3 and cycle_count(a)==3,'physical family')
  need(len(a)==4+(u-1).bit_length()+(v-1).bit_length(),'order')
  need(all(a[i]>>j&1 for i in range(2,len(a))for j in range(i+1,len(a))),'fixed port feedback bound')
  family+=1
examples=[]
for u,v,target in [(4,2,67),(8,7,533),(6,31,1861)]:
 a=physical(word(u,v));h,e,R,s,d=profile(a,0);need(h==target==E(u,v) and cycle_count(a)==3,'positive escape')
 examples.append({'u':u,'v':v,'H':h,'C':3,'word':word(u,v),'rows':a,'order':len(a)})
# The 533 example is literally the frozen positive witness in its birth labeling.
original=[2,268,905,881,999,967,903,779,513,3];birth=[1,3,4,5,6,2,7,8,9,0]
need(examples[1]['rows']==[sum(1<<j for j,v in enumerate(birth)if original[u]>>v&1)for u in birth],'literal 533 relabeling')
cap=5001;image={};trials=0
for u in range(2,(cap-1)//7+1):
 for v in range(1,(cap-1)//(7*u)+1):
  x=E(u,v);need(x>=7*u*v+1,'inverse lower bound');trials+=1
  if x<=cap:image.setdefault(x,[u,v])
need(len(image)==1581 and 821 not in image,'exact bounded family image')
stream=''.join(str(x)+'\n'for x in sorted(image)).encode()
out={'status':'PASS','general_labelled_outdegree_one_cases':general,'general_start_cases':start_checks,'module_partition_counts':partition_counts,'physical_family_graphs':family,'examples':examples,'image_ceiling':cap,'complete_inverse_candidate_pairs':trials,'distinct_family_values':len(image),'image_sha256':hashlib.sha256(stream).hexdigest(),'821_in_family':False,'scope':'constructive C3 feedback-at-most-two family, not full class'}
s=json.dumps(out,sort_keys=True,indent=2)+'\n';(P/'c3_triangle_gate_certificate.json').write_text(s);print(s,end='')
