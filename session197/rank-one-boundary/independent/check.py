from itertools import permutations,product
from math import factorial
import random,json

def need(x,msg):
 if not x:raise RuntimeError(msg)
def graph(n,code):
 a=[0]*n;k=0
 for u in range(n):
  for v in range(u+1,n):
   i,j=(u,v) if code>>k&1 else (v,u);a[i]|=1<<j;k+=1
 return tuple(a)
def endpoint(a,vs):
 s={v:0 for v in vs};e=s.copy()
 for p in permutations(vs):
  if all(a[u]>>v&1 for u,v in zip(p,p[1:])):
   if p:s[p[0]]+=1;e[p[-1]]+=1
 return s,e

def runs(a,r):
 D=[v for v in range(len(a)) if v!=r];out=[0]*(len(a)+1)
 for p in permutations(D):
  for mask in range(1<<(len(D)+1)):
   word=[]
   for i in range(len(D)+1):
    if mask>>i&1:word.append(r)
    if i<len(D):word.append(p[i])
   if all(a[u]>>v&1 for u,v in zip(word,word[1:])):out[mask.bit_count()]+=1
 return out

def check(a,r):
 O=[v for v in range(len(a)) if a[r]>>v&1];I=[v for v in range(len(a)) if a[v]>>r&1];D=O+I
 edges=[(u,v) for u in O for v in I if a[u]>>v&1]
 def match(k,used):
  if k==len(O):return 0
  return max([match(k+1,used)]+[1+match(k+1,used|{v}) for v in I if v not in used and (O[k],v) in edges])
 nu=match(0,set());R=runs(a,r);b=sum((-1)**j*v for j,v in enumerate(R));signed=0
 if D:
  for p in permutations(D):
   if p[0] not in I or p[-1] not in O:continue
   w=1
   for u,v in zip(p,p[1:]):w*=((a[u]>>v)&1)-int(u in I and v in O)
   signed+=w
 need(b==signed,'literal root gaps versus signed word sum')
 centers=0
 if nu==0:need(b==0,'rank0 boundary including singleton')
 if nu==1:
  so,eo=endpoint(a,O);si,ei=endpoint(a,I)
  need(R[2]>0 and not any(R[3:]),'exact maximum root runs')
  for c in D:
   if not all(c in edge for edge in edges):continue
   prod=so[c]*sum(ei[v] for v in I if a[c]>>v&1) if c in O else ei[c]*sum(so[u] for u in O if a[u]>>c&1)
   need(b==-prod,'literal endpoint factorization');centers+=1
  need(centers>0 and b<=0,'cover existence and sign')
 return nu,centers,b,R
stats={};roots=0;centers=0
for n in range(1,6):
 for code in range(1<<(n*(n-1)//2)):
  a=graph(n,code)
  for r in range(n):
   nu,c,b,R=check(a,r);stats[nu]=stats.get(nu,0)+1;centers+=c;roots+=1
rng=random.Random(1971213)
for n,rounds in ((6,12),(7,4)):
 for _ in range(rounds):
  a=graph(n,rng.getrandbits(n*(n-1)//2))
  for r in range(n):
   nu,c,b,R=check(a,r);stats[nu]=stats.get(nu,0)+1;centers+=c;roots+=1
hostiles=[]
for a,expect in [((16,5,9,3,14),[5,15,11,2,0,0]),((2,20,9,3,13),[1,11,17,4,0,0])]:
 nu,c,b,R=check(a,1);need(nu==2 and R==expect,'exact rank2 hostile')
 hostiles.append({'rows':a,'root':1,'runs':R,'boundary':b})
print(json.dumps({'status':'PASS','literal_rooted_graphs':roots,'all_cover_centers':centers,'rank_counts':stats,'hostiles':hostiles},sort_keys=True,indent=2))
