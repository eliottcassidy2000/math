#!/usr/bin/env python3
import itertools as it, math,json,functools,collections
from fractions import Fraction as F
from pathlib import Path
counts=collections.Counter()
def ck(x,s):
 counts[s]+=1
 if not x:raise RuntimeError(s)
# Geometry: integer-scaled signed sums, rather than producer pairwise Fraction points.
D=math.lcm(*(1+t*t for t in range(1,8)))
v=[((1-t*t)*(D//(1+t*t)),2*t*(D//(1+t*t))) for t in range(1,8)]
unit=0
for z in it.product([-1,0,1],repeat=7):
 x=sum(a*c for (a,b),c in zip(v,z));y=sum(b*c for (a,b),c in zip(v,z));sq=x*x+y*y
 support=sum(c!=0 for c in z)
 ck((sq==0)==(support==0),'signed_sum_injectivity')
 ck((sq==D*D)==(support==1),'signed_sum_induced_distance')
 unit+=sq==D*D
ck(unit==14,'unit_difference_types')
# Exact Gaussian elimination on K8 Gram matrix, independent of eigenvalue shortcut.
G=[[F(1 if i==j else 1,1 if i==j else 2) for j in range(7)] for i in range(7)]
det=F(1)
for k in range(7):
 p=G[k][k];ck(p!=0,'gram_pivot');det*=p
 for j in range(k,7):G[k][j]/=p
 for i in range(k+1,7):
  f=G[i][k]
  for j in range(k,7):G[i][j]-=f*G[k][j]
ck(det==F(1,16),'K8_Gram_determinant')
def tour(n,b):
 a=[[0]*n for _ in range(n)]
 for k,(u,v) in enumerate(it.combinations(range(n),2)):a[u][v]=(b>>k)&1;a[v][u]=1-a[u][v]
 return a
def profile(a):
 n=len(a)
 h={s:sum(all(a[u][v] for u,v in zip(p,p[1:])) for p in it.permutations([i for i in range(n) if s>>i&1])) for s in range(1,1<<n)}
 @functools.lru_cache(None)
 def part(s):
  if not s:return (1,)+(0,)*n
  out=[0]*(n+1);t=s;low=s&-s
  while t:
   if t&low:
    sub=part(s^t)
    for k in range(1,n+1):out[k]+=h[t]*sub[k-1]
   t=(t-1)&s
  return tuple(out)
 return [math.factorial(k)*x for k,x in enumerate(part((1<<n)-1))]
for n in range(3,6):
 for bits in range(1<<(n*(n-1)//2)):
  a=tour(n,bits);f=profile(a)
  # Census by degree complement, not producer explicit cyclic triples.
  c3=math.comb(n,3)-sum(math.comb(sum(row),2) for row in a)
  ck(f[n-2]==math.factorial(n-2)*(math.comb(n,3)+2*c3+3*math.comb(n,4)),'triangle_cover_A1')
  # General shift extraction rather than three-term pre-simplification.
  coeff=sum(F(f[d])*math.comb(d-1,n-3)*F(-1,2)**(d-1-(n-3)) for d in range(n-2,n+1))
  ck(coeff==2*math.factorial(n-2)*c3-F(math.factorial(n)*(n-2),12),'centered_coefficient_A2')
# Expand hypothetical polynomial after u=z+1/2: 240[(z+1/2)^2-1/4]^2 =240z²(1+z)².
shift=[F(0)]*5
for k,coef in [(4,F(240)),(2,F(-120)),(0,F(15))]:
 for j in range(k+1):shift[j]+=coef*math.comb(k,j)*F(1,2)**(k-j)
ck(shift==[0,0,240,480,240],'hypothetical_A3_A4_shift')
ck(4*shift[2]==960,'hypothetical_A5_probe')
def coeff(q):
 m=len(q)-1
 @functools.lru_cache(None)
 def go(used,last):
  out=[0]*(m+2)
  if used==(1<<m)-1:out[0]=1
  if last!=0 and (last==-1 or q[last][0]):
   z=go(used,0)
   for k in range(m+1):out[k+1]+=z[k]
  for j in range(1,m+1):
   if not used>>(j-1)&1 and (last==-1 or q[last][j]):
    z=go(used|1<<(j-1),j)
    for k in range(m+2):out[k]+=z[k]
  return tuple(out)
 return list(go(0,-1))
def build(B):
 r=len(B);q=[[0]*(2*r+1) for _ in range(2*r+1)]
 for u,w in it.combinations(range(2*r+1),2):
  if u==0:value=w<=r
  elif u<=r<w:value=B[u-1][w-r-1]
  else:value=1
  q[u][w]=int(value);q[w][u]=1-int(value)
 return q
def ryser(B):
 r=len(B)
 return sum((-1)**(r-z.bit_count())*math.prod(sum(B[i][j] for j in range(r) if z>>j&1) for i in range(r)) for z in range(1<<r))
Bs=[[[1,0,0],[1,0,0],[0,1,1]],[[1,0,0],[0,1,0],[1,0,1]]]
qs=list(map(build,Bs));cs=list(map(coeff,qs));scores=[list(map(sum,q)) for q in qs]
ck(scores==[[3,3,2,2,4,4,3]]*2,'same_labelled_scores')
ck([ryser(b) for b in Bs]==[0,1],'Ryser_Hall_hostile')
ck([c[4] for c in cs]==[0,6],'automaton_Hall_hostile')
# Producer arrays intentionally omit c0-containing external-only words by setting index0=0.
for c in cs:c[0]=0
saved=json.loads((Path(__file__).resolve().parent.parent/'tournament/addendum.json').read_text())
ck(cs==saved['word_coefficient_vectors'],'full_word_vectors_match')
for r in [1,2]:
 seen={}
 for b in range(1<<(r*r)):
  B=[[(b>>(r*i+j))&1 for j in range(r)] for i in range(r)]
  key=(tuple(map(sum,B)),tuple(map(sum,zip(*B))));p=ryser(B)
  ck(key not in seen or seen[key]==p,'small_margin_minimality')
  seen[key]=p
print(json.dumps({'status':'PASS','checks':dict(counts),'geometry_scale':D,'geometry_signed_sum_vectors':3**7,'K8_Gram_determinant':str(det),'word_coefficients':cs,'scores':scores,'hypothetical_unit_triangle_difference_F_prefix':[str(x) for x in shift]},indent=2,sort_keys=True))
