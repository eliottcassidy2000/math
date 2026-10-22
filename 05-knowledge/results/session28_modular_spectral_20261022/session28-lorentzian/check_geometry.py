from fractions import Fraction as F
from itertools import combinations, permutations, product
import json, math

def dot(a,b): return sum(x*y for x,y in zip(a[:-1],b[:-1]))-a[-1]*b[-1]
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def scale(a,k):return tuple(k*x for x in a)
def lattice(v):
    return (all(x.denominator==1 for x in v) or all(x.denominator==2 for x in v)) and (sum(v[:-1])-v[-1])%2==0
n=25
c=tuple(map(F,[12]*25+[70])); z=tuple(map(F,[14]*25+[60])); w=tuple(map(F,list(range(25))+[70]))
x=tuple(map(F,[1,1]+[0]*24))
assert dot(c,c)==-1300 and dot(z,z)==1300 and dot(c,z)==0 and lattice(x)
def rz(v):return add(v,scale(z,-2*dot(v,z)/dot(z,z)))
def cu(v):
    avg=sum(v[:-1])/25
    return tuple(2*avg-a for a in v[:-1])+(v[-1],)
def cz(v):return add(scale(v,-1),scale(c,2*dot(v,c)/dot(c,c)))
images={name:f(x) for name,f in [('normal_reflection',rz),('converse_normal_fixed',cu),('converse_normal_negated',cz)]}
assert all(not lattice(v) for v in images.values())
for f in (rz,cu,cz):
 for v in (x,c,z,w):assert dot(f(v),f(v))==dot(v,v) and f(f(v))==v
 assert f(c)==c
assert cu(w)==cz(w)==tuple(map(F,list(range(24,-1,-1))+[70]))
rstar=tuple([F(-1,2)]+[F(1,2)]*23+[F(3,2),F(5,2)])
nu=tuple([F(1,2)]*25+[F(5,2)])
assert lattice(rstar) and dot(rstar,rstar)==2 and dot(w,rstar)==-1 and dot(nu,rstar)==0
reflected=add(c,scale(rstar,-dot(c,rstar)))
assert lattice(reflected) and dot(reflected,reflected)==-1300 and reflected[-1]!=70
roots=[]
for i,j in combinations(range(25),2):
 r=tuple(F((k==i)-(k==j)) for k in range(26)); roots.append(r)
 assert dot(r,r)==2 and dot(r,c)==0 and dot(r,z)==0
 v=add(w,scale(r,-dot(w,r)))
 assert v[i]==w[j] and v[j]==w[i] and all(v[k]==w[k] for k in range(26) if k not in (i,j))
# Independent expected triangle and Hamiltonian quadratic coefficients from permutations.
counts=[]
for n in range(3,7):
 edges=list(combinations(range(n),2)); m=len(edges)
 idx={e:i for i,e in enumerate(edges)}
 B=[[int(v==i)-int(v==j) for i,j in edges] for v in range(n)]
 gram=[[sum(B[v][a]*B[v][b] for v in range(n)) for b in range(m)] for a in range(m)]
 H2=[[0]*m for _ in range(m)]
 for p in permutations(range(n)):
  arr=[(idx[tuple(sorted((a,b)))],1 if a<b else -1) for a,b in zip(p,p[1:])]
  for (a,sa),(b,sb) in combinations(arr,2):H2[a][b]+=sa*sb;H2[b][a]+=sa*sb
 factor=2*math.factorial(n-2)
 assert all(H2[a][b]==factor*(2*int(a==b)-gram[a][b]) for a in range(m) for b in range(m))
 # Expected triangles at non-Boolean rational x, exact probability sum.
 signs=[F((3*i+1)%7-3,4) for i in range(m)]
 def p(i,j):return (1+(signs[idx[(i,j)]] if i<j else -signs[idx[(j,i)]]))/2
 actual=sum(p(i,j)*p(j,k)*p(k,i)+p(j,i)*p(k,j)*p(i,k) for i,j,k in combinations(range(n),3))
 bx=[sum(b*t for b,t in zip(row,signs)) for row in B]
 formula=F(math.comb(n,3),4)+sum(t*t for t in signs)/4-sum(t*t for t in bx)/8
 assert actual==formula
 # Boolean score identity, exhaustive through5, deterministic small6 sample.
 universe=range(1<<m) if n<=5 else range(0,1<<m,101)
 tests=0
 for mask in universe:
  scores=[0]*n
  for e,(i,j) in enumerate(edges):scores[i if mask>>e&1 else j]+=1
  cyclic=sum(len({(i if (mask>>idx[(i,j)])&1 else j),(j if (mask>>idx[(j,k)])&1 else k),(i if (mask>>idx[(i,k)])&1 else k)})==3 for i,j,k in combinations(range(n),3))
  assert cyclic==math.comb(n,3)-sum(s*(s-1)//2 for s in scores)
  tests+=1
 counts.append({'n':n,'quadratic_matrix_entries':m*m,'boolean_triangle_tests':tests,'fractional_expected_triangles':str(actual)})
# metric numerator identity for random scores and explicit regular-center range.
import random
rng=random.Random(20261022); metric_tests=0
for _ in range(100):
 s=[0]*25;t=[0]*25
 for i,j in combinations(range(25),2):s[i if rng.randrange(2) else j]+=1;t[i if rng.randrange(2) else j]+=1
 a=2300-sum(y*(y-1)//2 for y in s);b=2300-sum(y*(y-1)//2 for y in t)
 vs=tuple(map(F,s+[70]));vt=tuple(map(F,t+[70]))
 assert -dot(vs,vt)==a+b+F(sum((i-j)**2 for i,j in zip(s,t)),2)
 assert dot(vs,c)==-1300 and dot(vs,vs)==-2*a and 0<a<=650
 metric_tests+=1
out={'status':'all exact assertions passed','root_swap_tests':len(roots),'forbidden_extensions_on_e0_plus_e1':{k:[str(y) for y in v] for k,v in images.items()},'rstar_regular_image':[str(y) for y in reflected],'hessian_and_triangle_checks':counts,'metric_tests':metric_tests}
print(json.dumps(out,indent=2,sort_keys=True))
