"""Independent literal signed matrices versus divisor-stratum formulas.
Standard library only; explicit failures survive Python -O.
"""
from math import gcd
from itertools import product
from collections import Counter
import json

def need(c,m):
 if not c: raise RuntimeError(m)
def divisors(n):return [d for d in range(1,n+1) if n%d==0]
def phi(n):return sum(gcd(x,n)==1 for x in range(n))
def primes(n):return [p for p in divisors(n) if p>1 and all(p%d for d in range(2,int(p**.5)+1))]
def order(u,n):
 x=u%n;k=1
 while x!=1:x=x*u%n;k+=1
 return k
def valuation(n):
 v=0
 while n%2==0:v+=1;n//=2
 return v
def arithmetic(n,u):
 out=[]
 for d in divisors(n):
  if d==1:continue
  o=order(u,d); neg=o%2==0 and pow(u,o//2,d)==d-1
  h=o//2 if neg else o;e=-1 if neg else 1;c=phi(d)//(2*h)
  out.extend([(h,e)]*c)
 return Counter(out)
def literal(n,u):
 r=(n-1)//2;dest=[];sign=[]
 for x in range(1,r+1):
  y=u*x%n;dest.append(min(y,n-y)-1);sign.append(1 if y<=r else -1)
 seen=set();out=[]
 for i in range(r):
  if i in seen:continue
  h=0;e=1;j=i
  while j not in seen:
   seen.add(j);h+=1;e*=sign[j];j=dest[j]
  need(j==i,'cycle closure');out.append((h,e))
 return dest,sign,Counter(out)
def factors(n,cycs,alpha):
 r=(n-1)//2;c=sum(cycs.values());z=t=0
 for (h,e),count in cycs.items():
  eta=e if alpha==-1 else ((-1)**h)*e
  if eta==1:z+=count
  else:t+=count
 return r-c,t,z
def matrix(dest,sgn,alpha):
 r=len(dest);a=[[int(i==j) for j in range(r)] for i in range(r)]
 for i,j in enumerate(dest):a[i][j]+=alpha*sgn[i]
 return a
def direct_snf(a):
 # Every signed-cycle matrix has only unit pivots followed by isolated 2 pivots.
 # This does not inspect cycles, divisor orders or proposed factor counts.
 a=[row[:] for row in a];n=len(a);k=0
 while k<n:
  pivot=next(((i,j) for i in range(k,n) for j in range(k,n) if abs(a[i][j])==1),None)
  if pivot is None:break
  i,j=pivot;a[k],a[i]=a[i],a[k]
  for row in a:row[k],row[j]=row[j],row[k]
  if a[k][k]==-1:a[k]=[-v for v in a[k]]
  for i in range(k+1,n):
   t=a[i][k]
   if t:a[i]=[v-t*w for v,w in zip(a[i],a[k])]
  for j in range(k+1,n):
   t=a[k][j]
   if t:
    for i in range(n):a[i][j]-=t*a[i][k]
  k+=1
 rem=[row[k:] for row in a[k:]]
 need(all(v in (-2,0,2) for row in rem for v in row),'remaining scalars')
 need(all(sum(v!=0 for v in row)<=1 for row in rem),'remaining rows')
 need(all(sum(row[j]!=0 for row in rem)<=1 for j in range(n-k)),'remaining columns')
 t=sum(v!=0 for row in rem for v in row)
 return k,t,n-k-t

def run():
 cases=snf_cases=modular_vectors=sign_words=0;profiles=Counter();examples={}
 for n in range(3,202,2):
  for u in range(1,n):
   if gcd(n,u)!=1:continue
   dest,sgn,cycles=literal(n,u);pred=arithmetic(n,u)
   need(cycles==pred,('cycle strata',n,u,cycles,pred))
   need(sum(h*c for (h,e),c in pred.items())==(n-1)//2,'dimension')
   ts={valuation(order(u,p)) for p in primes(n)}
   profile='mixed' if len(ts)>1 else str(next(iter(ts)))
   profiles[profile]+=1
   expected_minus=len(ts)==1 and next(iter(ts))>=1
   expected_plus=len(ts)==1 and next(iter(ts))!=1
   for alpha,expected in [(-1,expected_minus),(1,expected_plus)]:
    one,two,zero=factors(n,pred,alpha)
    need((zero==0)==expected,('profile criterion',n,u,alpha,ts))
    if n<=51:
     need(direct_snf(matrix(dest,sgn,alpha))==(one,two,zero),('SNF',n,u,alpha))
     snf_cases+=1
    if n<=11:
     for q in [2,3,4]:
      actual=0
      for values in product(range(q),repeat=len(dest)):
       actual+=all((values[i]+alpha*sgn[i]*values[j])%q==0 for i,j in enumerate(dest))
       modular_vectors+=1
      need(actual==q**zero*gcd(q,2)**two,('modular kernel',n,u,alpha,q))
   if n<=21:
    for a in [1,-1]:
     actual=0
     for values in product([-1,1],repeat=len(dest)):
      actual+=all(sgn[i]*values[j]==a*values[i] for i,j in enumerate(dest));sign_words+=1
     admissible=ts==({0} if a==1 else {1})
     need(actual==(2**sum(pred.values()) if admissible else 0),('sign count',n,u,a))
   if (n,u) in [(3,2),(5,2),(7,2),(9,2),(15,2),(17,2),(65,2),(73,2)]:
    examples[str(n)]={'profile':sorted(ts),'cycles':[[h,e,c] for (h,e),c in sorted(pred.items())],
      'minus_snf_counts':factors(n,pred,-1),'plus_snf_counts':factors(n,pred,1)}
   cases+=1
 # Concrete integral obstruction, not merely rational invertibility.
 dest,sgn,_=literal(5,2);a=matrix(dest,sgn,1)
 need(a==[[1,1],[-1,1]],'N5 literal matrix')
 need(all(not (x+y==1 and -x+y==0) for x in range(-10,11) for y in range(-10,11)),'N5 parity obstruction')
 # The exact equations imply x=y and 2x=1; the finite range is only a control.
 need(examples['15']['minus_snf_counts']==(4,2,1),'mixed minus hostile')
 need(examples['15']['plus_snf_counts']==(4,1,2),'mixed plus hostile')
 return {'status':'PASS','arithmetic_matrix_cases':cases,'odd_moduli_through':201,
  'direct_integer_snf_cases':snf_cases,'direct_snf_moduli_through':51,
  'exhaustive_modular_vectors':modular_vectors,'modular_kernel_moduli':[2,3,4],
  'modular_vector_input_moduli_through':11,'exhaustive_sign_words':sign_words,
  'sign_input_moduli_through':21,'profile_counts':dict(sorted(profiles.items())),
  'examples':examples,'scope':'normalized reflection-odd lattice; full even-coefficient odd module has separate fixed-point torsion'}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
