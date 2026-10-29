"""Exact source-mass moments; no external dependencies or assertion-only tests."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
import json
checks=0
def eq(a,b):
 global checks
 checks+=1
 if a!=b:raise RuntimeError((a,b))
def mobius(n):
 ans=1;p=2
 while p*p<=n:
  if n%p==0:
   n//=p;ans=-ans
   if n%p==0:return 0
   while n%p==0:n//=p
  p+=1
 return -ans if n>1 else ans
def divisors(n):return [d for d in range(1,n+1) if n%d==0]
def key(word):return min(word[i:]+word[:i] for i in range(len(word)))
def push(distributions):
 out=defaultdict(Q)
 for word in product(*(tuple(d) for d in distributions)):
  w=Q(1)
  for d,a in zip(distributions,word):w*=d[a]
  out[key(word)]+=w
 return dict(out)
def moment(d,r):return sum((x**r for x in d.values() if x),Q(0))
def mixed(a,b,r,s):return sum((a[x]**r*b[x]**s for x in a.keys()&b.keys() if a[x] and b[x]),Q(0))
def recursive_moment(p,k,r,d=3):
 @lru_cache(None)
 def f(k,r):
  if k==0:return sum((x**r for x in p if x),Q(0))
  return sum((Q(ell)**(r-1)*sum((mobius(e)*f(k-1,r*d*e//ell)**(ell//e) for e in divisors(ell)),Q(0)) for ell in divisors(d)))
 return f(k,r)
def recursive_mixed(p,q,k,r,s,d=3):
 @lru_cache(None)
 def f(k,r,s):
  if k==0:return sum((x**r*y**s for x,y in zip(p,q) if x and y),Q(0))
  return sum((Q(ell)**(r+s-1)*sum((mobius(e)*f(k-1,r*d*e//ell,s*d*e//ell)**(ell//e) for e in divisors(ell)),Q(0)) for ell in divisors(d)))
 return f(k,r,s)
def centralizer_moment(p,k,r):
 ans=sum((x**r for x in p if x),Q(0))
 pp=tuple(x*y for x in p for y in p)
 for j in range(k):ans*=recursive_moment(pp,j,r)
 return ans
def negindex(v,k):
 out=0
 for j in range(k):
  out+=((-((v//3**j)%3))%3)*3**j
 return out
def ka(f,k):
 if not k:return f[0]
 n=3**(k-1)
 return key(tuple(ka(f[i*n:(i+1)*n],k-1) for i in range(3)))
def kc(f,k):
 if not k:return f[0]
 n=3**(k-1)
 h=tuple((f[n+i],f[2*n+negindex(i,k-1)]) for i in range(n))
 return (kc(f[:n],k-1),ka(h,k-1))
def direct(p,k,canon):
 out=defaultdict(Q)
 for f in product(range(len(p)),repeat=3**k):
  w=Q(1)
  for a in f:w*=p[a]
  out[canon(f,k)]+=w
 return dict(out)
def entropy_loss(p,k,canon):
 sizes=defaultdict(int)
 for f in product(range(len(p)),repeat=3**k):sizes[canon(f,k)]+=1
 masses=direct(p,k,canon);out=Q(0)
 for key_,size in sizes.items():
  power=0
  while size>1:
   if size%3:raise RuntimeError("nonternary orbit size")
   size//=3;power+=1
  out+=masses[key_]*power
 return out
def predicted_loss(p,k):
 if not k:return Q(0)
 return 3*predicted_loss(p,k-1)+1-recursive_moment(p,k-1,3)
def main():
 for d in range(2,10):
  for p in [(Q(1,2),Q(1,2)),(Q(1,3),Q(2,3)),(Q(1,6),Q(1,3),Q(1,2))]:
   if len(p)**d>25000:continue
   src=dict(enumerate(p));out=push([src]*d)
   for r in range(5):eq(moment(out,r),recursive_moment(p,1,r,d))
   eq(sum(out.values()),1)
 for p in [(Q(1,2),Q(1,2)),(Q(1,3),Q(2,3))]:
  for k in range(3):
   da=direct(p,k,ka);dc=direct(p,k,kc)
   eq(entropy_loss(p,k,ka),predicted_loss(p,k))
   pp=tuple(x*y for x in p for y in p)
   eq(entropy_loss(p,k,kc),sum((predicted_loss(pp,j) for j in range(k)),Q(0)))
   for r in range(5):
    eq(moment(da,r),recursive_moment(p,k,r))
    eq(moment(dc,r),centralizer_moment(p,k,r))
   q=p[::-1];db=direct(q,k,ka)
   for r,s in [(1,1),(1,2),(2,1),(2,2)]:eq(mixed(da,db,r,s),recursive_mixed(p,q,k,r,s))
   # Global converse followed by reflection pairs C-orbits, but source mass need not be half.
   folded=defaultdict(Q)
   for f in product(range(2),repeat=3**k):
    w=Q(1)
    for a in f:w*=p[a]
    rev=tuple(1-f[negindex(i,k)] for i in range(3**k))
    folded[min(kc(f,k),kc(rev,k))]+=w
   dq=direct(q,k,kc)
   eq(moment(folded,2),moment(dc,2)+mixed(dc,dq,1,1))
   if p[0]==p[1]:eq(moment(folded,2),2*moment(dc,2))
 # Arbitrary noninvariant source: same 3x3 overlaps, different exact collision.
 sets=[]
 for j in range(5):sets.append({x:Q(1,16) for x in range(32) if (x>>j)&1==0})
 parity=lambda x:x.bit_count()%2
 laws=[]
 for last in [0,1]:
  ds=sets+[{x:Q(1,16) for x in range(32) if parity(x)==last}]
  a,b=ds[:3],ds[3:]
  matrix=[[mixed(x,y,1,1) for y in b] for x in a]
  eq(matrix,[[Q(1,32)]*3 for _ in range(3)])
  six=sum((__import__('math').prod(t.get(x,0) for t in ds) for x in range(32)),Q(0))
  collision=mixed(push(a),push(b),1,1)
  eq(collision,sum((__import__('math').prod(matrix[i][(i+s)%3] for i in range(3)) for s in range(3)),Q(0))-2*six)
  laws.append(collision)
 eq(laws,[Q(1534,16**6),Q(1536,16**6)])
 eq(recursive_moment((Q(1,2),Q(1,2)),1,2),Q(5,16))
 eq(recursive_moment((Q(1,2),Q(1,2)),2,3),Q(19775,2097152))
 print(json.dumps({'checks':checks,'status':'PASS','binary_A_collision':[str(recursive_moment((Q(1,2),Q(1,2)),k,2)) for k in range(4)],'binary_C_collision':[str(centralizer_moment((Q(1,2),Q(1,2)),k,2)) for k in range(4)],'heterogeneous_hostile':[str(x) for x in laws]},sort_keys=True))
if __name__=='__main__':main()
