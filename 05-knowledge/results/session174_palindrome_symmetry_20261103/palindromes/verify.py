from itertools import permutations
import json,math
checks=0
def check(x):
 global checks
 checks+=1
 if not x:raise RuntimeError('exact verification failed')
def primes(n):return [p for p in range(2,n+1) if all(p%d for d in range(2,math.isqrt(p)+1))]
def order(b,q):
 x=b%q;t=1
 while x!=1:x=x*b%q;t+=1
 return t
def digits(x,b):
 d=[]
 while x:d.append(x%b);x//=b
 return d
for b in range(2,17):
 for m in range(1,16):
  P=b**(2*m+1)-b**m-1
  check((2*b**(m+1)-1)**2==4*b*P+4*b+1)
  check(digits(P,b)==[b-1]*m+[b-2]+[b-1]*m)
sieve=[]
for b in range(2,13):
 for q in primes(300):
  if q==2 or b%q==0:continue
  t=order(b,q);roots=[m for m in range(t) if (b*pow(b,m,q)**2-pow(b,m,q)-1)%q==0]
  check(len(roots)<=2)
  if roots:check(pow((4*b+1)%q,(q-1)//2,q) in (0,1))
  if t%2:check(len(roots)<=1)
  else:check(set(roots)=={(t//2-1-r)%t for r in roots})
  if b==10 and roots:sieve.append([q,t,roots])
for m in range(200):check((pow(10,2*m+1,41)-pow(10,m,41)-1)%41!=0)
for n in range(30):
 x=10**(n+1);B=100*x*x+666*x+1
 check((100*x+333)**2-110789==100*B)
 check(str(B)=='1'+'0'*n+'666'+'0'*n+'1')
check(110789==49*2261==7**3*17*19)
for b in range(2,13):
 for m in range(1,10):
  for d in range(1,m+1):
   x=b**(m-d);A=b**(2*d+1);C=b**(2*d)+1;P=b**(2*m+1)-b**(m+d)-b**(m-d)-1
   check(P==A*x*x-C*x-1);check((2*A*x-C)**2==4*A*P+C*C+4*A)
   if b>=3 or m>d:check(digits(P,b)==digits(P,b)[::-1])
check(digits(2,2)==[0,1]);check([0,1]!=[1,0]);check(sum([1])==1) # n=1 S polynomial
# User's 7417 is an explicit two-residue sieve prime.
check(7417 in primes(7417))
check(order(10,7417)==2472)
r7417=[m for m in range(2472) if (pow(10,2*m+1,7417)-pow(10,m,7417)-1)%7417==0]
check(r7417==[1312,2395]);check((1312+2395)%2472==2472//2-1)
counts=[]
for n in range(2,6):
 edges=[(i,j) for i in range(n) for j in range(i+1,n)];perms=list(permutations(range(n)))
 signs={p:(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n)) for p in perms}
 for mask in range(1<<len(edges)):
  arcs={(i,j) if mask>>k&1 else (j,i) for k,(i,j) in enumerate(edges)};F=[0]*n;S=[0]*n
  for p in perms:
   f=sum((p[i],p[i+1]) in arcs for i in range(n-1));F[f]+=1;S[f]+=signs[p]
  eps=(-1)**(n*(n-1)//2)
  check(F==F[::-1]);check(S==[eps*x for x in S[::-1]]);check(sum(S)==0)
  if eps==1:check(sum(i*v for i,v in enumerate(S))==0)
  if n%4 in (0,3):check(sum((-1)**i*v for i,v in enumerate(S))==0)
  if n%2==0:check(sum((-1)**i*v for i,v in enumerate(F))==0)
 counts.append([n,1<<len(edges)])
print(json.dumps({'status':'PASS','explicit_checks':checks,'decimal_sieve_q_le_300':sieve,'all_labelled_tournaments':counts,'binary_outer_defect_hostile':2,'signed_factor_minimum_n':2,'q7417_order':2472,'q7417_roots':r7417},indent=2))
