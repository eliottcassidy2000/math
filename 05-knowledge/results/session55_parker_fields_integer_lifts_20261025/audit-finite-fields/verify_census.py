"""Independent certificate replay, sieve universe and Rabin irreducibility.
No producer functions are imported. Witness encodings use the supplied modulus.
"""
from field_arithmetic import K
import argparse,json,hashlib,math,time
from pathlib import Path
BASE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--census',type=Path,default=BASE.parent/'finite-fields/census.json')
parser.add_argument('--small-counts',type=Path)
args=parser.parse_args()
CENSUS=args.census
D=json.loads(CENSUS.read_text());start=time.time();limit=D['limit_exclusive']
# Independent Eratosthenes universe, not producer trial primality enumeration.
sieve=bytearray(b'\1')*limit;sieve[:2]=b'\0\0'
for p in range(2,math.isqrt(limit)+1):
 if sieve[p]:sieve[p*p:limit:p]=b'\0'*len(range(p*p,limit,p))
expected=[]
for p in range(3,limit,2):
 if sieve[p]:
  q=p;n=1
  while q<limit:expected.append((q,p,n));q*=p;n+=1
assert sorted(expected)==[(r['q'],r['p'],r['n']) for r in D['rows']]
def trim(a):
 while a and not a[-1]:a.pop()
 return a
def gcd(a,b,p):
 a=trim(a[:]);b=trim(b[:])
 while b:
  while len(a)>=len(b):
   c=a[-1]*pow(b[-1],-1,p)%p;off=len(a)-len(b)
   for j in range(len(b)):a[j+off]=(a[j+off]-c*b[j])%p
   trim(a)
  a,b=b,a
 return a
def field(r):
 F=K.__new__(K);F.p=r['p'];F.f=r['modulus_low_first'];F.n=r['n'];F.z=(0,)*F.n;F.one=(1,)+(0,)*(F.n-1)
 assert len(F.f)==F.n+1 and F.f[-1]==1 and all(0<=x<F.p for x in F.f)
 if F.n>1:
  x=(0,1)+(0,)*(F.n-2);z=x
  for k in range(1,F.n+1):
   z=F.pow(z,F.p)
   if k<=F.n//2:assert len(gcd(F.f,[ (a-b)%F.p for a,b in zip(z,x)],F.p))==1,(r,k)
  assert z==x
 return F
def decode(i,F):
 assert 0<=i<F.p**F.n
 out=[]
 for _ in range(F.n):out.append(i%F.p);i//=F.p
 return tuple(out)
def matrix(F,c,u,v):
 return [F.add(c,F.add(F.scale(s,u),F.scale(t,v))) for s,t in [(1,0),(-1,-1),(0,1),(-1,1),(0,0),(1,-1),(0,-1),(1,1),(-1,0)]]
def valid(F,m,sq=None):
 if len(set(m))!=9:return False
 if sq is not None:
  if any(x not in sq for x in m):return False
 elif any(x!=F.z and F.pow(x,(F.p**F.n-1)//2)!=F.one for x in m):return False
 lines=[(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
 return len({F.add(F.add(m[a],m[b]),m[c]) for a,b,c in lines})==1
small=json.loads((args.small_counts or CENSUS.parent/'small_counts.json').read_text());sm={r['q']:r for r in small};negs=[];counts=[];positive=0
for r in D['rows']:
 F=field(r);q=r['q']
 if r['witness'] is not None:
  assert valid(F,[decode(i,F) for i in r['witness']]),q;positive+=1
 if r['witness'] is None or q in sm:
  E=[decode(i,F) for i in range(q)];sq={F.mul(x,x) for x in E}
  assert len(sq)==(q+1)//2
  # Separate direct matrix path, exhaustive no early exits. Center-one
  # pair condition is mathematically necessary and checked independently.
  count0=sum(valid(F,matrix(F,F.z,F.one,v),sq) for v in E)
  C=[x for x in E if F.add(F.one,x) in sq and F.add(F.one,F.scale(-1,x)) in sq]
  count1=sum(valid(F,matrix(F,F.one,u,v),sq) for u in C for v in C)
  if r['witness'] is None:assert count0==count1==0;negs.append(q)
  if q in sm:assert (count0,count1)==(sm[q]['center_zero_u1'],sm[q]['center_one']),(q,count0,count1)
  counts.append({'q':q,'center_zero_u1':count0,'center_one':count1})
assert negs==D['parker_orders'];assert len(sm)==len(counts)
assert min(r['q'] for r in counts if r['center_one'])==59
out={'status':'PASS','census_sha256':hashlib.sha256(CENSUS.read_bytes()).hexdigest(),'field_count':len(expected),'positive_witnesses':positive,'negative_orders':negs,'all_moduli':'Rabin/Frobenius independently checked','small_counts_replayed':len(counts),'first_center_one':59,'counts':counts,'seconds':time.time()-start}
(BASE/'census-audit-results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
