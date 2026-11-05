from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from itertools import product
import random,json
H=Path(__file__).resolve().parent
checks=0

def need(c,m):
 global checks
 checks+=1
 if not c:raise RuntimeError(m)

def source(U):
 D,E,V,d=1,0,2,1
 for c in bin(U-1)[3:]:
  D,E,V,d=(D,0,2*V,d)if c=='1'else(D-E+V,d,2*V-d,d)
 need(V==U and E==U%2 and d==1,'literal O/Z response compiler')
 return D

def runs(n,k):
 out=[];i=0
 while i<k:
  if n>>i&1:i+=1;continue
  j=i
  while not(n>>i&1):i+=1
  out.append((j,i-j,n%2**j))
 return out

def control(n):
 k=n.bit_length();B=1<<k;R=int(f'{n:0{k}b}'[::-1],2);C=B-1-R
 D=source(n+1);rr=runs(n,k)
 correction=B*(1+len(rr))-sum(m*((B>>i)-(B>>(i+l)))for i,l,m in rr)
 need(B*D-n*C==correction,'exact run product')
 need(B<=correction<=B*(1+len(rr))<=B*(1+k//2),'positive correction bound')
 need((correction==B)==(len(rr)==0),'equality boundary')
 return D,C

def main():
 for n in range(1,100001):control(n)
 rng=random.Random(114891)
 for k in range(1,257):
  for _ in range(8):control((1<<(k-1))+rng.getrandbits(k-1))
 rectangles=0;max_ratio=F(0)
 for k in range(1,12):
  B=1<<k;pts=[(F(n,B),F(source(n+1),n))for n in range(B//2,B)]
  for xnum,ynum in product(range(8,17),range(17)):
   x,y=F(xnum,16),F(ynum,16)
   count=sum(a<x and b<y for a,b in pts)
   err=abs(count-B*(x-F(1,2))*y)
   need(err<=2*k+4,'anchored discrepancy');rectangles+=1;max_ratio=max(max_ratio,err/(2*k+4))
  for _ in range(200):
   a,b=sorted([F(rng.randrange(8,17),16),F(rng.randrange(8,17),16)])
   c,d=sorted([F(rng.randrange(17),16),F(rng.randrange(17),16)])
   count=sum(a<=x<b and c<=y<d for x,y in pts)
   need(abs(count-B*(b-a)*(d-c))<=8*k+16,'rectangle discrepancy');rectangles+=1
   if B*(b-a)*(d-c)>8*k+16:need(count>0,'physical rectangle witness')
 print(json.dumps(dict(status='PASS',checks=checks,exact_initial_indices=100000,large_integer_controls=2048,rectangle_checks=rectangles,largest_anchored_error_fraction_of_bound=str(max_ratio),proof_sha256=sha256((H/'PROOF.md').read_bytes()).hexdigest(),script_sha256=sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,sort_keys=True))
if __name__=='__main__':main()
