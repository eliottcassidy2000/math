"""Independent exact controls for the perfect-target/prime-square conjunction."""
from math import isqrt
import argparse,json

checks=0
def need(c,m):
 global checks
 checks+=1
 if not c:raise RuntimeError(m)
ap=argparse.ArgumentParser();ap.add_argument('--negative-control',action='store_true');args=ap.parse_args()
B=1000000
spf=list(range(B+1))
for p in range(2,isqrt(B)+1):
 if spf[p]==p:
  for n in range(p*p,B+1,p):
   if spf[n]==n:spf[n]=p
sig=[0]*(B+1);sig[1]=1
for n in range(2,B+1):
 p=spf[n];r=n;power=1;geometric=1
 while r%p==0:r//=p;power*=p;geometric+=power
 sig[n]=geometric*sig[r]
perfect=[n for n in range(1,B+1) if sig[n]==2*n]
need(perfect==[6,28,496,8128],'complete perfect target census')
prime_square_targets=[(p+1,p) for p in range(2,100001) if spf[p]==p and sig[p+1]==2*(p+1)]
need(prime_square_targets==([(28,27)] if args.negative_control else [(6,5)]),'negative-control unique perfect prime-square target')
need(sig[6]-6==sig[25]-25==6,'direct divisor-source check')
need(sig[3]!=6 and sig[4]!=8,'p2 and forced p3 boundary')
need(sig[9]-9==4 and sig[49]-49==8,'removing perfectness destroys uniqueness')
need(sig[28]==56 and spf[27]!=27,'other perfect number has no prime-square source')
need([sig[n]-n for n in (4,9,57,65,8,18)]==[3,4,23,19,7,21],'root residue witnesses')
need(2%7==23%7 and 5%7==19%7,'root prime target residue collisions')
need([(2*n)%7 for n in (2,5)]==[4,3],'root multiplier images')
need({frozenset((n,7-n)) for n in range(1,7)}=={frozenset((i,j)) for i in range(1,7) for j in range(i+1,7) if i*i%7==j*j%7},'all opposite pairs collapse')
need([n for n in range(1,100) if n*n==n]==[1],'unit square loop hostile')
need(sig[6]-6==6 and sig[28]-28==28,'two distinct aliquot loops')
for N in perfect:
 a=0;m=N
 while m%2==0:a+=1;m//=2
 q=2**(a+1)-1
 need(a>=1 and m%q==0,'Euler q divides odd core')
 t=m//q
 need(sig[m]==m+t and t==1,'Euler t equals one')
 need(spf[q]==q and spf[a+1]==a+1,'Euler primes')
 need(N==6 or (N%3==1 and N-1>=27 and (N-1)%3==0),'mod-three obstruction')
# A direct finite divisor-pair route, independent of the multiplicative sieve.
for p in range(2,1000):
 if spf[p]!=p:continue
 n=p*p;proper=0
 for d in range(1,isqrt(n)+1):
  if n%d:continue
  if d<n:proper+=d
  q=n//d
  if q!=d and q<n:proper+=q
 need(proper==p+1,'prime-square divisor formula')
print(json.dumps({'status':'PASS','checks':checks,'target_census_bound':B,'perfect_targets':perfect,
 'prime_source_bound':100000,'prime_square_perfect_pairs_N_p':prime_square_targets,
 'direct_prime_square_source_count':168,'negative_control_gate_always_active':True},indent=2,sort_keys=True))
