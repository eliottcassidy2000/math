from math import isqrt
import json

def need(c,m):
 if not c:raise RuntimeError(m)
def prime(n):
 if n<2:return False
 if n%2==0:return n==2
 if n%3==0:return n==3
 return all(n%d for d in range(5,isqrt(n)+1,2))
B=10000;s=[0]*(B+1)
for d in range(1,B//2+1):
 for n in range(2*d,B+1,d):s[n]+=d
perfect=[n for n in range(1,B+1) if s[n]==n]
need(perfect==[6,28,496,8128],'small perfect control')
companions=[n for n in perfect if prime(n-1)]
need(companions==[6],'prime square companion')
need(s[25]==6 and s[6]==6,'actual companion source')
rows=[]
for r in list(range(2,20))+[31]:
 q=(1<<r)-1
 if not prime(q):continue
 N=(1<<(r-1))*q
 need(((1<<r)-1)*(q+1)==2*N,'Euler constructed perfect')
 need(r==2 or N%3==1,'mod three exclusion')
 need(prime(N-1)==(r==2),'prime square branch unique')
 rows.append({'exponent':r,'N':N,'N_minus_1_prime':prime(N-1)})
need(not (prime(2) and s[3]==3),'p2 boundary')
print(json.dumps({'status':'PASS','sieve_bound':B,'perfect_targets':perfect,'prime_square_companion_targets':companions,'Mersenne_controls':rows},sort_keys=True,indent=2))
