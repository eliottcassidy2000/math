#!/usr/bin/env python3
"""Exact finite controls of the word-address algebra classification."""
from functools import lru_cache
from math import gcd
import json
COUNTS={}
def check(ok,kind,detail=None):
 COUNTS[kind]=COUNTS.get(kind,0)+1
 if not ok:raise RuntimeError((kind,detail))
@lru_cache(None)
def words(A):
 if A==0:return ((),)
 return tuple((a,)+w for a in range(1,A+1) for w in words(A-a))
def address(word,q):
 M=1<<sum(word);B=sum(pow(q,len(word)-i-1)*pow(2,sum(word[:i])) for i in range(len(word)))
 return -B*pow(pow(q,len(word),M),-1,M)%M
@lru_cache(None)
def address_map(A,q):
 M=1<<A
 return {(B,w):address(w,q)*(1<<(A-B)) for B in range(1,A+1) for w in words(B)}
def transfer(A,q,p):
 source=address_map(A,q);target=address_map(A,p)
 tau={0:0};tau.update({source[label]:target[label] for label in source})
 return [tau[k] for k in range(1<<A)]
@lru_cache(None)
def decoder(A,q,p,x):
 if A==0:return 0
 M=1<<A;x%=M
 if x%2==0:return 2*decoder(A-1,q,p,x//2)
 return pow(p,-1,M)*(2*decoder(A-1,q,p,(q*x+1)//2)-1)%M
def v2(n):
 n=abs(n);k=0
 while n and n%2==0:n//=2;k+=1
 return k
rows=[]
for A in range(1,9):
 M=1<<A;units=range(1,M,2);paircount=0
 for q in units:
  check(len(set(address_map(A,q).values()))==M-1,'address_partition',(A,q))
  for p in units:
   tau=transfer(A,q,p);c=q*pow(p,-1,M)%M
   criterion=(p-q)%(M//2)==0
   scalar=all(tau[k]==c*k%M for k in range(M))
   real=all(tau[-k%M]==-tau[k]%M for k in range(M))
   check(scalar==criterion,'scalar_iff',(A,q,p))
   check(real==criterion,'real_iff',(A,q,p))
   check(all(tau[(k+M//2)%M]==(tau[k]+M//2)%M for k in range(M)), 'odd_coset_pairing',(A,q,p))
   check(all(tau[k]==decoder(A,q,p,k) for k in range(M)),'independent_decoder',(A,q,p))
   decoder.cache_clear()
   if p!=q and A==v2(p-q)+2:
    check(all(tau[k]==(c*k+(M//2 if k%2 and k%4==q%4 else 0))%M for k in range(M)),'first_defect',(A,q,p))
    check((tau[1]+tau[-1])%M==M//2,'nonreal_cosine',(A,q,p))
   if A<=5:
    multiplicative=all(tau[(i+j)%M]==(tau[i]+tau[j])%M for i in range(M) for j in range(M))
    check(multiplicative==criterion,'full_product_table',(A,q,p))
   paircount+=1
 rows.append({'depth':A,'modulus':M,'multiplier_pairs':paircount})
# Arbitrarily high literal integers test the valuation-based first failure.
for q,p in [(3,5),(3,7),(3,11),(3,19),(3,35),(3,67),(5,133),(101,229)]:
 A=v2(p-q)+2;M=1<<A;tau=transfer(A,q,p);c=q*pow(p,-1,M)%M
 check(all(tau[k]==(c*k+(M//2 if k%2 and k%4==q%4 else 0))%M for k in range(M)),'large_first_failure',(q,p,A))
 check(all(transfer(A-1,q,p)[k]==(q*pow(p,-1,M//2)*k)%(M//2) for k in range(M//2)),'last_good_depth',(q,p,A))
check([transfer(3,3,5)[k] for k in (1,3,5,7)]==[7,1,3,5],'eighth_permutation')
print(json.dumps({'status':'PASS','counts':COUNTS,'total_gates':sum(COUNTS.values()),'universes':rows,'eighth_map':{'1/8':'7/8','3/8':'1/8','5/8':'3/8','7/8':'5/8'}},indent=2,sort_keys=True))
