from collections import defaultdict
import json
checks=0
def gate(x,m):
 global checks
 checks+=1
 if not x:raise RuntimeError(m)
def v2(n):return (n&-n).bit_length()-1
def decode(n,q,A):
 # Literal shortcut dynamics, stopping as soon as the prescribed valuation clock is reached.
 word=[];clock=0
 while clock<A:
  z=q*n+1;a=min(v2(z),A-clock);word.append(a);clock+=a;n=z>>a
 return tuple(word)
cache={}
def bank(q,A):
 key=(q%(1<<A),A)
 if key not in cache:
  M=1<<A;cache[key]={decode(n,q,A):n for n in range(1,M,2)}
  gate(len(cache[key])==M//2,('word bijection',q,A))
 return cache[key]
def transport(q,p,A):
 M=1<<A;out=[0]
 for k in range(1,M):
  v=v2(k);B=A-v;n=k>>v;word=decode(n,q,B)
  out.append(bank(p,B)[word]<<v)
 return out
rows=[]
for A in range(1,8):
 M=1<<A;count=0;good=0;first=0
 for p in range(1,M,2):
  for q in range(1,M,2):
   tau=transport(q,p,A);c=q*pow(p,-1,M)%M
   gate(sorted(tau)==list(range(M)),('permutation',A,p,q))
   multiplicative=all(tau[k]==k*tau[1]%M for k in range(M))
   real=all(tau[-k%M]==-tau[k]%M for k in range(M))
   condition=(p-q)%(1<<(A-1))==0
   gate(multiplicative==condition,('algebra iff',A,p,q))
   gate(real==condition,('real iff',A,p,q))
   gate(all(tau[k+M//2]==(tau[k]+M//2)%M for k in range(M//2)),('odd quotient descends',A,p,q))
   if condition:gate(all(tau[k]==c*k%M for k in range(M)),('explicit scalar',A,p,q))
   if p!=q and A==v2(p-q)+2:
    first+=1
    for k in range(M):
     expected=(c*k+(M//2 if k%2 and k%4==q%4 else 0))%M
     gate(tau[k]==expected,('first defect',A,p,q,k))
    gate(tau[M-1]==(-tau[1]+M//2)%M,('cosine imaginary on odds',A,p,q))
   # Literal full product tables for representative boundary pairs.
   if p==3 and q in (3,5,7) and M>q:
    product_ok=all(tau[(j+k)%M]==(tau[j]+tau[k])%M for j in range(M) for k in range(M))
    gate(product_ok==condition,('full product table',A,p,q))
   count+=1;good+=condition
 rows.append({'depth':A,'all_odd_residue_pairs':count,'multiplicative_pairs':good,'first_failure_pairs':first})
tau=transport(3,5,3);gate([tau[k] for k in (1,3,5,7)]==[7,1,3,5],'odd eighth map')
print(json.dumps({'status':'PASS','gates':checks,'decoder':'literal shortcut valuation stopping clock','rows':rows,'three_to_five_eighth_permutation':tau},indent=2))
