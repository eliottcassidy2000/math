from math import isqrt
from pathlib import Path
import json
checks=0
def ck(v,s):
 global checks
 checks+=1
 if not v:raise RuntimeError(s)
def components(adj,removed=set()):
 left=set(adj)-removed;out=[]
 while left:
  a=min(left);seen={a};todo=[a];left.remove(a)
  for v in todo:
   for w in adj[v]:
    if w in left:left.remove(w);seen.add(w);todo.append(w)
  out.append(sorted(seen))
 return out
cert=[]
for N in range(2,48,2):
 for d in range(1,N//2+1):
  if N%(2*d):continue
  side=lambda v:(v-1)//d%2
  adj={v:[w for w in range(1,N+1) if side(v)!=side(w) and isqrt(v+w)**2==v+w]for v in range(1,N+1)}
  isolated=[v for v in adj if not adj[v]];leaves=[v for v in adj if len(adj[v])==1];comps=components(adj)
  if isolated:reason='isolated';witness=isolated
  elif len(leaves)>2:reason='too many leaves';witness=leaves
  elif len(comps)>1:reason='disconnected';witness=comps
  else:
   A=[v for v in adj if side(v)==0]
   J={v:[w for w in A if w!=v and isqrt(v+w+d)**2==v+w+d]for v in A}
   marked=[v for v in A if isqrt(2*v+d)**2==2*v+d]
   if not marked:reason='no marked central endpoint';witness=[]
   else:
    ck((N,d)==(36,6),'unexpected residual')
    cut=components(J,{1,14});ck(len(cut)>3,'cut obstruction')
    reason='two vertex cut has four components';witness=cut
  cert.append(dict(N=N,d=d,reason=reason,witness=witness))
w=json.loads((Path(__file__).parent/'palindrome_witness_48_8.json').read_text());p=w['path'];N=w['N'];d=w['d']
tau=lambda v:v+d if (v-1)//d%2==0 else v-d
ck(sorted(p)==list(range(1,N+1)),'permutation')
s=[a+b for a,b in zip(p,p[1:])]
ck(all(isqrt(a)**2==a for a in s),'squares');ck(s==s[::-1],'palindrome')
ck(p[::-1]==[tau(v)for v in p],'twin reversal')
A=[min(v,tau(v))for v in p[:N//2]]
ck(len(set(A))==N//2,'pair transversal')
ck(all(isqrt(a+b+d)**2==a+b+d for a,b in zip(A,A[1:])),'quotient path')
ck(isqrt(2*A[-1]+d)**2==2*A[-1]+d,'central mark')
print(json.dumps(dict(status='PASS',checks=checks,minimal_order=48,lower_certificates=cert),indent=2))
