from math import isqrt
import json
checks=0
def ck(b):
 global checks
 checks+=1
 if not b:raise RuntimeError('check failed')
def neigh(N):return [k*k-N for k in range(isqrt(N)+1,isqrt(2*N-1)+1)]
for N in range(1,2001):
 a=neigh(N);b=[v for v in range(1,N) if isqrt(v+N)**2==v+N]
 ck(a==b);ck(len(a)==isqrt(2*N-1)-isqrt(N))
for n in range(3,31):
 v=1 if n<=7 else 8 if n<=16 else 16 if n==17 else 18
 ck(sum(u!=v and isqrt(u+v)**2==u+v for u in range(1,n+1))==1)
for origin in (31,32):
 add=lambda a,b:a+b-origin
 mul=lambda a,b:origin+(a-origin)*(b-origin)
 for a in range(origin,origin+20):
  ck(add(a,origin)==a);ck(mul(a,origin+1)==a)
  for b in range(origin,origin+20):
   ck(add(a,b)-origin==(a-origin)+(b-origin));ck(mul(a,b)-origin==(a-origin)*(b-origin))
   for c in range(origin,origin+10):
    ck(add(add(a,b),c)==add(a,add(b,c)));ck(mul(a,add(b,c))==add(mul(a,b),mul(a,c)))
for n in range(32,100001):
 q,r=divmod(n,32);terms=[32]*q if not r else [32]*(q-1)+[32+r]
 ck(sum(terms)==n);ck(all(32<=a<=63 for a in terms))
ck(neigh(32)==[4,17]);ck(neigh(33)==[3,16,31]);ck(all(isqrt(a+b)**2!=a+b for a in neigh(33) for b in neigh(33) if a<b))
print(json.dumps({'status':'passed','checks':checks,'degree_orders':[1,2000],'semiring_origins':[31,32],'apery_orders':[32,100000]},sort_keys=True))
