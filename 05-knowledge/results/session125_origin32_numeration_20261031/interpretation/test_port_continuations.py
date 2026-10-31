from math import isqrt
import json
checks=0
def ck(b):
 global checks
 checks+=1
 if not b:raise RuntimeError('failed')
for k in range(2,1002):
 a,b,r=1,k*k-1,k;c=a+b+2*r
 ck(c==(k+1)**2-1);ck(a*c+1==(a+r)**2);ck(b*c+1==(b+r)**2)
a,b=1,3
for i in range(1000):
 r=isqrt(a*b+1);c=a+b+2*r
 ck(r*r==a*b+1);ck(r==b-a);ck(c==3*b-a);ck(a*c+1==b*b);ck(b*c+1==(c-b)**2)
 a,b=b,c
for n in range(1,10001):
 expected=[k*k-1 for k in range(2,isqrt(n+1)+1)]
 actual=[v for v in range(1,n+1) if v!=1 and isqrt(v+1)**2==v+1]
 ck(actual==expected)
ck(1+15==16);ck(isqrt(1+21)**2!=1+21);ck(1-3*8+8*8==41)
print(json.dumps({'status':'PASS','checks':checks,'branches':1000,'neighbor_orders':[1,10000]},sort_keys=True))
