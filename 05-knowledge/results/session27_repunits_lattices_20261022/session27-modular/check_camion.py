import json
from math import gcd

def strong_mod(N):
 t=[pow(2,m*(m-1)//2,N) for m in range(N+1)]
 s=[0]*(N+1);row=[1]
 for m in range(1,N+1):
  row=[1]+[(row[k-1]+row[k])%N for k in range(1,m)]+[1]
  s[m]=(t[m]-sum(row[k]*s[k]*t[m-k] for k in range(1,m)))%N
 return s[N]
r={}
for n in [21,1729]:
 s=strong_mod(n); c=pow(2,(n-1)//2,n)
 assert s==c
 r[n]={'s(n) mod n':s,'c(n) mod n':c,'Camion':True}
assert pow(2,20,21)!=1
assert 1729==2**6*3**3+1
r['21_not_Carmichael_witness']={'a':2,'2^20 mod21':pow(2,20,21)}
print(json.dumps(r,indent=2))
