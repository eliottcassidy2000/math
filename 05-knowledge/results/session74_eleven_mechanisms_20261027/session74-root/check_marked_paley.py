from math import gcd,isqrt
import json

def need(b,m):
    if not b: raise RuntimeError(m)

def prime(p):return p>=2 and all(p%d for d in range(2,isqrt(p)+1))
rows={}
for p in range(3,100):
    if not prime(p):continue
    chi=lambda x:0 if x%p==0 else (1 if pow(x%p,(p-1)//2,p)==1 else -1)
    S=sum(a*chi(a) for a in range(p))
    R=sum(x*x%p for x in range(p))
    need(R==p*(p-1)//2+S,('moment',p))
    values=[sum(a*chi(a-v) for a in range(p)) for v in range(p)]
    need(values==[S-p*sum(chi(y) for y in range(p-v,p)) for v in range(p)],('carry',p))
    rows[p]=values
q=341;strata={}
for x in range(q):
    c=q//gcd(x,q); rec=strata.setdefault(c,[0,0,0]);rec[0]+=1;rec[1]+=x*x%q;rec[2]+=(-x*x)%q
K={1};todo=[1]
while todo:
    x=todo.pop()
    for g in (4,81):
        y=x*g%q
        if y not in K:K.add(y);todo.append(y)
qr={x*x%q for x in range(q) if gcd(x,q)==1}
need(K==qr and len(K)==75,'jointclock')
need(strata=={1:[1,0,0],11:[10,1364,2046],31:[30,6138,4092],341:[300,49104,53196]},'strata')
need(rows[11]==[-11,0,-11,0,11,22,11,0,-11,0,-11],'marked11')
print(json.dumps({'status':'PASS','odd_primes_checked':len(rows),'eleven_moments':rows[11],'q341_strata':strata,'joint_clock_order':len(K)},sort_keys=True,indent=2))
