import math, json
N=1000000
Q=[0]*(N+1)
for m in range(1,N+1):
 t=m//2
 Q[m]=Q[t]+(t+t%2 if m%2 else 0)
def step(s,a):
 h,r,e,d=s
 if a=='O':return h,2*r,0,d
 if a=='Z':return h+r-e,2*r-d,d,d
 return h+r-e,h+r,h-e,h
seeds={(1,2,0,1)}
for _ in range(4):seeds|={step(s,a) for s in tuple(seeds) for a in 'OZJ'}
rows=[]
for s in sorted(seeds):
 h,r,e,d=s; g=r-d
 coefficient=math.log(r/g)/d
 counts=[]
 for X in [10000,100000,1000000]:
  count=sum(h+r*m-d*Q[m]-e*(m%2)<=X for m in range(max(0,(X-h)//g)+1))
  counts.append([X,count,coefficient*X,(count-coefficient*X)/math.sqrt(X)])
 rows.append({'seed':s,'coefficient':coefficient,'counts':counts})
for m in range(10000):
 bits=[int(c) for c in bin(m)[:1:-1]] if m else []
 q=sum(b*(m//(2**(j+1))) for j,b in enumerate(bits))+sum(bits[j]*bits[j+1] for j in range(len(bits)-1))
 assert q==Q[m]
 u=sum(b/(2**(j+1)) for j,b in enumerate(bits))
 assert abs(q-m*u)<=len(bits)
print(json.dumps({'status':'PASS','digit_checks':10000,'physical_seed_count':len(seeds),'rows':rows},indent=2))
