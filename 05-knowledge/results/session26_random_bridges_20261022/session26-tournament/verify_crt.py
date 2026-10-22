import math,itertools,json

def states(m):
 return [(h,a,p) for h in range(m) for a in range(m) for p in [0,1] if m%2 or (h%2 and a%2==p)]
def check(x,msg):
 if not x:raise RuntimeError(msg)
rows=[]
for m in range(3,20):
 for n in range(m+1,20):
  if math.gcd(m,n)>1:continue
  target={(h%m,a%m,p,h%n,a%n,p) for h,a,p in states(m*n)}
  product={(h,a,p,H,A,P) for h,a,p in states(m) for H,A,P in states(n) if p==P}
  check(target==product,f'CRT {m} {n}');rows.append([m,n,len(target)])
for m in range(3,81):
 obs={( (h+a+p)%m,(h+2*a+2)%m,(h+4*a+4)%m) for h,a,p in states(m)}
 check(len(obs)==len(states(m)),f'strong {m}')
print(json.dumps({'coprime_pairs':len(rows),'crt_fiber_counts':rows,'strong_observer_moduli':'3..80','status':'PASS'},indent=2))
