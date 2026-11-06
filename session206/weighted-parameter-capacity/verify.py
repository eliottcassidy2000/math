import json,math
N=1000000
q=[0]*(N+1)
for m in range(1,N+1):
 t=m//2;q[m]=q[t]+(t+t%2 if m%2 else 0)
def step(s,a):
 h,r,e,d=s
 if a=='O':return h,2*r,0,d
 if a=='Z':return h+r-e,2*r-d,d,d
 return h+r-e,h+r,h-e,h
seeds={(1,2,0,1)}
for _ in range(3):seeds|={step(s,a) for s in tuple(seeds) for a in 'OZJ'}
fiber_checks=0
for s in sorted(seeds):
 h,r,e,d=s
 for p in range(1,16,2):
  g=math.gcd(r,d,p)
  for b in (0,1):
   counts=[0]*p
   for u in range(p):
    for v in range(p):counts[(h-e*b+r*u-d*v)%p]+=1
   for a in range(p):
    expected=p*g if (a-h+e*b)%g==0 else 0
    if counts[a]!=expected:raise RuntimeError('fiber count')
    fiber_checks+=1
rows=[]
for word in ('','JJ','OZJ'):
 s=(1,2,0,1)
 for a in word:s=step(s,a)
 h,r,e,d=s;c=[math.log(r/(r-d/2))/d,math.log((r-d/2)/(r-d))/d]
 counts=[[0]*3 for _ in (0,1)]
 for m in range(max(0,(N-h)//(r-d))+1):
  H=h+r*m-d*q[m]-e*(m%2)
  if H<=N:counts[m%2][H%3]+=1
 rows.append({'word':word or 'root','state':s,'parity_coefficients':c,'finite_counts_mod3':counts})
print(json.dumps({'status':'PASS','physical_states':len(seeds),'fiber_checks':fiber_checks,'cutoff':N,'rows':rows},indent=2))
