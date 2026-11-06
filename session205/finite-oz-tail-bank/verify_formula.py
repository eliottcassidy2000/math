import itertools,json

def step(s,a):
 h,r,e,d=s
 if a=='O':return h,2*r,0,d
 if a=='Z':return h+r-e,2*r-d,d,d
 return h+r-e,h+r,h-e,h

def Q(m):
 if not m:return 0
 return Q(m//2)+(m//2+(m//2)%2 if m%2 else 0)
seeds={(1,2,0,1)}
for n in range(4):
 for w in itertools.product('OZJ',repeat=n):
  s=(1,2,0,1)
  for a in w:s=step(s,a)
  seeds.add(s)
checks=0
for s in seeds:
 if step(step(s,'J'),'Z')!=step(step(s,'O'),'J'):raise RuntimeError('JZ=OJ')
 h,r,e,d=s
 for n in range(10):
  for w in itertools.product('OZ',repeat=n):
   z=s;m=0
   for i,a in enumerate(w):
    z=step(z,a);m+=(a=='Z')*2**i
   if z[0]!=h+r*m-d*Q(m)-e*(m%2):raise RuntimeError('formula')
   if m and not 0<=Q(m)<=m-1:raise RuntimeError('bound')
   if m and z[0]<h+(r-d)*m+d-e:raise RuntimeError('positive-input lower bound')
   if z[0]<h+(r-d)*m:raise RuntimeError('all-input lower bound')
   checks+=1
print(json.dumps({'status':'PASS','physical_seeds':len(seeds),'tail_checks':checks,'max_tail_length':9,'JZ_OJ_response_identity_checks':len(seeds)},indent=2))
