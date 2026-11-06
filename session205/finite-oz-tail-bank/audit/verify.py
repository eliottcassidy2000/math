import itertools,json

def step(s,a):
 h,r,e,d=s
 if a=='O':return h,2*r,0,d
 if a=='Z':return h+r-e,2*r-d,d,d
 return h+r-e,h+r,h-e,h

def q(n):
 if not n:return 0
 k=n//2
 return q(k)+(k+k%2 if n%2 else 0)
seeds={(1,2,0,1)};layer=seeds
for i in range(5):
 layer={step(s,a) for s in layer for a in 'OZJ'};seeds|=layer
cases=0
for s in seeds:
 h,r,e,d=s
 if step(step(s,'J'),'Z')!=step(step(s,'O'),'J'):raise RuntimeError('JZ identity')
 for m in range(512):
  b=[(m>>i)&1 for i in range(m.bit_length())]
  pair=sum(b[i]*b[j]*(2**(i-1-j)+(i==j+1)) for i in range(len(b)) for j in range(i))
  qm=q(m)
  if qm!=pair or qm%2 or not(0<=qm<=max(0,m-1)):raise RuntimeError('Q')
  actual=s
  for bit in b:actual=step(actual,'Z' if bit else 'O')
  out=h+r*m-d*qm-e*(m%2)
  if actual[0]!=out or out<h+(r-d)*m:raise RuntimeError('orbit/bound')
  if m and out<h+(r-d)*m+d-e:raise RuntimeError('positive-m bound')
  if m>(out-h)//(r-d):raise RuntimeError('inverse')
  cases+=1
print(json.dumps({'status':'PASS','physical_seed_states':len(seeds),'tail_cases':cases,'max_encoded_tail':511,'JZ_OJ_identity':True},indent=2))
