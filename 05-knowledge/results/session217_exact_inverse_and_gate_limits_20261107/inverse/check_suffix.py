import random,json

def f(z):
 bits=[i for i in range(z.bit_length()) if z>>i&1]
 return sum((1<<(j-i-1))+(j==i+1) for i in bits for j in bits if i<j)
def rev(z,b):
 return int(f'{z:0{b}b}'[::-1],2) if b else 0
def solve(h,A,E,tau,b,N):
 X=(2<<b)*A+h-tau-N; B=2*h
 t=2
 while A*(4**t)<24*h*(1<<b):t+=1
 t=min(t,b)
 states=[(A,E,0,0)]
 for i in range(t):
  ns=[]
  for K,S,C,l in states:
   ns.append((2*K,0,C,l))
   ns.append((2*K+h,h,C+K+S-B*(1<<(b-i-1)),l+(1<<i)))
  states=ns
 out=[];maxwidth=0
 for K,S,C,l in states:
  m=b-t;x=X-C
  lo=max(0,-(-(x-h)//(K+h)))
  hi=min((1<<m)-1,(x+B*((1<<m)-1))//K)
  assert t==b or hi-lo<=0,(h,A,E,tau,b,N,t,lo,hi)
  for u in range(lo,hi+1):
   z=l+(u<<t)
   if A*z+h*f(z)-B*rev(z,b)+E*(z&1)==X:out.append(z)
 return sorted(out)
r=random.Random(217)
checks=0
for h in [5,9]:
 for b in range(11):
  for case in range(20):
   A=r.randrange(h,500);E=r.randrange(h+1);tau=r.randrange(A+1)
   vals={}
   for z in range(1<<b):
    N=(2<<b)*A+h-tau-(A*z+h*f(z)-2*h*rev(z,b)+E*(z&1))
    vals.setdefault(N,[]).append(z)
   for N in r.sample(list(vals),min(12,len(vals)))+[r.randrange(1,10000)]:
    assert solve(h,A,E,tau,b,N)==vals.get(N,[])
    checks+=1
print(json.dumps({'status':'PASS','fixed_suffix_inverse_tests':checks,'b_max':10,'h':[5,9],'note':'Synthetic coordinate controls, not graph census.'},indent=2))
