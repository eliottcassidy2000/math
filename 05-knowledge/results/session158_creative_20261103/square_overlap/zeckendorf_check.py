from pathlib import Path
from importlib.util import spec_from_file_location, module_from_spec
s=spec_from_file_location('check',str(Path(__file__).with_name('check.py'))); m=module_from_spec(s);s.loader.exec_module(m)

def H(T):return m.dp(T)
def triangle_gadget(n):
 edges={(0,7),(5,9),(9,14)}
 T=[[i!=j and ((i<j)^((min(i,j),max(i,j)) in edges)) for j in range(n)] for i in range(n)]
 return T

def endpoint(T,v):
 n=len(T); table=[[0]*n for i in range(1<<n)]
 for i in range(n):table[1<<i][i]=1
 for mask in range(1,1<<n):
  for x in range(n):
   c=table[mask][x]
   if c:
    for y in range(n):
     if not mask>>y&1 and T[x][y]: table[mask|1<<y][y]+=c
 return table[-1][v]

def zeck(N):
 F=[0,0,1,2]
 while F[-1]<=N:F.append(F[-1]+F[-2])
 out=[]
 for i in range(len(F)-1,1,-1):
  if F[i]<=N:out.append(i);N-=F[i]
 if N:raise Exception('greedy remainder')
 if any(abs(i-j)==1 for i in out for j in out):raise Exception('consecutive')
 return out

def pack(indices):
 p=1
 for i in indices:p*=1+2**(i+2)
 return p

def unpack(p):
 result=[]
 while p>1:
  v=p-1;k=0
  while v%2==0:k+=1;v//=2
  if k<4 or p%(1+2**k):raise ValueError('not an encoding')
  result.append(k-2);p//=1+2**k
 return result

T=triangle_gadget(15)
core=[row[:10] for row in T[:10]]
without=[row[:9] for row in core[:9]]
print('Core H:',H(core),'core-minus-10 ending at 6:',endpoint(without,5))
print('15-vertex square-sum core + first attachment H:',H(T))
if H(core)!=265 or endpoint(without,5)!=0 or H(T)!=4505:raise Exception('attachment failed')
for N in range(10000):
 I=zeck(N)
 if sorted(unpack(pack(I)))!=sorted(I):raise Exception(('roundtrip',N,I))
print('Zeckendorf 2-adic round-trip verified for N=0,...,9999')
for N in range(5):print(N,zeck(N),265*pack(zeck(N)))
