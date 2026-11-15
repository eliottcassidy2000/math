from itertools import combinations
from model import host, prof

def trans(n):
 return [((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]

def join(A,B):
 l=len(A);r=len(B)
 return [a|(((1<<r)-1)<<l) for a in A]+[b<<l for b in B]

def packet(B,mask):
 h,q,start,w=prof(B);v=[i for i in range(len(B)) if (mask>>i)&1]
 return h,q,sum(start[i] for i in v),sum(w[i][j] for i in v for j in v)

def closer(B,A,mask):
 n=len(B);m=len(A);out=join(B,A)
 for i in range(n+m):
  if i>=n or not (mask>>i)&1:out[i]|=1<<(n+m)
 return out+[mask]

def contraction(B,x,y):
 R=[i for i in range(len(B)) if i not in (x,y)]
 out=[]
 z=0
 for i,u in enumerate(R):
  a=sum(((B[u]>>v)&1)<<j for j,v in enumerate(R))
  if ((B[u]>>x)&1) and ((B[u]>>y)&1):a|=1<<len(R)
  else:z|=1<<i
  out.append(a)
 return out+[z]

pairs=0
for n in range(2,6):
 for code in range(1<<(n*(n-1)//2)):
  B=host(n,code);h,q,start,w=prof(B)
  for x,y in combinations(range(n),2):
   W=w[x][y]+w[y][x]
   assert W>=prof(contraction(B,x,y))[0] and W%2==0 and W>=2
   pairs+=1
print('contraction and parity checks:',pairs)

C=host(3,2); P=host(3,2); A=trans(1)
hA,qA,_,_=packet(A,0)
hC,qC,MC,WC=packet(C,5);hP,qP,_,_=packet(P,0)
hS,qS,_,_=packet(C,0)
for t in range(4):
 S=join(C,trans(t));B=join(join(P,C),S);mask=5<<len(P)
 h,q,M,W=packet(B,mask)
 assert (h,M,W)==(hP*hC*hS,hP*MC*(qS<<t),hP*WC*hS)
 assert packet(closer(B,A,mask),0)[0]==hP*(hC-WC)*hS*hA+(hP*MC*qS*qA<<t)

for d in range(1,6):
 for t in range(4):
  B=trans(d+t+1);marks=1|(1<<d)
  h,q,M,W=packet(B,marks)
  assert (M,W)==((1+(1<<d))<<t,(1<<(d-1))+(d==1))
  assert packet(closer(B,A,marks),0)[0]==1-W+2*M
print('dyadic rays and literal closer checks: pass')
