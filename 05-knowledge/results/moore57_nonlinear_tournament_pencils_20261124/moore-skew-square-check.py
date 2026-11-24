from itertools import combinations
from math import comb, ceil

def add(edges,u,v):
 edges.add(tuple(sorted((u,v))))

def pentagon():
 E=set()
 for i in range(5):add(E,i,(i+1)%5)
 return 5,E

def petersen():
 V=list(combinations(range(5),2))
 return 10,{(i,j) for i,j in combinations(range(10),2) if not set(V[i])&set(V[j])}

def HS():
 E=set()
 for a in range(5):
  for b in range(5):
   add(E,5*a+b,5*a+(b+1)%5)
   add(E,25+5*a+b,25+5*a+(b+2)%5)
   for c in range(5):add(E,5*a+b,25+5*c+(a*c+b)%5)
 return 50,E

def minimum(n,E):
 E=sorted(E);d=2*len(E)//n;r=n*d*d
 for w in range(1<<len(E)):
  b=[0]*n
  for i,(u,v) in enumerate(E):
   s=1 if not (w>>i)&1 else -1
   b[u]+=s;b[v]-=s
  r=min(r,(n*d*d-sum(x*x for x in b))//4)
 return r

def inspect(n,E):
 A=[set() for _ in range(n)]
 for u,v in E:A[u].add(v);A[v].add(u)
 d=len(A[0]);assert n==d*d+1 and all(len(x)==d for x in A)
 Q={};mid={}
 def sign(a,b):return 1 if a<b else -1
 for u,v in combinations(range(n),2):
  m=A[u]&A[v]
  assert len(m)==(0 if v in A[u] else 1)
  if m:
   w,=m;mid[u,v]=w;Q[u,v]=sign(u,w)*sign(v,w)
 negative=sum(v==-1 for v in Q.values())
 assert negative==sum(sum(u<w for u in A[w])*sum(u>w for u in A[w]) for w in range(n))
 st=nt=pt=0
 for u,v,x in combinations(range(n),3):
  pairs=((u,v),(v,x),(u,x))
  if not all(s in Q for s in pairs):continue
  a,b,c=(mid[s] for s in pairs)
  z=Q[u,v]*Q[v,x]*Q[u,x]
  if a==b==c:assert z==1;st+=1
  else:
   assert len({u,v,x,a,b,c})==6
   zc=sign(u,a)*sign(a,v)*sign(v,b)*sign(b,x)*sign(x,c)*sign(c,u)
   assert z==-zc
   nt+=z==-1;pt+=z==1
 assert st==n*comb(d,3) and (nt+pt+st)*6==n*(n-d-1)*(n-2*d-1)
 print(d,n,negative,nt//2,pt//2)

if __name__=='__main__':
 for graph in (pentagon,petersen,HS):inspect(*graph())
 assert minimum(*pentagon())==1
 assert minimum(*petersen())==6
 n,d=3250,57
 f=29*ceil(n*(d-8)/4)
 bound=n*(f-(d*(n-d-1)-comb(d,2))//2)//6
 assert (f,bound)==(1154577,576551625)
 print('minima 1, 6; Moore57 lower bounds',f,bound)
