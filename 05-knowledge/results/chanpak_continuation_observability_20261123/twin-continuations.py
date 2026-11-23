from collections import defaultdict
from itertools import combinations,permutations
from functools import cache
A8=(16,1,3,7,14)
A10=(20,1,2,7,14)
def sub(A,I):
 return tuple(sum((A[i]>>k&1)<<j for j,k in enumerate(I)) for i in I)
@cache
def h(A):
 if not A:return 1
 return sum(all(A[i]>>j&1 for i,j in zip(P,P[1:])) for P in permutations(range(len(A))))
def signature(A,S):
 U=[0,0];V=[0,0];W=[0,0];K=[0,0]
 for k in range(1,len(A)+1):
  for P in permutations(range(len(A)),k):
   if not all(A[i]>>j&1 for i,j in zip(P,P[1:])):continue
   a,b=((S>>P[0])&1),1-((S>>P[-1])&1)
   w=h(sub(A,tuple(i for i in range(len(A)) if i not in P)))
   U[k%2]+=a*w;V[k%2]+=b*w;W[k%2]+=a*b*w;K[k%2]+=w
 return (W[0],U[1],U[0],V[1],V[0]),tuple(K)

def apex(A,S):
 n=len(A)
 return tuple(x+(0 if S>>i&1 else 1<<n) for i,x in enumerate(A))+(S,)

def tt(n):
 return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))

def join(A,B):
 n,m=len(A),len(B)
 return tuple(a|(((1<<m)-1)<<n) for a in A)+tuple(b<<n for b in B)

def cyclic(A,k):
 n=0
 for I in combinations(range(len(A)),k):
  for W in permutations(I[1:]):
   P=(I[0],)+W
   n+=all(A[i]>>j&1 for i,j in zip(P,P[1:]+P[:1]))
 return n

def autom(A):
 return sum(all(A[i]>>j&1==A[P[i]]>>P[j]&1
                for i in range(len(A)) for j in range(len(A))) for P in permutations(range(len(A))))

def probes(A,S):
 z=h(apex(A,S))
 l=[h(apex(join(tt(j),A),(S<<j)+(1<<j)-1)) for j in (1,2)]
 r=[h(apex(join(A,tt(j)),S)) for j in (1,2)]
 U=(r[0]-z)//2;V=(l[0]-z)//2
 return (z-9)//2,U,(r[1]-z)//2-2*U,V,(l[1]-z)//2-2*V

if __name__=="__main__":
 D=defaultdict(list)
 for k,A in ((8,A8),(10,A10)):
  assert h(A)==9 and autom(A)==1 and (cyclic(A,3),cyclic(A,5))==(3,1)
  for S in range(32):
   v,K=signature(A,S)
   assert K==(33,42) and probes(A,S)==v
   D[v].append((k,S))
  atlas={r:sorted((h(apex(A,S)),S) for S in range(32) if S.bit_count()==r) for r in range(6)}
  print("core",k,": (height, mask) by apex outdegree",atlas)
 print("distinct signatures",len(D))
 assert len(D)==60
 print("collisions",[(v,vals) for v,vals in D.items() if len(vals)>1])

 expected={((8,0),(10,0)):(27,31),((8,31),(10,31)):(17,19),
           ((10,10),(10,20)):(65,77),((10,11),(10,21)):(59,65)}
 observed={}
 for vals in D.values():
  if len(vals)==1:continue
  observed[tuple(vals)]=tuple(h(apex(apex({8:A8,10:A10}[a],S),1)) for a,S in vals)
 assert observed==expected
 print("sixth probe",observed)
