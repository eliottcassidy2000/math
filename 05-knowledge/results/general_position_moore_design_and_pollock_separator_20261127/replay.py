from itertools import combinations,product
def petersen():
 v=list(combinations(range(5),2));return [[int(not set(a)&set(b)) for b in v] for a in v]
def hs():
 A=[[0]*50 for x in range(50)]
 def join(i,j):A[i][j]=A[j][i]=1
 for i,x in product(range(5),repeat=2):
  join(5*i+x,5*i+(x+1)%5);join(25+5*i+x,25+5*i+(x+2)%5)
  for j in range(5):join(5*i+x,25+5*j+(x+i*j)%5)
 return A
cases=[
 ('Petersen',petersen(),[0,2,3,5,6,9],2,3),
 ('Hoffman-Singleton',hs(),[1,4,6,9,11,13,15,17,20,24,27,28,31,32,35,39,40,43,47,48],3,4),
]
for name,A,S,h,dout in cases:
 n=len(A);d=sum(A[0]);s=len(S);T=[i for i in range(n) if i not in S]
 assert all(sum(row)==d for row in A)
 assert all(sum(A[i][k]*A[k][j] for k in range(n))+A[i][j]-(d-1)*(i==j)==1
            for i in range(n) for j in range(n))
 assert s*(d+h)==n*(h+1)
 assert all(sum(A[i][j] for j in S)==(1 if i in S else dout) for i in range(n))
 assert all(sum(A[i][k]*A[j][k] for k in T)==1-A[i][j] for i,j in combinations(S,2))
 print(name,'sharp GP',s,'blocks',len(T),'block size',dout)
from sympy import Poly, symbols, oo
from math import isqrt
x=symbols('x')
def selector(c,Q,B,p4):
 P=Poly(x**5-Q*x**3-2*B*x**2+(2*Q*Q-p4)*x/4,x)
 D=P.diff()
 if D.gcd(D.diff()).degree()!=0 or D.count_roots(-oo,oo)!=4:return []
 boxes=[{} for _ in range(5)]
 R=isqrt((8*Q)//5)
 for y in range(max(-R,1-c),R+1):
  if D.eval(y)==0:continue
  i=D.count_roots(-oo,y)
  boxes[i][P.eval(y)]=y
 L=sorted(set.intersection(*(set(a) for a in boxes)))
 return [(z,tuple(a[z] for a in boxes)) for z in L]
for c in [10,17,18,20,100]:
 print(c,selector(c,335,-105,143554))
def transitive_max(n,S):
 U=(1<<n)-1
 M=[sum(1<<((i+j)%n) for j in S) for i in range(n)]
 best=();good=bytearray(U+1);good[0]=1
 for a in range(1,U+1):
  v=a
  while v:
   b=v&-v;i=b.bit_length()-1;t=a^b
   if t & ~M[i]==0 and good[t]:
    good[a]=1
    if a.bit_count()>len(best):best=tuple(k for k in range(n) if a>>k&1)
    break
   v^=b
 return best
for S in ([1,2,3,4,5,6,7,8],[16,2,3,4,5,6,7,8],[16,2,3,13,5,11,7,8]):
 print('Fermat tournament',S,transitive_max(17,S))
