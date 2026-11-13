from itertools import permutations
from functools import cache
BANK=[
 ((20,1,2,7,14),(9,10,2),(-5,4,2),(1,13)),
 ((44,17,2,6,13,30),(25,28,6),(-13,10,6),(5,9)),
 ((36,1,2,7,15,30),(15,22,6),(-11,4,6),(7,11)),
 ((34,48,3,7,13,28),(41,41,10),(-11,11,10),(3,15)),
]
def H(R,L):
 n=len(R)
 V=[(R[0]<<(L-1))|sum(1<<j for j in range(i+1,L)) for i in range(L)]
 V+=[((R[i]&~1)<<(L-1))|(((1<<L)-1) if R[i]&1 else 0) for i in range(1,n)]
 @cache
 def f(S,j):
  if not S:return 1
  return sum(f(S^(1<<i),i) for i in range(len(V)) if S&V[j]>>0&(1<<i))
 return sum(f(((1<<len(V))-1)^(1<<j),j) for j in range(len(V)))

def reach(R):
 S=1
 for _ in R:
  for i in range(len(R)):
   if S>>i&1:S|=R[i]
 return S==(1<<len(R))-1
for R,a,coeff,classes in BANK:
 n=len(R); O=[i for i in range(1,n) if R[0]>>i&1]; I=[i for i in range(1,n) if R[i]&1]
 assert min(len(O),len(I))==2
 assert any(R[i]>>b&1 and R[j]>>c&1 for i,j in permutations(O,2) for b,c in permutations(I,2))
 assert reach(R) and reach([sum(1<<j for j,v in enumerate(R) if v>>i&1) for i in range(n)])
 h1,h2,h3=[H(R,L) for L in (1,2,3)]
 assert (h1,h2,h3)==(a[0],a[0]+2*a[1],a[0]+6*a[1]+6*a[2])
 A,B,C=coeff
 assert (A,B,C)==(a[0]-2*a[1]+3*a[2],a[1]-3*a[2],a[2])
 for L in range(1,8):assert H(R,L)==A+B*2**L+C*3**L
 assert tuple(sorted(((A+C)%16,(A+3*C)%16)))==classes
 print(R, coeff, classes)
for k in range(4,15):
 mod=2**k
 values={(A+C*pow(3,L,mod))%mod for _,_,(A,B,C),_ in BANK for L in range(k,k+2**(k-3))}
 assert values==set(range(1,mod,2))
print('all odd residues through precision 14')
def join(X,Y):
 n=len(X);m=len(Y)
 return tuple(x|(((1<<m)-1)<<n) for x in X)+tuple(y<<n for y in Y)
def cywrap(A):
 n=len(A)
 return tuple(a|(1<<n) for a in A)+(1<<(n+1),(1<<n)-1)
def rev(A):
 return tuple(sum(1<<j for j,a in enumerate(A) if a>>i&1) for i in range(len(A)))
def dp(A):
 n=len(A)
 @cache
 def p(S,v):
  if S==1<<v:return 1
  return sum(p(S^(1<<v),j) for j in range(n) if S>>j&1 and A[j]>>v&1)
 @cache
 def h(S):return sum(p(S,v) for v in range(n) if S>>v&1) if S else 1
 F=(1<<n)-1
 return h(F),sum(h(S)*h(F^S) for S in range(F+1))
C3=(2,4,1); C4=(10,4,1,6); C5=tuple((1<<((i+1)%5))|(1<<((i+2)%5)) for i in range(5))
for X,ans in zip([(0,),C3,C4,C5],[1,3,5,15]):
 assert dp(X)[0]==ans
 T=C3
 S=join(join(join(join(T,(0,)),X),(0,)),rev(T))
 h,q=dp(S)
 assert h==9*ans and q==4*dp(T)[1]**2*dp(X)[1] and dp(cywrap(S))[0]==h+q
print('mirror sandwich and cyclic wrap controls passed')
