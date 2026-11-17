from collections import Counter,defaultdict
from functools import cache
def G(rows):
 n=len(rows);f=(1<<n)-1;B=defaultdict(list)
 for u in range(n):
  F=f^((1<<(u+1))-1); b=1<<u
  @cache
  def d(s,v):
   z=s^(1<<v)
   if not z:return rows[u]>>v&1
   return sum(d(z,w) for w in range(u+1,n) if z>>w&1 and rows[w]>>v&1)
  k=F
  while k:
   if k.bit_count()>=2 and k.bit_count()%2==0:
    c=sum(d(k,w) for w in range(u+1,n) if k>>w&1 and rows[w]>>u&1)
    if c:B[b].append((k|b,c))
   k=(k-1)&F
 @cache
 def fcn(S):
  if not S:return {(0,0):1}
  b=S&-S;P=Counter({(c,s+1):v for (c,s),v in fcn(S^b).items()})
  for u,k in B[b]:
   if S&u==u:
    for (c,s),v in fcn(S^u).items():P[c+1,s]+=k*v
  return P
 return fcn(f)
def cyclic(n,S):return tuple(sum(1<<((i+j)%n) for j in S) for i in range(n))
def energy(rows):
 n=len(rows); A=[[0 if i==j else 2*(r>>j&1)-1 for j in range(n)] for i,r in enumerate(rows)]
 return sum((sum(a*b for a,b in zip(A[i],A[j]))+1)**2 for i in range(n) for j in range(n) if i!=j)
def P(g):
 s=Counter()
 for (i,j),v in g.items():s[i+j]+=v*2**i
 return tuple(sorted(s.items()))
from itertools import product,combinations,permutations
from fractions import Fraction
from math import comb,prod
from random import Random

def need(ok):
 if not ok:raise ArithmeticError("packing identity mismatch")

def ev(g,x,y):
 return sum(a*x**i*y**j for (i,j),a in g.items())

def induced(t,I):
 return tuple(sum((t[i]>>j&1)<<k for k,j in enumerate(I)) for i in I)

def ordinal(A,B):
 return tuple(a|(((1<<len(B))-1)<<len(A)) for a in A)+tuple(b<<len(A) for b in B)

def substitute(Q,i,A):
 k=len(A); sets=[(j,) if j<i else (j+k-1,) for j in range(len(Q))]
 sets[i]=tuple(range(i,i+k)); T=[]
 for j,r in enumerate(Q):
  e=sum(1<<z for h,I in enumerate(sets) if r>>h&1 for z in I)
  T.extend(e+(a<<i) for a in A) if j==i else T.append(e)
 return tuple(T)

def normal(T):
 n=len(T)
 return all(t.bit_count()==(n-1)//2 for t in T)

def variance_check(T):
 n=len(T);g=G(T);e=energy(T)
 t=n*(n*n-1)//24; tv=3*t//n; m=Fraction(n+1,4)
 a=comb(t,2)-n*comb(tv,2)+comb(n,2)*m*(m-1)/2
 need(g.get((2,n-6),0)==a+Fraction(e,64))
 need(g.get((1,n-5),0)==Fraction(n*(n*n-1)*(n-2)*(n-3),160)-Fraction(e,32))
 # A deletion sum is the partial derivative with x fixed.
 d=Counter()
 for (c,s),v in g.items():
  if s:d[c+s-1]+=s*v*2**c
 D=Counter()
 for v in range(n):
  D.update(dict(P(G(induced(T,[i for i in range(n) if i!=v])))))
 need(d==D)

def orbit_tournament(l,m,sets):
 n=l*m; T=[0]*n
 for i in range(m):
  for j in range(l):
   k=l*i+j
   for d in range(1,(l+1)//2):T[k]|=1<<(l*i+(j+d)%l)
   for b in range(i+1,m):
    for z in range(l):
     x=b*l+z
     if (z-j)%l in sets[i,b]:T[k]|=1<<x
     else:T[x]|=1<<k
 return tuple(T)

def modular_prediction(l,m,sets):
 # Disjoint unoriented odd quotient cycles, grouped with their reverse.
 D={(i,j):len(s)%l for (i,j),s in sets.items()}
 D.update({(j,i):(-v)%l for (i,j),v in list(D.items())})
 bank=[]
 for r in range(3,m+1,2):
  for S in combinations(range(m),r):
   for Q in permutations(S[1:]):
    if Q[0]>Q[-1]:continue
    W=(S[0],)+Q
    V=[D[i,j] for i,j in zip(W,W[1:]+W[:1])]
    b=sum(prod(V[:j]+V[j+1:]) for j in range(r))%l
    bank.append((sum(1<<i for i in S),b))
 def predict(x,y):
  U=y**l-x*pow(2,-1,l); X=x**l-x
  @cache
  def calc(s):
   if not s:return 1
   bit=s&-s
   return (U*calc(s^bit)+X*sum(c*calc(s^z) for z,c in bank if z&bit and z&s==z))%l
  return calc((1<<m)-1)
 return predict

def run():
 flat=cyclic(7,(1,2,4));rough=cyclic(7,(1,2,3))
 need(energy(flat)==0 and energy(rough)==448)
 gp,gr=G(flat),G(rough)
 need((gp[1,2],gr[1,2],gp[2,1],gr[2,1])==(42,28,7,14))
 need(P(gp)==((1,48),(3,112),(5,28),(7,1)))
 need(P(gr)==((1,34),(3,112),(5,28),(7,1)))
 for n in (5,7,9):
  for signs in product((-1,1),repeat=(n-1)//2):
   variance_check(cyclic(n,tuple(j*s for j,s in enumerate(signs,1))))
 # Non-circulant regular controls by triangle flips, preserving degrees.
 t=flat; r=Random(767)
 for _ in range(10):
  C=[(a,b,c) for a,b,c in permutations(range(7),3)
     if a<b and a<c and t[a]>>b&1 and t[b]>>c&1 and t[c]>>a&1]
  a,b,c=r.choice(C);u=list(t)
  for i,j in ((a,b),(b,c),(c,a)):u[i]^=1<<j;u[j]^=1<<i
  t=tuple(u);need(normal(t));variance_check(t)
 A=(0,65,3,7,15,31,61);B=(126,108,120,48,98,64,8)
 g,d=G(A),G(B)
 need(P(g)==P(d)==((3,8),(5,8),(7,1)))
 for x,y in ((2,3),(-1,4),(4,2)):
  need(ev(d,x,y)-ev(g,x,y)==x*(x-2*y)*y)
 need(sum(ev(G(induced(B,[i for i in range(7) if i!=v])),2,1) for v in range(7))-
      sum(ev(G(induced(A,[i for i in range(7) if i!=v])),2,1) for v in range(7))==-4)
 for Q in ((2,4,1),(2,4,9,3),(2,12,9,17,39,15)):
  # The last scaffold is only used if it is a tournament.
  n=len(Q)
  if not all((Q[i]>>j&1)+(Q[j]>>i&1)==1 for i in range(n) for j in range(i)):continue
  need(P(G(substitute(Q,0,A)))==P(G(substitute(Q,0,B))))
 triple=(2,4,1)
 for t in (A,B,flat,(0,),(2,0)):
  g=G(t);h=G(substitute(triple,0,t))
  for x,y in ((2,1),(8,1),(3,4),(4,2),(0,1),(-1,2)):
   need(2*ev(h,x,y)==2*y*y*ev(g,x,y)+x*(ev(g,x+2,y+1)-ev(g,x-2,y-1)))
 for l,m in ((3,2),(3,3),(3,4),(5,2),(7,1)):
  for _ in range(2):
   S={(i,j):{k for k in range(l) if r.randrange(2)} for i,j in combinations(range(m),2)}
   g=G(orbit_tournament(l,m,S))
   pred=modular_prediction(l,m,S)
   for x,y in ((2,1),(3,2),(0,1),(1,0),(-3,-1),(2,2)):
    need(ev(g,x,y)%l==pred(x,y))
 for z,N in enumerate((2697,2625,2721)):
  t=orbit_tournament(3,3,{(0,1):{0},(1,2):{0},(0,2):{z}})
  need(ev(G(t),2,1)==N)
 # A local two-arm cover reader consumes exactly one colour of agreement.
 for U,V in (((0,),(2,0)),(A,rough)):
  gu,gv=G(U),G(V)
  Q=(2,4,1)
  T=substitute(substitute(Q,0,U),len(U),V)
  h=G(T)
  for q in range(1,5):
   a=lambda i:ev(gu,2*i,i)
   b=lambda i:ev(gv,2*i,i)
   want=q*(a(q)*b(q)+(a(q+1)-a(q))*(b(q+1)-b(q))+(a(q-1)-a(q))*(b(q-1)-b(q)))
   need(ev(h,2*q,q)==want)
 K=(38,8,2,5,15,30); L=(58,4,1,6,14,30)
 a,b=G(K),G(L)
 need(a=={(0,6):1,(1,3):4,(1,1):2,(2,0):1})
 need(b=={(0,6):1,(1,3):4,(1,1):4} and P(a)==P(b))
 pieces=[
 [(0,1,3,71,47,79,151,63),(2,32,3,23,39,77,159,63),(0,17,35,39,77,211,143,31)],
 [(2,4,9,19,103,79,143,63),(2,16,3,7,45,79,159,63),(256,1,11,275,327,287,431,319,6),(0,)]
 ]
 L=[tuple(prod(ev(G(t),2*q,q) for t in Z) for q in (1,2,3)) for Z in pieces]
 need(L==[(56551,3390208000,6312366804375),(395857,23731456000,67837683612825)])
 need([b-7*a for a,b in zip(*L)]==[0,0,23651115982200])
 print('Hadamard defect, bivariate collision, diagonal transfer and cyclic symmetry: PASS')

if __name__=='__main__':run()
