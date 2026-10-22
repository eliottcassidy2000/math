from itertools import combinations, permutations
from fractions import Fraction as F
from math import comb, factorial

def tournament(n,mask):
 a=[0]*n
 for b,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(i,j) if mask>>b&1 else (j,i)
  a[u]|=1<<v
 return a

def profile(a):
 n=len(a); p=[0]*(n+1)
 for w in permutations(range(n)):
  b=sum(not(a[u]>>v&1) for u,v in zip(w,w[1:])); g=n-1-b
  for s in range(g+1):p[1+b+s]+=comb(g,s)
 return p

def ham(a):
 n=len(a); d=[[0]*n for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 for s in range(1,1<<n):
  for v in range(n):
   if d[s][v]:
    t=a[v]&~s
    while t:
     b=t&-t;w=b.bit_length()-1;t-=b;d[s|b][w]+=d[s][v]
 return sum(d[-1])

def Q(m,q):
 return sum((-1)**(q-j)*comb(q,j)*j**m for j in range(q+1))
def expect(p,m):
 return sum(F(p[r],4**r)*(Q(m,r+1)+4*Q(m,r)+4*Q(m,r-1)) for r in range(1,len(p)))
def channels(p):
 k=len(p)-1;c=[F(0)]*(k+2)
 for r in range(1,k+1):
  for q,w in [(r+1,1),(r,4),(r-1,4)]:
   for j in range(q+1):c[j]+=F(p[r]*w,4**r)*(-1)**(q-j)*comb(q,j)
 return c

def direct(a,m):
 k=len(a);n=k+m;total=0;lo=None;hi=0
 for mask in range(1<<(k*m)):
  b=a.copy()+[sum(1<<v for v in range(k+i+1,n)) for i in range(m)]
  for i in range(k):
   for j in range(m):
    u,v=(i,k+j) if mask>>(i*m+j)&1 else (k+j,i)
    b[u]|=1<<v
  h=ham(b);total+=h;lo=h if lo is None else min(lo,h);hi=max(hi,h)
 return F(total,1<<(k*m)),(lo,hi)

def check(x,msg):
 if not x:raise RuntimeError(msg)
count=0; masks=0;profiles=0
for k in range(1,6):
 for mask in range(1<<(k*(k-1)//2)):
  a=tournament(k,mask);p=profile(a);profiles+=1
  c3=sum(all(sum(bool(a[u]>>v&1) for v in S if v!=u)==1 for u in S) for S in combinations(range(k),3))
  check(p[k]==factorial(k),'top profile')
  if k>=2:check(p[k-1]==comb(k,2)*factorial(k-1),'second profile')
  if k>=3:check(p[k-2]==factorial(k-2)*(comb(k,3)+2*c3+3*comb(k,4)),'third profile')
  c=channels(p)
  check(c[k+1]==F(factorial(k),4**k),'top channel')
  if k>=2:check(c[k]==F((k+1)*factorial(k),4**k),'second channel')
  for m in range(9):check(expect(p,m)==sum(c[j]*j**m for j in range(k+2)),'channel transform')
  if k<=4:
   for m in range(4 if k<=3 else 3):
    v,span=direct(a,m);check(v==expect(p,m),f'direct {k} {mask} {m}');count+=1;masks+=1<<(k*m)
print('PASS profiles',profiles,'all labelled cores of orders 1..5')
print('PASS independent attachment-mask Hamiltonian subset DP',count,'core/base pairs',masks,'attachments')
for mask,name in [(7,'T3'),(5,'C3')]:
 a=tournament(3,mask);p=profile(a)
 print(name,'ordered profile',p[1:],'channels',list(map(str,channels(p))),'expectations m0..6',list(map(str,[expect(p,m) for m in range(7)])))
a=tournament(3,5);b=tournament(3,7)
for m in range(20):check(expect(profile(a),m)-expect(profile(b),m)==F(2**m+2+0**m,2),'C3 difference')
print('PASS C3-T3 exact finite-head formula m0..19')
# Nonisomorphic SCC order collision, distinguished by sorted outdegrees.
cycle=tournament(3,5)
left=[14]+[x<<1 for x in cycle]
right=[x|8 for x in cycle]+[0]
check(profile(left)==profile(right),'collision')
print('PASS nonisomorphic same-profile hostile', sorted(x.bit_count() for x in left),sorted(x.bit_count() for x in right),'profile',profile(left)[1:])
print('Hostile actual-H versus expectation C3+m1:',direct(a,1))
print('PASS exact rational arithmetic; no probabilistic sampling; no inherited strong/class filter')
recovered=0
for k in range(1,6):
 for mask in range(1<<(k*(k-1)//2)):
  p=profile(tournament(k,mask));h=(k+1)//2;low=[0]*(k+1)
  for m in range(h):
   low[m+1]=(expect(p,m)-expect(low,m))*F(4**m,factorial(m))
  B=[0]*k
  for b in range(h):
   B[b]=low[b+1]-sum(comb(k-1-j,b-j)*B[j] for j in range(b))
   B[k-1-b]=B[b]
  got=[0]+[sum(comb(k-1-b,r-1-b)*B[b] for b in range(r)) for r in range(1,k+1)]
  check(got==p,'half-prefix tomography');recovered+=1
print('PASS first ceil(k/2) expectation tomography',recovered,'labelled cores')
from fractions import Fraction

def det(M):
 M=[list(map(F,row)) for row in M];d=F(1);n=len(M)
 for i in range(n):
  q=next((j for j in range(i,n) if M[j][i]),None)
  if q is None:return F(0)
  if q!=i:M[i],M[q]=M[q],M[i];d=-d
  pivot=M[i][i];d*=pivot
  for j in range(i+1,n):
   w=M[j][i]/pivot
   for h in range(i+1,n):M[j][h]-=w*M[i][h]
 return d
reflection_count=0
for k in range(1,6):
 for mask in range(1<<(k*(k-1)//2)):
  c=channels(profile(tournament(k,mask)))
  reflected=[sum(c[j]*comb(j,i)*(-2)**(j-i)*(-1)**i for j in range(i,len(c))) for i in range(len(c))]
  check(reflected==[(-1)**(k+1)*a for a in c],'spectral reflection')
  for order in range(3 if k%2==0 else 2):
   check(sum(c[j]*factorial(j)/factorial(j-order)*(-1)**(j-order) for j in range(order,len(c)))==0,'root multiplicity')
  reflection_count+=1
print('PASS polynomial reflection on',reflection_count,'cores')
p=profile(tournament(4,63));c=channels(p)
check(p[1:]==[1,14,36,24],'T4 profile')
check(c==[F(j,32) for j in [-5,-9,6,22,15,3]],'T4 channels')
H=[[expect(p,i+j+1) for j in range(5)] for i in range(5)]
d=det(H);vand=1
for i in range(1,6):
 for j in range(i+1,6):vand*= (j-i)**2
product=F(vand)
for j in range(1,6):product*=c[j]*j
check(d==product and d<0,'negative shifted Hankel')
print('PASS T4 negative tail atom a1=',c[1],'; shifted 5x5 Hankel det=',d)
