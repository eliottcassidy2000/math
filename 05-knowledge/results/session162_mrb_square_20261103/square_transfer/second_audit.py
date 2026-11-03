from itertools import permutations,combinations
from functools import lru_cache
from math import isqrt,prod

def covers_brute(d,r):
 n=d+1
 # all paths are directed in ascending order except d->0
 A=lambda a,b: (a<b and (a,b)!=(0,d)) or (a,b)==(d,0)
 total=0
 for p in permutations(range(n)):
  for cuts in combinations(range(n+r-1),r-1):
   lens=[];pre=-1
   for t in cuts+(n+r-1,):lens.append(t-pre-1);pre=t
   offset=0;ok=True
   for l in lens:
    for i in range(offset,offset+l-1):
     if not A(p[i],p[i+1]):ok=False;break
    offset+=l
   total+=ok
 return total
for d in range(1,6):
 for r in range(1,5):
  got=covers_brute(d,r);want=r**(d+1)+r*((r+1)**(d-1)-(r-1)**(d-1))
  if got!=want:raise ValueError(('cover profile',d,r,got,want))
print('ordered path covers verified by permutations/weak compositions, d=1..5,r=1..4')

def path_data(p,start):
 n=len(p)
 if sorted(p)!=list(range(1,n+1)):raise ValueError('not permutation')
 if any(isqrt(a+b)**2!=a+b for a,b in zip(p,p[1:])):raise ValueError('not square')
 prev=start;es=[];ds=[]
 for a,b in zip(p,p[1:]):
  d=abs(a-b);s=1
  while True:
   s+=1
   if (s*s-d)%2==0 and (s*s-d)//2>prev:break
  x=(s*s-d)//2;y=(s*s+d)//2
  if not (x>prev and y-x==d and isqrt(x+y)**2==x+y):raise ValueError('embedding')
  es.append((x,y));ds.append(d);prev=y
 return es,ds

def factors(ds,g):
 return prod(1 if d==1 else 1+2**(d-1) for d in ds), 2**g*prod(2**(d+1)+2*(3**(d-1)-1) for d in ds)

p=[8,1,15,10,6,3,13,12,4,5,11,14,2,7,9]
es,ds=path_data(p,5)
expected=[(9,16),(25,39),(58,63),(70,74),(83,86),(93,103),(112,113),(124,132),(144,145),(159,165),(179,182),(194,206),(218,223),(241,243)]
if es!=expected:raise ValueError(('intervals',es))
g=es[-1][1]-5-sum(d+1 for d in ds)
P,R=factors(ds,g)
if (g,P,R)!=(143,464858633025922867875,2286430711049262759319343471689019920235263576776074977652877326066121533571465216):raise ValueError(('15 path',g,P,R))
print('15 path: all 14 embeddings, g, P and R match exact printed claims')

A=[3,22,14,2,23,13,12,4,21,15,10,6,19,17,8,1,24,25,11,5,20,16,9,7,18]
B=[18,7,9,16,20,5,11,25,24,1,8,17,19,6,10,15,21,4,12,13,3,22,14,2,23]
ea,da=path_data(A,5);eb,db=path_data(B,5)
if sorted(da)!=sorted(db) or set(ea)==set(eb):raise ValueError('collision')
if (ea[-1][1],eb[-1][1])!=(846,623):raise ValueError('embeddings end')
ga=846-5-sum(d+1 for d in da);gb=846-5-sum(d+1 for d in db)
PA,RA=factors(da,ga);PB,RB=factors(db,gb)
if (ga,gb)!=(592,592) or (PA,RA)!=(PB,RB):raise ValueError('padded collision')
if PA!=10918691908903704879820836527116367337011525929450009908515625:raise ValueError('collision P')
print('25 paths: square edges, distinct embedded edges, identical g, P, R after padding')

def bitmatrix(n,back):
 A=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   if (j,i) in back:A[j]|=1<<i
   else:A[i]|=1<<j
 return A

def path_recursive(A):
 n=len(A);allv=(1<<n)-1
 @lru_cache(None)
 def go(available,v):
  if not available:return 1
  z=0;c=A[v]&available
  while c:
   b=c&-c;c-=b;z+=go(available^b,b.bit_length()-1)
  return z
 return sum(go(allv^(1<<i),i) for i in range(n))

def strong(A):
 n=len(A)
 for start in range(n):
  seen=1<<start;cur=seen
  while cur:
   b=cur&-cur;cur-=b;j=b.bit_length()-1
   y=A[j]&~seen;seen|=y;cur|=y
  if seen!=(1<<n)-1:return False
 return True
for j,want in enumerate((6435,9945,13455)):
 back={(5,i) for i in range(5)}|{(6,5),(8,6),(15,8)}
 if j>=1:back.add((2,0))
 if j>=2:back.add((3,0))
 M=bitmatrix(16,back);got=path_recursive(M)
 if got!=want or not strong(M):raise ValueError(('16',j,got))
print('16 vertex graphs: strong, recursive deletion counts 6435 9945 13455')
