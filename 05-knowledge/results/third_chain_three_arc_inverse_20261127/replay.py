from itertools import permutations, combinations
from collections import Counter
def patt(m):
 for B in combinations(combinations(range(m),2),3):
  if len(set(sum(B,())))==m: yield B
def comp(m,B):
 v=Counter()
 for p in permutations(range(m)):
  if any(i>j and (j,i) not in B for i,j in zip(p,p[1:])): continue
  edges=list(zip((-1,)+p,p+(m,)))
  req=[k for k,e in enumerate(edges) if e in B]
  for k in range(1<<len(req)):
   drop={req[j] for j in range(len(req)) if k>>j&1}
   h=tuple(sum(i<=g<j for n,(i,j) in enumerate(edges) if n not in drop) for g in range(m-1))
   v[h]+=(-1)**len(drop)
 return {h:c for h,c in v.items() if c}
from itertools import product
from math import prod

def emit():
 print(62)
 for m in range(3,7):
  for B in patt(m):
   u=comp(m,B)
   print(m,len(u),*sum(B,()))
   for h,c in sorted(u.items()): print(c,*h)

def gaps(m,B,cap,prev=()):
 g=len(prev)
 if g==m-1:
  yield prev
  return
 arcs=[(i,j) for i,j in B if i<=g<j]
 top=min((cap-sum(prev[i:g]) for i,j in arcs),default=0)
 for a in range(int((g,g+1) in B),top+1):
  yield from gaps(m,B,cap,prev+(a,))

def exact(C,L):
 return sum(c*prod(b**t for b,t in zip(v,L)) for v,c in C.items())

def verify189():
 totals=[];near=Counter()
 for m in range(3,7):
  count=0
  for B in patt(m):
   C=comp(m,B)
   for L in gaps(m,B,7):
    h=exact(C,L);count+=1
    assert h!=189
    if 181<=h<=197:near[h]+=1
  totals.append(count)
 assert totals==[21,2227,13902,19885]
 print("189 exclusion:",totals,"total",sum(totals),"nearby",sorted(near.items()))

def ham(n,B):
 N=1<<n
 A=[sum(1<<j for j in range(n) if i!=j and ((i<j)^((min(i,j),max(i,j)) in B))) for i in range(n)]
 for v in range(n):
  u=1<<v
  while True:
   w=u
   for i in range(n):
    if u>>i&1:w|=A[i]
   if w==u:break
   u=w
  assert u==N-1
 D=[[0]*n for _ in range(N)]
 for i in range(n):D[1<<i][i]=1
 for s in range(1,N):
  for i,h in enumerate(D[s]):
   if not h:continue
   a=A[i]&~s
   while a:
    v=a&-a;a-=v;D[s|v][v.bit_length()-1]+=h
 return sum(D[-1])

def form(a,b,c,d,e):
 return 1+2**(b+c+d+2)+2**(a+b+1)+2**(a+c+d+2)*3**b+2**(d+e+1)+2**(b+c+e+2)*3**d+2**(a+b+d+e+2)+2**(a+c+e+2)*3**(b+d)

def decode(h):
 def val2(x):
  if x<=0:raise ValueError
  return (x&-x).bit_length()-1
 try:
  z=h-1
  S=val2(z)-2;z-=2**(S+2)
  Q=val2(z)-1;z-=2**(Q+1)
  R=val2(z)-2
  if (S+Q-R)%2:return None
  b=(S+Q-R)//2;a=Q-b;cd=S-b
  if b<1:return None
  z-=2**(R+2)*3**b
  de=val2(z)-1;z-=2**(de+1)
  ce=val2(z)-b-2
  if (cd+ce-de)%2:return None
  c=(cd+ce-de)//2;d=cd-c;e=ce-c
  L=(a,b,c,d,e)
  slope=(8,1,4,1,16)
  lower=max([2]+[(10*x+10*y)//(10*y+1) for x,y in zip(L,slope)])
  upper=min(x//y for x,y in zip(L,slope))
  return L if lower<=upper and form(*L)==h else None
 except ValueError:return None

def checks():
 verify189()
 assert [max(len(comp(m,B)) for B in patt(m)) for m in range(3,7)]==[5,17,21,27]
 for n,B,h in [(8,[(0,3),(0,4),(0,7)],63),(8,[(0,4),(0,7),(1,4),(1,7)],189),(9,[(0,4),(1,7),(2,8)],567),(11,[(0,6),(1,7),(1,10)],1701),(13,[(0,12),(1,4),(2,8)],15309)]:
  assert ham(n,B)==h
  print("Strong witness:",h,n,B)
 B=[(0,2),(1,4),(3,5)]
 C=comp(6,B)
 for L in product(range(1,4),repeat=5):assert exact(C,L)==form(*L)
 num=0
 for m in (2,10,20,34):
  for L in product(*(range(m*t,m*t+m//10+1) for t in (8,1,4,1,16))):
   assert decode(form(*L))==L;num+=1
 assert decode(189) is None
 print("Seven-term compiler checked on 243 tuples; exact box inversions:",num)

if __name__=='__main__':
 import sys
 if len(sys.argv)>1 and sys.argv[1]=='--terms':emit()
 else:checks()
