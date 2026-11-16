from functools import lru_cache
from itertools import product, combinations

def dp(r,root=None):
 n=len(r);N=1<<n
 d=[[0]*n for _ in range(N)]
 for j in range(n) if root is None else (root,):d[1<<j][j]=1
 for i,row in enumerate(d):
  for j,v in enumerate(row):
   if not v:continue
   z=r[j]&~i
   while z:
    b=z&-z;z^=b
    d[i|b][b.bit_length()-1]+=v
 return [1]+[sum(t) for t in d[1:]]

def prof(t,x):
 d=dp(t);z=dp(t,x);N=len(d)-1
 ev,od=(sum(z[s]*d[N^s] for s in range(1,N+1) if s.bit_count()%2==k) for k in (0,1))
 return d[-1],d[N^(1<<x)],ev,od

def tail(t,x,w):
 t=list(t)
 for b in map(int,w):
  n=len(t);m=0 if b else 1<<x
  t=[a|((not(m>>j&1))<<n) for j,a in enumerate(t)]+[m]
 return tuple(t)

def rc(w):return w[::-1].translate(str.maketrans('01','10'))

@lru_cache(None)
def f(z):
 if not z:return 1
 u,p=divmod(z,2)
 return f(u)+(0 if p else u+(u&1))

def response(profile,w):
 assert w and w[-1]=='0'
 h,g,r,s=profile
 z=int(rc(w),2)
 x=int(w,2)
 b=s-g
 hr=h+g*(f(z)-1)+(r+b)*z+(r-b)*(w[0]=='0')
 rr=((1<<len(w))*(r+b)+g*x)//2
 return hr,g,rr,g+rr

@lru_cache(None)
def inverse(N,a,b,c=0):
 if N<b:return ()
 res=[0] if N==b else []
 if N>=2*a+b and (N-b)%2==0:
  res += [2*u for u in inverse(N,2*a+b,b,b) if u]
 if N>=a+b+c and (N-a-c-b)%2==0:
  res += [2*u+1 for u in inverse(N-a-c,2*a,b,0)]
 return tuple(sorted(res))

def select(C,C0,g,r):
 return tuple(z for z in inverse(C-C0+g,2*r,g) if z)

def from_slope(S,L,g,r):
 I=S-(1<<L)*r
 if I%g or not 0<=I//g<1<<(L-1):return None
 return format(2*(I//g),f'0{L}b')

def diff(a,b):
 return tuple((i,j) for i,j in combinations(range(len(a)),2) if (a[i]^b[i])>>j&1)

def main():
 p=(2,12,9,17,39,79,31)
 q=(2,20,17,7,41,79,31)
 assert prof(p,1)==(41,9,20,29) and prof(q,1)==(43,9,20,29)
 for L in range(1,6):
  rec=set()
  for bits in product('01',repeat=L-1):
   w=''.join(bits)+'0'
   a,b=(prof(tail(t,1,w),1) for t in (p,q))
   assert a==response((41,9,20,29),w)
   assert b==response((43,9,20,29),w)
   assert b[0]-a[0]==2 and a[1:]==b[1:]
   assert from_slope(a[2],L,9,20)==w
   rec.add(a[2])
  assert rec=={20*(1<<L)+9*j for j in range(1<<(L-1))}
 # A non-stable seed is covered as well: remove two newest near-sinks.
 for L in range(1,5):
  t=(18,8,3,21,38,15);x=0
  assert prof(t,x)==(31,9,17,22)
  for bits in product('01',repeat=L-1):
   w=''.join(bits)+'0'
   assert prof(tail(t,x,w),x)==response(prof(t,x),w)
 examples=(
 (378,4,1,(2,12,17,37,43,71,31),(2,12,17,69,43,15,55),295),
 (10008,347,1,(2,12,33,21,7,27,63,127),(2,8,51,85,3,27,55,127),7353),
 (30090,410,3,(2,8,3,21,71,31,175,63),(2,8,67,21,71,31,43,127),18750),
 )
 for C,z,x,a,b,R in examples:
  h,g,r,s=prof(a,x)
  assert s-r==g and prof(b,x)==(h+2,g,r,s)
  assert z in select(C,h+1,g,r)
  w=rc(format(z,'b'))
  A,B=(tail(t,x,w) for t in (a,b))
  assert prof(A,x)==(C-1,g,R,g+R) and prof(B,x)==(C+1,g,R,g+R)
  assert diff(A,B)==diff(a,b)
  print('centre, reverse address, tail word, half-slope, edits:',C,z,w,R,diff(a,b))
if __name__=='__main__':main()
