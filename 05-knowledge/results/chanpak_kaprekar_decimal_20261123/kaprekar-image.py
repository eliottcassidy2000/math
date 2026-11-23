from itertools import combinations_with_replacement
from math import isqrt
def W(n,k):return 'O'+''.join('I' if k>>j&1 else 'O' for j in range(n-2))+'I'
def hh(w):
 h,x,prev=1,0,0
 for t in w:
  if t=='I':h+=x+prev
  x=2*x+(t=='O');prev=(t=='O')
 return h
def is_prime(n):
 return n>1 and all(n%d for d in range(2,isqrt(n)+1))
vals={1:''}
R={111*d+10*e for d in range(1,10,2) for e in range(d+1)}
for n in range(2,13):
 for k in range(1<<(n-2)):
  w=W(n,k);h=hh(w)
  if h<=max(R)*9:vals.setdefault(h,w)
def inverse(x):
 if x in vals:return [vals[x]]
 for d in range(3,isqrt(x)+1,2):
  if x%d:continue
  a,b=inverse(d),inverse(x//d)
  if a and b:return a+b
 return None

# Hard prime or exceptional-factor cofactors, guaranteed by literal endpoint DP.
FIXED={
 105:(56,405,401,406,0,158,447,273,49),
 343:(280,1,411,2,330,31,47,379,106),
 827:(70,108,256,53,263,21,444,319,43),
 1069:(234,48,331,290,237,4,42,110,243),
}

def count(A):
 n=len(A);N=1<<n
 assert all(not row>>i&1 for i,row in enumerate(A))
 assert all((row>>j&1)+(A[j]>>i&1)==1 for i,row in enumerate(A) for j in range(i))
 d=[[0]*n for _ in range(N)]
 for i in range(n):d[1<<i][i]=1
 for S in range(1,N):
  for i,x in enumerate(d[S]):
   if not x:continue
   z=A[i]&~S
   while z:
    b=z&-z;z-=b;d[S|b][b.bit_length()-1]+=x
 return sum(d[-1])

def rows(w):
 n=len(w);a=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]+[0]
 for j,t in enumerate(w):
  if t=='I':a[j]|=1<<n
  else:a[n]|=1<<j
 return a

def find(target):
 if target in FIXED:return [FIXED[target]]
 if target in vals:return [vals[target]]
 for d in range(3,isqrt(target)+1,2):
  if target%d:continue
  a,b=find(d),find(target//d)
  if a and b:return a+b

def fac(s):
 a=sorted(s)
 return int(''.join(a[::-1]))-int(''.join(a))

if __name__=='__main__':
 from math import prod
 for width in (2,3,4,5):
  m=width//2
  images={sum((10**(width-i-1)-10**i)*d for i,d in enumerate(ds[::-1]))
          for ds in combinations_with_replacement(range(10),m)}
  if width<=4: assert images=={fac(f"{j:0{width}d}") for j in range(10**width)}
  assert len({x%10**m for x in images})==len(images)
  odd={x for x in images if x%2}
  for t in sorted(odd):
   fs=find(t)
   assert fs, t
   assert prod(count(f if not isinstance(f,str) else rows(f)) for f in fs)==t
  print("all odd image heights",width,len(odd))

 for p in (3,5,11,17,41):
  base=p-1
  pool=set(range(1,base,2))
  lengths=[]
  for x in sorted(pool):
   if x not in pool:continue
   y=x;length=0
   while y in pool:
    pool.remove(y);length+=1
    t=y*(base-1)
    fs=find(t)
    assert fs and prod(count(f if not isinstance(f,str) else rows(f)) for f in fs)==t
    y=abs(2*y-p)
   assert y==x
   lengths.append(length)
  print("lucky base",base,"and cycle lengths",lengths)
