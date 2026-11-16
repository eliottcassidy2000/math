from itertools import combinations
from fractions import Fraction

def require(t):
 if not t: raise AssertionError

def dp(a):
 d=[[0]*len(a) for _ in range(1<<len(a))]
 for i in range(len(a)): d[1<<i][i]=1
 for m,row in enumerate(d):
  for i,c in enumerate(row):
   if not c:continue
   z=a[i]&~m
   while z:
    b=z&-z;z-=b;d[m|b][b.bit_length()-1]+=c
 return d

def new_mask(k):
 return sum(1<<(5-j) for j in range(3) if k>>j&1)

def main():
 n=6;window=range(3,6)
 edges=[(i,j) for i in range(n) for j in range(i+1,n) if j-i<=3]
 fixed=[sum(1<<j for j in range(i+4,n)) for i in range(n)]
 maxima=[Fraction(0)]*8
 raw=[0]*8
 for code in range(1<<len(edges)):
  a=fixed[:]
  for k,(i,j) in enumerate(edges):
   if code>>k&1:a[j]|=1<<i
   else:a[i]|=1<<j
  d=dp(a);h=[sum(row) for row in d]
  for k in range(8):
   A=new_mask(k);bound=Fraction(1);count=1
   for b in window:
    if not A>>b&1:continue
    bound+=1;count+=1;S=a[b]
    while S:
     e=sum(c for i,c in enumerate(d[S]) if not A>>i&1)
     if e:
      bound+=Fraction(e,h[S]);count+=1
     S=(S-1)&a[b]
   maxima[k]=max(maxima[k],bound)
   raw[k]=max(raw[k],count)
 expected=[1,9,17,19,33,33,39,31]
 require(len(edges)==12 and maxima==expected)
 require(raw==[1,9,17,19,33,33,39,32])
 for k,p in enumerate(("012345","015234","041235","041523","301245","301524","304125","304152")):
  p=list(map(int,p));ix={v:i for i,v in enumerate(p)}
  require(all(ix[i]<ix[j] for i,j in combinations(range(n),2) if j-i>3))
  a=[sum(1<<j for j in range(n) if ix[i]<ix[j]) for i in range(n)]
  A=new_mask(k);b=[z|((1<<n) if not A>>i&1 else 0) for i,z in enumerate(a)]+[A]
  require(sum(dp(a)[-1])==1 and sum(dp(b)[-1])==expected[k])
 print("32768 exact fractional bounds and eight independent equality witnesses:",expected)
if __name__=="__main__":main()
