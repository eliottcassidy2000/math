from math import isqrt
from model import prof

def trans(n):
 return [((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]

def reverse_complement(a,r):
 return int(''.join('1' if b=='0' else '0' for b in f'{a:0{r}b}'[::-1]),2)

def h_word(w,m):
 # H = 1+2I with a uniform O(m)-step recurrence.
 i=0;prefix=0
 for j,b in enumerate(f'{w:0{m}b}'):
  if j and b=='0':i+=(prefix>>1)+(prefix&1)
  prefix=2*prefix+int(b)
 return 1+2*i

def placement(N):
 m=max(1,N.bit_length()-1)
 while 4*(1<<m)<5*N:m+=1
 assert 2*(1<<m)<5*N<=4*(1<<m)
 r,e=divmod(m,2);L=1<<r;s=r//4;g=1<<s
 # ceil(sqrt(N/2^e)) exactly
 a0=isqrt(N>>e)
 if (a0*a0)<<e < N:a0+=1
 v=reverse_complement(a0%g,s) if s else 0
 den=4*(g<<e); M=den*a0+4*v+e*(den//4)
 delta=M*M-((den*den*N)>>e)
 t=(isqrt(delta)+e*(den//4))//(g*den)
 a,b=a0+t*g,a0-t*g
 if not (L//2<=a<L and L//2<=b<L):return None
 A=(a<<(r+e))|reverse_complement(b,r)
 B=reverse_complement(A,m)
 assert B==(b<<(r+e))+(e<<r)+reverse_complement(a,r)
 assert A>>(m-1)==1 and A%2==0
 H=h_word(A,m)
 assert abs((H<<m)-A*B)<=(m//2+1)*(1<<m)
 assert abs(H-N)**8 <= 20**8*N**3
 return m,A,H

if __name__=='__main__':
 from random import Random
 rng=Random(635671)
 direct=0
 for m in range(6,11):
  for N in range(2*(1<<m)//5+1,4*(1<<m)//5+1,max(1,(1<<m)//40)):
   value=placement(N)
   if value is None:continue
   k,A,H=value
   O=int(f'{A:0{k}b}'[::-1],2)
   # Construct the apex directly; this does not use the closer formula.
   B=trans(k)
   T=[v|((1<<k) if not (O>>i)&1 else 0) for i,v in enumerate(B)]+[O]
   assert prof(T)[0]==H
   direct+=1
 samples=0
 for m in range(20,250):
  minN=2*(1<<m)//5+1;maxN=4*(1<<m)//5
  for N in [minN,maxN]+[rng.randrange(minN,maxN+1) for j in range(32)]:
   assert placement(N) is not None
   samples+=1
 print('direct tournament replays:',direct,'large placement controls:',samples)
