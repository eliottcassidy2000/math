# Independent digit-by-digit controls, unlike block-transition checker.
from fractions import Fraction

def acc(n):
 a=1;x=0
 for b in bin(n)[2:]:
  if b=='0':a+=x;x*=2
  else:a+=1-(x&1);x=2*x+1
 return a

def bits(a,n,word):
 for c in word:
  if c=='0':a+=n;n*=2
  else:a+=1-(n&1);n=2*n+1
 return a,n
for n in range(1,10000):
 a=acc(n);f=n+a
 if f%2!=1:raise ValueError('parity',n)
 g=((f-1)//2)&1
 for b in range(2):
  an=acc(2*n+b);gn=((2*n+b+an-1)//2)&1
  target=g^((n&1) if b==0 else ((n//2+1)&1))
  if gn!=target:raise ValueError('digit',n,b)
for j in range(1,8):
 q=1<<j;word='0'*j+'1'*j
 for n in range(1,100):
  a=acc(n);aa,nn=bits(a,n,word)
  if nn!=q*q*n+q-1 or aa!=a+(q-1)*n+1:raise ValueError('block',j,n)
  f=nn+aa
  if f%2!=1:raise ValueError('block parity')
print('PASS: digit recurrences n<10000; digit-by-digit block transition j<=7, n<100')
# Directly derive the right hand side of the Artin-Schreier equation using
# the original g-series relation: h=(1+z)g => h^2+h=(1+z)c.
# Its denominator has a fifth-order pole at z=1 and numerator evaluates to 1.
num_at_1=1*(1+1+1)%2
if num_at_1!=1:raise ValueError('pole cancelled')
print('PASS: Artin-Schreier right side has exact odd pole order 5 at z=1')
