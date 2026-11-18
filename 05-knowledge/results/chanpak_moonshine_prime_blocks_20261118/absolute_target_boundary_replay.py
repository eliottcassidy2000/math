from unit_gain_replay import K,forward,H
from math import gcd

def find(p,N):
 h,e,t,d=p
 if h==N:return 0
 if N<=d:return None
 ratio=(2*N-2*d)//(t-d)
 if ratio<1:return None
 k=ratio.bit_length()-1
 for w in range(1,1<<k):
  if h+t*w-d*K(w)-e*(w&1)==N:return w
 return None
for N in (32561,32563):
 for p in ((15,3,19,5),(15,3,21,5)):
  assert find(p,N) is None
 print('all-word target exclusion',N,'PASS')
A=(126,108,56,112,98,64,4);B=list(A);B[2]^=1<<5;B[5]^=1<<2;B=tuple(B)
assert H(forward(A,3,list(map(int,'1000101101'))))==32561
assert H(forward(B,3,list(map(int,'1000101101'))))==32563
print('literal positive controls elsewhere: PASS')
for d in range(3,40,2):
 orbit=set();x=1
 while x not in orbit:orbit.add(x);x=x*2%d
 signs=orbit|{(-x)%d for x in orbit}
 for D in range(-d+1,d):
  if not D:continue
  predicted=D%d in signs
  observed=False
  for n in range(1,4*d+10):
   for A in (-1,1):
    num=A-(1<<n)*D
    if num%d:continue
    S=num//d
    if abs(S)<(1<<n):
     assert S%2
     observed=True
  assert observed==predicted
print('source response-orbit classification: PASS')
