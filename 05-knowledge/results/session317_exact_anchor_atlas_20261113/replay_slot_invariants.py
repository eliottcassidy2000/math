from functools import cache
from fractions import Fraction as F
from itertools import combinations
from math import lcm
A=(112,81,35,7,12,26,60); B=(96,81,35,7,13,26,60); C=(104,81,35,6,13,26,60)
D=(100,1,2,7,15,30,62); E=(72,45,1,4,15,29,62)
def flip(R,i,j):
 S=list(R);S[i]^=1<<j;S[j]^=1<<i;return tuple(S)
def vals(R,K=3):
 n=len(R);N=(1<<n)-1
 assert all((R[i]>>j&1)+(R[j]>>i&1)==1 for i,j in combinations(range(n),2))
 for i in range(n):
  u=1<<i
  for _ in range(n):
   for j in range(n):
    if u>>j&1:u|=R[j]
  assert u==N
 @cache
 def h(S,i):
  if S==1<<i:return 1
  return sum(h(S^(1<<i),j) for j in range(n) if S&R[i]&(1<<j))
 H=[1]+[sum(h(S,i) for i in range(n) if S>>i&1) for S in range(1,N+1)]
 @cache
 def q(S,k):
  if not k:return int(not S)
  r=S;v=0
  while True:
   v+=q(S^r,k-1)*H[r]
   if not r:return v
   r=(r-1)&S
 return tuple(q(N,k) for k in range(1,K+1))
assert flip(A,0,4)==B and flip(B,0,3)==C
for r,expected in [(A,(153,1812,11283)),(B,(159,1836,11349)),(C,(159,1812,11253)),(D,(37,740,6135)),(E,(45,732,5823))]:
 assert vals(r)==expected
 print(expected)
def W(R):
 n=len(R);return tuple(x|1<<n for x in R)+(1<<(n+1),(1<<n)-1)
assert vals(W(D),2)==(777,15156) and vals(W(E),2)==(777,14484)
c=[F(0),F(1),F(0)]
for k in range(2,8):c.append(c[k-1]+F(1-k*k,k)*c[k])
assert c[3:7]==[1,F(-8,3),11,F(-832,15)]
for delta in (1,2):
 nodes=[1,9,16];values=[c[k]/(k**delta*(k*k-4)) for k in (1,3,4)]
 print('Newton inverse constants',delta)
 for l in (2,3):
  U=values[:l];V=nodes[:l];co=[]
  for j in range(l):
   co.append(U[0]);U=[(u-U[0])/(v-V[0]) for v,u in zip(V[1:],U[1:])];V=V[1:]
  k=lcm(*(t.denominator for t in co))
  print(l-1,co,lcm(2,k),2*k)
G=(98,92,17,5,9,30,60);J=(82,36,17,7,10,29,62)
assert vals(G)==(99,1332,9153) and vals(J)==(103,1328,8877)
assert vals(W(G),2)==(1431,23436) and vals(W(J),2)==(1431,22860)
print('neutral arrivals',vals(W(G),2),vals(W(J),2))
