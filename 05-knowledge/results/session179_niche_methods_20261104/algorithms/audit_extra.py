from itertools import product
from pathlib import Path
import json

def ck(x,s):
 if not x:raise RuntimeError(s)
def matrix(w):
 # Independent homogeneous 2x2 left multiplication, no author import.
 m=((1,0),(0,1))
 for b in w:
  q=((3,1),(0,2)) if b else ((1,0),(0,2))
  m=tuple(tuple(sum(q[i][k]*m[k][j] for k in range(2)) for j in range(2)) for i in range(2))
 return m[0][0],m[1][1],m[0][1]
def itinerary(n,L):
 w=[]
 for j in range(L):
  b=n&1;w.append(b);n=(3*n+1)//2 if b else n//2
 return tuple(w),n
checks=0;illegal=0;zero=0
for L in range(1,10):
 B=1<<L; words=list(product((0,1),repeat=L)); traces={n:itinerary(n,L) for n in range(B)}
 for w in words:
  A,BB,C=matrix(w);ck(BB==B,'matrix denominator')
  legal=[]
  for n in range(B):
   division=(A*n+C)%B==0;native=traces[n][0]==w
   ck(division==native,'illegal residue rejection');checks+=1
   if division:legal.append(n)
   else:illegal+=1
  ck(len(legal)==1,'unique residue');checks+=1
  if A<B:
   # Independent direct walk along the arithmetic progression.
   r=legal[0];n=r
   while n<=0 or itinerary(n,L)[1]>=n:n+=B
   formula=r+B*max(0,(C-(B-A)*r)//((B-A)*B)+1)
   ck(n==formula,'minimum source independent search');checks+=1
  elif A==B:zero+=1
# Explicit conservative PDA counter semantics through weighted +/- unit updates.
for L in range(1,13):
 for w in product((0,1),repeat=L):
  sign=0;stack=0
  for bit in w:
   step=1 if bit else -1
   for _ in range(2 if bit else 3):
    if stack==0:sign=step;stack=1
    elif sign==step:stack+=1
    else:stack-=1
  ck((sign*stack<=0)==(2*sum(w)<=3*(L-sum(w))),'PDA signed counter');checks+=1
  if sign*stack<=0:
   a=sum(w);ck(3**(5*a)*32**L<=2**(5*L)*27**L,'conservative bound');checks+=1
print(json.dumps({'status':'PASS','checks':checks,'incompatible_residue_controls':illegal,'positive_length_equal_slope_words':zero},indent=2))
