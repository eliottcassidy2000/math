from collections import deque
def mul(A,B):
 a,b,c,d=A; e,f,g,h=B
 return a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h
L=(1,1,0,1);R=(1,0,1,1)
def red(M):
 a,b,c,d=M;w=''
 while True:
  if a>=c and b>=d:
   a,b=a-c,b-d;w+='L'
  elif c>=a and d>=b:
   c,d=c-a,d-b;w+='R'
  else:return w,(a,b,c,d)
M=(11,4,23,9)
print(red(M))
allcores={(c+r,b,c,b+s) for r in range(1,8) for s in range(1,8)
          for b in range(8) for c in range(8) if c*s+b*r+r*s==7}
seen={red(M)[1]}; todo=deque(seen); maxout=0
while todo:
 A=todo.popleft()
 for letter,B in [('L',L),('R',R)]:
  w,N=red(mul(A,B));maxout=max(maxout,len(w))
  if N not in seen:seen.add(N);todo.append(N)
print(len(allcores),len(seen),maxout,sorted(allcores-seen))
assert len(allcores)==len(seen)==25 and maxout==7 and allcores==seen
def wordval(w):
 a,b=1,1
 for s in reversed(w): a,b=(a+b,b) if s=='L' else (a,a+b)
 return a,b
def addr(a,b):
 w=''
 while a!=b:
  if a>b:a-=b;w+='L'
  else:b-=a;w+='R'
 assert a==1
 return w
from math import gcd
from itertools import product
def emit(w):
 out,A=red(M)
 for s in w:
  v,A=red(mul(A,L if s=='L' else R));out+=v
 a,b,c,d=A
 g=gcd(a+b,c+d)
 return out+addr((a+b)//g,(c+d)//g)
for n in range(10):
 for t in product('LR',repeat=n):
  w=''.join(t);x,y=wordval(w)
  a,b=11*x+4*y,23*x+9*y
  g=gcd(a,b)
  assert g==(7 if (x+y)%7==0 else 1)
  assert wordval(emit(w))==(a//g,b//g)
print('1023 primitive address checks passed')
