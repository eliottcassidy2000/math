from functools import cache
from collections import Counter
B=(2,80,3,7,173,143,189,15);M=255
def counts(B):
 @cache
 def e(S,i):
  if S==1<<i:return 1
  return sum(e(S^(1<<i),j) for j in range(8) if S>>j&1 and B[j]>>i&1)
 return e
e=counts(B);s=counts(tuple(M^(1<<i)^B[i] for i in range(8)))
@cache
def H(S):return sum(e(S,i) for i in range(8) if S>>i&1)
@cache
def f(S,r):
 if r==0:return int(S==0)
 if S.bit_count()<r:return 0
 if r==1:return H(S)
 A=S;z=0
 while A:
  z+=H(A)*f(S^A,r-1);A=(A-1)&S
 return z
def marked(r,v,h):
 A=M;z=0
 while A:
  if A>>v&1:z+=h(A,v)*f(M^A,r-1)
  A=(A-1)&M
 return z
from itertools import product
W=[tuple((i+j)%3 for j in range(k)) for i in range(3) for k in range(1,25)]
all1=Counter(tuple(w.count(i) for i in range(3)) for w in W)
st1=Counter(tuple(w.count(i) for i in range(3)) for w in W if w[0]==0)
en1=Counter(tuple(w.count(i) for i in range(3)) for w in W if w[-1]==0)
en2=Counter(tuple(w.count(i)+z.count(i) for i in range(3)) for w in W if w[-1]==0 for z in W)
st2=Counter(tuple(w.count(i)+z.count(i) for i in range(3)) for w in W if w[0]==0 for z in W)
def val(p,tab):return sum(mult*p(a)*f(M,b)*f(M,c) for (a,b,c),mult in tab.items() if a<=8 and b<=8 and c<=8)
HT=val(lambda a:f(M,a),all1)
print('H',HT,HT%27)
for v in range(8):
 S=val(lambda a:marked(a,v,s),st1);E=val(lambda a:marked(a,v,e),en1);D=val(lambda a:marked(a,v,e),en2);R=val(lambda a:marked(a,v,s),st2)
 print(v,(HT-S+E+D)%27,(HT-E+S+R)%27,HT-S+E+D,HT-E+S+R)
