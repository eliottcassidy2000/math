from functools import cache
from itertools import permutations,product
D=(2,4,1);C=(2,12,1,5)
def pathcount(a):
 @cache
 def f(S,v):
  if not S:return 1
  w=S&a[v];h=0
  while w:
   z=w&-w;w-=z
   h+=f(S^z,z.bit_length()-1)
  return h
 F=(1<<len(a))-1
 return sum(f(F^(1<<v),v)for v in range(len(a)))
def host(P,Q,flip):
 l=P.bit_length();r=Q.bit_length()
 x=l+3;a=l+4;n=l+5+r;row=[0]*n
 L=list(range(l));DC=list(range(l,l+4));R=list(range(l+5,n))
 def add(i,j):row[i]|=1<<j
 for V in (L,R):
  for i,u in enumerate(V):
   for v in V[i+1:]:add(u,v)
 for i,j in product(L,DC+R):add(i,j)
 for i,j in product(DC,R):add(i,j)
 for i,j in product(range(4),repeat=2):
  if C[i]>>j&1:add(l+i,l+j)
 for j,i in enumerate(L):
  if P>>(l-j-1)&1:add(a,i)
  else:add(i,a)
 for j,i in enumerate(R):
  if Q>>j&1:add(i,a)
  else:add(a,i)
 add(a,l);add(a,l+1);add(l+2,a)
 if flip:add(a,x)
 else:add(x,a)
 assert all((row[i]>>j&1)+(row[j]>>i&1)==1 for i in range(n)for j in range(i))
 return row
def quotient(incoming,outgoing):
 a=list(D)+[outgoing]
 for i in range(3):
  if incoming>>i&1:a[i]|=8
 return a
K=[pathcount(quotient(2,3)),pathcount(quotient(4,5))]
assert K[1]-K[0]==2
def menu(G,first=None,last=None):
 n=len(G);F=(1<<n)-1
 @cache
 def h(S):
  if not S:return 1
  return sum(all(G[i]>>j&1 for i,j in zip(p,p[1:]))for p in permutations(i for i in range(n)if S>>i&1))
 total=0
 for S in range(1,F+1):
  for p in permutations(i for i in range(n)if S>>i&1):
   if first is not None and p[0] not in first:continue
   if last is not None and p[-1] not in last:continue
   if all(G[i]>>j&1 for i,j in zip(p,p[1:])):total+=h(F^S)
 return total,sum(h(S)*h(F^S)for S in range(F+1))
assert menu(D,last={1})[0]==3 and menu(D,first={0,2})[0]==6
assert menu(C,last={2,3})==(15,32) and menu(C,first={0,1})==(12,32)
@cache
def p(v):
 if not v:return 0
 z,b=divmod(v,2)
 return p(z)+(z+z%2 if not b else 0)
for P,Q in product(range(8),repeat=2):
 hs=[pathcount(host(P,Q,b))for b in(0,1)]
 assert hs[1]-hs[0]==2+9*Q-6*P+Q%2-2*(P%2)
for t in range(9):
 h=11+768*t*t+138*t+5*(p(6*t)+p(4*t))
 hs=[pathcount(host(6*t,4*t,b))for b in(0,1)]
 assert hs==[h,h+2]
 print(t,6*t,4*t,hs)
print('64 literal address pairs and 9 exact families checked.')
