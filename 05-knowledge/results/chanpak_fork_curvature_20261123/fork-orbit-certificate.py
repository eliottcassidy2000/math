"""Exact certificate for modular continuation orbits from positive physical words."""
from collections import deque
from math import isqrt, prod

C = [
 # prime, cyclic word, index, prime factorization of p^d-1
 (3,'W',1,{2:1}),
 (5,'W',1,{2:3,3:1}),
 (7,'WPWWPPPWPWPPWWPP',6,{2:1,3:2,19:1}),
 (11,'WP',1,{2:1,5:2,3221:1}),
 (13,'WPWWPPPWPWPPWWPP',4,{2:3,3:2,7:1,61:1,157:1}),
 (17,'WPPWPPPWPWWP',1,{2:7,3:2,5:1,29:1,41761:1}),
 (19,'WWPWPPWWWPPPWWWPWWWWPWW',54,{2:1,3:4,127:1,523:1,29989:1}),
 (23,'VPVVVVVVVVVVPPV',2,{2:1,11:2,3937230404603:1}),
 (29,'WPPWWWWWWPPWWPWPWWWPWPWP',12,{2:3,3:1,5:1,7:2,574995877:1,88009573:1}),
 (31,'WPWWPWWWWPWPWWWWWWWWWWWPWWPPWWWWPWPWWWWWWW',3,{2:1,3:2,5:2,11:1,331:1,2521:1,17351:1,327412201:1}),
]
def prime(q):
 return q==2 or (q>2 and q%2 and all(q%a for a in range(3,isqrt(q)+1,2)))

def verify(p,word,index,factors):
 d=(p-1)//2; N=p**d-1; o=N//index
 assert N==prod(l**e for l,e in factors.items()) and all(map(prime,factors))
 I=[[int(i==j) for j in range(d)] for i in range(d)]
 D=[[s*a for a in row] for s,row in enumerate(I,1)]
 def mul(A,B):
  return [[sum(a*b for a,b in zip(row,col))%p for col in zip(*B)] for row in A]
 def act(A,v):
  return tuple(sum(a*b for a,b in zip(row,v))%p for row in A)
 def power(A,n):
  R=I
  while n:
   if n&1:R=mul(R,A)
   A=mul(A,A);n//=2
  return R
 def wrap(sigma):
  W=[[s*s*a for a in row] for s,row in enumerate(I,1)]
  W[-1][-1]=(W[-1][-1]+sigma*d)%p
  for s in range(1,d):W[s-1][s]=s;W[s][s-1]=(-s-1)%p
  return W
 P=mul(D,D); W=wrap(1); V=mul(mul(D,wrap(-1)),D)
 def inv(A):
  T=[row+iden for row,iden in zip(A,I)]
  for i in range(d):
   j=next(j for j in range(i,d) if T[j][i])
   T[i],T[j]=T[j],T[i]
   T[i]=[x*pow(T[i][i],-1,p)%p for x in T[i]]
   for j in range(d):
    if i!=j:T[j]=[(x-T[j][i]*y)%p for x,y in zip(T[j],T[i])]
  return [r[d:] for r in T]
 if p==23:
  # Explicit nonzero kernel of even-order W via first d-1 rows.
  ker=[1]
  for s in range(1,d):
   ker.append(((ker[s-2] if s>1 else 0)-s*ker[s-1])%p)
  assert act(W,ker)==(0,)*d
  maps={'P':P,'V':V}
 else:
  maps={'P':P,'W':W}
 for M in maps.values():assert mul(M, inv(M))==I
 S=I
 for a in word:S=mul(maps[a],S)
 assert power(S,o)==I
 assert all(power(S,o//l)!=I for l in factors if o%l==0)
 # A primitive degree-d factor is compulsory in the characteristic polynomial.
 l=next(l for l in factors if o%l==0 and pow(p,d,l)==1 and all(pow(p,j,l)!=1 for j in range(1,d)))
 Ms=[I]
 for _ in range(d-1):Ms.append(mul(S,Ms[-1]))
 K_inv=inv([[M[i][0] for M in Ms] for i in range(d)])
 def key(v):
  c=act(K_inv,v)
  a=[[sum(cj*M[i][j] for cj,M in zip(c,Ms))%p for j in range(d)] for i in range(d)]
  return tuple(row[0] for row in power(a,o))
 e=tuple(row[0] for row in I); seen={key(e)};todo=deque([e])
 while len(seen)<index and todo:
  v=todo.popleft()
  for j in range(100):
   for M in maps.values():
    u=act(M,v);k=key(u)
    if k not in seen:seen.add(k);todo.append(u)
   if len(seen)==index:break
   v=act(S,v)
 assert len(seen)==index
 print(p,'orbit',N,'cyclic index',index,'primitive prime',l,'word',word,flush=True)

if __name__=='__main__':
 for case in C:verify(*case)
