def require(x):
 if not x: raise AssertionError("audit failed")
from collections import deque
exec(open('lift.py').read().split('C={}')[0])
m=1729; T=19152
require(power(Q,T,m)==I)
p=31;A=power(Q,T,p); C=mul(mul(B,A,p),inv(B,p),p)
require(A==(27,5,7,22) and C==(12,15,13,6))
start=tuple(x%p for x in z);seen={start};todo=deque([start])
while todo:
 v=todo.popleft()
 for X in (A,C):
  w=act(X,v,p)
  if w not in seen: seen.add(w);todo.append(w)
require(seen=={(x,y) for x in range(p) for y in range(p) if x or y})
print('Q^19152 = identity modulo 1729; A,C mod 31:',A,C)
print('Kernel fibre orbit modulo 31:',len(seen),'all nonzero vectors')
