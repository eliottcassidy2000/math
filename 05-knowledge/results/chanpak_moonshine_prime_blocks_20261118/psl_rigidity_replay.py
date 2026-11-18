from itertools import product
from unit_gain_replay import state,forward,H

def psl(p):
 def norm(A):return min(A,tuple(-x%p for x in A))
 def mul(A,B):
  a,b,c,d=A;e,f,g,h=B
  return norm(((a*e+b*g)%p,(a*f+b*h)%p,(c*e+d*g)%p,(c*f+d*h)%p))
 I=norm((1,0,0,1))
 G={norm(A) for A in product(range(p),repeat=4) if (A[0]*A[3]-A[1]*A[2])%p==1}
 def order(a):
  x=I
  for k in range(1,len(G)+1):
   x=mul(x,a)
   if x==I:return k
 def close(a,b):
  S={I};todo=[I]
  while todo:
   x=todo.pop()
   for y in(a,b):
    z=mul(x,y)
    if z not in S:S.add(z);todo.append(z)
  return S
 target={5:(12,3),7:(24,4),11:(60,5)}[p]
 found=None
 for a in sorted(G):
  if order(a)!=2:continue
  for b in sorted(G):
   if order(b)!=3 or order(mul(a,b))!=target[1]:continue
   S=close(a,b)
   if len(S)==target[0]:found=S;break
  if found is not None:break
 assert found is not None and len(G)==p*len(found)
 cosets=[];unseen=set(G)
 while unseen:
  rep=min(unseen);C={mul(rep,s) for s in found};cosets.append(C);unseen-=C
 assert len(cosets)==p
 a=next(a for a in G if order(a)==2)
 perm=[]
 for C in cosets:
  x=mul(a,next(iter(C)));perm.append(next(i for i,D in enumerate(cosets) if x in D))
 assert any(i!=j for i,j in enumerate(perm)) and all(perm[perm[i]]==i for i in range(p))
 print('PSL action',p,'group',len(G),'stabilizer',len(found),'involution',perm)
for p in(5,7,11):psl(p)
for p in(5,7,11):
 D={1,2} if p==5 else {x*x%p for x in range(1,p)}
 A=tuple(sum(1<<j for j in range(p) if (j-i)%p in D) for i in range(p))
 h,e,t,d=state(A,0)
 multiplier={5:1,7:2,11:3}[p]
 perm=[multiplier*i%p for i in range(p)]
 assert all(((A[i]>>j)&1)==((A[perm[i]]>>perm[j])&1) for i in range(p) for j in range(p))
 Q=forward(A,0,[1]);r=p
 hh,ee,tt,dd=state(Q,r);eta=tt-hh
 assert (hh,ee,tt,dd)==(h+t-e,h-e,h+t,h)
 defect=hh+2*eta-3*dd
 assert H(forward(Q,r,[0,1]))==3*hh+2*eta
 assert H(forward(Q,r,[1,1]))==3*hh+2*eta+defect
 assert defect%{5:1,7:3,11:5}[p]==0
 print('literal affine source',p,h,eta,hh,defect)
print('projective action and rooted-rigidity controls: PASS')
