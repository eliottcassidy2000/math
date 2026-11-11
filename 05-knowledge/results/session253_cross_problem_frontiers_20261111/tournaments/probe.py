from itertools import combinations,product
from math import factorial,comb
import json

def tour(n,mask):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a[i][j]=bool(mask>>k&1);a[j][i]=not a[i][j]
 return a

def profile(a):
 n=len(a);N=1<<n
 end=[[0]*n for _ in range(N)];h=[0]*N
 for v in range(n):end[1<<v][v]=1
 for s in range(1,N):
  for v in range(n):
   if s>>v&1:
    r=s^(1<<v)
    if r:end[s][v]=sum(end[r][u]*a[u][v] for u in range(n) if r>>u&1)
  h[s]=sum(end[s])
 pc=[[0]*(n+1) for _ in range(N)];pc[0][0]=1
 for s in range(1,N):
  root=s&-s;t=s
  while t:
   if t&root:
    for k in range(1,n+1):pc[s][k]+=h[t]*pc[s^t][k-1]
   t=(t-1)&s
 return [pc[-1][k]*factorial(k) for k in range(n+1)]

def sub(q,blocks):
 own=[i for i,b in enumerate(blocks) for _ in b];local=[j for b in blocks for j in range(len(b))]
 return [[blocks[i][local[u]][local[v]] if i==j else q[i][j] for v,j in enumerate(own)] for u,i in enumerate(own)]

def trans(n):return [[int(i<j) for j in range(n)] for i in range(n)]
C3=[[0,1,0],[0,0,1],[1,0,0]]
Q4=[[0,1,0,0],[0,0,1,0],[1,0,0,1],[1,1,0,0]]
small=[tour(n,m) for n in range(1,4) for m in range(1<<(n*(n-1)//2))]
profiles={str(a):profile(a) for a in small}
checks=0
for q in [C3,trans(3)]:
 for bs in product(small,repeat=3):
  if sum(map(len,bs))>7:continue
  H=profile(sub(q,bs))[1]
  eps=[trans(1 if len(b)%2 else 2) for b in bs]
  rhs=profile(sub(q,eps))[1]
  for b in bs:rhs*=profiles[str(b)][1]
  if (H-rhs)%3:raise RuntimeError(('mod3',bs,H,rhs))
  checks+=1
hostile=None
for n in range(1,6):
 for mask in range(1<<(n*(n-1)//2)):
  a=tour(n,mask);F=profile(a)
  if F[1]%5==0 and F[2]%5!=0:
   H=profile(sub(C3,[a,trans(1),trans(1)]))[1]
   hostile={'n':n,'mask':mask,'rows':a,'ordered_profile':F,'context_H':H};break
 if hostile:break

def jet_join(F,n,m,K):
 # Q_S(j), using F(k)*binomial(j,k)
 vals=[sum(F[k]*comb(j,k) for k in range(1,min(j,len(F)-1)+1)) for j in range(K)]
 return [0]+[sum((-1)**(k-j)*comb(k,j)*pow(vals[j],n,m) for j in range(k+1))%m for k in range(1,K)]
period_checks=0
for p,a,period in [(2,1,1),(2,2,2),(2,3,2),(3,1,2),(3,2,6),(5,1,4)]:
 m=p**a;K=1
 while factorial(K)%m:K+=1
 for b in [trans(1),C3,Q4]:
  F=profile(b)
  for n in range(a,a+9):
   if jet_join(F,n,m,K)!=jet_join(F,n+period,m,K):raise RuntimeError('period')
   period_checks+=1
print(json.dumps({'mod3_direct_substitution_checks':checks,'prime5_absorption_hostile':hostile,'join_jet_period_checks':period_checks},indent=2))
extra={}
for name,q in [('C3',C3),('Q4',Q4)]:
 extra[name+'_parity_germs']={''.join(map(str,e)):profile(sub(q,[trans(i) for i in e]))[1]%3 for e in product((1,2),repeat=len(q))}
extra['Q4_ordered_profile']=profile(Q4)
cover=[]
for a in range(1,8):
 m=3**a;P=2*3**(a-1)
 x={ (1+pow(2,r+1,m))%m for r in range(P)}
 y={ (5+pow(2,r+5,m))%m for r in range(P)}
 if len(x|y)!=m:raise RuntimeError('cover')
 cover.append({'a':a,'period':P,'first':len(x),'second':len(y),'union':len(x|y)})
extra['two_seed_modular_cover']=cover
padding_checks=0
for S in [trans(1),C3,Q4]:
 F=profile(S);Q2=2*F[1]+(F[2] if len(F)>2 else 0)
 for r in range(4):
  block=sub(trans(2),[S,trans(r)]) if r else S
  H=profile(sub(C3,[block,trans(1),trans(1)]))[1]
  if H!=F[1]+2**r*Q2:raise RuntimeError('padding')
  padding_checks+=1
extra['padding_checks']=padding_checks
print(json.dumps(extra,indent=2))
