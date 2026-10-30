"""Independent permutation-order and literal-centralizer checks. Standard library only."""
from itertools import permutations,product
from math import gcd
from fractions import Fraction
from collections import Counter
import json
checks=0
def gate(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)
def mul(a,b):return tuple(a[x] for x in b)
def inv(a):
 b=[0]*len(a)
 for i,x in enumerate(a):b[x]=i
 return tuple(b)
def order(a):
 seen=set();o=1
 for i in range(len(a)):
  if i in seen:continue
  j=i;n=0
  while j not in seen:seen.add(j);n+=1;j=a[j]
  o=o*n//gcd(o,n)
 return o
def close(gens):
 I=tuple(range(len(gens[0])));S={I};todo=[I]
 while todo:
  a=todo.pop()
  for b in gens:
   c=mul(a,b)
   if c not in S:S.add(c);todo.append(c)
 return S
def sign(a):return (-1)**sum(a[i]>a[j] for i in range(len(a)) for j in range(i+1,len(a)))
def parity_group(n):return set(a for a in permutations(range(n)) if sign(a)==1)
def direct(A,B):
 n=len(next(iter(A)))
 return {a+tuple(n+x for x in b) for a in A for b in B}
def matrix(A,x):
 out=0
 for i in range(3):
  if (x>>i)&1:out^=A[i]
 return out
M=[A for A in product(range(1,8),repeat=3) if len({matrix(A,x) for x in range(8)})==8]
GL={tuple(matrix(A,x)-1 for x in range(1,8)) for A in M}
gate(len(GL)==168,'GL order')
D8=close([(1,2,3,0),(0,3,2,1)])
S4=set(permutations(range(4)));C3=close([(1,2,0)])
groups=[('D8',D8),('S4',S4),('GL3(2)',GL),('A6',parity_group(6)),('A7',parity_group(7)),('D8xC3',direct(D8,C3)),('S4xC3',direct(S4,C3))]
rows=[]
for name,G in groups:
 orders={a:order(a) for a in G};seen=set();weights=[]
 for x in sorted(G):
  if x in seen or orders[x] not in (2,4):continue
  cls={mul(mul(g,x),inv(g)) for g in G};seen|=cls
  C={g for g in G if mul(g,x)==mul(x,g)}
  rho=Fraction(sum(orders[g]%2==1 for g in C),len(C))
  gate(len(cls)*len(C)==len(G),(name,'centralizer index'))
  weights.append((orders[x],rho,len(cls),len(C)))
 c=sum(o==2 for o,*_ in weights)
 mass=Fraction(sum(o%2==0 for o in orders.values()),len(G))
 gate(mass==sum(w for o,w,*_ in weights),(name,'decoration sum'))
 gate(mass==Fraction(2*c+1,8),(name,'D8 formula'))
 gate(sorted(w for o,w,*_ in weights if o==2)==[Fraction(1,8)]+[Fraction(1,4)]*(c-1),(name,'involution weights'))
 gate([w for o,w,*_ in weights if o==4]==[Fraction(1,4)],(name,'order4 weight'))
 rows.append({'group':name,'size':len(G),'orders':dict(sorted(Counter(orders.values()).items())),'involution_classes':c,'even_mass':str(mass),'decorations':[{'order':o,'odd_centralizer_mass':str(w),'class_size':cl,'centralizer_size':ce} for o,w,cl,ce in weights]})
# Enumerate odd subgroups through actual cyclic groups and Sylow-7 normalizers.
C3s={frozenset(close([g])) for g in GL if order(g)==3}
C7s={frozenset(close([g])) for g in GL if order(g)==7}
N7s=set()
for H in C7s:
 a=next(g for g in H if order(g)==7)
 N=frozenset(g for g in GL if mul(mul(g,a),inv(g)) in H)
 gate(len(N)==21,'Singer normalizer');N7s.add(N)
gate((len(C3s),len(C7s),len(N7s))==(28,8,8),'odd subgroup census')
affine=[];affine_decorations=[];audited_sizes=set()
I=tuple(range(7))
for H in [frozenset({I}),*C3s,*C7s,*N7s]:
 # Convert linear action on nonzero vectors into permutations of eight affine points.
 A=[]
 for h in H:
  linear=(0,)+tuple(x+1 for x in h)
  for v in range(8):A.append(tuple(v^linear[x] for x in range(8)))
 mass=Fraction(sum(order(a)%2==0 for a in A),len(A))
 predicted=1-sum(Fraction(1,1+sum(h[x]==x for x in range(7))) for h in H)/len(H)
 gate(mass==predicted,('norm kernel',len(H)))
 gate(mass=={1:Fraction(7,8),3:Fraction(5,8),7:Fraction(1,8),21:Fraction(3,8)}[len(H)],('cube mass',len(H)))
 if len(H) not in audited_sizes:
  audited_sizes.add(len(H));A=set(A);seen=set();dec=[];O={g:order(g) for g in A}
  for x in sorted(A):
   if O[x]!=2 or x in seen:continue
   seen|={mul(mul(g,x),inv(g)) for g in A}
   C={g for g in A if mul(g,x)==mul(x,g)}
   dec.append(Fraction(sum(O[g]%2 for g in C),len(C)))
  expected={1:[Fraction(1,8)]*7,3:[Fraction(1,8)]*2+[Fraction(3,8)],7:[Fraction(1,8)],21:[Fraction(3,8)]}[len(H)]
  gate(sorted(dec)==sorted(expected),('affine decorations',len(H)))
  affine_decorations.append({'H_order':len(H),'weights':list(map(str,sorted(dec)))})
 affine.append((len(H),str(mass)))
print(json.dumps({'status':'PASS','gates':checks,'D8_controls':rows,'odd_subgroup_census':{'1':1,'3':len(C3s),'7':len(C7s),'21':len(N7s)},'affine_groups_checked':len(affine),'affine_decorations':affine_decorations},indent=2))
