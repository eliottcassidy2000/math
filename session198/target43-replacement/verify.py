"""Exact scalar-local obstruction and constructive replacement branches."""
from pathlib import Path
from itertools import permutations
import hashlib,json,math
import positive_model as m
from positive_model import P,GEN,BASE,I,mv,mm,tr,dot,old,new,source_bank,build,profile,require
from core import inventories,cycles,glue
import bilinear,c3
HERE=Path(__file__).resolve().parent

def count_perm(T):
 n=len(T)
 H=sum(all(T[x]>>y&1 for x,y in zip(p,p[1:]))for p in permutations(range(n)))
 C=sum(all(T[x]>>y&1 for x,y in zip((0,)+p,p+(0,)))for p in permutations(range(1,n)))
 return H,C

def hp_dp(T):
 n=len(T);A=[[0]*n for _ in range(1<<n)]
 for i in range(n):A[1<<i][i]=1
 for S in range(1,1<<n):
  for v in range(n):
   if A[S][v]:
    rest=T[v]&~S
    while rest:
     b=rest&-rest;rest-=b;A[S|b][b.bit_length()-1]+=A[S][v]
 return sum(A[-1])

def value(f,y,z):return dot(y,mv(P[f],z))
def raw_value(f,s,t):return bilinear.value(f,s,t)if f in('N','M')else c3.value(f,s,t)
def graph(f,yl,zr):
 V,c=build(yl);W,a=build(zr)
 return bilinear.construct(f,V,c,W,a)[0]if f in('N','M')else c3.construct(f,V,c,W,a)[0]
def iso(A,B):
 for p in permutations(range(len(A))):
  if all(bool(A[i]>>j&1)==bool(B[p[i]]>>p[j]&1)for i in range(len(A))for j in range(len(A))):return p
 return None

# Check inherited positive readouts directly against raw four-coordinate laws.
basis=[tuple(int(i==j)for i in range(4))for j in range(4)]
for f in P:
 for y in basis:
  for z in basis:require(value(f,y,z)==raw_value(f,old(y),old(z)),'matrix formula')
for A in P.values():
 require(all(A[i][j]>=P['N'][i][j]for i in range(4)for j in range(4)),'all forms dominate N')
for G in GEN.values():
 A=mm(tr(G),P['N'])
 require(all(3*A[i][j]>=4*P['N'][i][j]for i in range(4)for j in range(4)),'one-letter growth')
require(64*43>27*87,'three-letter window cutoff')

b7=source_bank(7);S=sorted(b7,key=lambda y:(old(y)[2],old(y)));c2=[];hits2=[]
for f in ('N','M'):
 for i,y in enumerate(S):
  for z in S[i:]:
   u,v=old(y)[2],old(z)[2]
   low=3*u*v+1 if f=='N'else 6*u*v+u+v+1
   if low>43:continue
   h=value(f,y,z);row=dict(family=f,left=old(y),right=old(z),H=h)
   c2.append(row)
   if h==43:hits2.append(row)
require(len(b7)==7 and len(c2)==11 and len(hits2)==1,'complete C2 inverse counts')
require(hits2[0]['left']==(1,0,2,1)and hits2[0]['right']==(3,0,6,1),'unique pair')

# Independent complete C3 scan with the inherited conservative U<=floor(43/3).
b14=source_bank(14);oldbank={old(y):''.join(w)for y,w in b14.items()};oldbank[(1,1,1,1)]='singleton'
c3support=set();admissible=0
for f in ('A','B','C','T','J0','J1','J2','J3'):
 for s in oldbank:
  for t in oldbank:
   if not c3.admissible(f,s,t):continue
   admissible+=1;h=c3.value(f,s,t)
   if h<=43:c3support.add(h)
require(len(oldbank)==20 and c3support=={15,27,29,31,33},'complete C3 inverse')

# Every extension of total length<=2 and all ten final families.
words={'':I}
for k in(1,2):
 for w,A in list(words.items()):
  if len(w)==k-1:
   for g,G in GEN.items():words[w+g]=mm(G,A)
y=BASE;z=mv(GEN['O'],mv(GEN['Z'],BASE));reach={};readouts=set();trials=0
for l,L in words.items():
 for r,R in words.items():
  if len(l)+len(r)>2:continue
  for f,A in P.items():
   T=mm(mm(tr(L),A),R);readouts.add(T);h=dot(y,mv(T,z));trials+=1
   if h<=87:reach.setdefault(h,[]).append((f,l,r))
require(trials==340 and len(readouts)==298,'whole finite suffix universe')
require(set(reach)=={43,69,77,79,81,83},'exact descendant spectrum')

# Literal parent classes, both C2 children, and the C3 upper child.
parent=graph('N','','ZO');reverseparent=graph('N','ZO','')
converse_map=iso(bilinear.converse(parent),reverseparent)
require(converse_map is not None,'ordered parent presentations are converses')
literal=[]
for f,a,b,h,C in [('N','','ZO',43,2),('N','ZO','',43,2),('M','Z','Z',85,2),('N','','OZJ',87,2),('J0','','ZJ',87,3)]:
 T=graph(f,a,b);require(count_perm(T)==(h,C),'literal physical special graph')
 require(hp_dp(T)==h,'independent subset check')
 literal.append(dict(family=f,left_word=a,right_word=b,H=h,C=C,rows=T))

births=[]
for mask in range(1<<len(parent)):
 T=parent+[0];n=len(parent)
 for i in range(n):
  if mask>>i&1:T[n]|=1<<i
  else:T[i]|=1<<n
 h=hp_dp(T);births.append(h)
require(85 not in births and 87 not in births,'every arbitrary one-vertex extension')

# Two exact target-dependent replacement identities, on all four-vectors.
qtri=mv(GEN['Z'],BASE)
for q in basis:
 D,E,U,d=old(q);n=value('N',BASE,mv(GEN['O'],q))
 require(n==3*D+12*U-2*d,'parent row')
 require(value('M',qtri,q)-2*n==3*D-4*E-2*U,'lower replacement identity')
 require(value('J0',BASE,mv(GEN['J'],q))-2*n==9*D-3*E-9*U+4*d,'upper replacement identity')

rays=[]
for k in range(1,101):
 word='O'+'Z'*k+'JOZ';q=BASE
 for letter in word:q=mv(GEN[letter],q)
 t=1<<(k-1);D,E,U,d=old(q)
 require((D,E,U,d)==(36*t-1,6*t-1,42*t+1,6*t-1),'lower ray profile')
 require(math.gcd(D,E,U,d)==1,'primitive lower source')
 n=value('N',BASE,mv(GEN['O'],q));child=value('M',qtri,q)
 require(n==600*t+11 and child==2*n-1,'lower ray exact child')
 wordp='OJO'+'Z'*k+'JOZ';qp=BASE
 for letter in wordp:qp=mv(GEN[letter],qp)
 Dp,Ep,Up,dp=old(qp)
 require((Dp,Ep,Up,dp)==(108*t-7,18*t-3,126*t-5,18*t-3),'positive C2 ray profile')
 require(math.gcd(Dp,Ep,Up,dp)==1,'primitive positive C2 source')
 np=value('N',BASE,mv(GEN['O'],qp));cp=value('M',qtri,qp)
 require(np==1800*t-75 and cp==2*np+1,'positive C2 ray exact child')
 word2='O'*(k-1)+'J';u=BASE
 for letter in word2:u=mv(GEN[letter],u)
 x=(1<<k)+1
 require(old(u)==(x,1,x,1),'upper ray profile')
 n2=value('N',BASE,mv(GEN['O'],u));child2=value('J0',BASE,mv(GEN['J'],u))
 require(n2==15*(1<<k)+13 and child2==2*n2+1,'upper ray exact child')
 if k<=3:
  T,r=build(word);require(profile(T,r)==(D,E,U,d),'literal lower ray source')
  for fam,l,r,h,C in [('N','',word+'O',n,2),('M','Z',word,child,2),('N','',word2+'O',n2,2),('J0','',word2+'J',child2,3),('N','',wordp+'O',np,2),('M','Z',wordp,cp,2)]:
   T=graph(fam,l,r);require(inventories(T)['H']==h and cycles(T)==C,'literal ray outputs')
  rays.append(dict(k=k,lower_source=(D,E,U,d),lower_parent=n,lower_child=child,upper_source=old(u),upper_parent=n2,upper_child=child2,positive_C2_source=old(qp),positive_C2_parent=np,positive_C2_child=cp))
receipt=dict(status='PASS',C2_source_cap=7,C2_sources=[old(y)for y in S],C2_candidates=c2,C2_hits=hits2,C3_source_cap=14,C3_sources=len(oldbank),C3_admissible_pairs=admissible,C3_support=sorted(c3support),extension_tests=trials,distinct_extension_readouts=len(readouts),descendant_spectrum_through87=reach,literal_graphs=literal,parent_converse_map=converse_map,arbitrary_births=128,birth_counts=births,lower_upper_ray_algebra_tests=100,literal_ray_parameters=rays,dependency_hashes={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted((HERE/'dependencies').glob('*'))if p.is_file()},scope='Exact target43 C2/C3 inverse and source-extension failure, plus actual restricted target-dependent replacement branches; not general odd coverage.')
out=json.dumps(receipt,sort_keys=True,indent=2)+'\n';(HERE/'VALIDATION.json').write_text(out);print(out,end='')
