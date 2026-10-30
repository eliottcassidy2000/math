#!/usr/bin/env python3
"""Exact controls for the two cross-lane maps in SYNTHESIS.md."""
from fractions import Fraction as Q
from itertools import product
from collections import Counter
import json
checks=0
def gate(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)
# Q(zeta_8), represented modulo z^4+1: no floating roots of unity.
def zpow(n):
 n%=8;v=[Q(0)]*4;v[n%4]=Q(1 if n<4 else -1);return tuple(v)
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def scale(c,a):return tuple(c*x for x in a)
def mult(a,b):
 out=[Q(0)]*4
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   k=i+j;out[k%4]+=x*y*(1 if k<4 else -1)
 return tuple(out)
root2=add(zpow(1),scale(-1,zpow(3)))
rows=[]
for j,chi2,chiminus2 in [(1,1,1),(3,-1,1),(5,-1,-1),(7,1,-1)]:
 original=scale(Q(1,2),mult(root2,add(zpow(j),zpow(-j))))
 image=scale(Q(1,2),mult(root2,add(zpow(7*j),zpow(5*j))))
 gate(original==scale(chi2,zpow(0)),('chi2',j))
 gate(image==scale(-chiminus2,zpow(2)),('imaginary transport',j))
 weight=add(scale(Q(88,105),zpow(0)),scale(Q(32,105),original))
 imageweight=add(scale(Q(88,105),zpow(0)),scale(Q(32,105),image))
 gate(weight[0] in (Q(8,7),Q(8,15)) and not any(weight[1:]),('physical weight',j))
 gate(imageweight[2]!=0,('nonreal image',j))
 rows.append({'j':j,'chi2':chi2,'chi_minus2':chiminus2,'weight':str(weight[0]),'image_constant':str(imageweight[0]),'image_i_coefficient':str(imageweight[2])})
# All 16 residue pairs; residue1 denotes a legal multiplier congruent to1.
def address(word,q):
 prefix=0;value=0
 for i,a in enumerate(word,1):
  value-=2**prefix*pow(q,-i,8);prefix+=a
 return value%8
words=[(3,),(1,2),(2,1),(1,1,1)]
pairrows=[]
for q,p in product((1,3,5,7),repeat=2):
 tau={address(w,q):address(w,p) for w in words}
 gate(len(tau)==4,('complete primitive addresses',q,p))
 real_positive=True;fixed=True
 for j in (1,3,5,7):
  original=scale(Q(1,2),mult(root2,add(zpow(j),zpow(-j))))
  transported=scale(Q(1,2),mult(root2,add(zpow(tau[1]*j),zpow(tau[7]*j))))
  weight=add(scale(Q(88,105),zpow(0)),scale(Q(32,105),transported))
  real_positive &= not any(weight[1:]) and weight[0]>0
  fixed &= transported==original
  if p%4==q%4:
   c=q*pow(p,-1,8)%8
   sign=1 if c==1 else -1
   gate(transported==scale(sign,original),('specific weight sign',q,p,j))
  else:
   gate(weight[2]!=0,('bad congruence nonreal',q,p,j))
 gate(real_positive==(p%4==q%4),('positivity gate',q,p))
 gate(fixed==(p==q),('weight invariance gate',q,p))
 pairrows.append({'q_mod8':q,'p_mod8':p,'positive':real_positive,'fixed':fixed})
I=(1,2,4)
def act(A,x):
 r=0
 for j in range(3):
  if x>>j&1:r^=A[j]
 return r
def mm(A,B):return tuple(act(A,x) for x in B)
def order(A):
 B=I
 for k in range(1,169):
  B=mm(B,A)
  if B==I:return k
 raise RuntimeError('matrix order')
def affine_order(v,A):
 w=0;B=I
 for k in range(1,337):
  w^=act(B,v);B=mm(B,A)
  if w==0 and B==I:return k
 raise RuntimeError('affine order')
# A=B=C=D=1 normalized gauge: OR plus the linear character111.
raw_or=[int(v!=0) for v in range(8)]
normalized=[raw_or[v]^(v.bit_count()%2) for v in range(8)]
gate(sum(raw_or)==7 and sum(normalized)==3,'raw/normalized zero-reading distinction')
gate(sum(normalized[v]^(v.bit_count()%2) for v in range(8))==7,'normalized named111 reading')
profiles=Counter()
for A in product(range(1,8),repeat=3):
 if len({act(A,x) for x in range(8)})!=8:continue
 ell=order(A)
 if ell%2==0:continue
 values=[]
 for v in range(8):
  norm=0;x=v
  for _ in range(ell):norm^=x;x=act(A,x)
  f=int(norm!=0)
  gate(f==int(affine_order(v,A)%2==0),('order predicate',A,v))
  values.append(f)
 if ell==1:gate(values==[0]+[1]*7,'identity OR phase')
 elif ell==3:
  gate(sum(values)==4,'linear phase mass')
  gate(all(values[x^y]==values[x]^values[y] for x in range(8) for y in range(8)),'linear phase law')
 elif ell==7:gate(not any(values),'zero phase')
 else:raise RuntimeError('unexpected odd order')
 profiles[(ell,sum(values))]+=1
gate(dict(profiles)=={(1,7):1,(3,4):56,(7,0):48},'complete odd matrix profiles')
means=[]
for h,n3,n7,expected in [(1,0,0,Q(7,8)),(3,2,0,Q(5,8)),(7,0,6,Q(1,8)),(21,14,6,Q(3,8))]:
 mass=(Q(7,8)+n3*Q(1,2))/h
 gate(mass==expected and 1+n3+n7==h,('mixture',h))
 means.append({'complement_order':h,'even_mass':str(mass)})
print(json.dumps({'status':'PASS','checks':checks,'all16_weight_transports':pairrows,'cyclotomic_transport':rows,'odd_matrix_phase_profiles':[{'matrix_order':ell,'negative_reading_numerator':count,'matrices':multiplicity} for (ell,count),multiplicity in sorted(profiles.items())],'group_mixtures':means},indent=2))
