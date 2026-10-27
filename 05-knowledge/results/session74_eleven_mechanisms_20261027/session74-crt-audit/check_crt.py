#!/usr/bin/env python3
"""Independent standard-library golden15/341 audit, exhaustive on1364 source cells."""
from collections import Counter
from fractions import Fraction
from math import gcd
import json

checks=0
def eq(a,b,name):
 global checks
 checks+=1
 if a!=b:raise RuntimeError(f'{name}: {a!r} != {b!r}')

def rawmul(x,y):
 a,b=x;c,d=y
 return a*c+b*d,a*d+b*c+b*d

def rawpow(x,n):
 z=(1,0)
 for _ in range(n):z=rawmul(z,x)
 return z

# Quotient ideal lattice basis (682,0),(-162,2).
def red(x):
 a,b=x;k,b=divmod(b,2)
 return (a+162*k)%682,b

def mul(x,y):return red(rawmul(x,y))
def power(x,n):
 z=(1,0)
 for _ in range(n):z=mul(z,x)
 return z

def f4mul(x,y):
 a,b=x&1,x>>1;c,d=y&1,y>>1
 return ((a*c+b*d)%2)+2*((a*d+b*c+b*d)%2)

def coords(x):
 a,b=x
 return (a%2)+2*(b%2),(a+4*b)%11,(a+19*b)%31

def project(x):
 a,b=x;return (a+81*b)%341

def tax(x):return (-x*x)%341

phi=(0,1)
eq(rawpow(phi,15),(377,610),'phi15 Fibonacci')
eq(376*986-610*610,-1364,'ideal determinant')
eq((-493*376+305*610,-493*610+305*986),(682,0),'first lattice generator')
eq((118*376-73*610,118*610-73*986),(-162,2),'second lattice generator')
eq(73*493-118*305,-1,'unimodular lattice change')
eq(red((376,610)),(0,0),'generator zero')
eq(red((610,986)),(0,0),'phi generator zero')
eq(power(phi,15),(1,0),'period15')
eq(coords(phi),(2,4,19),'phi coordinates')
eq((81*81-81-1)%341,0,'quotient polynomial')
eq(pow(81,15,341),1,'quotient period')

states=[(a,b) for a in range(682) for b in range(2)]
eq(len(set(coords(x) for x in states)),1364,'coordinate bijection')
fibers=Counter(project(x) for x in states)
eq(set(fibers.values()),{4},'projection fibers')
eq(len(fibers),341,'projection onto')
for x in states:
 c=coords(x);px=mul(phi,x)
 eq(coords(px),(f4mul(2,c[0]),4*c[1]%11,19*c[2]%31),'phi coordinate action')
 eq(project(px),81*project(x)%341,'golden quotient equivariance')
 eq(coords(mul((4,0),x)),(0,4*c[1]%11,4*c[2]%31),'scalar4 coordinate action')
 eq(project(mul((4,0),x)),4*project(x)%341,'scalar4 quotient equivariance')
 eq(project(x)%11,c[1],'CRT11')
 eq(project(x)%31,c[2],'CRT31')
kernel=[x for x in states if project(x)==0]
eq(len(kernel),4,'kernel size')
eq({coords(x) for x in kernel},{(z,0,0) for z in range(4)},'kernel F4')
scalar4image={mul((4,0),x) for x in states}
eq(len(scalar4image),341,'scalar4 not invertible at source')

# Golden period5 -> period15 repetition injection.
eq(rawpow(phi,5),(3,5),'phi5')
eq(rawpow(phi,10),(34,55),'phi10')
Q=(38,60)
eq(rawmul((2,5),Q),(376,610),'factorization')
eq(red(Q),(124,0),'Q quotient representative')
eq(coords(Q),(0,3,0),'Q coordinates')
eq(project(Q),124,'Q projection')
J=[mul(Q,(a,0)) for a in range(11)]
eq(len(set(J)),11,'source injection')
eq({project(x) for x in J},{31*a for a in range(11)},'source conductor11 axis')
eq(mul(Q,Q),mul(Q,(3,0)),'J is not multiplicative')
eq(mul(Q,Q)!=red(Q),True,'non-ring hostile')
e11=mul(Q,(4,0));eq(mul(e11,e11),e11,'local idempotent')
for a in range(11):
 eq(coords(J[a]),(0,3*a%11,0),'source exact coordinates')
 eq(mul(phi,J[a]),J[4*a%11],'source golden equivariance')
 eq(mul((4,0),J[a]),J[4*a%11],'two clocks coincide on source')
 for b in range(11):
  eq(mul(J[a],J[b]),J[3*a*b%11],'twisted source multiplication')

K={1};pending=[1]
while pending:
 x=pending.pop()
 for g in (4,81):
  y=x*g%341
  if y not in K:K.add(y);pending.append(y)
eq(len(K),75,'K order')
eq(K,{x*x%341 for x in range(341) if gcd(x,341)==1},'K unit square image')
action_kernel=[k for k in K if all(k*project(x)%341==project(x) for x in J)]
eq(len(action_kernel),15,'source action kernel')
eq(len({k%11 for k in K}),5,'source action image')
eq(len({k*project(J[1])%341 for k in K}),5,'source nonzero orbit')

source_tax=sum(tax(project(x)) for x in J)
ambient_tax=sum(tax(x) for x in range(341))
full_ring_tax=sum(tax(project(x)) for x in states)
eq(source_tax,2046,'eleven source tax')
eq(ambient_tax,59334,'ambient tax')
eq(full_ring_tax,4*ambient_tax,'full ring source multiplicity')
eq(Fraction(source_tax,11),186,'source normalized tax')
eq(Fraction(ambient_tax,341),174,'ambient normalized tax')
# A genuinely nonuniform source measure; no replacement by count support.
weights=list(range(1,12))
weighted_direct=sum(weights[a]*tax(project(J[a])) for a in range(11))
push={project(J[a]):weights[a] for a in range(11)}
eq(weighted_direct,sum(w*tax(x) for x,w in push.items()),'weighted pushforward')

# Period3 kernel: all1364 legal length15 rooted words, with an independent affine reader.
Q35=(0,0)
for j in range(5):
 z=rawpow(phi,3*j);Q35=(Q35[0]+z[0],Q35[1]+z[1])
eq(coords(Q35),(1,0,0),'period3 repetition idempotent')
eq(mul(Q35,Q35),red(Q35),'period3 kernel idempotent')
words=[]
def grow(prefix):
 if len(prefix)==15:
  if not(prefix[0] and prefix[-1]):words.append(tuple(prefix))
  return
 grow(prefix+[0])
 if not prefix or prefix[-1]==0:grow(prefix+[1])
grow([])
eq(len(words),1364,'rooted legal word count15')
word_images=[];kernel_words=[]
for w in words:
 image=(0,0);A=Fraction(1);B=Fraction(0)
 for bit in w:
  y=mul(image,phi);image=red((y[0]+bit,y[1]))
  if bit:A=3*A;B=3*B+1
  else:A/=2;B/=2
 word_images.append(image)
 if project(image)==0:
  source=B/(1-A)
  eq(all(w[i]==w[i%3] for i in range(15)),True,'kernel word period3')
  kernel_words.append({'word':''.join(map(str,w)),'source':str(source),'ring_label':image})
eq(len(set(word_images)),1364,'golden word bijection15 independent finite check')
eq(len(kernel_words),4,'kernel word count')
eq({r['source'] for r in kernel_words},{'0','1','2','4'},'rational source kernel')

print(json.dumps({
 'status':'PASS','checks':checks,'source_cells':1364,
 'quotient':{'order':341,'fiber_size':4,'kernel_representatives':kernel,
             'phi_image':81,'scalar4_image_size_at_source':341},
 'lift':{'Q_coefficients':[38,60],'Q_reduced':[124,0],'Q_CRT':[0,3,0],
         'projected_sources':[project(x) for x in J],
         'source_cardinality':11,'K_order':75,'K_source_action_order':5,
         'K_source_action_kernel_order':15,'local_idempotent':e11},
 'period3_kernel':{'Q_coefficients':Q35,'Q_reduced':red(Q35),'words':kernel_words},
 'weights':{'source_tax':source_tax,'ambient_tax':ambient_tax,'full_ring_tax':full_ring_tax,
            'source_normalized_tax':186,'ambient_normalized_tax':174,
            'nonuniform_source_weights':weights,'nonuniform_source_tax':weighted_direct},
 'hostiles':['scalar4 is nonunit on full R15',
             'J(a)J(b)=J(3ab), not J(ab)',
             'original eleven-source mass is not ambient uniform mass',
             'the K action on the selected source forgets its order15 factor']
},indent=2,sort_keys=True))
