#!/usr/bin/env python3
"""Exact bounded CRT / golden-source / bracket-axis bridge checks."""
from fractions import Fraction
from collections import Counter
import json
CHECKS=0
def check(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise RuntimeError(msg)
def mul(x,y):
 a,b=x;c,d=y
 return a*c+b*d,a*d+b*c+b*d
def add(x,y):return x[0]+y[0],x[1]+y[1]
def scale(c,x):return c*x[0],c*x[1]
def power(n):
 x=(1,0)
 for _ in range(n):x=mul(x,(0,1))
 return x
def norm(x):a,b=x;return a*a+a*b-b*b
def read(w):
 x=(0,0)
 for bit in w:x=add(mul(x,(0,1)),(bit,0))
 return x
def words(n):
 def rec(w):
  if len(w)==n:
   if not(w[0] and w[-1]):yield tuple(w)
   return
  yield from rec(w+[0])
  if not w or not w[-1]:yield from rec(w+[1])
 yield from rec([])
def fixed(w):
 a,b,h=1,0,0
 for bit in w:
  if bit:a,b=3*a,3*b+2**h
  else:h+=1
 return Fraction(b,2**h-a)
def standard(x):return 3*x+1 if x.numerator%2 else x/2
def crt(x):
 a,b=x
 return (a%2,b%2),(a+4*b)%11,(a+19*b)%31
def mod341(x):return (x[0]+81*x[1])%341
def order(x,p):
 y=1
 for n in range(1,p+1):
  y=y*x%p
  if y==1:return n
 raise RuntimeError('order missing')
g=add(power(15),(-1,0));A=(188,305)
check(g==(376,610)==scale(2,A),'g15 factor2')
check(norm(A)==-341,'A norm')
check((A[0]+81*A[1])==73*341,'eval kernel')
check((81*81-81-1)==19*341,'phi eval relation')
check(mul(scale(-2,power(9)),mul((-4,1),(7,-2)))==g,'full exact factorization')
check(norm((-4,1))==11 and norm((7,-2))==31,'split norms')
check(crt((0,1))==((0,1),4,19),'phi CRT')
check(order(4,11)==5 and order(19,31)==15,'split clocks')
check(mul(mul((0,1),(0,1)),(0,1))[0]%2==1 and mul(mul((0,1),(0,1)),(0,1))[1]%2==0,'F4 clock3')
# All actual rooted sources at the standard-map clock15.
rows=[];seen=set();fibres=Counter();zero_sources=[]
for w in words(15):
 v=read(w);c=crt(v);r=mod341(v);x=fixed(w);y=x
 check(c not in seen,'all golden source labels distinct');seen.add(c)
 check(r%11==c[1] and r%31==c[2],'CRT341 projection')
 for bit in w:
  check(y.denominator%2==1 and y.numerator%2==bit,'actual source parity')
  y=standard(y)
 check(y==x,'actual source return')
 check(crt(read(w[1:]+w[:1]))==crt(mul(v,(0,1))),'one-step equivariance')
 fibres[r]+=1
 if r==0:zero_sources.append(str(x))
 rows.append((w,v,c,r,x))
check(len(rows)==1364 and len(seen)==4*11*31,'full source/CRT census')
check(set(fibres)==set(range(341)) and set(fibres.values())=={4},'projection uniform four-to-one')
check(sorted(zero_sources)==['0','1','2','4'],'kernel actual sources')
# Same-source period5 and period3 lifts, with their integer costs retained.
q5=add(add((1,0),power(5)),power(10));q3=(0,0)
for j in range(5):q3=add(q3,power(3*j))
check(q5==(38,60) and crt(q5)==((0,0),3,0) and mod341(q5)==124,'period5 lift phase')
check(crt(q3)==((1,0),0,0),'period3 kernel lift')
lift5=[]
for w in words(5):
 v=read(w);base=(v[0]+4*v[1])%11;lift=read(w*3);x=fixed(w)
 check(lift==mul(q5,v),'period5 exact read')
 check(fixed(w*3)==x,'same period5 source')
 check(crt(lift)==((0,0),3*base%11,0),'axis phase3')
 check(mod341(lift)==124*base%341,'axis residue map')
 lift5.append({'source':str(x),'base_label':base,'mod341':mod341(lift),'odd_count5':sum(w),'odd_count15':sum(w*3)})
check({r['mod341'] for r in lift5}=={31*j for j in range(11)},'full conductor11 axis')
for w in words(3):
 check(crt(read(w*5))[1:]==(0,0),'period3 source kernel')
 check(fixed(w*5)==fixed(w),'same period3 source')
# Projected K carrier, not a full-R15 scalar4 group action.
K={pow(4,a,341)*pow(81,b,341)%341 for a in range(5) for b in range(15)}
H11={x*x%11 for x in range(1,11)};H31={x*x%31 for x in range(1,31)}
expected={x for x in range(341) if x%11 in H11 and x%31 in H31}
check(K==expected and len(K)==75,'K is QR11 times QR31')
check(order(4,341)==5 and order(81,341)==15,'K generator orders')
check(pow(19,12,31)==4,'generator CRT exponent')
check(sum(r in K for _,_,_,r,_ in rows)==300,'four golden cells over each K element')
check(all(r['mod341'] not in K for r in lift5),'K unit filtering loses all period5 sources')
check(crt((4,0))[0]==(0,0),'integer4 not full-ring unit')
check(crt((81,0))[0]!=(crt((0,1))[0]),'phi unequal integer81 before quotient')
# Exact inherited conductor11 square weights, uniform on the eleven retained sources.
z=[r['mod341'] for r in lift5];squares=Counter(t*t%341 for t in z)
positive_sum=sum(squares.values())
R=sum(c*r for r,c in squares.items());T=sum(c*((-r)%341) for r,c in squares.items())
check(positive_sum==11 and squares[0]==1,'square source mass')
check(len(squares)==6 and all(c==2 for r,c in squares.items() if r),'axis square fibres2')
check((R,T)==(1364,2046),'conductor11 weighted residue sums')
check(Fraction(2*T,341)==12,'conductor11 exact correction12')
axis31=[11*x for x in range(31)]
T31=sum((-(z*z))%341 for z in axis31)
check(T31==4092 and Fraction(2*T31,341)==24,'conductor31 correction24')
# Nonuniform weights must retain the source-specific phase124x, not replace it by31x.
hostile=[]
for x in range(11):
 a=(-(124*x)**2)%341;b=(-(31*x)**2)%341
 if a!=b:hostile.append([x,a,b])
check(bool(hostile),'source phase changes nonuniform observables')
receipt={'status':'PASS','checks':CHECKS,'g15':g,'A':A,'norm_A':norm(A),'CRT_factors':[4,11,31],'sources':len(rows),'projection_fibres':dict(sorted(fibres.items())),'kernel_sources':sorted(zero_sources),'q5':q5,'q5_mod341':124,'q3':q3,'lift5_sources':sorted(lift5,key=lambda z:z['base_label']),'K_order':len(K),'golden_sources_over_K':300,'period5_sources_retained_by_K':0,'axis11_square_multiplicities':dict(sorted(squares.items())),'axis11_R':R,'axis11_T':T,'axis11_bracket_correction':'12','axis31_T':T31,'axis31_bracket_correction':'24','phase_hostiles':hostile}
print(json.dumps(receipt,sort_keys=True,indent=2))
