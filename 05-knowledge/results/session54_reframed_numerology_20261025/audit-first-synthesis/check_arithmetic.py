from math import gcd,comb,isqrt
from fractions import Fraction
import json
def need(c,m):
 if not c:raise RuntimeError(m)
counts=[]
for r in range(1,6):
 q=9**r;M=6*9**(r-1);x=1;cs=[0,0,0]
 for k in range(M):
  n=x;b=0;good=True
  for _ in range(r):
   n,d=divmod(n,9);b+=d>=3;good&=d<=5
  for budget in range(3):cs[budget]+=good and b<=budget
  x=x*2%q
 need(x==1,'clock');need(cs==[2*3**(r-1)*sum(comb(r,j) for j in range(min(b,r)+1)) for b in range(3)],'digit counts');counts.append(cs)
sets=[[],[]]
for k in range(2001):
 n=1<<k;b=0;good=True
 while n:n,d=divmod(n,9);b+=d>=3;good&=d<=5
 for budget in range(2):
  if good and b<=budget:sets[budget].append(k)
need(sets==[[0,1,13],[0,1,2,7,13]],'full-digit guard')
ns=[]
for a in range(6):
 for b in range(3):
  for c in range(2):
   for d in range(2):
    for e in range(2):
     n=2**a*3**b*5**c*7**d*13**e
     if n<3:continue
     U={i*i%n for i in range(n) if gcd(i,n)==1};G={min(i,(-i)%n) for i in U}
     if len(G)<=3:ns.append(n)
need(sorted(ns)==[3,4,5,6,7,8,9,10,12,13,14,15,16,18,20,21,24,26,28,30,36,40,42,48,56,60,72,84,120,168],'complete conductor screen')
forms={}
for D in range(-31,-2):
 if D%4 not in (0,1):continue
 f=[]
 for a in range(1,isqrt((-D)//3)+1):
  for b in range(-a,a+1):
   if (b*b-D)%(4*a):continue
   c=(b*b-D)//(4*a)
   if c<a or gcd(gcd(a,b),c)>1 or ((abs(b)==a or a==c) and b<0):continue
   f.append([a,b,c])
 forms[D]=f
need(all(len(f)==1 for D,f in forms.items() if -D<15),'degree2 threshold')
need([D for D,f in forms.items() if -D<19 and len(f)==2]==[-15],'sole degree2 discriminant')
need(all(len(f)<4 for f in forms.values()),'degree4 threshold')
for T,s in [(159616,Fraction(387,100)),(840064,Fraction(435,100)),(13839104,Fraction(565,100))]:
 x=Fraction(157,50)*s;term=Fraction(1);v=term
 for i in range(1,40):term*=x/i;v+=term
 need(v>T+2079,'strict exponential threshold')
print(json.dumps({'status':'PASS','counts':counts,'finite_exponent_sets':sets,'conductors':sorted(ns),'small_discriminant_class_numbers':{str(D):len(f) for D,f in forms.items()},'rational_trace_bounds':True},indent=2))
