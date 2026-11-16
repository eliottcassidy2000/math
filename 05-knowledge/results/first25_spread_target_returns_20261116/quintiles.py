from math import gcd,lcm,isqrt
from fractions import Fraction
from itertools import combinations
pr=[p for p in range(2,100) if all(p%d for d in range(2,isqrt(p)+1))]
f={2,3,5,7,11,13};b={17,19,23,31,47};mm={29,41,59,71}
assert [tuple(sum(p in z for p in pr[i:i+5])for z in [f,b,mm])for i in range(0,25,5)]==[(5,0,0),(1,3,1),(0,2,1),(0,0,2),(0,0,0)]
def square(x):return isqrt(x)**2==x
rails=[(i,16-i,20+i) for i in range(1,6)]
assert all(x+y==16 and y+z==36 for x,y,z in rails)
assert [(u,v) for a in rails for u,v in zip(a,a[1:])if square(pr[u-1]+pr[v-1])]==[(1,15)]
a=b0=1; ws=[]
for _ in range(4):
 a,b0=(5*a+3*b0)//2,(7*a+5*b0)//2;m=(a*a-1)//3
 assert (3*m+1,7*m+1)==(a*a,b0*b0)
 ws.append(m)
assert ws==[5,120,2760,63365]
def add(q,old):
 a,m=old;b,n=q;g=gcd(m,n)
 if (b-a)%g:return None
 z=(b-a)//g*pow(m//g,-1,n//g)%(n//g)
 return (a+m*z)%lcm(m,n),lcm(m,n)
D={(0,1):1}; tags={}; ratios=[]
for i,p in enumerate(pr):
 if p>2:
  d=next(v for v in range(1,p) if pow(2,v,p)==1)
  a=next((v for v in range(d) if (99*pow(2,v,p)+23)%p==0),None)
  if a is not None:
   tags[p]=(a,d); E=D.copy()
   for x,c in D.items():
    t=add((a,d),x)
    if t is not None:E[t]=E.get(t,0)-c
   D={k:v for k,v in E.items() if v}
 if i%5==4:ratios.append(sum(Fraction(v,x[1]) for x,v in D.items()))
assert ratios==[Fraction(3,4),Fraction(41,72),Fraction(41,72),Fraction(131893,248820),Fraction(5392213,10201620)]
assert tags[29][0]%4==tags[37][0]%4==3 and tags[29][1]%4==tags[37][1]%4==0
c=[(a,d) for p,(a,d) in tags.items() if p not in (29,37)]
L=lcm(*(d for a,d in c))
def hit(k,ps):return any(k%d==a for a,d in ps)
assert L==61209720 and all(hit(k,c)!=hit(k+L//p,c)for p,k in [(2,5),(3,21)])
a,d=tags[83]
assert next(k for k in range(a,1000,d) if not hit(k,[r for q,r in tags.items()if q!=83]))==917
print(ws,rails,ratios,L)
kept={p:c for p,c in tags.items() if p not in (29,37)}
for p,k in zip((5,13,17,19,53,59,61,67,83),(3,13,17,45,77,93,60,56,917)):
 assert [q for q,(a,d) in kept.items() if k%d==a]==[p]
