#!/usr/bin/env python3
import json,math

def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def F(n,t):
 c=n&-n;return n-n//c+2*t*c-1

def classified(t):
 out=[];c=2
 while c<=2*t-2:
  d=2*c
  while d<=8*(t-1)**2:
   D=c+d-1;g=math.gcd(2*c-1,2*d-1)
   if (2*t-1)%(D//g)==0:
    j=(2*t-1)*(2*d-1)//D; jp=4*t-2-j
    n=c*(2*t*(d-1)+j);m=d*(2*t*(c-1)+jp)
    require(n&-n==c and m&-m==d,('cores',t,c,d))
    require(n!=m and F(n,t)==m and F(m,t)==n,('cycle',t,c,d))
    out.append(tuple(sorted((n,m))))
   d*=2
  c*=2
 return sorted(out)
count=0
for t in range(1,301):
 cs=classified(t);count+=len(cs)
 for a in range(1,15):
  c=1<<a
  for b in range(a+1,20):
   d=1<<b;D=c+d-1;g=(1<<math.gcd(a+1,b+1))-1
   require(math.gcd(D,2*d-1)==g,('gcd',a,b))
   integral=((2*t-1)*(2*d-1))%D==0
   require(integral==((2*t-1)%(D//g)==0),('criterion',t,a,b))
   if integral:require(c<=2*t-2 and d<=8*(t-1)**2,('bound',t,a,b))
for t in range(1,21):
 want=set(classified(t));seen=set()
 for n in range(2,100001,2):
  m=F(n,t)
  if m!=n and F(m,t)==n:seen.add(tuple(sorted((n,m))))
 require(seen=={cy for cy in want if min(cy)<=100000},('direct census',t))
for a in range(1,101):
 c=1<<a;t=c//2+1;d=2*c*c
 n=c*(2*c**3+4*c*c+c-1);m=d*(c*c+c-1)
 require(F(n,t)==m and F(m,t)==n and d==8*(t-1)**2,('sharp',a))
print(json.dumps({'status':'PASS','classified_parameters':300,'two_cycles_counted':count,'literal_parameters':20,'literal_even_starts_each':50000,'sharp_family_exponents':100},sort_keys=True,indent=2))
