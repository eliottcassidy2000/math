#!/usr/bin/env python3
"""Independent small-universe audit: recurrence words and permutation cycles, no producer imports."""
import json
from math import gcd,lcm
from pathlib import Path
checks=0
def check(ok,label):
 global checks
 checks+=1
 if not ok: raise RuntimeError(label)
def per(w):
 # KMP border algorithm, not exhaustive shifts.
 p=[0]*len(w)
 for i in range(1,len(w)):
  j=p[i-1]
  while j and w[i]!=w[j]:j=p[j-1]
  if w[i]==w[j]:j+=1
  p[i]=j
 d=len(w)-p[-1]
 return d if len(w)%d==0 else len(w)
def orbit(n):
 a,b=0,1; out=[]
 while not out or (a,b)!=(0,1):
  out.append((a,b));a,b=b,(a+b)%n
  check(len(out)<=10*n*n,'orbit bound')
 return out
def obs(n,stride):
 o=orbit(n); states=[o[j%len(o)] for j in range(0,lcm(stride,len(o)),stride)]
 return len(states),per([(2*b-a)%n for a,b in states]),states
summary={}
for n in range(2,101):
 for stride,c in [(1,4),(2,2)]:
  m,t,states=obs(n,stride);u=n;e=0
  while u%5==0:u//=5;e+=1
  mu=obs(u,stride)[0] if u>1 else 1
  expected=mu if not e else lcm(mu,c*5**(e-1))
  check(t==expected,('trace',n,stride,t,expected))
  scalar=next((k for k,(a,b) in enumerate(states[1:],1) if a==0),m)
  check(m//scalar in ((1,2,4) if stride==1 else (1,2)), 'scalar loss')
  check(m//t in (1,5),'trace loss')
  # Coefficient of B for Q is F_j; for A is F_(2j), because A=Q+I.
  if e:
   tags={}
   for a,b in states:
    aa=(b-a)%n if stride==1 else (b-2*a)%n
    tau=1 if stride==1 else 3; det=-1 if stride==1 else 1
    T=(2*aa+tau*a)%n; U=(tau*aa+(tau*tau-2*det)*a)%n
    tag=(T,U,(a%5**e)//5**(e-1))
    check(tag not in tags,'cyclic sidecar collision');tags[tag]=(a,b)
# Entire algebra, including noncyclic states, and explicit radical square-zero.
for n in range(2,41):
 for tau,det in [(1,-1),(3,1)]:
  e=0;u=n
  while u%5==0:e+=1;u//=5
  tags=set();rad=[]
  for a in range(n):
   for b in range(n):
    T=(2*a+tau*b)%n;U=(tau*a+(tau*tau-2*det)*b)%n
    tag=(T,U, (b%5**e)//5**(e-1) if e else 0)
    check(tag not in tags,('algebra sidecar',n,tau));tags.add(tag)
    if T==U==0:rad.append((a,b))
  check(len(rad)==gcd(n,5),'radical size')
  for a,b in rad:
   for c,d in rad:
    check((a*c-det*b*d)%n==0 and (a*d+b*c+tau*b*d)%n==0,'radical square zero')
# Actual permutations on all board points: cycle-lcm order, distinct from affine iteration.
def perm_order(p):
 seen=set();order=1
 for start in range(len(p)):
  if start in seen:continue
  j=start;k=0
  while j not in seen:seen.add(j);k+=1;j=p[j]
  check(j==start,'bijection cycle')
  order=lcm(order,k)
 return order
affine_cases=0;hostile={}
for n in range(2,12):
 for r in range(1,11):
  f,g=0,1
  for _ in range(r):f,g=g,f+g
  A=(g,f,f,g-f)
  def v(x,y):return ((A[0]*x+A[1]*y)%n,(A[2]*x+A[3]*y)%n)
  linear=[n*v(x,y)[0]+v(x,y)[1] for x in range(n) for y in range(n)]
  t=perm_order(linear)
  image={((x-v(x,y)[0])%n,(y-v(x,y)[1])%n) for x in range(n) for y in range(n)}
  d=gcd(g-1,f,g-f-1);D=abs((g-1)*(g-f-1)-f*f)
  check(n*n//len(image)==gcd(n,d)*gcd(n,D//d),'SNF cokernel')
  for bx in range(n):
   for by in range(n):
    p=[n*((v(x,y)[0]+bx)%n)+(v(x,y)[1]+by)%n for x in range(n) for y in range(n)]
    actual=perm_order(p)
    # Read return translation from permutation power at the origin.
    j=0
    for _ in range(t):j=p[j]
    cx,cy=divmod(j,n)
    check(actual==t*n//gcd(n,cx,cy),'affine order')
    if (bx,by) in image:check(actual==t,'removable order')
    affine_cases+=1
    if (n,r,bx,by) in [(5,4,1,0),(11,5,1,3)]:
     hostile[f'N{n}_r{r}_b{bx},{by}']={'linear':t,'affine':actual,'removable':(bx,by) in image}
m,t,states=obs(15,2)
check(states[10]==(0,11),'mixed scalar 11I')
check(not any(a==0 and b==14 for a,b in states),'no negative identity')
# Ramified metallic tests use norm valuation (not explicit ideal enumeration).
metallic=[]
for prime,metal in [(5,1),(5,4),(13,3),(17,8)]:
 for k in (1,2):
  n=prime**k;seen=set();orbits=0
  gc,gb=metal*metal+1,metal*(metal*metal+2)
  def advance(x):
   c,b=x
   return ((gc*c+gb*b)%n,(gc*b+gb*c+metal*gb*b)%n)
  for c in range(n):
   for b in range(n):
    x=(c,b)
    if x in seen:continue
    word=[];y=x
    while y not in seen:
     seen.add(y);word.append((2*y[0]+metal*y[1])%n);y=advance(y)
    check(y==x,'metallic orbit closure')
    if x==(0,0):check(len(word)==per(word)==1,'zero clocks');continue
    norm=c*c+metal*c*b-b*b
    v=0
    while norm%prime==0 and norm:
     v+=1;norm//=prime
    check(norm!=0 and v<2*k,'norm valuation defined')
    r=2*k-v
    check(len(word)==prime**(r//2),'metallic full clock')
    check(per(word)==prime**((r-1)//2),'metallic trace clock')
    orbits+=1
  metallic.append({'p':prime,'a':metal,'k':k,'states':n*n,'nonzero_orbits':orbits})
# Independent SNF closed factors against recurrence entries, through stride300.
fibs=[0,1]
for _ in range(301):fibs.append(fibs[-1]+fibs[-2])
for r in range(1,301):
 d=gcd(fibs[r+1]-1,fibs[r],fibs[r-1]-1)
 D=abs((fibs[r+1]-1)*(fibs[r-1]-1)-fibs[r]**2)
 if r%2:expected_d=2 if r%3==0 else 1;expected_e=D//expected_d
 elif r%4==2:expected_d=expected_e=fibs[r//2-1]+fibs[r//2+1]
 else:expected_d=fibs[r//2];expected_e=5*expected_d
 check((d,D//d)==(expected_d,expected_e),'closed Smith factors')
summary={'status':'PASS','checks':checks,'metallic_exhaustive':metallic,'clock_moduli':'2..100, Q and A','full_algebra_sidecar_moduli':'2..40, Q and A','affine_moduli':'2..11, strides 1..10, every carry','smith_strides':300,'affine_cases':affine_cases,'hostiles':hostile,'algorithms':['Fibonacci pair recurrence','KMP trace word period','full permutation cycle-lcm affine order','exhaustive algebra trace-sidecar injectivity and radical products'],'optimization_safe':True}
print(json.dumps(summary,indent=2))
Path(__file__).with_name('results.json').write_text(json.dumps(summary,indent=2)+'\n')
