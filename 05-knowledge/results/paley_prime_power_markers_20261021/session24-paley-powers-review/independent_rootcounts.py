"""Root counts via gcd with X^(p^d)-X, without building extension fields."""
import json,math

def trim(a):
 while a and a[-1]==0:a.pop()
 return a

def rem(a,b,p):
 a=a[:]; inv=pow(b[-1],p-2,p)
 while len(a)>=len(b):
  shift=len(a)-len(b);c=a[-1]*inv%p
  for i,v in enumerate(b):a[i+shift]=(a[i+shift]-c*v)%p
  trim(a)
 return a

def prod(a,b,m,p):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%p
 return rem(trim(c),m,p)

def power(a,n,m,p):
 r=[1]
 while n:
  if n%2:r=prod(r,a,m,p)
  a=prod(a,a,m,p);n//=2
 return r

def gcd(a,b,p):
 while b:a,b=b,rem(a,b,p)
 return a
out=[]
for p,e in [(3,9),(7,9),(3,15)]:
 u=e//3;r=p**u;f=[p-1]+[0]*(r-1)+[1,p-1]
 D=1 if p==3 else (2 if p%3==1 else 0)
 for d in range(1,e+1):
  if e%d:continue
  a=power([0,1],p**d,f,p)
  a+= [0]*max(0,2-len(a));a[1]=(a[1]-1)%p
  actual=len(gcd(f,trim(a),p))-1
  g=math.gcd(u,d);expected=D if d//g==1 else p**g+1
  if actual!=expected:raise RuntimeError((p,e,d,actual,expected))
  out.append({'p':p,'e':e,'d':d,'g':g,'actual_roots':actual,'expected_M_d':expected})
print(json.dumps(out,indent=2,sort_keys=True))
