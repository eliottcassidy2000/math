"""Literal small oracle examples for the generalized Collatz wrapper."""
from itertools import product
from fractions import Fraction
from math import isqrt

def primes():
 p=3
 while True:
  if not any(p%d==0 for d in range(2,isqrt(p)+1)):yield p
  p+=2

def wrapper(D):
 # D uses only input/work counter I (and can be replaced by a counter program).
 N,B,I,F,S=range(5)
 code={}
 def put(name,*op):code[name]=op
 def clear(name,c,nxt):put(name,'dec',c,name,nxt)
 def copy(name,c,dst,nxt):
  put(name,'dec',c,name+'s',name+'r')
  put(name+'s','inc',S,name+'d')
  put(name+'d','inc',dst,name)
  put(name+'r','dec',S,name+'c',nxt)
  put(name+'c','inc',c,name+'r')
 clear('start',I,'fuel0')
 clear('fuel0',F,'scratch0')
 clear('scratch0',S,'copyN')
 copy('copyN',N,I,'copyB')
 copy('copyB',B,F,'g0')
 put('restart','inc',B,'start')
 put('end','halt')
 for i,op in enumerate(D):
  name='g'+str(i)
  put(name,'dec',F,name+'p','restart')
  def dest(d):return {'OK':'end','BAD':'restart'}.get(d,'g'+str(d))
  put(name+'p',op[0],op[1],*(dest(x)for x in op[2:]))
 P=primes(); ps=[next(P)for _ in range(5)]; state={k:[next(P),next(P)] for k in code}
 rules={}
 for k,op in code.items():
  for z in (0,1):
   f=state[k][z]
   if op[0]=='halt':rules[f]=();continue
   typ,j,t,*e=op
   if typ=='inc':rules[f]=(Fraction(ps[j]*state[t][1-z],f),)
   else:rules[f]=(Fraction(state[t][1-z],f*ps[j]),Fraction(state[e[0]][1-z],f))
 def step(n):
  if n%2:return (n+1)//2 if n%4==1 else (n-1)//2
  if n%4!=2:return n+1
  hit=[f for f in rules if n%f==0]
  if len(hit)!=1 or n%(hit[0]**2)==0:return n+1
  for q in rules[hit[0]]:
   if n%q.denominator==0:return (n//q.denominator)*q.numerator
  return n+1
 def enc(name,a,z=0):
  n=2*state[name][z]
  for p,k in zip(ps,a):n*=p**k
  return n
 return step,enc,list(code),ps[B]
def finishes(f,n,limit):
 for _ in range(limit):
  if n==1:return True
  n=f(n)
 return n==1

# Total machine whose verdict is always OK, after draining I.
f,enc,names,_=wrapper([('dec',2,0,'OK')])
for k in names:
 for z in (0,1):
  for counters in product(range(2),repeat=5):
   if not finishes(f,enc(k,counters,z),2000):raise ValueError((k,z,counters))
# Total machine whose verdict is BAD iff protected input >=2 after a clean reset.
g,enc,names,bp=wrapper([('dec',2,1,'OK'),('dec',2,'BAD','OK')])
for x in (0,1,2,3):
 n=enc('start',(x,0,0,0,0))
 if finishes(g,n,10000)!=(x<2):raise ValueError(x)
print('literal wrapper: arbitrary toy control starts normalize; BAD epochs do not converge')
