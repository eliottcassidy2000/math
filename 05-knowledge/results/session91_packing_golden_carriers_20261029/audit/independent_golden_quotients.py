#!/usr/bin/env python3
"""Independent exact golden quotient, cyclic word fibres and modulus lift audit.
Uses integer pairs and adjugate congruences, not producer normal forms.
"""
from itertools import product
import json
CHECKS=0
def need(p,m):
 global CHECKS
 CHECKS+=1
 if not p:raise RuntimeError(m)
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def mul(x,y):
 a,b=x;c,d=y
 return(a*c+b*d,a*d+b*c+b*d)
def power(x,n):
 y=(1,0)
 for _ in range(n):y=mul(x,y)
 return y
def modulus(n):
 a,b=power((0,1),n);return(a-1,b)
def norm(x):a,b=x;return a*a+a*b-b*b
def key(x,g):
 a,b=g;u,v=x;D=abs(norm(g))
 return(((a+b)*u-b*v)%D,(-b*u+a*v)%D)
def words(n):
 def rec(w):
  if len(w)==n:
   if not(w[0] and w[-1]):yield tuple(w)
  else:
   yield from rec(w+[0])
   if not w or not w[-1]:yield from rec(w+[1])
 yield from rec([])
def read(w):
 x=(0,0)
 for bit in w:x=add(mul(x,(0,1)),(bit,0))
 return x
def lucas(n):
 a,b=2,1
 for _ in range(n):a,b=b,a+b
 return a

def main():
 rows=[];banks={}
 for n in range(1,21):
  g=modulus(n);fib={};bank=list(words(n));banks[n]=bank
  for w in bank:
   x=read(w);q=key(x,g);fib.setdefault(q,[]).append(w)
   need(key(read(w[1:]+w[:1]),g)==key(mul(x,(0,1)),g),'shift equivariance')
  need(len(bank)==lucas(n),'Lucas cyclic count')
  need(len(fib)==abs(norm(g))==lucas(n)-1-(-1)**n,'quotient index')
  fat=[ws for ws in fib.values() if len(ws)>1]
  if n%2:need(not fat,'odd bijection')
  else:
   need(len(fat)==1 and len(fat[0])==3,'one triple fibre')
   need(set(fat[0])=={(0,)*n,tuple(i%2 for i in range(n)),tuple(1-i%2 for i in range(n))},'three boundary words')
   need(key(read(fat[0][0]),g)==(0,0),'zero fibre')
  rows.append(dict(period=n,rooted_words=len(bank),quotient_order=len(fib),largest_fibre=max(map(len,fib.values()))))
 # Exact period lift versus quotient reduction, independent of words in even fibre.
 lifts=0
 for m in range(1,11):
  gm=modulus(m)
  reps={key(read(w),gm):read(w) for w in banks[m]}
  for r in range(1,8):
   gn=modulus(m*r);Q=(0,0)
   for j in range(r):Q=add(Q,power((0,1),j*m))
   image={}
   for q,v in reps.items():
    z=mul(Q,v);image[key(z,gn)]=q
    need(key(z,gm)==key((r*v[0],r*v[1]),gm),'reduction after lift is r')
    need(key(mul(Q,mul((0,1),v)),gn)==key(mul((0,1),z),gn),'phi-equivariant lift')
   need(len(image)==len(reps),'lift injective')
   for w in banks[m]:need(key(read(w*r),gn)==key(mul(Q,read(w)),gn),'word repetition lift')
   lifts+=1
 # Unmarked quotient cannot distinguish the three rational sources at even periods.
 need(key(read((0,0)),modulus(2))==key(read((1,0)),modulus(2))==key(read((0,1)),modulus(2)),'even source hostile')
 print(json.dumps(dict(status='PASS',checks=CHECKS,periods=rows,lift_pairs=lifts,scope='rooted words through20 and exact period lifts m<=10,r<=7'),indent=2))
if __name__=='__main__':main()
