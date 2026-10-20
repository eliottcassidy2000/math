#!/usr/bin/env python3
"""Exact Gaussian-integer controls for phase-centered commutator energy."""
import json,math,random
rng=random.Random(904210)
checks=0
def req(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def sub(a,b):return[(x-u,y-v) for (x,y),(u,v) in zip(a,b)]
def norm(a):return sum(x*x+y*y for x,y in a)
def phase(a,j):
 for _ in range(j):a=[(-y,x) for x,y in a]
 return a
mods=[]
for m in range(5,81):
 if math.gcd(m,6)!=1:continue
 i2=pow(2,-1,m);i3=pow(3,-1,m);i4=pow(4,-1,m)
 def U(f):return[f[2*x%m] for x in range(m)]
 def V(f):return[f[(2*x-1)*i3%m] for x in range(m)]
 def Ui(f):return[f[i2*x%m] for x in range(m)]
 def Vi(f):return[f[(3*x+1)*i2%m] for x in range(m)]
 for rep in range(20):
  f=[(rng.randrange(-5,6),rng.randrange(-5,6)) for _ in range(m)]
  sx=sum(x for x,y in f);sy=sum(y for x,y in f);f=[(m*x-sx,m*y-sy) for x,y in f]
  req(sum(x for x,y in f)==0 and sum(y for x,y in f)==0,'mean-zero')
  w=U(V(Ui(Vi(f))))
  req(w==[f[(x+i4)%m] for x in range(m)],'homomorphic permutation convention')
  for j in range(4):
   # Replacing generators by zeta^-1 U,V leaves W unchanged.
   left=norm(sub(w,f));right=8*(norm(sub(U(f),phase(f,j)))+norm(sub(V(f),phase(f,j))))
   req(left<=right,('commutator energy',m,rep,j))
 mods.append(m)
print(json.dumps({'status':'PASS','checks':checks,'moduli':mods,'vectors_per_modulus':20,'phases':['1','i','-1','-i'],'arithmetic':'Gaussian integers, no floating point','scope':'exact controls of algebraic comparison; sine spectral lower bound is proved by cyclic Fourier diagonalization'},indent=2))
