#!/usr/bin/env python3
"""Exact actual-content formulas for the two legal endpoint chains."""
from math import gcd,isqrt
import json
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise ArithmeticError(msg)
def carry(w):
 b=0
 for i,t in enumerate(w):b=(3 if t else 1)*b+(t<<i)
 return b
def edited(source,edges):
 w=source[:]
 for i in edges:
  j=(i+1)%len(w);need(w[i]!=w[j],'nonactual edge');w[i],w[j]=w[j],w[i]
 return w
def egcd(a,b):
 oldr,r=a,b;olds,s=1,0;oldt,t=0,1
 while r:
  q=oldr//r;oldr,r=r,oldr-q*r;olds,s=s,olds-q*s;oldt,t=t,oldt-q*t
 return oldr,olds,oldt
def short_relation(p,r):
 N=isqrt(p)+(isqrt(p)**2<p);bins={}
 for j in range(N+1):
  residue=j*r%p;key=N*residue//p
  if key in bins:
   k=bins[key];a=-(j*r//p-k*r//p);b=j-k;m=a*p+b*r
   need(1<=b<=N and -N<=a<=0 and 0<abs(m)<N,'short relation bounds')
   return N,a,b,m
  bins[key]=j
 raise ArithmeticError('pigeonhole failed')
def small_integer(name,N,a,b,m):
 num,den=(10,13) if name=='plus' else(9,5)
 E=(2**(b-a)*num**m-3**b*den**m) if m>=0 else (2**(b-a)*den**(-m)-3**b*num**(-m))
 bound=(52 if name=='plus' else 36)**N
 need(0<abs(E)<bound,'nonzero short integer bound')
 return E
rows=[];maxima={};nonunits=[];shapes=0;crossings=0;signs={'negative':0,'positive':0};squaredepth=[]
for p in range(3,161):
 for r in range(p//2+1,p):
  if gcd(p,r)!=1:continue
  q=p+r;D=2**q-3**p;source=[((i+1)*p)//q-(i*p)//q for i in range(q)]
  Nshort,aa,bb,mm=short_relation(p,r)
  shapes+=1;signs['negative' if D<0 else 'positive']+=1
  labels={(-j*r)%p:(j+(j*r+p-1)//p-1)%q for j in range(p)}
  plus_edges=[labels[p-1-r],labels[p-1]]
  minus_edges=[(labels[r]-1)%q,(labels[0]-1)%q]
  for name,edges,N,num,den in [
   ('plus',plus_edges,gcd(2*13**p-10**p,3*13**r-2*10**r),10,13),
   ('minus',minus_edges,gcd(9**p-2*5**p,2*9**r-3*5**r),9,5)]:
   w=edited(source,edges);B=carry(w);g=gcd(abs(D),B);crossings+=q-1 in edges
   need(g==N,'exact content formula')
   E=small_integer(name,Nshort,aa,bb,mm);need(E%g==0,'short integer not divisible bycontent')
   need(gcd(N,6*num*den)==1,'bad forbidden prime')
   need(D%N==0 and B%N==0,'sufficiency fails')
   if N>1:
    z=num*pow(den,-1,N)%N
    need(pow(z,p,N)==2%N and 2*pow(z,r,N)%N==3%N,'clock equations')
    _,u,v=egcd(p,r);zc=pow(2,u,N)*pow(3*pow(2,-1,N)%N,v,N)%N
    need(z==zc and gcd(z-1,N)==1,'unique clock root')
    if name=='plus':need((13*z-10)%N==0,'linear plus root')
    else:need((5*z-9)%N==0,'linear minus root')
    row=dict(type=name,p=p,q=q,r=r,D=D,content=g,edges=edges)
    nonunits.append(row)
    if any(N%(ell*ell)==0 for ell in range(2,100)) and len(squaredepth)<10:squaredepth.append(row)
   if name not in maxima or N>maxima[name]['content']:
    maxima[name]=dict(type=name,p=p,q=q,r=r,content=N,edges=edges)
# Explicit prime-power-depth controls beyond the main batch.
primepowers=[]
for name,p,r,expected in [('plus',418,237,11849),('minus',203,165,169)]:
 q=p+r;D=2**q-3**p;source=[((i+1)*p)//q-(i*p)//q for i in range(q)]
 labels={(-j*r)%p:(j+(j*r+p-1)//p-1)%q for j in range(p)}
 edges=[labels[p-1-r],labels[p-1]] if name=='plus' else[(labels[r]-1)%q,(labels[0]-1)%q]
 w=edited(source,edges);g=gcd(abs(D),carry(w));N=(gcd(2*13**p-10**p,3*13**r-2*10**r) if name=='plus' else gcd(9**p-2*5**p,2*9**r-3*5**r))
 need(g==N==expected,'primepower depth control')
 n,a,b,m=short_relation(p,r);E=small_integer(name,n,a,b,m);need(E%g==0,'primepower short relation')
 primepowers.append(dict(type=name,p=p,q=q,r=r,content=g,edges=edges,short_relation=[a,b,m],short_integer=E))
# The known norm barriers have tiny actual content; test independently.
need(gcd(2*13**28-10**28,3*13**17-2*10**17)==1,'p28 content')
need(gcd(2*13**18-10**18,3*13**13-2*10**13)==17,'near-criticalp18')
print(json.dumps(dict(status='PASS_EXACT',checks=checks,universe='3<=p<=160, p/2<r<p, gcd(p,r)=1, both endpoint chains, actual cyclic edits',shapes=shapes,paths=2*shapes,shape_clock_signs=signs,crossing_cut_paths=crossings,maximum_contents=maxima,nonunit_count=len(nonunits),primepower_controls=primepowers,nonunit_rows=nonunits),indent=2))
