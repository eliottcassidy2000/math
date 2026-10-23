#!/usr/bin/env python3
from fractions import Fraction as F
from pathlib import Path
import json

def T(m,x):
 if m==0:return F(1)
 a,b=F(1),x
 for _ in range(1,m):a,b=b,2*x*b-a
 return b

def U(m,x):
 if m==0:return F(1)
 a,b=F(1),2*x
 for _ in range(1,m):a,b=b,2*x*b-a
 return b
P=lambda m,c:13-12*T(m,F(13-c,12))
D=lambda c:(c-1)*(c-25)
C=lambda t:13-6*(t+1/t)
checks=0
for den in range(1,8):
 for num in range(-30,31):
  c=F(num,den);x=F(13-c,12)
  for m in range(1,13):
   assert D(P(m,c))==D(c)*U(m-1,x)**2
   if c<1:assert U(m-1,x)>0 and P(m,c)<1
   checks+=1
  for m in range(1,5):
   for n in range(1,5):assert P(m,P(n,c))==P(m*n,c)
assert P(2,F(13))==25 and D(F(13))==-144
assert F(-22,5)+F(-3,5)==-5
A=((2,1),(1,1))
mm=lambda a,b:tuple(tuple(sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2))for i in range(2))
M=((1,0),(0,1));fib=[0,1]
for i in range(2,61):fib.append(fib[-1]+fib[-2])
rows=[]
for m in range(1,31):
 M=mm(M,A);tr=M[0][0]+M[1][1];cm=13-6*tr
 assert cm==P(m,F(-5))
 assert D(cm)==180*fib[2*m]**2
 rows.append({'m':m,'trace':tr,'path_forests':fib[2*m],'c':cm,'discriminant':D(cm)})
assert [r['c']for r in rows[:5]]==[-5,-29,-95,-269,-725]
for i in range(2,len(rows)):assert rows[i]['c']==3*rows[i-1]['c']-rows[i-2]['c']-13
result={'rational_grid_polynomial_checks':checks,'golden_orbit_through_m30':rows,'hostile':{'c':13,'D':-144,'P2c':25,'D_after':0},'scope':'Exact finite controls; all-m proof in addendum.'}
Path(__file__).with_name('powering_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'checks':checks,'first_five':rows[:5],'scope':result['scope']},indent=2))
