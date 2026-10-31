#!/usr/bin/env python3
"""Exact arithmetic controls for affine rays; no Hamiltonian-extension claims."""
import json,math
from fractions import Fraction
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise ValueError(msg)
def blocks(q,c,H,count):
 A=q*q;out=[]
 for t in range(count):
  Q=A**t;need(c*(Q-1)%(A-1)==0,'integral offset');C=c*(Q-1)//(A-1)
  out.append([Q*x+(-C if j%2==0 else C) for j,x in enumerate(H)])
 return out
pairs=0
for q in range(3,202,2):
 A=q*q;h=(A-1)//2
 for r in range(2,2*q+2,2):
  c=A+3-r*r
  eligible=abs(c)<=h
  need(eligible==((A+7)<=2*r*r<=(3*A+5)),'complete parameter interval')
  if not eligible:continue
  pairs+=1;need(c%4==0,'two-prefix parity refinement');B=blocks(q,c,[1,3],8);P=sum(B,[])
  need(all(x>0 for x in P) and all(x<y for x,y in zip(P,P[1:])),'two-prefix positive increasing')
  for t,H in enumerate(B):
   need(sum(H)==(2*q**t)**2,'two-prefix internal square')
   if t+1<len(B):need(H[-1]+B[t+1][0]==(r*q**t)**2,'two-prefix boundary square')
  need(sum(B[1:],[])==[A*x+(-c if j%2==0 else c) for j,x in enumerate(sum(B[:-1],[]))],'two-prefix nested lift')
limit=1000000
actual_plus=set();actual_minus=set()
for q in range(1,limit+1,2):
 for sign,out in [(1,actual_plus),(-1,actual_minus)]:
  square=(q*q+7)//2 if sign==1 else (3*q*q+5)//2;r=math.isqrt(square)
  if r*r==square and r%2==0:out.add((q,r))
def orbit(seeds,matrix):
 out=set();a,b,c,d=matrix
 for q,r in seeds:
  while q<=limit:out.add((q,r));q,r=a*q+b*r,c*q+d*r
 return out
need(orbit([(1,2),(5,4)],(3,4,2,3))==actual_plus,'positive-extreme complete Pell controls')
need(orbit([(1,2),(3,4)],(5,4,6,5))==actual_minus,'negative-extreme complete Pell controls')
specific=[]
for q,s,r,c,expected_descents in [(7,5,8,4,0),(5,6,7,6,20)]:
 need(c==q*q+s*s-6-r*r and (s-r)%2!=0,'four-prefix parameters')
 H=[1,3,6,s*s-6];B=blocks(q,c,H,21);P=sum(B,[])
 need(len(set(P))==len(P) and min(P)>0,'specific four-prefix simplicity')
 for t,block in enumerate(B):
  for j,m in enumerate([2,3,s]):need(block[j]+block[j+1]==(m*q**t)**2,'four-prefix inner square')
  if t+1<len(B):need(block[-1]+B[t+1][0]==(r*q**t)**2,'four-prefix boundary square')
 need(sum(x>y for x,y in zip(P,P[1:]))==expected_descents,'specific monotonicity pattern')
 # Prepending disjointness: no low literal is an image of any positive integer.
 for x in H:
  for sign in [-1,1]:need(not((x-sign*c)%(q*q)==0 and (x-sign*c)//(q*q)>=1),'positive-image residue exclusion')
 if q==5:
  for t,block in enumerate(B):
   Q=25**t;need(block==[(3*Q+1)//4,(13*Q-1)//4,(23*Q+1)//4,(121*Q-1)//4],'Q25 closed form')
   if t+1<len(B):need(B[t+1][0]-block[-1]==(1-23*Q)//2,'Q25 exact descent')
 else:
  for t,block in enumerate(B):
   Q=49**t;need(block==[(11*Q+1)//12,(37*Q-1)//12,(71*Q+1)//12,(229*Q-1)//12],'Q49 closed form')
 specific.append({'q':q,'s':s,'r':r,'c':c,'checked_blocks':21,'descents':expected_descents})
C=blocks(5,0,[1,3,6,75],2);need(C[0][3]==C[1][1]==75,'generic four-prefix counterexample')
# Exact unrestricted-integer c positivity / monotonicity controls for several q.
for q in [3,5,7,9]:
 A=q*q;d=A-1
 for c in range(-4*d,2*d+1):
  beta=Fraction(c,d)
  need((-3*d<=c<=d)==(1-beta>=0 and 3+beta>=0),'exact positivity coefficient classification')
  need((-d<=c<=A-5)==(1+beta>=0 and A-3-(A+1)*beta>=0),'exact increasing coefficient classification')
print(json.dumps({'status':'passed','checks':checks,'two_prefix_parameter_pairs':pairs,'odd_q_tested':[3,201],'two_prefix_blocks_per_pair':8,'Pell_q_bound':limit,'positive_extreme_solutions':sorted(actual_plus),'negative_extreme_solutions':sorted(actual_minus),'specific_four_prefix_controls':specific,'generic_collision_counterexample':C},indent=2))
