#!/usr/bin/env python3
"""Exact outer packet signs, power-of-three pair masses, and clock controls.
Every observed ROOT sign is obtained by a terminating integer computation.
No unobserved source is assigned a sign.
"""
import argparse,json,math
from fractions import Fraction as Q
from pathlib import Path
from collections import defaultdict
from collatz_pair_mass_20261013b import children,parents,step

def need(x,msg):
 if not x:raise RuntimeError(msg)
def T(n):return (3*n+1)//2 if n%2 else n//2
def primitive(n):return n%2==1 and n%8!=5
def w(n):return 2 if n%3==2 else -1

def pair_masses(L=12,B=12,F=12):
 back={(0,0):1}
 for _ in range(B):back=step(back,parents)
 result=[]
 for ell in range(2,L+1):
  forward={(0,3**ell):1};rows=[]
  for k in range(F+1):
   need(sum(forward.values())==2**k,'forward multiplicity')
   count=sum(v*back.get(s,0) for s,v in forward.items())
   rows.append(dict(depth=k+B,unmerged_numerator=2**(k+B)-count,denominator=2**(k+B)))
   if k<F:forward=step(forward,children)
  result.append(dict(ell=ell,translation=3**ell,rows=rows))
 # Independent literal integer-cylinder replay at 12 bits. Equal odd counts
 # exclude isolated source coincidences not shared by the cylinder.
 K=B;checks=0
 for record in result:
  e=record['translation'];absorbed=0
  for r in range(2**K):
   a=2**(K+40)+r;b=a+e;oa=ob=0
   for _ in range(K):
    oa+=a%2;ob+=b%2;a=T(a);b=T(b)
   absorbed+=a==b and oa==ob
  need(2**K-absorbed==record['rows'][0]['unmerged_numerator'],'direct literal pair replay')
  checks+=2**K
 return result,checks

def root_table(limit):
 signs=bytearray(limit+1);signs[1]=1 # 1 -> +1; 2 -> -1
 total=0
 for n in range(2,limit+1):
  a=n;k=0
  while a>=n:
   a=T(a);k+=1
   need(k<100000,'unresolved source: do not invent a ROOT sign')
  signs[n]=signs[a] if k%2==0 else 3-signs[a];total+=k
 return signs,total

def signs_census(max_k=20,max_d=16):
 signs,steps=root_table(4*2**max_k);rows=[];literal_checks=0
 sg=lambda n:1 if signs[n]==1 else -1
 for k in range(4,max_k+1):
  N=2**k
  for d in range(max_d+1):
   q=3**d;b=2**d;lo=(2*N+1)//b+1;hi=(4*N+1)//b
   signed=mass=squared=count=0;by_class={r:[0,0,0] for r in range(8)}
   for h in range(lo,hi+1):
    t=q*h-1
    if not primitive(t):continue
    s=b*h-1;v=w(s)*sg(s)
    signed+=v;mass+=abs(w(s));squared+=w(s)**2;count+=1
    by_class[h%8][0]+=v;by_class[h%8][1]+=abs(w(s));by_class[h%8][2]+=1
    need(w(s)==(2 if h%3==0 else -1),'depth-independent weight')
    if h==lo or h==hi or h%2053==0:
     x=s
     for unused in range(d):
      need(x%2==1,'literal all-odd outer prefix');x=T(x)
     need(x==t,'exact terminal after d shortcut steps')
     a=t;clock=0
     while a!=1:a=T(a);clock+=1;need(clock<100000,'unresolved terminal')
     need(sg(s)==(-1)**(d+clock),'full source and terminal first-hit clocks')
     literal_checks+=1
   need(abs(Q(mass)-Q(N,b))<=12,'isolated absolute packet law')
   need(abs(Q(squared)-Q(3*N,2*b))<=18,'isolated squared packet law')
   rows.append(dict(log2_N=k,depth=d,signed=signed,absolute=mass,squared=squared,count=count,by_h_mod8=by_class))
 # Three forced steps on each primitive class: exactly 2,2,3 odd steps.
 controls=[]
 for r in (1,3,7):
  for d in range(20):
   for z in range(20):
    a=8*3**d*z+r;b=a+8*3**d;oa=0
    for unused in range(3):oa+=a%2;need(a%2==b%2,'shared three-bit prefix');a=T(a);b=T(b)
    need(oa==(3 if r==7 else 2) and b-a==3**(d+oa),'translation after three literal steps')
    controls.append((r,d,z))
 return dict(rows=rows,certified_source_limit=len(signs)-1,integer_steps=steps,literal_outer_paths=literal_checks,three_step_controls=len(controls))

def fibonacci_controls():
 fib=[(0,1)];a,b=0,1
 while True:
  a,b=b,(a+b)%9
  if (a,b)==(0,1):break
  fib.append((a,b));need(len(fib)<100,'finite period bound')
 need(len(fib)==24 and len(set(fib))==24,'Fibonacci mod9 exact 24-cycle')
 # Direct matrix-power tests, not inferred from sequence matching.
 mul=lambda x,y:tuple(sum(x[2*i+k]*y[2*k+j] for k in range(2))%9 for i in range(2) for j in range(2))
 powers=[(1,0,0,1)];F=(1,1,1,0)
 for k in range(24):powers.append(mul(powers[-1],F))
 need(powers[8]==(7,3,3,4) and powers[24]==(1,0,0,1),'matrix lift identity')
 need(powers[12]==(8,0,0,8),'antipodal half-return')
 # h -> h+8 gives three points in each clock fibre. F^8=I+3B.
 need(all(powers[8][i]==((1,0,0,1)[i]+3*(2,1,1,1)[i])%9 for i in range(4)),'unipotent order-three lift')
 seen=set();cycles=[]
 for aa in range(9):
  for bb in range(9):
   if (aa%3,bb%3)==(0,0) or (aa,bb) in seen:continue
   cycle=[];x,y=aa,bb
   while (x,y) not in seen:
    seen.add((x,y));cycle.append((x,y));x,y=(x+y)%9,x
   cycles.append(dict(length=len(cycle),norms=sorted({(x*x-x*y-y*y)%9 for x,y in cycle})))
 need(len(seen)==72 and all(c['length']==24 for c in cycles),'three norm-labelled primitive orbits')
 orbit={(pow(2,j,9)*p[0]%9,pow(2,j,9)*p[2]%9) for j in range(6) for p in powers[:24]}
 need(len(orbit)==72,'Fibonacci plus doubling covers every primitive pair')
 periods={}
 for mod in (24,248,365,613):
  x,y=0,1;k=0
  while True:
   x,y=y,(x+y)%mod;k+=1
   if (x,y)==(0,1):break
  periods[mod]=k
 kathara=[0,1]
 for unused in range(10):kathara.append(sum(kathara))
 need(kathara==[0,1]+[2**i for i in range(10)],'user-supplied cumulative recurrence')
 for cost in range(1,30):
  need(sum(math.comb(cost-1,d-1) for d in range(1,cost+1))==2**(cost-1),'all positive compositions')
 return dict(primitive_cycles=cycles,doubling_joint_orbit=len(orbit),period_controls=periods,kathara=kathara,period=24,pairs=fib,matrix_power8=powers[8],matrix_power12=powers[12],matrix_power24=powers[24])

def remainder_controls():
 z=Q(11,12);F=z/(1-z);G=3*z/(2-z)
 need(F==11 and G/F==Q(3,13),'carry average ratio')
 need(Q(2)+Q(1,2)/(1-G/F)==Q(53,20) and Q(1,3)<Q(1,2)/(1-G/F),'uniform generating coefficient')
 need(36*Q(53,20)/(1-z*z)+Q(8,3)<600,'uniform numerical constant')
 need(Q(12**24,11**22)<2**10,'new uniform exponent below one-half')
 need(Q(8,7)**10*Q(16,7)<Q(12**12,11**11),'Chernoff tail below new main exponent')
 # Exhaust all small words to compare actual normalized offsets with
 # the positive generating sum and the previous worst-word estimate.
 def words(d,K):
  if d==0:yield ()
  else:
   for a in range(1,K-d+2):
    for tail in words(d-1,K-a):yield (a,)+tail
 rec=[]
 for d in range(5):
  for K in range(max(d,1),13):
   total=Q(0);count=0
   for k in range((K-d)//2+1):
    for word in words(d,K-2*k):
     alpha=Q(1);beta=Q(0)
     for a in reversed(word):
      aa=1 if a%2 else 2;c=(3+2**aa-2**a)//3
      alpha*=Q(2**a,3);beta=(2**a*beta-c)/3
     alpha*=4**k;beta=4**k*beta+Q(4**k-1,3)
     total+=abs(beta/alpha)+2;count+=1
   upper=Q(53,20)*z**(-K)*F**d/(1-z*z)
   need(total<=upper,'sum of actual carries bounded by generating series')
   rec.append(dict(d=d,K=K,words=count,carry_sum=str(total),generating_upper=str(upper)))
 return dict(exponent=(12*math.log2(12)-11*math.log2(11))/10,constant=600,checks=rec)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--max-k',type=int,default=20);a=p.parse_args()
 pairs,checks=pair_masses();print('pair laws complete',checks,flush=True)
 r=dict(status='FINITE-EXACT; no universal ROOT claim',pairs=pairs,pair_literal_cylinders=checks,outer=signs_census(a.max_k),fibonacci=fibonacci_controls(),remainder=remainder_controls())
 Path(a.output).write_text(json.dumps(r,sort_keys=True,separators=(',',':'))+'\n')
 print('outer rows',len(r['outer']['rows']),'source limit',r['outer']['certified_source_limit'],'new exponent',r['remainder']['exponent'])
