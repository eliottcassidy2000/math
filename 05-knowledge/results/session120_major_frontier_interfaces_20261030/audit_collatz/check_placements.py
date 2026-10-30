"""Independent exhaustive retained-pair and all-placement audit."""
from fractions import Fraction as F
import json

def need(p,msg):
    if not p: raise RuntimeError(msg)
q=5**10
M=25*q
T=2*5**9
m=pow(3*pow(2,-1,M)%M,T,M)
u=(m,(m-1)%M)
U=((u[0]-1)//q,u[1]//q)
need(U==(4,4),'two-digit allones tangent')
score=sum(T//(2*5**(e-1)) for e in range(1,11))
need(F(5*T,8)-score==F(1,4),'allones defect')
need(F(12207031,19531250)*T-score==F(1,5),'allones price')
checks=0
primitive_checks=0
for alpha in range(25):
 for beta in range(25):
  c=(1+q*alpha,q*beta)
  for a in range(5):
   predicted=((1+q*(alpha+a*U[0]))%M,q*(beta+a*U[1])%M)
   for l in range(a+1):
    blocks=[u]*l+[c]+[u]*(a-l)
    state=(1,0)
    returns=[]
    for index,(mb,tb) in enumerate(blocks,1):
     state=(mb*state[0]%M,(mb*state[1]+tb)%M)
     if state[0]%(5*q)==1 and state[1]%(5*q)==0:
      returns.append(index)
    need(state==predicted,'all-placement full two-digit invariant')
    if (alpha+4*a)%5==0 and (beta+4*a)%5==0:
     need(returns==[a+1],'only terminal block boundary returns at level11')
     tangent=(((state[0]-1)//(5*q))%5,(state[1]//(5*q))%5)
     need(tangent==(((alpha+4*a)//5)%5,((beta+4*a)//5)%5),'quotient tangent')
     primitive_checks+=1
    checks+=1
# Algebraic j-transfer, without assuming existence of any core score.
for a in range(5):
 for core_clock in range(1,20):
  for core_score in range(20):
   A=core_clock+a*T
   S11=core_score+a*score+1
   need(5*A-8*S11==5*core_clock-8*core_score+2*a-8,'defect transfer')
print(json.dumps({'status':'PASS','allones_two_digit_pair':U,
 'all_pair_placement_checks':checks,'admissible_pair_placement_checks':primitive_checks,
 'native_internal_prefix_semantics':'proved in AUDIT.md','actual_core_existence_claimed':False},indent=2))
