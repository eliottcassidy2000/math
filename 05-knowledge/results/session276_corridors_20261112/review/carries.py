"""Independent integer arithmetic controls; physical consumer is already audited."""
from math import gcd
import json

def g(x):
 value=0
 while x:
  half=x>>1
  if not(x&1):value+=half+(half&1)
  x=half
 return value

def F(A,B):return 15381+3357*(g(A)+g(B))+50004*A*B+25915*A+20732*B+1865*(A&1)+1492*(B&1)
checks=0
for k in range(1,8):
 for d in range(64):
  x=(1<<k)*(2*d+1)
  assert g(x)-g(x-1)==((1<<(k+1))-3)*d+(1<<k)-(d&1);checks+=1
for k in range(1,6):
 for l in range(1,6):
  ak=(3357-50004)*(1<<(k+1))-3*3357;al=(3357-50004)*(1<<(l+1))-3*3357
  for d in range(16):
   for e in range(16):
    A=(1<<k)*(2*d+1)-1;B=(1<<l)*(2*e+1)
    C=(50004-3357)*((1<<l)-(1<<k))+4810+ak*(d&1)-al*(e&1)+3357*((e&1)-(d&1))
    pred=2*ak*(d//2)-2*al*(e//2)+C
    assert F(A+1,B-1)-F(A,B)==pred;checks+=1
examples=[]
for j in [0,1,2,7,100,10**6,10**30]:
 A=29527+1362656*j;B=29140+1344752*j
 assert (A+1)%16==8 and B%8==4
 assert ((A+1)//8-1)//2%2==1 and (B//4-1)//2%2==0
 value=F(A,B);assert F(A+1,B-1)==value+4
 examples.append({'j':j,'N':value,'interval':[value,value+2,value+4,value+6]});checks+=1
print(json.dumps({'status':'PASS','integer_identity_checks':checks,'family_examples':examples,'scope':'exact integer compiler identities, not direct enumeration of exponentially many paths'},indent=2))
