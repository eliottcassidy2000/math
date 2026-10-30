from fractions import Fraction as F
import json
checks=0
def gate(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)
rows=[]
for a in range(101):
 c=F(a,20)
 for z in range(261):
  b=F(z,20);i=b-c
  lo=max(F(0),b-c-2,(b-3)/2,c-3);hi=b-2*c+2
  actual=i>=0 and lo<=hi
  predicted=c<=4 and b>=max(c,2*c-2,4*c-7)
  gate(actual==predicted,('projection',c,b))
  if actual:
   for p in (lo,hi):
    gate(i>=0 and p>=0 and i<=2+p and b<=2+2*i-p and b<=3+2*p and b<=3+i+p,('original rows',c,b,p))
for a in range(81):
 c=F(a,20)
 if c<=2:i,p,t=F(0),F(0),c
 elif c<=F(5,2):i,p,t=c-2,F(0),2*c-2
 else:i,p,t=3*c-7,2*c-5,4*c-7
 gate(t-i==c and t==max(c,2*c-2,4*c-7),('boundary',c))
 gate(i<=2+p and t<=2+2*i-p and t<=3+2*p and t<=3+i+p and i>=0 and p>=0,('boundary rows',c))
gate(max(F(0),9-4-2,F(9-3,2),4-3)==9-2*4+2==3,'unique four')
print(json.dumps({'status':'PASS','exact_gates':checks,'projection_grid':'c=0..5 and b=0..13, both in steps of1/20','boundary':'81 rational points, c0..4','scope':'rowwise potential cone, not certification of circuit theorem'},indent=2))
