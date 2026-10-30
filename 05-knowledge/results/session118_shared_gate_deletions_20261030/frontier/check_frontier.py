from fractions import Fraction as F
import json

def require(ok,why):
 if not ok:raise RuntimeError(why)
def rows(i,p,t):
 return i>=0 and p>=0 and i<=2+p and t<=2+2*i-p and t<=4+2*p and t<=3+i+p
count=feasible=0
for ci in range(101):
 c=F(ci,20)
 for bi in range(261):
  b=F(bi,20);i=b-c
  lo=max(F(0),b-c-2,(b-4)/2,c-3);hi=b-2*c+2
  predicted=0<=c<=4 and b>=max(c,2*c-2,4*c-8)
  actual=i>=0 and lo<=hi
  require(actual==predicted,('projection',c,b))
  if actual:
   for p in [lo,hi,(lo+hi)/2]:require(rows(i,p,b),('witness',c,b,p))
   feasible+=1
  count+=1
boundary=0
for ci in range(81):
 c=F(ci,20)
 i,p,t=(F(0),F(0),c) if c<=2 else ((c-2,F(0),2*c-2) if c<=3 else (3*c-8,2*c-6,4*c-8))
 require(rows(i,p,t) and t-i==c,'boundary')
 require(t==max(c,2*c-2,4*c-8),'sharp boundary')
 boundary+=1
# Coefficient-four weights must satisfy p=b-6 and i=b-4; b>=8.
for bi in range(161,261):
 b=F(bi,20);lo=max(F(0),b-6,(b-4)/2,F(1));hi=b-6
 require(lo==hi,'coefficient-four fiber')
print(json.dumps({'status':'PASS','projection_points':count,'feasible_pairs':feasible,'boundary_points':boundary,'coefficient_four_minimum_penalty':8,'optimal_weights':[1,4,2,8],'scope':'strengthened row constraints, not external theorem certification'},indent=2,sort_keys=True))
