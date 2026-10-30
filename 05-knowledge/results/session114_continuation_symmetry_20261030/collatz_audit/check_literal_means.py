from fractions import Fraction as F
from math import gcd,ceil,log
import ast,json
ns={}
source='check_periodic_independent.py'
t=ast.parse(open(source).read())
# Import only definitions/constants without running the separate enumeration.
t.body=[node for node in t.body if isinstance(node,(ast.Import,ast.ImportFrom,ast.FunctionDef)) or (isinstance(node,ast.Assign) and any(isinstance(x,ast.Name) and x.id=='INF' for x in node.targets))]
exec(compile(t,source,'exec'),ns)
v5,depth,params=ns['v5'],ns['depth'],ns['parameters']
INF=float('inf')
def affine(word):
 A=0;B=0;d=len(word)
 for i,a in enumerate(word):B+=3**(d-i-1)*2**A;A+=a
 return F(3**d,2**A),F(B,2**A)
def mean(w,u):
 m,b=affine(w);t,r=params(m,b);answer=F(0);Ds=[]
 for j in range(1,len(w)+1):
  n,c=affine(w[:j]+u);D=depth(m,b,n,c);Ds.append(D)
  if D<=r:z=F(int(D),t)
  else:z=F(4*r+1,4*t)-(F(0) if D==INF else F(1,4*t*5**int(D-r)))
  answer+=z
 return answer,Ds
rows=[];count=0
for w in [(1,),(2,),(5,),(1,2),(2,1),(1,)*5,(1,)*10+(2,),(5,5),(3,1,3),(1,2,1,2)]:
 for u in [(),(1,),(2,),(3,1)]:
  rho,Ds=mean(w,u);P=Q=1;B=0
  for a in u:P,Q,B=3*P,2**a*Q,3*B+Q
  total=0;maxheight=0
  for k in range(1,201):
   for a in w:
    P,Q,B=3*P,2**a*Q,3*B+Q
    c=int(v5(gcd(P-Q,B)));total+=c;maxheight=max(maxheight,c)
   if k in (1,10,50,200):
    # Direct finite sum must match its exact finite-E expected value within one per layer.
    # Take E at least every observed height; tail error is computed analytically.
    E=max(maxheight,8);m,b=affine(w);t,r=params(m,b)
    truncated=F(0)
    for D in Ds:
     for e in range(1,E+1):
      if e<=D:truncated+=F(1,t*5**max(0,e-r))
    if abs(total-k*truncated)>len(w)*E:raise RuntimeError((w,u,k,total,truncated,E))
    count+=1
  rows.append({'w':w,'u':u,'mean':str(rho),'sum200':total,'maximum_observed_height':maxheight,'error200':str(total-200*rho)})
print(json.dumps({'status':'PASS','literal_prefix_checks':count,'rows':rows},indent=2))
