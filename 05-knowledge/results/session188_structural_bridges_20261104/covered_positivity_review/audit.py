"""Independent positive-square audit; no producer imports. Standard library."""
from pathlib import Path
import hashlib,itertools,json,subprocess,sys
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parent/'covered_positivity'
def require(ok,why):
 if not ok:raise RuntimeError(why)
def plus(*ps):
 d={}
 for p in ps:
  for e,c in p.items():d[e]=d.get(e,0)+c
 return {e:c for e,c in d.items() if c}
def times(p,q):
 d={}
 for e,c in p.items():
  for f,b in q.items():
   g=tuple(x+y for x,y in zip(e,f));d[g]=d.get(g,0)+c*b
 return {e:c for e,c in d.items() if c}
def scale(p,c):return {e:b*c for e,b in p.items() if b*c}
def formal(m):
 q=[{}]+[{tuple(int(i==j) for i in range(m)):1} for j in range(m)]+[{},{}]
 lhs={};rhs=times(times(q[1],q[1]),q[1])
 for j in range(1,m+1):
  a,b,c=q[j-1:j+2]
  lhs=plus(lhs,scale(times(plus(c,scale(b,2)),plus(times(b,b),scale(times(a,c),-1))),3*j))
  z=plus(b,scale(c,-1));rhs=plus(rhs,times(times(z,z),plus(scale(b,j+2),scale(c,2*j+1))))
 for j in range(1,m+2):
  a,b,c=q[j-1:j+2];z=plus(a,scale(c,-1));rhs=plus(rhs,scale(times(times(z,z),plus(a,scale(b,3),scale(c,2))),j))
 require(lhs==rhs,'independent cubic coefficient identity')
 return [m,len(lhs)]
def check(v):
 q=[0]+list(v)+[0,0];m=len(v)
 p=sum(j*(q[j+1]+2*q[j])*(q[j]**2-q[j-1]*q[j+1]) for j in range(1,m+1))
 r=q[1]**3+sum(j*(q[j-1]-q[j+1])**2*(q[j-1]+3*q[j]+2*q[j+1]) for j in range(1,m+2))+sum((q[j]-q[j+1])**2*((j+2)*q[j]+(2*j+1)*q[j+1]) for j in range(1,m+1))
 require(3*p==r,'numeric identity')
 if all(x>=0 for x in v) and any(v):
  last=max(j+1 for j,x in enumerate(v) if x)
  require(3*p>=v[0]**3+(2*last+3)*v[last-1]**3>0,'boundary strictness')
 return p
formals=[formal(m) for m in range(1,31)]
numeric=0
for m in range(1,7):
 for v in itertools.product(range(4),repeat=m):check(v);numeric+=1
signed=[check(v) for v in ([1,-3,2],[-1,0,4],[0,0,-2],[0,0,3],[1])]
normal=subprocess.check_output([sys.executable,str(SOURCE/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(SOURCE/'verify.py')]);require(normal==optimized,'producer normal versus optimized')
require(normal==(SOURCE/'receipt.json').read_bytes()==(SOURCE/'receipt_optimized.json').read_bytes(),'producer receipts reproduce')
print(json.dumps({'status':'PASS','formal_cutoffs':formals,'numeric_sequences':numeric,'signed_and_boundary_controls':signed,'producer_replayed_normal_and_optimized':True,'source_sha256':{p:hashlib.sha256((SOURCE/p).read_bytes()).hexdigest() for p in ('PROOF.md','verify.py','receipt.json','receipt_optimized.json')},'scope':'Exact positive-square identity and lower bound. Tournament application depends separately on the physical response-gap identity; no unequal-profile or strong-ambient extension.'},sort_keys=True,indent=2))
