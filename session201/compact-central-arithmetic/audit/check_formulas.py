"""Compare every coefficient of the twelve written formulas with verified matrices."""
import ast,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
if not(ROOT/'MENUS.json').is_file():ROOT=ROOT/'session201-compact-central'
vars=('A','p','u','B','q','v');zero=(0,)*6

def add(x,y,sign=1):
 z=dict(x)
 for k,v in y.items():z[k]=z.get(k,0)+sign*v
 return{k:v for k,v in z.items()if v}
def multiply(x,y):
 z={}
 for a,b in x.items():
  for c,d in y.items():
   k=tuple(u+v for u,v in zip(a,c));z[k]=z.get(k,0)+b*d
 return{k:v for k,v in z.items()if v}
def parse(x):
 if isinstance(x,ast.Constant)and isinstance(x.value,int):return{zero:x.value}
 if isinstance(x,ast.Name)and x.id in vars:
  k=tuple(int(v==x.id)for v in vars);return{k:1}
 if isinstance(x,ast.UnaryOp)and isinstance(x.op,ast.USub):return{k:-v for k,v in parse(x.operand).items()}
 if isinstance(x,ast.BinOp):
  l,r=parse(x.left),parse(x.right)
  if isinstance(x.op,ast.Add):return add(l,r)
  if isinstance(x.op,ast.Sub):return add(l,r,-1)
  if isinstance(x.op,ast.Mult):return multiply(l,r)
 raise ValueError('unsupported formula syntax')
menus={r['word']:r['matrix']for r in json.loads((ROOT/'MENUS.json').read_text())['5']}
lines=[line for line in(ROOT/'FORMULAS.md').read_text().splitlines()if line.startswith('- ')]
seen=[]
for line in lines:
 name,formula=line[2:].split(': ',1);w=''if name=='empty'else name;matrix=menus[w];expected={}
 for i,left in enumerate(('A','p','u','1')):
  for j,right in enumerate(('B','q','v','1')):
   term=parse(ast.parse(left+'*'+right,mode='eval').body)
   expected=add(expected,{k:v*matrix[i][j]for k,v in term.items()})
 actual=parse(ast.parse(formula,mode='eval').body)
 if actual!=expected:raise RuntimeError('coefficient mismatch '+name)
 seen.append(name)
if len(seen)!=12 or len(set(seen))!=12:raise RuntimeError('formula count')
print(json.dumps(dict(status='PASS',formal_polynomial_coefficients_checked=12,words=seen,
 formula_sha256=hashlib.sha256((ROOT/'FORMULAS.md').read_bytes()).hexdigest(),
 verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,sort_keys=True))
