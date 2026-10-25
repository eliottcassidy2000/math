import sympy as S
import ast,json,runpy
from pathlib import Path
base=Path(__file__).parent
checks=[]
def require(ok,label):
 if not ok:raise RuntimeError(label)
 checks.append(label)
def zero(expr,label):require(S.cancel(expr)==0,label)
p,q,r,s,m=S.symbols('p q r s m',nonzero=True)
a=(p*p-q*q+r*r-s*s)**2
b=(m*m*(p*p-q*q)+r*r-s*s)**2/m**2
c=4*(p*p-s*s)*(q*q-r*r)
e=4*(m*m*p*p-s*s)*(m*m*q*q-r*r)/m**2
H=lambda a,b:a**4-6*a*a*b*b+b**4
zero(a+e-b-c-(m*m-1)*(H(r,s)-m*m*H(p,q))/m**2,'cross harmonic')
zero(a+c-b-e-(m*m-1)*((r*r+s*s)**2-m*m*(p*p+q*q)**2)/m**2,'diagonal norm')
# Elliptic duplication identity using the defining equation to eliminate y^2.
x,y,N=S.symbols('x y N',nonzero=True)
h=(x*x+N*N)/y
xd=(3*x*x-N*N)**2/(4*y*y)-2*x
num=S.together(h*h-4*xd).as_numer_denom()[0]
require(S.rem(S.Poly(num,y),S.Poly(y*y-x**3+N*N*x,y)).is_zero,'hypotenuse squared is 4 doubled x')
# General rational triangle identities and inverse map.
a,b,h,N=S.symbols('a b h N',nonzero=True)
X=N*(a+h)/b;Y=2*N*N*(a+h)/b**2
num=S.together((Y*Y-X**3+N*N*X).subs(N,a*b/2)).as_numer_denom()[0]
require(S.rem(S.Poly(num,h),S.Poly(h*h-a*a-b*b,h)).is_zero,'triangle to curve')
# Substitution above can leave h^3; use Groebner reduction if necessary.
# For centers c=h^2, u^2=c-d,w^2=c+d: legs (w+u)/2,(w-u)/2.
u,w=S.symbols('u w');aa=(w+u)/2;bb=(w-u)/2
zero(aa*aa+bb*bb-(w*w+u*u)/2,'center triangle hypotenuse')
zero(aa*bb/2-(w*w-u*u)/8,'center triangle area d/4')
# Compile independent audit of every existing assert into non-optimizable explicit raises.
class Convert(ast.NodeTransformer):
 def visit_Assert(self,n):
  return ast.copy_location(ast.If(test=ast.UnaryOp(op=ast.Not(),operand=n.test),body=[ast.Raise(exc=ast.Call(func=ast.Name(id='RuntimeError',ctx=ast.Load()),args=[n.msg or ast.Constant('verification failed at source line '+str(n.lineno))],keywords=[]),cause=None)],orelse=[]),n)
for source in [base.parent/'norm-quartic-repair/certificate.py',base.parent/'norm-quartic-repair/polynomial_certificate.py',base.parent/'congruent-curve/check.py']:
 src=Path(source);tree=Convert().visit(ast.parse(src.read_text()));ast.fix_missing_locations(tree)
 dest=base/('session57-'+src.parent.name+'_'+src.name)
 dest.write_text(ast.unparse(tree)+'\n')
 require(not any(isinstance(n,ast.Assert) for n in ast.walk(tree)),'exception safe '+src.name)
print(json.dumps({'independent_symbolic_checks':checks},indent=2))
