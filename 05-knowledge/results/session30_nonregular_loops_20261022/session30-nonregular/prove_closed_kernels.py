"""Exact scalar rational identities certify all path-length coefficients."""
import sympy as s
from decorated_engine import symbolic_kernel
x,y,z=s.symbols('x y z')
def fs(t):
 return [t/(1-t),-t/(1+t),-t*t*(1-t*t)/((1-t)**2*(1+t*t)),-t*t*(1-t*t)/((1+t)**2*(1+t*t)),t*t/(1+t*t)]
def total(ed,aut,vars):
 sy=symbolic_kernel(ed,aut);f=[fs(v)for v in vars]
 return s.factor(sum(s.Rational(c.numerator,c.denominator)*s.prod(f[j][a]for j,a in enumerate(key))for key,c in sy.items()))
p=total(((0,1),),2,[x]);l=total(((0,0),(0,1)),2,[x,y]);Y=total(((0,1),(0,2),(0,3)),6,[x,y,z])
expected_l=4*x*x*y*y/((1-x**4)*(1+y*y))
F=lambda t:2*t*t/(1+t*t)
V=lambda a,b:-2*a*b*(1+a*b)/((1+a*a)*(1+b*b)) # sum (-1)^a F_(a+b) x^a y^b
# Direct summation of +/-i resolvents yields V: -[R_-i(a)R_i(b)+R_i(a)R_-i(b)]
I=s.I
R=lambda t,w:w*t/(1-w*t)
V=lambda a,b:-(R(a,-I)*R(b,I)+R(a,I)*R(b,-I))
expected_Y=-(F(x)*F(y)*F(z)+V(x,y)*F(z)+V(y,z)*F(x)+V(z,x)*F(y))/6
for name,actual,expected in [('path',p,x*x/(1+x*x)),('lollipop',l,expected_l),('Y',Y,expected_Y)]:
 diff=s.factor(actual-expected)
 print(name,'actual=',actual,'exact_difference=',diff)
 if diff!=0:raise RuntimeError('closed kernel identity failed')
