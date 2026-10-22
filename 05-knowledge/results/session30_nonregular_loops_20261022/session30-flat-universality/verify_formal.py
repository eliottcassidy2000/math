"""Exact symbolic conversion from the proved flat-class expansions."""
from pathlib import Path
import sympy as s
import json
u,z=s.symbols('u z');cut=lambda f:s.series(f,u,0,3).removeO().expand()
g=z*z*u+(z**4-z*z)*u*u;h=-z*z*u+(z*z+z**4)*u*u
figure=cut((2*g*h-s.Rational(3,2)*g*g)/u)
theta_diag=z**6*u*u/3;theta_off=z**4*u*u
vac=cut(figure+theta_diag+theta_off)
A1=cut(vac-z*z*u+4*z**4*u*u)
assert s.expand(A1+(z*z+s.Rational(7,2)*z**4)*u-(12*z**4-s.Rational(8,3)*z**6)*u*u)==0
A2=(7*z**4+26*z**6)*u*u
X0=z*z*(1-u)+z**6*u*u/3
D=lambda f:s.expand(z*s.diff(f,z))
c1=D(X0);c2=D(c1);c3=D(c2);c4=D(c3)
b1=c1;b2=c2+c1*c1;b3=c3+3*c1*c2+c1**3;b4=c4+4*c1*c3+3*c2*c2+6*c1*c1*c2+c1**4
p1=(b2-b1)/2;p2=(3*b4-2*b3-3*b2+2*b1)/24
shift=(D(D(A1))+(2*c1-1)*D(A1))/2
Y1=cut(A1+p1*u);Y2=cut(A2+A1*A1/2+(A1*p1+shift)*u+p2*u*u)
k1=cut(s.exp(X0.subs(z,1)-1)*(1+Y1.subs(z,1)))
k2=cut(s.exp(X0.subs(z,1)-1)*(1+Y1.subs(z,1)+Y2.subs(z,1)))
assert k1==1-s.Rational(5,2)*u+s.Rational(20,3)*u*u
assert k2==1-s.Rational(5,2)*u+s.Rational(65,8)*u*u
out={k:str(v) for k,v in [('figure_eight',figure),('vacuum_A1',vac),('A1',A1),('Y1',Y1),('Y2',Y2),('kappa1_over_e',k1),('kappa2_over_e',k2)]}
Path(__file__).with_name('symbolic-conversion.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
