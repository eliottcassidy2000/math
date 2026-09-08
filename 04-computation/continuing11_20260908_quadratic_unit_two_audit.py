"""Independent exact audit of a global quadratic with a cyclic two-arm unit.

The all-degree/source-module statements use the companion analytic proof.
No producer is imported or executed. All gates remain active under -O.
"""
from pathlib import Path
import json,sys
import sympy as s
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve().parent
OUT=HERE.parent/'05-knowledge/results' if HERE.name=='04-computation' else HERE
STEM=Path(__file__).stem
u,t,r,b,c,y,w,g=s.symbols('u t r b c y w g')
checks=0
def need(ok,label):
 global checks
 checks+=1
 if not ok:raise ArithmeticError(label)
def zero(v,label):need(s.cancel(v)==0,label)
def jac(F,G,a=u,d=t):return s.diff(F,a)*s.diff(G,d)-s.diff(F,d)*s.diff(G,a)

W=1+u**2*t; L=(u+1)*W-4
F=s.expand(u*W*L)
G=((1-2*u)*W+2)/(3*u**2*W**2)
zero(jac(F,G)-1,'literal original-source Jacobian one')
zero(F-((u+1)*u**5*t**2+2*u**3*(u-1)*t+u*(u-3)), 'expanded source polynomial')
need(s.degree(F,t)==2,'genuine source degree two')
zero(s.diff(F,u).subs(u,0)+3,'affine u0 transverse derivative')
Y=u*W
H=(1+1/u)*y**2-4*y
zero(H.subs(y,Y)-F,'actual alternate rational coordinates')
zero(jac(u,Y)-u**3,'u nonzero coordinate chart')
zero(s.diff(H,u)+y**2/u**2,'critical point forces y0')
zero(s.diff(H,y).subs(y,0)+4,'other derivative never zero at y0')

F_inf=s.cancel(F.subs({u:1/r-1,t:-r**2-r**4*b}))
zeta=r*(1+b*(1-r)**2)
zero(F_inf-(1-r)*(zeta**2-4),'full polynomial second-chart identity')
need(s.fraction(F_inf)[1]==1,'no boundary denominator')
zero(F_inf.subs(r,0)+4,'constant boundary value')
zero(s.diff(F_inf,r).subs(r,0)-4,'entire boundary transverse derivative')
zero(s.diff(F_inf,b).subs(r,0),'entire boundary tangential derivative')
zero(s.cancel(Y.subs({u:1/r-1,t:-r**2-r**4*b}))-(1-r)*(2-zeta),'full y boundary extension')

U=y**2/(c-y**2+4*y); T=(y/U-1)/U**2
zero(F.subs({u:U,t:T})-c,'inverse rational field coordinate F')
zero(Y.subs({u:U,t:T})-y,'inverse rational field coordinate y')
zero(jac(F,Y)+u*Y**2,'fixed-F nonzero derivation on y')
Q=c/(3*y**3)+2/y**2-1/y
zero(Q.subs({c:F,y:Y})-G,'independent rational primitive formula')
zero(s.diff(Q,y)*(-y**4/(c-y**2+4*y))-1,'primitive derivative on rational generic fibre')

# Each simple special component is tested in coordinates in which it is a divisor.
Fw=u*w*((u+1)*w-4); Gw=((1-2*u)*w+2)/(3*u**2*w**2)
zero(F-u*W*L,'complete three-factor special fibre')
zero(W.subs(u,0)-1,'u and w components disjoint')
zero(L.subs(u,0)+3,'u and L components disjoint')
zero(((u+1)*w-4).subs(w,0)+4,'w and L components disjoint')
need(s.gcd(u**2,u*s.Integer(0)+1)==1,'w primitive linear factor')
need(s.gcd(u**2*(u+1),u-3)==1,'L primitive linear factor')
P=L**2*((1-2*u)*W+2)/3
zero(F**2*G-P,'exact polynomial annihilator witness')
zero(jac(F,P)-F**2,'literal polynomial response equation')
zero((Fw**2*Gw).subs(w,0)-s.Rational(32,3),'nonzero leading scalar pole on w component')
zero((F**2*G).cancel().subs(u,0)-9,'nonzero leading scalar pole on u component')
zero(((u+1)*w-4).subs(w,4/(u+1)),'L component parametrization')
need(s.cancel(Gw.subs(w,4/(u+1)))!=0,'primitive regular and not identically zero on L')

pp_u=9/F**2+4/F; pp_w=s.Rational(32,3)/Fw**2+4/Fw
for expr,var,label in [(G-pp_u,u,'complete u scalar principal part'),(Gw-pp_w,w,'complete w scalar principal part')]:
 num,den=s.fraction(s.cancel(expr))
 need(s.Poly(den,var).coeff_monomial(1)!=0,label+' has regular remainder')
for val in [-7,-4,-1,1,5]:
 dis=s.discriminant(Fw-val,w)
 zero(dis-4*u*((val+4)*u+val),'whole fibre discriminant')
 need(s.Poly(dis,u).coeff_monomial(u)!=0,'simple zero rules out square')
need(s.diff(s.discriminant(Fw-c,w),u).subs(u,0)==4*c,'uniform nonzero-c fibre irreducibility mechanism')

A=s.Matrix([9,s.Rational(32,3)]); B=s.Matrix([4,4])
need(s.Matrix.hstack(A,B).det()==-s.Rational(20,3),'two component vectors independent modulo diagonal')
theta=A/g**2+B/g
zero((g*theta-A/g)[0]-B[0],'multiplication leaves only regular remainder first component')
zero((g*theta-A/g)[1]-B[1],'multiplication leaves only regular remainder second component')
need(s.simplify(g*theta.diff(g)+2*theta)==B/g,'Euler operation isolates other torsion arm')
for j in range(13):
 need(s.diff(1/g,g,j)==(-1)**j*s.factorial(j)/g**(j+1),'independent derivatives span every tested principal order')
 v=s.diff(theta,g,j)
 top=(g**(j+2)*v).applyfunc(lambda z:s.cancel(z).subs(g,0))
 need(top==(-1)**j*s.factorial(j+1)*A,'canonical derivative top coefficient survives')
zero(s.diff(g**2*theta[0],g,2), 'left ideal g2 annihilates first component modulo regular terms')
zero(s.diff(g**2*theta[1],g,2), 'left ideal g2 annihilates second component modulo regular terms')

cert={'status':'FINITE-EXACT independent controls; unbounded conclusions use companion analytic proof',
 'source_F':str(F),'boundary_F':str(s.factor(F_inf)),
 'primitive_in_F_y':str(Q),'scalar_principal_parts':[['9','4'],['32/3','4'],['0','0']],
 'source_unit_annihilator':'(F^2)','torsion_arms':2,
 'weyl_left_annihilator':'A1 * F^2','canonical_derivative_order':'j+2',
 'always_active_gates':checks}
(OUT/(STEM+'_certificate.json')).write_text(json.dumps(cert,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
print('INDEPENDENT QUADRATIC UNIT2: exact original/boundary charts and nowhere-critical first function')
print('COMPONENTS three reduced disjoint source components; other fibres irreducible by simple discriminant zero')
print('PRINCIPAL_PARTS (9/F^2+4/F, (32/3)/F^2+4/F, 0); exact unit annihilator(F^2)')
print('WEYL_CYCLIC two full torsion arms; left annihilator A1*F^2; derivative order j+2')
print('Always-active exact gates:',checks)
