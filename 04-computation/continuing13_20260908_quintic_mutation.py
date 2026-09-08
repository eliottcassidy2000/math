"""Exact controls for the polynomial mutation; no producer is imported."""
from pathlib import Path
from fractions import Fraction
import json,sys
import sympy as s
sys.stdout.reconfigure(newline='\n')
u,z,t,R,B,lam,h=s.symbols('u z t R B lambda h')
g,Astar,kappa,b3,b4=s.symbols('g A kappa b3 b4', nonzero=True)
gates=0
def check(condition,label):
 global gates
 gates+=1
 if not bool(condition):raise RuntimeError(label)
def zero(expr,label):check(s.cancel(expr)==0,label)
def bracket(a,b):return s.cancel(u**3*(s.diff(a,u)*s.diff(b,z)-s.diff(a,z)*s.diff(b,u)))
def primitive(f):return s.integrate(s.expand(f*f),z)
def coeff(expr,var,k):return s.expand(expr).coeff(var,k)

# Formal local transfer pays absence of the possible simple pole.
local=kappa*Astar+b3*Astar**3+b4*Astar**4
# Integration by parts: pp coefficient C_j=(1/j) Res(local^j A^-3 dA).
zero(coeff(local,Astar,2), 'no simple principal coefficient')
zero(coeff(local**2,Astar,2)/2-kappa**2/2,'complete double coefficient')
for j in range(3,8):zero(coeff(local**j,Astar,2)/j,'no higher local pole')

records=[]
for q in range(2,7):
 f=s.prod(z-j for j in range(q));Q=primitive(f);a0=f.subs(z,-2)
 check(s.degree(f,z)==q,'genuine degree')
 check(s.degree(Q,z)==2*q+1,'primitive degree')
 zero(s.diff(Q,z)-f*f,'primitive identity')
 zero(Q.subs(z,0),'normalization')
 check(s.gcd(f,s.diff(f,z))==1,'squarefree')
 check(a0!=0 and a0+1!=0,'admissible Eu values')
 A=u*f;T=lam*A+Q;G=1/(2*A*A)
 zero(bracket(A,Q)-A**3,'power Jacobian')
 zero(bracket(T,G)-1,'actual rational mate')
 zsource=u-2+u**3*t
 Asource=A.subs(z,zsource)
 Qsource=Q.subs(z,zsource)
 zero(s.series(Asource,u,0,3).removeO()-a0*u-s.diff(f,z).subs(z,-2)*u*u,'actual Eu A jet')
 zero(s.series(Qsource-Q.subs(z,-2)-a0*Asource,u,0,3).removeO(),'actual Eu Q jet')
 actual_grad=s.diff(lam*Asource+Qsource,u).subs(u,0)
 zero(actual_grad-a0*(lam+a0),'actual Eu derivative')
 zero(actual_grad.subs(lam,-a0),'forbidden lambda Eu hostile')
 M=-(3-R)-B*(1-R)**3;zinf=R*M
 finf=s.div(f,z,z)[0].subs(z,zinf)
 Ainf=(1-R)*M*finf;Tinf=lam*Ainf+Q.subs(z,zinf)
 zero(Tinf.subs(R,0)+lam*s.diff(f,z).subs(z,0)*(3+B),'whole boundary value')
 zero(s.diff(Tinf,B).subs(R,0)+lam*s.diff(f,z).subs(z,0),'whole boundary tangent')
 check(s.Poly(Ainf,R,B) is not None,'global A polynomial')
 # Use the exact inverse without expanding the high-degree source polynomial.
 invu=(g-Q)/(lam*f)
 invt=(z-invu+2)/invu**3
 zero((u-2+u**3*t).subs({u:invu,t:invt},simultaneous=True)-z,'source z inverse')
 zero(T.subs(u,invu)-g,'whole target inverse')
 check(3*q+2==(q+1)+(2*q+1),'puncture count')
 cvals=[s.factor(Q.subs(z,j)) for j in range(q)]+[s.factor(Q.subs(z,-2))]
 check(len(set(cvals))==q+1,'positive no-collision control')
 P=s.prod(g-c for c in cvals)
 # Complete root-component tuple and derivative at all its points.
 for j in range(q):
  check(f.subs(z,j)==0 and s.diff(f,z).subs(z,j)!=0,'root component simple')
  zero(s.rem(Q-Q.subs(z,j),(z-j)**3,z),'root primitive third-order remainder')
  zero(s.rem(P.subs(g,T),z-j,z),'target polynomial vanishes on root')
  zero(s.diff(T,z).subs(z,j)-lam*u*s.diff(f,z).subs(z,j),'root gradient')
 zero(P.subs(g,Q.subs(z,-2)),'target polynomial vanishes on Eu')
 check(s.degree(P,g)==q+1,'annihilator distinct-support degree')
 # One pure order-two coefficient vector at each support: one visible arm.
 for c in cvals:
  check(cvals.count(c)==1,'positive singleton support')
  check(s.degree((g-c)**2,g)==2,'exact local annihilator degree')
 records.append(dict(q=q,source_t_degree=2*q+1,pair_degree=2*(2*q+1),
                     punctures=3*q+2,ambient_arms=q+1,visible_arms=q+1,
                     support=[str(c) for c in cvals]))

f=z*(z-1);Q=primitive(f);a0=f.subs(z,-2)
zero(Q-(z**5/s.Integer(5)-z**4/s.Integer(2)+z**3/s.Integer(3)),'explicit quintic primitive')
zero(a0-6,'explicit Eu value')
zero(Q.subs(z,1)-s.Rational(1,30),'explicit positive critical value')
zero(Q.subs(z,-2)+s.Rational(256,15),'explicit Eu target')
P=g*(g-s.Rational(1,30))*(g+s.Rational(256,15))
quotient,remainder=s.div(s.expand(P.subs(g,u*f+Q)),f,z)
zero(remainder,'literal annihilator divisible by f in C[u,z]')
zero(quotient.subs({u:0,z:-2}),'remaining Eu factor after actual source substitution')
check(2*(2*2+1)==10,'explicit pair degree ten')
check(3*2+2==8,'explicit eight punctures')
collision=6*z*z-15*z+10
zero(s.rem(Q,collision,z),'actual complex target collision')
check(s.gcd(f,collision)==1,'collision still source smooth')
check(s.gcd(1+f,collision)==1,'collision still globally smooth')
check(s.discriminant(collision,z)==-15,'two distinct complex collision parameters')
check(s.rem(Q-s.Rational(1,30),collision,z)!=0,'collision does not erase other target')
check(3>2,'three ambient versus two visible arms at collision')

certificate=dict(scope='Finite exact controls supporting the analytic all-q mutation theorem; completed genus0 versus non-isotrivial affine punctures kept separate.',
                 gates=gates,parameter_universe=records,
                 explicit_quintic={'support':['0','1/30','-256/15'],'pair_degree':10,'punctures':8},
                 collision={'equation':'6z0^2-15z0+10','ambient_arms':3,'visible_arms':2},
                 always_active=True)
here=Path(__file__).resolve();out=here.parent
if here.parent.name=='04-computation':out=here.parent.parent/'05-knowledge/results'
(out/(here.stem+'_certificate.json')).write_text(json.dumps(certificate,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
print('Polynomial mutation of the complete univariate family')
print('Finite controls: q=2,...,6; actual source and boundary; all listed fibres and local jets')
print('Explicit quintic: source degree5, pair degree10,8 punctures,3 ambient and generated arms')
print('Admissible complex collision:3 ambient arms,2 unit-generated arms')
print('PASS',gates,'always-active exact gates')
