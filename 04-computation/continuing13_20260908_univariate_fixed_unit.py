"""Exact controls for the whole univariate W2 family and its fixed unit module.

Analytic all-degree proofs are in the paired report. No producer is imported.
"""
from pathlib import Path
from hashlib import sha256
import json,sys
import sympy as S
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve()
DEST=HERE.parent.parent/'05-knowledge/results' if HERE.parent.name=='04-computation' else HERE.parent
u,z,t,h,r,b,g,c,d,ell=S.symbols('u z t h r b g c d ell')
gates=0
def need(ok,msg):
 global gates
 gates+=1
 if not ok:raise ArithmeticError(msg)
def zero(e,msg):need(S.cancel(e)==0,msg)
def jac(F,G,X,T):return S.diff(F,X)*S.diff(G,T)-S.diff(F,T)*S.diff(G,X)

alpha=-2*h
zs=u-2*h+u**3*t
J=h*h*(3-h*r)+b*(1-h*r)**3
zero(zs.subs({u:1/r-h,t:-r*r-r**4*b},simultaneous=True)+r*J,'actual complete z boundary chart')
zero(jac(u,zs,u,t)-u**3,'actual source coordinate Jacobian')

# Formal first three jets prove the complete principal parts at the retained source line.
f0,f1,f2,f3,Q0=S.symbols('f0 f1 f2 f3 Q0')
qjet=Q0+f0*f0*(u+u**3*t)+f0*f1*u*u
fjet=f0+f1*u+f2*u*u/2
zero(S.expand(qjet-Q0-f0*u*fjet).coeff(u,1),'Eu exact missing first-order remainder')
zero(S.expand(qjet-Q0-f0*u*fjet).coeff(u,2),'Eu exact missing second-order remainder')
need(S.expand(qjet-Q0-f0*u*fjet).coeff(u,3)!=0,'Eu omitted higher coordinate jet is generally nonzero')

records=[]
for degree in range(2,9):
 f=S.prod(z-j for j in range(degree)).expand()
 Q=S.integrate(f*f,z);Q=S.expand(Q-Q.subs(z,0))
 need(S.degree(f,z)==degree,'declared exact source-degree universe')
 zero(Q.diff(z)-f*f,'exact primitive polynomial differential')
 need(S.degree(Q,z)==2*degree+1,'exact rational-pair degree')
 need(S.gcd(f,f.diff(z))==1,'declared f is squarefree')
 zero(f.subs(z,0),'globality zero at boundary-origin label')
 need(f.diff(z).subs(z,0)!=0,'boundary tangential derivative nonzero')
 need(f.subs(z,-2)!=0,'retained source line submersion')
 need(Q.subs(z,1)!=0,'at least one nonzero component principal coefficient in this control')
 F=u*f
 zero(u**3*jac(F,Q,u,z)-F**3,'actual Hamiltonian polynomial identity before division')
 zero(u**3*jac(F,Q/F**3,u,z)-1,'rational mate exact Jacobian')
 zero(F.subs(u,c/f)-c,'complete nonzero fibre inverse')
 ts=(z-u+2*h)/u**3
 zero(zs.subs(t,ts)-z,'complete original-source t reconstruction')
 fi=S.Poly(f,z)
 Finf=(1-h*r)*sum(fi.nth(i)*(-1)**i*r**(i-1)*J**i for i in range(1,degree+1))
 zero(Finf.subs(r,0)+f.diff(z).subs(z,0)*(3*h*h+b),'complete boundary restriction')
 zero(S.diff(Finf,b).subs(r,0)+f.diff(z).subs(z,0),'complete boundary submersion derivative')
 for root in range(degree):
  need(f.diff(z).subs(z,root)!=0,'every source component is reduced')
  need(S.rem(Q-Q.subs(z,root),(z-root)**3,z)==0,'complete root principal part has no lower poles')
 eu_f=f.subs(z,-2)
 need(eu_f!=0,'Eu second principal coefficient is nonzero')
 need(3*degree>2*degree+1,'all-root vanishing degree contradiction in the declared universe')
 if degree<=3:
  fs=f.subs(z,zs);Qs=Q.subs(z,zs);Fs=u*fs
  zero(jac(Fs,Qs,u,t)-Fs**3,'independent full original polynomial Jacobian control')
  numerator=S.Poly(S.expand(Qs-Q.subs(z,-2*h)-f.subs(z,-2*h)*Fs),u)
  for i in range(3):zero(numerator.nth(i),'complete Eu remainder divisible by u cubed')
 records.append(dict(degree=degree,f=str(f),Q=str(Q),ambient_arms=degree,unit_arms=2,pair_degree=2*degree+1))

# First cubic gain: three ambient arms, two unit arms.
f=z*(c+d*z*z)
Q=d*d*z**7/7+2*c*d*z**5/5+c*c*z**3/3
zero(Q.diff(z)-f*f,'symmetric cubic complete primitive')
zero(f.subs(z,-2*h)+2*h*(c+4*d*h*h),'symmetric cubic actual source submersion wall')
rho=S.symbols('rho')
zero(Q.subs({z:rho,c:-d*rho*rho})-8*d*d*rho**7/105,'symmetric cubic nonzero root coefficient')

# Hostiles: each required hypothesis excludes an actual failure.
bad=z*z
zero(bad.subs(z,0),'multiple-root hostile source value')
zero(bad.diff(z).subs(z,0),'multiple-root hostile source gradient')
need((z+1).subs(z,0)!=0,'nonzero constant f creates a true R^-1 boundary pole')
bad=z*(z+2)
zero(bad.subs(z,-2),'Eu collision hostile creates actual source critical line')
linear=z
Qlinear=z**3/3
zero(Qlinear.subs(z,0),'degree-one boundary all root coefficients vanish')
need(3*S.degree(linear,z)==S.degree(Qlinear,z),'degree-one boundary defeats strict degree contradiction')

# A genuine critical-value collision must not be mistaken for rank collapse.
a=(7+S.sqrt(-7))/14
f=z*(z-1)*(z-a)
Q=S.integrate(S.expand(f*f),z)
zero(S.simplify(Q.subs(z,1)),'one nonzero source-root principal coefficient can vanish')
need(S.simplify(Q.subs(z,a))!=0,'another root retains the required coefficient direction')
need(S.simplify(f.subs(z,-2))!=0,'collision control remains source-submersive at Eu')

# Full Weyl relation compiler: formal action, separate from geometric identities.
E=S.Matrix([1,0]);B=S.Matrix([0,1])
theta=E/g**3+B/g**2
def proj(v):
 return v.applyfunc(lambda q:S.Add(*[term for term in S.expand(q).as_ordered_terms() if term.as_powers_dict().get(g,0)<0]))
def mul(v):return proj(g*v)
def der(v):return v.diff(g)
need(mul(mul(mul(theta)))==S.zeros(2,1),'g cubed kills the unit')
need((theta+der(mul(theta))+der(der(mul(mul(theta)))))/2 != S.zeros(2,1),'hostile half applies only to the second derivative term')
need(theta+der(mul(theta))+der(der(mul(mul(theta))))/2==S.zeros(2,1),'complete exact second Weyl generator')
need(mul(theta)==-der(E/g)+B/g,'first C[derivative] basis column')
need(mul(mul(theta))==E/g,'second C[derivative] basis column')
zero(S.det(S.Matrix([[-ell,1],[1,0]]))+1,'basis determinant is minus one')
for j in range(9):
 exact=(-1)**j*S.factorial(j+2)*E/(2*g**(j+3))+(-1)**j*S.factorial(j+1)*B/g**(j+2)
 need((theta.diff(g,j)-exact).applyfunc(S.cancel)==S.zeros(2,1),'canonical derivatives exact all displayed pole orders')

cert=dict(status='Analytic all-degree theorem with FINITE-EXACT independent controls',
 family='u=x-h,z=u-2h+u^3*t,F=u*f(z)',
 globality='f(0)=0',submersion='f squarefree and f(-2*h)!=0',
 main_universe='all complex polynomials f of degree q>=2 satisfying these conditions',
 primitive='Q(0)=0,Qprime=f^2,G=Q(z)/F^3',
 constants='C(F)',nonzero_fibre='C[z,f(z)^-1]',pair_degree='2*q+1',
 full_source_arms='q',unit_generated_arms=2,unit_order=3,
 exact_unit_annihilator=['D*g^3','D*(1+derivative*g+derivative^2*g^2/2)'],
 principal_parts={'Eu':'Q(-2h)/g^3+f(-2h)/g^2','root_rho':'Q(rho)/g^3'},
 controls=records,gates=gates,
 scope='Whole univariate ansatz on fixed W2; original C[x,t] response ring; no classification of arbitrary cubics or global regular-mate claim')
raw=json.dumps(cert,sort_keys=True,separators=(',',':')).encode()+b'\n'
(DEST/(HERE.stem+'_certificate.json')).write_bytes(raw)
print('WHOLE FAMILY: global iff f(0)=0; submersive iff f is squarefree and f(-2h)!=0.')
print('ALL q>=2: exact unit order3; fixed pointed two-arm Weyl module; ambient torsion has q arms.')
print('PAIR DEGREE: exactly2q+1. Degree3 already has three ambient arms; unit-order1 remains open outside this ansatz.')
print('CERTIFICATE_SHA256',sha256(raw).hexdigest())
print('Always-active exact gates:',gates)
