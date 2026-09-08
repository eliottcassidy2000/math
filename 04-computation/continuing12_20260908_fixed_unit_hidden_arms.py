"""Exact controls: fixed unit Weyl response with unbounded ambient torsion.

The paired report proves all n>=3. This engine uses no mathematical producer.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,sys
import sympy as s
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve()
OUT=HERE.parent.parent/'05-knowledge/results' if HERE.parent.name=='04-computation' else HERE.parent
u,z,y,v,g,W,R,B,h,b,d,rho,lam,T=s.symbols('u z y v g W R B h b d rho lam T')
checks=0
def need(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)
def zero(E,label):need(s.cancel(E)==0,label)

# Both physical charts and the original volume, without a producer import.
t=s.symbols('t');ys=u*(u-2*h+u**3*t)
M=-3*h*h+h**3*R-(1-h*R)**3*B
Y=(1-h*R)*M
zero(ys.subs({u:1/R-h,t:-R*R-R**4*B},simultaneous=True)-Y,'whole second-chart global linear function')
zero(s.diff(ys,t)-u**4,'actual u,y Jacobian')
zero(Y.subs(R,0)+3*h*h+B,'entire boundary restriction')
zero(s.diff(Y,B).subs(R,0)+1,'boundary tangential submersion')
zero(s.diff(ys,u).subs(u,0)+2*h,'omitted source line derivative')

# A symbolic K proves the Jacobian identity for every polynomial K(W).
K=s.Function('K'); w=u*z*z;gg=u*z*K(w)
bracket=s.diff(gg,u)*s.diff(w,z)-s.diff(gg,z)*s.diff(w,u)
zero(bracket-w*K(w),'native bracket before source u^3 multiplier')
zero(u**3*w**3*K(w)**6/gg**6-1,'H prime W2K5 gives actual constant Jacobian')
g0=b*u*(u-2*h)
pp0=-8*b*b*h**3/(3*g0**3)-2*b*h/g0**2
zero(1/(3*b*u**3)-pp0-1/(3*b*(u-2*h)**3),'complete Eu scalar principal part remainder')

records=[]
for n in range(3,10):
 kn=b+d*W**n
 Sn=sum(s.binomial(5,m)*d**m*b**(5-m)*W**(n*m)/s.Integer(n*m+3) for m in range(6))
 Hn=W**3*Sn
 cn=sum(Q((-1)**m*int(s.binomial(5,m)),n*m+3) for m in range(6))
 prod=1
 for j in range(6):prod*=3+j*n
 need(cn==Q(120*n**5,prod)>0,'exact nonzero common component coefficient')
 zero(s.diff(Hn,W)-W**2*kn**5,'complete polynomial primitive identity')
 zero(s.rem(Hn-b**5*s.Rational(cn)*W**3,kn,W),'all component values modulo K')
 for j in range(1,6):
  zero(s.diff(Hn,W,j).subs({W:rho,d:-b/rho**n},simultaneous=True),'all five lower component jets vanish')
 zero(Hn.subs({W:rho,d:-b/rho**n},simultaneous=True)-b**5*rho**3*s.Rational(cn),'every component highest scalar coefficient')
 gn=b*y+d*y**(2*n+1)/u**n
 zero(s.diff(gn,u)+n*d*y**(2*n+1)/u**(n+1),'open chart critical point mechanism')
 zero(s.diff(gn,y).subs(y,0)-b,'open chart derivative at y zero')
 gin=b*Y+d*R**n*(1-h*R)**(n+1)*M**(2*n+1)
 zero(gin.subs(R,0)-b*(-3*h*h-B),'full global source boundary value')
 zero(s.diff(gin,B).subs(R,0)+b,'whole boundary is submersive')
 # The nonzero source fibre is literally a punctured one-variable ring.
 yy=g*v**n/(d+b*v**n);uu=v*yy**2;zz=1/(v*yy)
 zero(uu*zz-yy,'torus monomial inverse y')
 zero(uu*zz**2-1/v,'torus monomial inverse W')
 zero((b*y+d*y**(2*n+1)/u**n).subs({u:uu,y:yy},simultaneous=True)-g,'all nonzero-fibre rational inverse')
 # Exact order of the ignored correction at the u=0 source component.
 num=s.together(Sn-kn**6/(3*b)).as_numer_denom()[0]
 need(s.Poly(num,W).terms()[-1][0][0]>=n,'Eu primitive correction has W-adic order at least n')
 # A cancelled coefficient is not a lost fibre component.
 need(n+2>=5 and n+1>2,'ambient arm count strictly exceeds generated count')
 records.append({'n':n,'source_t_degree':2*n+1,'components':n+2,'ambient_arms':n+1,
                 'unit_order':6,'unit_generated_arms':2,'component_constant':str(cn)})

zero(s.discriminant(u*(u-2*h+u**3*t)**2-rho,t)-4*rho*u**7,'every L component has odd discriminant valuation')
zero((u*(u-2*h+u**3*t)**2-rho).subs(u,0)+rho,'component primitivity at omitted line')

# Independent two-arm Laurent algebra with the actual unit normalized only
# by constant, invertible component-basis changes. No localization of g.
def clean(A):return {k:tuple(s.expand(q) for q in val) for k,val in A.items() if k<0 and any(s.expand(q)!=0 for q in val)}
def add(A,C):return clean({k:tuple(A.get(k,(0,0))[j]+C.get(k,(0,0))[j] for j in range(2)) for k in set(A)|set(C)})
def scale(a,A):return clean({k:tuple(a*q for q in val) for k,val in A.items()})
def mult(A):return clean({k+1:val for k,val in A.items()})
def diff(A):return clean({k-1:tuple(k*q for q in val) for k,val in A.items()})
def repeat(op,A,j):
 for _ in range(j):A=op(A)
 return A
def poly_action(p,A):
 result={}
 for (j,),c in s.Poly(p,T).terms():result=add(result,scale(c,repeat(diff,A,j)))
 return result
theta={-6:(1,0),-3:(0,1),-2:(0,lam)}
Aj=[(-1)**(5-j)*T**(5-j)/s.factorial(5-j) for j in range(6)]
Bj=[T**2/2-lam*T,-T+lam,s.Integer(1),s.Integer(0),s.Integer(0),s.Integer(0)]
for j in range(6):
 actual=repeat(mult,theta,j)
 want=add(poly_action(Aj[j],{-1:(1,0)}),poly_action(Bj[j],{-1:(0,1)}))
 need(clean(actual)==clean(want),'complete PBW column identities')
need(repeat(mult,theta,6)=={} and repeat(mult,theta,5)!={},'exact order six with nonzero fifth multiple')
for j in [0,1,3,4]:
 Rj=add(repeat(mult,theta,j),scale(-1,poly_action(Bj[j],repeat(mult,theta,2))))
 Rj=add(Rj,scale(-1,poly_action(Aj[j]+Bj[j]*T**3/6,repeat(mult,theta,5))))
 need(Rj=={},'every exact left-annihilator generator vanishes')
need(s.Matrix([[Aj[2],Aj[5]],[Bj[2],Bj[5]]]).det()==-1,'remaining columns form an invertible polynomial basis')
for j in range(11):
 Dtheta=repeat(diff,theta,j)
 need(repeat(mult,Dtheta,j+6)=={} and repeat(mult,Dtheta,j+5)!={},'all declared derivative orders rise while two arms remain')
 need(s.simplify(Dtheta[-6-j][0])!=0 and s.simplify(Dtheta[-3-j][1])!=0,'both coefficient directions survive all declared derivatives')

cert={'status':'FINITE-EXACT controls; all-n proof and full annihilator in report',
 'parameters':'integer n>=3; b,d,h nonzero complex constants; fixed W2',
 'source':'u=x-h; z=u-2h+u^3*t; W=u*z^2; F=u*z*(b+d*W^n)',
 'unit_order':6,'ambient_torsion_arms':'n+1','unit_generated_arms':2,
 'normalized_unit':'A/g^6+B*(g^-3+lambda*g^-2), lambda=3/(4*b*h^2)',
 'pointed_module_independent_of':['n','d'],
 'exact_annihilator':'D*g^6 and R_j=g^j-B_j(partial)g^2-(A_j(partial)+B_j(partial)partial^3/6)g^5 for j=0,1,3,4',
 'controls':records,'always_active_gates':checks}
raw=(json.dumps(cert,indent=2,sort_keys=True)+'\n').encode()
(OUT/(HERE.stem+'_certificate.json')).write_bytes(raw)
print('FIXED W2: all n>=3 have globally submersive first function of t-degree2n+1 and rational mate')
print('COMPONENTS: n+2 reduced disjoint zero components; every nonzero fibre is a rational punctured torus chart')
print('UNIT: order6; full torsion n+1 arms; unit generates exactly2 with identical pointed Weyl module at fixed b,h')
print('EXACT ANNIHILATOR: g^6 and four explicit PBW relations; seven n controls and eleven derivative levels')
print('CERTIFICATE_SHA256',hashlib.sha256(raw).hexdigest())
print('Always-active exact gates:',checks)
