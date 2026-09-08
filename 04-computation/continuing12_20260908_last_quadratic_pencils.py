"""Exact controls for the final two independent global W2 quadratic pencils.

All-parameter arguments are in the paired report; no mathematical producer
is imported.  Every gate remains active under Python optimization.
"""
from pathlib import Path
from hashlib import sha256
import json, sys
import sympy as S
sys.stdout.reconfigure(encoding='utf-8', newline='\n')
HERE = Path(__file__).resolve()
DEST = HERE.parent.parent/'05-knowledge/results' if HERE.parent.name == '04-computation' else HERE.parent
u,x,t,w,v,a,d,h,p,k,s,c,r,b,C = S.symbols('u x t w v a d h p k s c r b C')
gates=0
def need(ok,msg):
 global gates
 gates+=1
 if not ok: raise ArithmeticError(msg)
def zero(e,msg): need(S.cancel(e)==0,msg)
def jac(f,g,X,T): return S.diff(f,X)*S.diff(g,T)-S.diff(f,T)*S.diff(g,X)

# The first pencil: complete successive leading-term exclusions.
p0,p1,p2,p3,p4 = S.symbols('p0 p1 p2 p3 p4')
N=a*x+d
P=p4*x**4+p3*x**3+p2*x*x+p1*x+p0
Q=p4*x*x+p3*x+s
Delta=S.expand(P*P-4*N*Q)
for j,var,previous in [(8,p4,{}),(6,p3,{p4:0}),(4,p2,{p4:0,p3:0}),(2,p1,{p4:0,p3:0,p2:0})]:
 zero(S.expand(Delta.subs(previous)).coeff(x,j)-var*var,'linear pencil successive top coefficient')
Fl=(a*x+d)*t*t+p*t+s
zero(S.diff(Fl,x)-a*t*t,'linear pencil source critical condition')
zero(S.diff(Fl,t).subs(t,0)-p,'linear pencil source noncritical surviving derivative')
zero(S.det(S.Matrix([[d,p*p-4*d*s],[a,-4*a*s]]))+a*p*p,'linear pencil exact rank')
Fli=S.cancel(Fl.subs({x:1/r,t:-r*r-r**4*b},simultaneous=True))
need(S.denom(Fli)==1,'linear pencil globally polynomial')
zero(S.diff(Fli,r).subs(r,0),'linear pencil whole boundary critical normal')
zero(S.diff(Fli,b).subs(r,0),'linear pencil whole boundary critical tangential')

# Elliptic pencil: exact whole coefficient and rank matching.
N=u**5*(a*u**3+d)
P=u**3*(2*a*u**3-4*a*h*u*u+p*u+2*d+k)
Q=s+a*u**4-4*a*h*u**3+(4*a*h*h+p)*u*u+(d+k-2*h*p)*u
Nx=S.Poly(S.expand(N.subs(u,x-h)),x)
Px=S.Poly(S.expand(P.subs(u,x-h)),x)
Qx=S.Poly(S.expand(Q.subs(u,x-h)),x)
zero(Px.nth(6)-2*Nx.nth(8),'elliptic global degree-six coefficient')
zero(Px.nth(5)-2*Nx.nth(7),'elliptic global degree-five coefficient')
for j,value in [(4,Nx.nth(8)),(3,Nx.nth(7)),(2,Px.nth(4)-Nx.nth(6)),(1,Px.nth(3)-Nx.nth(5))]:
 zero(Qx.nth(j)-value,'elliptic complete global lower matching')
D=S.expand((P*P-4*N*Q)/u**5)
zero(D-(-4*d*s+(k*k+8*d*h*p)*u+2*(p*k-8*a*d*h*h)*u*u+(p*p-8*a*h*k-4*a*s)*u**3),'elliptic complete discriminant')
zero(S.det(S.Matrix([[d,-4*d*s],[a,p*p-8*a*h*k-4*a*s]]))-d*(p*p-8*a*h*k),'elliptic exact pencil rank')
Y=u*(u*w-2*h)
Fc=s+a*Y*Y+p*Y+d*u*w*w+k*u*w
zero(Fc.subs(w,1+u*u*t)-(N*t*t+P*t+Q),'elliptic source factorization')
M=2*h-r*(h*h+b*(1-h*r)**2)
Xi=(1-h*r)*M
Yi=-(1-h*r)*(h*h*(3-h*r)+b*(1-h*r)**3)
Fci=s+a*Yi*Yi+p*Yi+d*r*(1-h*r)*M*M+k*Xi
zero(Fc.subs({u:1/r-h,w:r*M},simultaneous=True)-Fci,'elliptic full boundary polynomial')
zero(Fci.subs({r:0,b:C-3*h*h})-(s+a*C*C-p*C+2*k*h),'elliptic boundary value')
zero(S.diff(Fci,b).subs({r:0,b:C-3*h*h})-(2*a*C-p),'elliptic boundary tangential derivative')
normal=(-2*a*C+p)*4*h*(C-2*h*h)+4*d*h*h-k*C
zero(S.diff(Fci,r).subs({r:0,b:C-3*h*h})-normal,'elliptic boundary normal derivative')
sol={p:-k*k/(8*d*h),a:-k**3/(64*d*d*h**3)}
zero((k*k+8*d*h*p).subs(sol),'nonzero translation first exact equation')
zero((p*k-8*a*d*h*h).subs(sol),'nonzero translation second exact equation')
zero((p*p-8*a*h*k).subs(sol)-9*k**4/(64*d*d*h*h),'nonzero translation independence forces nonzero k')
criticalC=4*d*h*h/k
zero((2*a*C-p).subs(sol).subs(C,criticalC),'nonzero translation actual boundary critical tangential')
zero(normal.subs(sol).subs(C,criticalC),'nonzero translation actual boundary critical normal')
zero(Fci.subs({h:0,k:0})-(s+a*b*b-p*b+d*r**3*b*b),'zero translation complete boundary chart')
zero(S.diff(Fci,b).subs({h:0,k:0,r:0,b:p/(2*a)}),'nonzero quadratic boundary coefficient critical tangential')
zero(S.diff(Fci,r).subs({h:0,k:0,r:0,b:p/(2*a)}),'nonzero quadratic boundary coefficient critical normal')

# The surviving elliptic family, in the actual source and added chart.
L=d*w+p*u
g=u*w*L
G=1/(3*u*u*L)
ga=S.expand(g.subs(w,1+u*u*t))
Ga=G.subs(w,1+u*u*t)
zero(S.diff(ga,u).subs(u,0)-d,'elliptic source u=0 actual noncritical derivative')
zero(S.diff(ga,t).subs(u,0),'elliptic source u=0 tangential derivative')
zero(S.diff(g,u)-w*(d*w+2*p*u),'elliptic open gradient normal')
zero(S.diff(g,w)-u*(2*d*w+p*u),'elliptic open gradient tangential')
zero(S.diff(g,w).subs(w,0)-p*u*u,'elliptic open zero-w hostile excluded')
zero(S.diff(g,w).subs(w,-2*p*u/d)+3*p*u*u,'elliptic open nonzero-w hostile excluded')
ginf=S.cancel(ga.subs({u:1/r,t:-r*r-r**4*b},simultaneous=True))
zero(ginf-(-p*b+d*r**3*b*b),'elliptic surviving complete boundary chart')
zero(S.diff(ginf,b).subs(r,0)+p,'elliptic surviving boundary gradient never zero')
zero(jac(ga,Ga,u,t)-1,'elliptic rational mate source Jacobian')
zero(jac(u,(1+u*u*t)/u,u,t)-u,'elliptic u,v source Jacobian')
zero(jac(ga,(1+u*u*t)/u,u,t)-3*ga,'elliptic exact derivation on v')
zero(g.subs(w,u*v)-u**3*v*(d*v+p),'elliptic cubic field equation')
zero(g*g*G-w*w*L/3,'elliptic polynomial repair of exact order two')
zero(g*G-w/(3*u),'elliptic first repair has genuine source pole')
zero(ga.subs(u,0),'source special Eu component')
zero(L.subs(w,0)-p*u,'source Ew and EL disjoint when u invertible')
zero(L.subs({u:0,w:1})-d,'source Eu and EL disjoint')

# Scalar principal parts retain the original coordinate jet w=1+u²t.
ppEu=d/(3*ga*ga)+p/(3*d*ga)
rem=S.cancel(Ga-ppEu)
need(S.denom(rem).subs(u,0)!=0,'complete Eu scalar principal-part remainder regular')
zero(S.limit(u*u*Ga,u,0)-1/(3*d),'Eu leading Laurent coefficient')
zero(S.limit(u*(Ga-1/(3*d*u*u)),u,0)+p/(3*d*d),'Eu second Laurent coefficient')
need(S.denom(G.subs(w,0))!=0,'Ew primitive is regular in generic component local ring')
zero((g*G).subs(w,-p*u/d)+p/(3*d),'EL complete simple scalar coefficient')
A2=S.Matrix([d/S.Integer(3),0])
A1=S.Matrix([p/(3*d),-p/(3*d)])
zero(S.det(S.Matrix.hstack(A2,A1))+p/9,'two independent component coefficient vectors')

# Irreducible nonzero fibres and genus-one squarefree quartic, no rational-field shortcut.
DeltaC=S.discriminant(ga-c,t)
zero(DeltaC-u**5*(p*p*u**3+4*d*c),'complete nonzero source fibre discriminant')
zero((ga-c).subs(u,0)+c,'nonzero source fibre primitive in t')
quartic=u*(p*p*u**3+4*d*c)
zero(S.discriminant(quartic,u)+27*p**4*(4*d*c)**4,'elliptic squarefree quartic discriminant')
need(S.Poly(quartic,u).degree()==4,'elliptic smooth completion genus formula uses degree four')

# Exact Weyl action in the inherited component-principal-part model.
z=S.symbols('z')
theta=A2/z**2+A1/z
def proj(e):
 return e.applyfunc(lambda q:S.Add(*[term for term in S.expand(q).as_ordered_terms() if term.as_powers_dict().get(z,0)<0]))
def m(e):return proj(z*e)
def diff(e):return e.diff(z)
need(m(m(theta))==S.zeros(2,1),'g squared annihilates unit tuple')
need(m(theta)==A2/z,'g unit recovers complete first direction')
need(m(diff(theta))+2*theta==A1/z,'Euler polynomial recovers complete second direction')
for j in range(7):
 expected=(-1)**j*S.factorial(j)*A1/z**(j+1)+(-1)**j*S.factorial(j+1)*A2/z**(j+2)
 need((theta.diff(z,j)-expected).applyfunc(S.cancel)==S.zeros(2,1),'canonical connection exact pole-order growth')
controls=[]
for dd,pp in [(1,1),(1,-2),(-3,2),(S.Rational(2,3),S.Rational(-5,7))]:
 need(S.det(S.Matrix.hstack(A2,A1).subs({d:dd,p:pp}))!=0,'declared parameter control two-arm rank')
 zero((jac(ga,Ga,u,t)-1).subs({d:dd,p:pp}),'declared exact rational mate parameter control')
 controls.append(dict(d=str(dd),p=str(pp),rank=2))

cert=dict(status='FINITE-EXACT controls; all-parameter analytic proofs in paired report',
 pencils=['span{1,x}','u^5 span{1,u^3}'],
 first_pencil_global_submersions=False,
 elliptic_global_submersion='h=a=k=0, d*p!=0; target shift arbitrary',
 surviving_first_function='s+x*(1+x^2*t)*(d*(1+x^2*t)+p*x)',
 rational_mate='1/(3*x^2*(d*(1+x^2*t)+p*x))',
 field='C(g)(v,u), u^3=g/(v*(d*v+p)); geometrically integral cubic extension',
 generic_completed_fibre_genus=1,source_special_components=3,full_torsion_arms=2,
 exact_unit_order=2,unit_weyl_annihilator='D*g^2',
 principal_parts=['d/(3*g^2)+p/(3*d*g)','0','-p/(3*d*g)'],
 controls=controls,gates=gates,
 scope='Independent two-dimensional exact quadratic pencils on fixed W2; original C[x,t] response ring; no dependent-pencil or global regular-mate claim')
raw=json.dumps(cert,sort_keys=True,separators=(',',':')).encode()+b'\n'
(DEST/(HERE.stem+'_certificate.json')).write_bytes(raw)
print('LINEAR PENCIL: every source-submersive member has the entire added divisor critical.')
print('ELLIPTIC PENCIL: whole global-submersion locus h=a=k=0, d*p!=0; all charts paid.')
print('ELLIPTIC SURVIVOR: rational mate, generic completed genus one, exact source unit order two.')
print('POINTED WEYL MODULE: exactly D/Dg^2, same unit module as the earlier genus-zero family.')
print('CERTIFICATE_SHA256',sha256(raw).hexdigest())
print('Always-active exact gates:',gates)
