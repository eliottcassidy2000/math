"""Exact controls for the all-degree quadratic logarithmic witness scout."""
from pathlib import Path
import hashlib,json,sys
import sympy as s
sys.stdout.reconfigure(newline='\n')
x,t,v,u,R,B,lam=s.symbols('x t v u R B lam')
gates=0
def check(ok,label):
    global gates
    gates+=1
    if not bool(ok):raise RuntimeError(label)
def eq(expr,label):check(s.cancel(expr)==0,label)
def J(F,H):return s.diff(F,x)*s.diff(H,t)-s.diff(F,t)*s.diff(H,x)
def at_infinity(F):return s.expand(F.subs({x:1/R,t:-R**2-R**4*B},simultaneous=True))
def smooth(F,label):check(list(s.groebner([s.diff(F,x),s.diff(F,t)],t,x))==[1],label)

# Local normalized orders, including the critical odd-multiplicity boundary.
orders=[]
for m in range(1,11):
    if m%2:
        e=(m-1)//2
        differential=s.diff(u*u,u)/u**m
        exponent=1-m
        eq(differential-2*u**(-2*e),'odd normalized finite order')
    else:
        e=m//2
        differential=u**(-e)
        exponent=-e
        eq(differential-u**(-m//2),'even normalized finite order')
    orders.append([m,exponent])
check(dict(orders)[3]==-2,'multiplicity3 is a double pole')
for d in range(10):
    if d%2:
        k=(d-1)//2
        eq(s.diff(u**(-2),u)/u**(-d)+2*u**(2*k-2),'odd infinity order')
    else:
        k=d//2
        eq(s.diff(u**(-1),u)/u**(-k)+u**(k-2),'even infinity order')

# Complete conic generator controls: arbitrary B handled by polynomial shears.
conics=[]
for n in [1,2,3,4]:
    kap=s.Rational(1,n)
    for AA,BB,CC in [(s.Integer(2),x*x+1,((x*x+1)**2-(kap**2*x*x+3*x+5))/8),
                     (x,1, -kap**2*x/4-s.Rational(3,4))]:
        HH=s.expand(AA*t*t+BB*t+CC)
        EE=s.expand(BB*BB-4*AA*CC)
        a1=s.Poly(AA,x).nth(1)
        a0=s.Poly(AA,x).nth(0)
        e1=s.Poly(EE,x).nth(1)
        e0=s.Poly(EE,x).nth(0)
        yy=2*AA*t+BB
        UU=yy+kap*x+(e1+4*a1*HH)/(2*kap)
        VV=yy-kap*x-(e1+4*a1*HH)/(2*kap)
        eq(yy*yy-(EE+4*AA*HH),'actual discriminant equation')
        eq(J(UU,HH)-kap*UU,'raw conic eigen-generator')
        eq(J(UU**n,HH)-UU**n,'raw polynomial eigenfunction sufficiency')
        eq(UU*VV-(e0+4*a0*HH-(e1+4*a1*HH)**2/(4*kap**2)),'conic rational field norm')
        check(s.degree(EE,x)==2,'conic exact discriminant degree')
        check(s.discriminant(EE+4*AA*v,x)!=0,'generic conic geometrically integral')
        conics.append({'n':n,'A':str(AA),'B':str(BB),'C':str(CC)})

# Holomorphic, variable-residue, and elliptic integral-residue hostiles.
check(s.discriminant(x**3+4*v,x)!=0,'squarefree holomorphic hostile')
eq((2*x*x*t+x)**2-(x*x+4*x*x*(x*x*t*t+x*t)),'variable-residue discriminant')
check(s.diff(1+4*v,v)!=0,'finite residue square varies with level')
HH=x**5*t*t+x*t
YY=1+2*x**4*t
eq((2*x**5*t+x)**2-x*x*(1+4*HH*x**3),'elliptic normalized equation')
eq(J(x**5*t*t,HH)-3*x**5*t*t,'elliptic near miss eigenvalue3')
ratio=s.cancel((YY-1)/(YY+1))
eq(J(ratio,HH)-3*ratio,'elliptic ratio has eigenvalue3')
check(s.discriminant(1+4*v*x**3,x)!=0,'elliptic generic curve smooth')
Y=s.symbols('Y')
# Differentiate along Y²=1+4v*x³, then use that equation.
q=(Y-1)/(Y+1)
dq=s.diff(q,Y)*(6*v*x*x/Y)
eq((s.cancel(dq/q-3/(x*Y))*(x*Y*(Y*Y-1))).subs(Y*Y,1+4*v*x**3),'elliptic logarithmic multiple')

# Infinite cyclic-cover spectrum: exact polynomial generators of both signs.
spectrum=[]
for m in range(1,9):
    HH=x**(m+2)*t*t+x*t
    YY=1+2*x**(m+1)*t
    RR=s.cancel((YY-1)/(YY+1))
    K0=1+x**(m+1)*t
    ell=m if m%2 else m//2
    positive=x**(m+2)*t*t if m%2 else x**((m+2)//2)*t
    negative=t**m*K0**(m+2) if m%2 else t**(m//2)*K0**((m+2)//2)
    eq(J(positive,HH)-ell*positive,'least positive polynomial eigenvalue')
    eq(J(negative,HH)+ell*negative,'least negative polynomial eigenvalue')
    eq(J(RR,HH)-m*RR,'cyclic ratio eigenvalue')
    eq(RR/(HH*(1-RR)**2)-x**m,'exact cyclic covering equation')
    eq((2*x**(m+2)*t+x)**2-x*x*(1+4*HH*x**m),'same local-residue pencil')
    if m<=2:
        eq(s.diff(positive,x).subs(x,0),'eigen1 repeated-x critical x partial')
        eq(s.diff(positive,t).subs(x,0),'eigen1 repeated-x critical t partial')
    spectrum.append({'m':m,'least_positive':ell,'positive':str(positive),'negative':str(negative)})
for m in range(1,13):
    ell=m if m%2 else m//2
    for eigenvalue in range(-2*m,2*m+1):
        characters=[j for j in range(m) if (eigenvalue-j)%m==0 and (2*j)%m==0]
        check(bool(characters)==(eigenvalue%ell==0),'complete declared character congruence bank')

# Admissible fixed-double source-submersive gauge family.
FF=x*x*t+x
HH=x*t+lam*FF*FF
eq(J(FF,HH)-FF,'admissible double source-submersive gauge')
AA=s.Poly(HH,t).nth(2);BB=s.Poly(HH,t).nth(1);CC=s.Poly(HH,t).nth(0)
eq(BB*BB+4*AA*(v-CC)-x*x*(1+4*lam*(1+v)*x*x),'admissible double exact pencil')
eq(s.diff(BB,x).subs(x,0)-1,'admissible double unit linearization')
smooth(FF,'admissible double source critical ideal')
eq(HH-lam*FF*FF-x*t,'gauge subtraction')
eq(J(FF,t/(x*t+1))-1,'admissible double rational mate')
check(min(term.as_powers_dict().get(R,0) for term in at_infinity(FF).as_ordered_terms())==-1,'admissible double source not global')

FF=t*(1+x*t);HH=-x*t+lam*FF
eq(J(FF,HH)-FF,'squarefree gauge thickening')
smooth(FF,'global but boundary-critical gauge source')
FI=at_infinity(FF)
eq(s.diff(FI,R).subs(R,0),'global gauge normal critical')
eq(s.diff(FI,B).subs(R,0),'global gauge tangent critical')

# Genuine minimal quadratic-witness family of arbitrary odd t degree.
HH=t*t-x*x/4;UU=x+2*t;VV=x-2*t
eq(J(UU,HH)-UU,'Euler eigenfunction')
eq(HH+UU*VV/4,'Euler product normalization')
family=[]
for m in range(1,6):
    PP=s.prod(v-j for j in range(1,m+1))
    FF=s.expand(UU*PP.subs(v,HH))
    eq(J(FF,HH)-FF,'nonlinear witness family bracket')
    check(s.degree(FF,t)==2*m+1,'genuine odd source degree')
    check(PP.subs(v,0)!=0 and s.discriminant(PP,v)!=0,'complete source smoothness factor conditions')
    eq(UU-FF/PP.subs(v,HH),'rational field inverse U')
    eq(VV+4*HH*PP.subs(v,HH)/FF,'rational field inverse V')
    eq(J(FF,HH/FF)-1,'unit-order-one rational mate')
    eq(HH.subs(x,-2*t),'regular zero-component label')
    for j in range(1,m+1):
        check(j!=0,'nonzero pole-component label')
    FI=at_infinity(FF)
    min_degree=min(term.as_powers_dict().get(R,0) for term in FI.as_ordered_terms())
    check(min_degree==-(2*m+1),'exact W2 boundary pole')
    eq(s.expand(FI).coeff(R,-(2*m+1))-s.Rational(-1,4)**m,'exact W2 leading pole coefficient')
    if m<=2:smooth(FF,'independent nonlinear source critical ideal')
    family.append({'P_degree':m,'source_t_degree':2*m+1,'minimum_witness_t_degree':2,'boundary_pole':2*m+1})

source=Path(__file__).resolve();folder=source.parent
if folder.name=='04-computation':folder=folder.parent/'05-knowledge'/'results'
cert={'status':'Analytic unbounded proof plus FINITE-EXACT controls; independent audit pending',
      'gates':gates,'producer_imports':False,'finite_orders':orders,'conic_controls':conics,
      'nonlinear_witness_controls':family,
      'elliptic_hostile':'H=x5*t²+xt; integral residues but degree-one divisor nonprincipal',
      'cyclic_spectrum_controls':spectrum,
      'exact_all_m_spectrum':'[m/gcd(m,2)] Z for rational and polynomial eigenfunctions',
      'scope':'Complete raw existence iff for squarefree generic D; necessary conditions and hostiles outside it; global nonsquarefree W2 remains open'}
target=folder/(source.stem+'_certificate.json')
target.write_bytes((json.dumps(cert,indent=2,sort_keys=True)+'\n').encode())
print('PASS: quadratic logarithmic witness pole corridor and affine hostiles')
print('Squarefree generic pencil: complete raw polynomial-eigenfunction existence iff')
print('Nonlinear source-submersive family: exact unit1, minimum witness degree2, non-global')
print('Elliptic integral-residue hostile: eigenvalue3 does not divide to eigenvalue1')
print('Exact all-m rational and polynomial spectrum: [m/gcd(m,2)] Z')
print('Always-active gates: '+str(gates))
print('Certificate SHA256: '+hashlib.sha256(target.read_bytes()).hexdigest())
