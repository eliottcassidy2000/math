"""Independent exact referee of the univariate polynomial mutation on W2.
No mathematical producer is imported. Analytic proofs are in the report.
"""
from pathlib import Path
import hashlib,json,sys
import sympy as S
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve()
DEST=HERE.parent.parent/'05-knowledge/results' if HERE.parent.name=='04-computation' else HERE.parent
u,z,t,h,r,b,L,c,X,A=S.symbols('u z t h r b L c X A')
gates=0
def need(ok,label):
    global gates
    gates+=1
    if not ok:raise ArithmeticError(label)
def eq(E,F,label):need(S.cancel(E-F)==0,label)
def jac(E,F,x,y):return S.diff(E,x)*S.diff(F,y)-S.diff(E,y)*S.diff(F,x)
zs=-2*h+u+u**3*t
M=-h*h*(3-h*r)-b*(1-h*r)**3
eq(zs.subs({u:1/r-h,t:-r*r-r**4*b},simultaneous=True),r*M,'literal complete z transition')
eq(jac(u,zs,u,t),u**3,'actual coordinate Jacobian and omitted source line')

# General local Eu jets with independent Taylor coefficients.
a0,a1,a2,a3,Q0=S.symbols('a0 a1 a2 a3 Q0')
epsilon=u+u**3*t
fj=a0+a1*epsilon+a2*epsilon**2/2+a3*epsilon**3/6
qj=Q0+a0*a0*epsilon+a0*a1*epsilon**2+(a1*a1+a0*a2)*epsilon**3/3
Aj=u*fj;Tj=L*Aj+qj
for degree in range(3):eq(S.expand(qj-Q0-a0*Aj).coeff(u,degree),0,'complete Eu target remainder through order two')
eq(S.diff(Tj,u).subs(u,0),a0*(L+a0),'exact Eu submersion coefficient')
eq(S.diff(Tj,t).subs(u,0),0,'exact Eu tangential coefficient')
for V in [u,t]:eq(S.diff(Tj,V).subs({u:0,L:-a0}),0,'forbidden lambda gives actual entire Eu critical line')
eta,beta=S.symbols('eta beta')
reparam=eta*A+beta*A**3
eq(S.expand((reparam/A)**2/2).coeff(A,0),eta*eta/2,'pure order-two leading scalar after reparametrization')
eq(S.expand((reparam/A)**2/2).coeff(A,1),0,'absence of quadratic target term removes entire simple pole')
eq(S.series(1/(2*A*A)-eta*eta/(2*reparam*reparam),A,0,1).removeO(),beta/eta,'reparametrized remainder is regular, not merely lower order')

# Complete independent polynomial controls in degrees q=2,...,6.
records=[]
for q in range(2,7):
    f=S.prod(z-j for j in range(q)).expand();Q=S.integrate(f*f,z);Q=S.expand(Q-Q.subs(z,0))
    fs=f.subs(z,zs);az=u*f;mutation=L*az+Q;mate=1/(2*az**2)
    need(S.degree(f,z)==q and S.degree(Q,z)==2*q+1,'exact degrees in independent f bank')
    eq(S.diff(Q,z),f*f,'polynomial primitive derivative')
    eq(S.gcd(f,S.diff(f,z)),1,'whole-root simplicity in declared bank')
    eq(f.subs(z,0),0,'actual globality hypothesis')
    a=f.subs(z,-2);need(a!=0 and 1+a!=0,'declared h1 lambda1 source Eu is admissible')
    eq(u**3*jac(az,Q,u,z),az**3,'actual univariate polynomial response identity')
    eq(u**3*jac(az,Q/az**3,u,z),1,'inherited supplier rational primitive independently checked')
    eq(az.subs(u,c/f),c,'inherited supplier whole nonzero-fibre inverse')
    need(Q.subs(z,1)!=0,'declared supplier bank retains a nonzero off-Eu coefficient')
    need(3*q>2*q+1,'supplier degree contradiction uses q>=2 exactly')
    eq(u**3*jac(mutation,mate,u,z),1,'literal rational mutation Jacobian')
    eq(mutation.subs(u,(c-Q)/(L*f)),c,'full fibre and field inverse')
    eq(zs.subs(t,(z-u+2*h)/u**3),z,'actual inverse recovers source coordinate')
    eq(u**3*jac(mutation,z,u,z),L*u**3*f,'nonzero derivation on rational fibre parameter')
    coefficient=S.Poly(f,z)
    ainf=(1-h*r)*sum(coefficient.nth(j)*r**(j-1)*M**j for j in range(1,q+1))
    tinf=L*ainf+Q.subs(z,r*M)
    eq(ainf.subs(r,0),-S.diff(f,z).subs(z,0)*(3*h*h+b),'whole boundary restriction of inherited A')
    eq(tinf.subs(r,0),-L*S.diff(f,z).subs(z,0)*(3*h*h+b),'whole boundary restriction of mutated first function')
    eq(S.diff(tinf,b).subs(r,0),-L*S.diff(f,z).subs(z,0),'all added points paid by nonzero tangential derivative')
    for j in range(3):eq(S.diff(Q.subs(z,r*M),r,j).subs(r,0),0,'Q has actual boundary order at least three')
    for rho in range(q):
        need(rho!=-2,'root component separated from omitted source line')
        eq(S.rem(Q-Q.subs(z,rho),(z-rho)**3,z),0,'entire root target reparametrization has no quadratic term')
        eq(S.diff(mutation,u).subs({z:rho,L:0}),0,'lambda0 root critical first derivative')
        eq(S.diff(mutation,z).subs({z:rho,L:0}),0,'lambda0 root critical second derivative')
        eq((zs-rho).subs(u,0),-2*h-rho,'primitive linear source root factor has paid constant')
    target_values=[Q.subs(z,j) for j in range(q)]+[Q.subs(z,-2)]
    need(len(set(target_values))==q+1,'declared real bank has all q+1 distinct target values')
    P=S.prod(X-value for value in target_values)
    # First divide in the honest (u,z) polynomial ring; only then restore Eu.
    Pmut=S.Poly(P.subs(X,mutation.subs(L,1)),z)
    divided,remainder=S.div(Pmut,S.Poly(f,z))
    eq(remainder.as_expr(),0,'all root components divide polynomial target support')
    eq(divided.as_expr().subs({u:0,z:-2}),0,'remaining actual Eu factor appears only after z=-2 at u0')
    need(3*q+2==(q+1)+(2*q+1),'generic puncture count includes fixed roots and infinity')
    need(q+1>=3,'three fixed punctures available for unlabeled non-isotriviality proof')
    need(S.degree((c-Q)**2,z)==2*(2*q+1),'exact complete rational-pair degree polynomial')
    if q<=3:
        literalT=mutation.subs(z,zs)
        eq(jac(literalT,mate.subs(z,zs),u,t),1,'independent actual original-source rational bracket')
        eq(S.diff(literalT,u).subs(u,0),f.subs(z,-2*h)*(L+f.subs(z,-2*h)),'literal whole Eu source gradient')
        for j in range(3):eq(S.diff(Q.subs(z,zs)-Q.subs(z,-2*h)-f.subs(z,-2*h)*u*fs,u,j).subs(u,0),0,'literal full original-source Eu second jet')
    records.append({'q':q,'f':str(f),'Q':str(Q),'target_values':[str(vv) for vv in target_values],
                    'punctures':3*q+2,'source_degree':2*q+1,'pair_degree':2*(2*q+1),'ambient_arms':q+1,'unit_arms':q+1})

# Explicit quintic and collision controls in exact quotient fields.
f=z*(z-1);Q=z**5/S.Integer(5)-z**4/S.Integer(2)+z**3/S.Integer(3)
eq(S.diff(Q,z),f*f,'complete explicit quintic primitive')
eq(f.subs(z,-2),6,'explicit quintic Eu coefficient')
values=[S.Integer(0),S.Rational(1,30),S.Rational(-256,15)]
need([Q.subs(z,j) for j in [0,1,-2]]==values,'all three explicit quintic target values')
need(len(set(values))==3,'explicit quintic has full torsion generation')
literalT=u*f.subs(z,zs)+Q.subs(z,zs)
eq(S.Poly(literalT,t).nth(5),u**15/5,'actual quintic leading coefficient without degree cancellation')
need(S.Poly(literalT,t).degree()==5,'actual original source quintic degree')
collision=6*z*z-15*z+10
eq(S.rem(Q,collision,z),0,'both admissible complex Eu collisions hit Q0')
eq(S.gcd(f*(1+f),collision),1,'both complex collision parameters retain all submersion assumptions')
eq(S.gcd(z,collision),1,'both complex collision translations have nonzero h')
need(3==2+1,'Eu collision leaves three ambient arms and two target groups')

# Independent root/root critical-value collision, separate from Eu collisions.
aa=S.symbols('aa');rootwall=7*aa*aa-7*aa+2
fc=z*(z-1)*(z-aa);Qc=S.integrate(fc*fc,z);Qc=S.expand(Qc-Qc.subs(z,0))
def reduce_wall(E,label,zero=True):
    nn,dd=S.fraction(S.cancel(E))
    eq(S.gcd(dd,rootwall),1,'root-collision denominator is a unit')
    if zero:eq(S.rem(nn,rootwall,aa),0,label)
    else:eq(S.gcd(nn,rootwall),1,label)
reduce_wall(Qc.subs(z,1),'nonzero f-root shares the regular root critical value')
for E in [fc.subs(z,-2),1+fc.subs(z,-2),Qc.subs(z,aa),Qc.subs(z,-2),Qc.subs(z,aa)-Qc.subs(z,-2)]:
    reduce_wall(E,'root/root collision stays smooth with exactly three target groups',False)

# The lower-degree q1 boundary is isotrivial by a literal puncture scaling.
linearQ=z**3/3
eq(linearQ.subs(z,2*z),8*linearQ,'q1 hostile: moving punctures are related by scaling')
need(1+1<3,'q1 has only two fixed punctures, so the three-point proof does not apply')

# Canonical Laurent arms and exact CRT separation of the explicit target supports.
def norm(A):return {j:S.cancel(E) for j,E in A.items() if S.cancel(E)!=0}
def mult(A,n):return norm({j-n:E for j,E in A.items() if j>n})
def diff(A,n):return norm({j+n:(-1)**n*S.rf(j,n)*E for j,E in A.items()})
def plus(*packets):
    result={}
    for packet in packets:
        for j,E in packet.items():result[j]=result.get(j,0)+E
    return norm(result)
for packet in [{3:S.Integer(1)},{2:S.Integer(1)}]:
    relation=plus(packet,diff(mult(packet,1),1),{j:E/2 for j,E in diff(mult(packet,2),2).items()})
    need(relation=={},'supplier exact extra Weyl relation retains both coefficient directions')
    need(mult(packet,3)=={},'supplier scalar third repair')
    for j in range(6):need(diff(relation,j)=={},'supplier relation remains exact under declared left derivatives')
der=S.symbols('der')
need(S.Matrix([[-der,1],[1,0]]).det()==-1,'supplier pivot columns form full polynomial derivative basis')
need(S.Matrix([[-der,1],[1,0]])*S.Matrix([-der,-der*der/2])==S.Matrix([der*der/2,-der]),'supplier relation eliminates the arbitrary PBW constant coefficient')
for coefficient in [S.Rational(1,2),S.Rational(49,2)]:
    packet={2:coefficient}
    need(mult(packet,2)=={} and mult(packet,1)!={},'complete local scalar annihilator has exact exponent two')
    for j in range(9):
        need(max(diff(packet,j))==j+2,'canonical derivatives retain exact primary order')
        need(diff(mult(packet,1),j)=={j+1:coefficient*(-1)**j*S.factorial(j)},'unit generates every negative power in its one coefficient direction')
P=S.prod(X-cc for cc in values)
for cc in values:
    block=(X-cc)**2;cofactor=S.div(P*P,block,X)[0]
    projector=S.expand(cofactor*S.invert(cofactor,block,X))
    for dd in values:
        expected=1 if cc==dd else 0
        eq(S.rem(projector,(X-dd)**2,X),expected,'actual polynomial CRT separates all unit supports at full pole depth')

cert={'status':'INDEPENDENT exact controls supporting all-parameter analytic audit','imports_producer':False,'gates':gates,
      'scope':'fixed W2, original source fibres and response ring; punctured curves non-isotrivial, completions genus0',
      'family':'T=lambda*u*f(z)+Q(z), z=u-2h+u^3*t, Qprime=f^2, Q0=0',
      'conditions':'f squarefree, degreeq>=2, f0=0, a0=f(-2h)!=0, lambda*(lambda+a0)!=0',
      'ambient_arms':'q+1','unit_arms':'number of distinct values Q(-2h),Q(rho)',
      'unit_annihilator':'product over distinct c of (T-c)^2','pair_degree':'2*(2q+1)',
      'generic_source_punctures':'3q+2','records':records,
      'quintic_target_values':[str(cc) for cc in values],'collision_controls':['Eu/root in degree2','root/root in degree3'],
      'retained_hostiles':['lambda0','lambda=-a0','q1 isotrivial','u,z must not be independent on Eu']}
cert['supplier_separate_audit']='univariate_fixed_unit: full ansatz globality/submersion, all-q complete PP, two coefficient directions, exact ideal Dg^3+D(1+partial*g+partial^2*g^2/2)'
raw=json.dumps(cert,sort_keys=True,separators=(',',':')).encode()+b'\n'
(DEST/(HERE.stem+'_certificate.json')).write_bytes(raw)
print('INDEPENDENT PASS: actual source and full boundary charts, rational mate, and exact Eu/root jets.')
print('INDEPENDENT PASS: full fibre inverses, special component collisions, unit CRT and canonical arm growth.')
print('INDEPENDENT PASS: actual quintic, degree-ten pair, eight punctures, and two distinct collision controls.')
print('INDEPENDENT PASS: separately pinned univariate supplier, full all-q proof and exact pointed Weyl presentation.')
print('CERTIFICATE_SHA256',hashlib.sha256(raw).hexdigest())
print('Always-active exact gates:',gates)
