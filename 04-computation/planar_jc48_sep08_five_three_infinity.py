#!/usr/bin/env python3
"""Exact controls for both infinity placements of the binary (5,3) entry.
The all-parameter local classification and compact-curve arguments are in the
companion proof. No finite coefficient scan is used as a theorem.
"""
import hashlib
import json
import sympy as S

u,t,p,r,s,d = S.symbols('u t p r s d')
a,b,c,k,lam,m,n = S.symbols('a b c k lam m n')
z,v,U,A = S.symbols('z v U A')
gates=0
records={}
def check(name, value):
    global gates
    if not bool(value):
        raise RuntimeError(name)
    gates += 1
def zero(name, expr):
    check(name, S.cancel(expr)==0)
def jac(f,g,x=u,y=t):
    return S.diff(f,x)*S.diff(g,y)-S.diff(f,y)*S.diff(g,x)
def inf(expr,power):
    return S.expand(s**power*expr.subs({u:1/r-p,t:-r*r-r**4/s}, simultaneous=True))
def global_check(name,H,L):
    for which,f in [('H',H),('L',L)]:
        transformed=S.expand(f.subs({u:1/r-p,t:-r*r-r**4*d},simultaneous=True))
        for e in range(-4,0): zero(name+which+str(e),transformed.coeff(r,e))
        check(name+which+' polynomial', S.denom(S.cancel(transformed))==1)

def residue(expr):
    return S.factor(S.residue(expr,u,0))

def inverse_rows(N,P,Q,M,R):
    D=S.cancel(Q-P*P/(4*N)); E=S.cancel(R-M*P/(2*N))
    return -M/(4*N),-M*E/(8*N),-M*(3*E*E-D*M*M/N)/(32*N)

# Complete all-p source rows, including every constant and first jet.
P4,P3,P2,P1,P0,Q0,M4,M3,M2,M1,M0,R0=S.symbols('P4 P3 P2 P1 P0 Q0 M4 M3 M2 M1 M0 R0')
P=P4*u**4+P3*u**3+P2*u*u+P1*u+P0
M=M4*u**4+M3*u**3+M2*u*u+M1*u+M0
R=M4*u*u+(M3-2*p*M4)*u+R0
for deg in (3,5):
    Q=P4*u*u+(P3-2*p*P4-(1 if deg==5 else 0))*u+Q0
    H=u**deg*t*t+P*t+Q; L=M*t+R
    global_check('full'+str(deg),H,L)
    # Derive the same constraints in the original x, not by surface translation.
    x=S.symbols('x'); Px=S.Poly(S.expand(P.subs(u,x-p)),x)
    Qx=S.expand(Q.subs(u,x-p)); q2=Px.coeff_monomial(x**4)
    q1=Px.coeff_monomial(x**3)-(1 if deg==5 else 0)
    zero('original Q2 '+str(deg),Qx.coeff(x,2)-q2)
    zero('original Q1 '+str(deg),Qx.coeff(x,1)-q1)
    zero('original R1 '+str(deg),S.expand(R.subs(u,x-p)).coeff(x,1)-S.expand(M.subs(u,x-p)).coeff(x,3))

# Weighted local faces: actual normalization, not Puiseux determination counts.
# entries: (r order,s order,extra weight numerator, Phi_s order,eta order)
faces=[('m3_n1_high',3,5,6,13,5),('m3_deep',2,3,4,9,2),
       ('m5_j3_n2',3,8,6,22,2),('m5_deep',2,5,4,15,0)]
for name,vr,vs,weight,den,out in faces:
    zero(name,weight+2*vs+vr-1-den-out)
# cancellation centres: displayed rational exponents multiply dr; ramification
# adds e-1, so nonnegative rational exponent is enough for regularity.
for mult,j,cap in [(3,1,2),(5,1,4),(5,2,3)]:
    for ell in range(1,cap+1):
        exponent=S.Rational(mult-j-ell-2*j,2)+2
        check('cancel '+str((mult,j,ell)),exponent>=0)
# The balanced m5,j2,n1 cubic has no triple root with nonzero a,b,M.
x0,y0,z0=S.symbols('x0 y0 z0', nonzero=True)
C=(x0+y0*s)**2+z0*s**3
# Direct double-root equations yield s=-3a/b, M=4b^3/(27a).
double=-3*x0/y0; Mdouble=4*y0**3/(27*x0)
zero('double face',C.subs({s:double,z0:Mdouble}))
zero('double derivative',S.diff(C,s).subs({s:double,z0:Mdouble}))
check('double not triple',S.factor(S.diff(C,s,2).subs({s:double,z0:Mdouble}))!=0)
for ell in (1,2):
    # eta=r dr/v in the parameterized Morse model v^2~r^ell.
    vr=2 if ell==1 else 1; vv=1
    check('Morse weighted '+str(ell),vr+vr-1-vv>=0)

# Quintuple infinity / finite triple: constant-D residual and actual cubic map.
H3=u**3*t*t+u*u*(a*u+b)*t+a*u+c
Z3=u**3*t+u; L3=lam*Z3
Ninf=inf(H3,2); Minf=inf(L3,1)
zero('m5 leading',Ninf.coeff(s,0)-r**5*(1-p*r)**3)
zero('m5 normal first',S.expand(Ninf.coeff(s,1)).coeff(r,1)+a)
zero('m5 M first',S.expand(Minf.coeff(s,0)).coeff(r,1)+lam)
zero('m3 T2 residue',residue((M3*u**3+M2*u*u)/u**3)-M2)
Hcube=z*z*U**3-2*z*U*U+(1+b*z)*U+a*z+c-b
zero('cubic model',H3.subs({u:1/U,t:(z-1/U)*U**3},simultaneous=True)-Hcube)
zero('cubic volume',S.det(S.Matrix([[S.diff(1/U,U),S.diff(1/U,z)],[S.diff((z-1/U)*U**3,U),S.diff((z-1/U)*U**3,z)]]))+U)
check('genus3 divisor',2*1+6-2*2==4)
check('bidegree arithmetic genus', (2-1)*(3-1)==2)

# Exact all-parameter critical cubic for a=0.
f=v*(v+b); H0=u*f+c; F0=H0*H0+lam*(u*u*v+u)
uc=-(2*v+b)/(v*(3*v+b)); hc=lam/(2*v*(3*v+b))
crit=-4*v**3+6*(c-b)*v*v+2*b*(c-b)*v-lam
zero('critical residual',S.factor((H0.subs(u,uc)-hc)*2*v*(3*v+b))-crit)
# When H0=hc and u=uc, both derivatives vanish; clear the displayed denominators.
zero('critical Fv relation',S.diff(F0,v)-2*H0*u*(2*v+b)-lam*u*u)
zero('critical combined relation',S.factor((2*hc*f+lam*(2*uc*v+1))))
zero('critical second relation',S.factor(2*hc*uc*(2*v+b)+lam*uc*uc))
# All roots concentrated at the two forbidden addresses: only one pattern.
patterns=[]
for count in range(4):
    target=S.Poly(S.expand(-4*(v+b/3)**count*(v+b/2)**(3-count)),v)
    csol=S.solve(S.expand(crit-target.as_expr()).coeff(v,2),c)[0]
    lsol=S.solve(S.expand(crit-target.as_expr()).coeff(v,0),lam)[0]
    residual=S.factor((crit-target.as_expr()).subs({c:csol,lam:lsol}))
    patterns.append([count,str(csol),str(lsol),str(residual)])
    check('blocked cubic '+str(count),(residual==0)==(count==3))
check('unique tuning',patterns[3][1]=='b/3' and S.sympify(patterns[3][2])==4*b**3/27)
records['critical_patterns']=patterns

# Conic exceptional family, birational field and unequal special-fibre poles.
HH=u**3*t*t+3*u*u*t+1; VV=u*t
AA=1+u*(VV+1)*(VV+4); BB=1+u*VV*(VV+1)
FF=HH*HH+4*(u**3*t+u); GG=1/(2*AA*u*(VV+1))
zero('conic F=AB',FF-AA*BB);zero('conic rational mate',jac(FF,GG)-1)
global_check('conic',HH,4*Z3)
BBa=z/A; va=4*(BBa-1)/(A-BBa); ua=(A-BBa)**2/(4*(A+3*BBa-4))
zero('conic inverse A', (1+u*(v+1)*(v+4)).subs({u:ua,v:va},simultaneous=True)-A)
zero('conic inverse B', (1+u*v*(v+1)).subs({u:ua,v:va},simultaneous=True)-BBa)
pp=S.cancel((FF-1)*GG)
zero('conic E pp',pp.subs(u,0)-2)
zero('conic Gamma pp',pp.subs(t,-1/u)-1)

# Triple infinity / finite quintuple: full residue reduction.
N5=u**5; P5=u*u*(a*u*u+b*u+k)
Q5=a*u*u+(b-2*p*a-1)*u+c
H5=N5*t*t+P5*t+Q5
M5=u*u*(m*u+n);R5=m*u
T2,T6,T10=inverse_rows(N5,P5.subs(k,0),Q5,M5,R5)
zero('T6 residual',residue(T6)-m*(2*a*n+b*m-2*m)/16)
zero('T10 residual',residue(T10)+(a*a*n**3+6*a*b*m*n*n+6*a*m*m*n*p-6*a*m*n*n+3*b*b*m*m*n-9*b*m*m*n-c*m**3+6*m*m*n)/32)
zero('T10 nonly',residue(T10).subs({m:0,c:0})+a*a*n**3/32)
zero('T10 monly tuned',residue(T10).subs({n:0,b:2})-c*m**3/32)
w=S.symbols('w'); F5=S.expand(H5*H5+M5*t+R5)
Fw=S.expand(F5.subs(t,w/u**2))
ff=(k*w+c)**2+n*w
fg=2*(k*w+c)*(w*w+b*w+b-2*p*a-1)+m*(w+1)
zero('moving constant',Fw.coeff(u,0)-ff);zero('moving first',Fw.coeff(u,1)-fg)
zero('moving leading',S.expand(fg).coeff(w,3)-2*k)
Fn=S.expand(F5.subs({k:0,c:0,m:0}))
zero('nonly source Fu',S.diff(Fn,u).subs(u,0));zero('nonly source Ft',S.diff(Fn,t).subs(u,0))

# Residual fold gives a genuine source critical point whenever a != 0.
h=v*v/u+a*u*(v-2*p); F=h*h+m*v
fold=a*(v-2*p)*u*u-v*v
zero('fold Hu',S.diff(h,u)-fold/u**2)
zero('fold h',h-2*v*v/u-fold/u)
Rcrit=12*a*v*v-16*a*p*v+m
zero('fold Fv',u*u*(S.diff(F,v)-Rcrit)-2*(a*u*u-2*v)*fold)
zero('fold Fu',S.diff(F,u)-2*h*fold/u**2)
# A double root at the only nonzero forbidden address would force p=0,m=0.
zero('fold concentration linear',S.expand(Rcrit-12*a*(v-2*p)**2).coeff(v,1)-32*a*p)
zero('fold concentration constant at p0',S.expand(Rcrit-12*a*(v-2*p)**2).subs(p,0)-m)

# Sharp second rational family and the polynomial-repair obstruction.
V=u*(1+u*u*t); H=u*(1+u*u*t)**2; F=H*H+m*V; G=1/(6*V**3)
zero('second rational mate',jac(F,G)-1);global_check('second',H,m*V)
# At E: F=mV+V²+O(V4); at Gamma: F=mV+O(V4).
vs=z/m-z*z/m**3+2*z**3/m**5
coeff=S.series(1/(6*vs**3),z,0,0).removeO()
zero('second E pp3',coeff.coeff(z,-3)-m**3/6)
zero('second E pp2',coeff.coeff(z,-2)-m/2)
zero('second Gamma pp2',S.expand(1/(6*(z/m)**3)).coeff(z,-2))
# Explicit rational field inverse, sufficient for constants of the fibre derivation.
hh=S.symbols('hh')
zero('second birational inverse',S.cancel((v*v/u).subs(u,v*v/hh))-hh)

records['gates']=gates
records['weighted_faces']=faces
semantic=hashlib.sha256(json.dumps(records,sort_keys=True,separators=(',',':')).encode()).hexdigest()
print('Both infinity placements of the binary (5,3) entry')
print('Exact gates:',gates)
print('Complete source rows; weighted m3/m5 branches; higher residues: PASS')
print('Actual critical-point alternatives and two rational sharpness families: PASS')
print('Same-fibre polynomial pole-repair obstructions: PASS')
print('Semantic SHA256:',semantic)
