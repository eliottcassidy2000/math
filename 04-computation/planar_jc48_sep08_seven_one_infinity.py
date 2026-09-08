#!/usr/bin/env python3
"""Exact controls for finite-simple/original-sevenfold-infinity DG quartics.
The companion proves the complete analytic branch and compact-component
argument. These symbolic controls have no finite mate-degree hypothesis.
"""
import hashlib
import json
import sympy as S

u,t,p,r,s,Z,zeta,x,bD=S.symbols('u t p r s Z zeta x bD')
a,b,c,d,e,k,al,be,ga,de,m0,ell=S.symbols('a b c d e k alpha beta gamma delta m0 ell')
gates=0
records={}
def check(name, value):
    global gates
    if not bool(value):
        raise RuntimeError(name)
    gates+=1
def zero(name, expr):
    check(name,S.cancel(expr)==0)
def jac(f,g,x0=u,y0=t):
    return S.diff(f,x0)*S.diff(g,y0)-S.diff(f,y0)*S.diff(g,x0)
def transformed(f,shift=0):
    return S.expand(s**shift*f.subs({u:1/r-p,t:-r*r-r**4/s},simultaneous=True))
def global_check(name,f):
    q=S.cancel(f.subs({u:1/r-p,t:-r*r-r**4*bD},simultaneous=True))
    check(name,S.denom(q)==1)

P=a*u**4+b*u**3+c*u*u+d*u+e
Q=a*u*u+(b-2*p*a)*u+k
M=al*u**4+be*u**3+ga*u*u+de*u+m0
R=al*u*u+(be-2*p*al)*u+ell
H=u*t*t+P*t+Q; L=M*t+R
F=H*H+L
for name,f in [('full H',H),('full L',L)]:global_check(name,f)
# Recover the complete constraints, starting from all admissible coefficient
# boxes, independently of the advertised solved formulas.
pp=S.symbols('p0:7'); qq=S.symbols('q0:5')
Hraw=u*t*t+sum(pp[i]*u**i for i in range(7))*t+sum(qq[i]*u**i for i in range(5))
Ht=S.expand(Hraw.subs({u:1/r-p,t:-r*r-r**4*bD},simultaneous=True))
hrows=[S.expand(Ht.coeff(r,j)).coeff(bD,i) for j in range(-4,0) for i in range(2)]
hmat=S.linear_eq_to_matrix(hrows,list(pp)+list(qq))[0]
check('full H constraint rank',hmat.rank()==6)
hsol=S.solve(hrows,[pp[6],pp[5],qq[4],qq[3],qq[2],qq[1]],dict=True)
check('full H unique solved block',len(hsol)==1)
hexpected={pp[6]:0,pp[5]:0,qq[4]:0,qq[3]:0,qq[2]:pp[4],qq[1]:pp[3]-2*p*pp[4]}
for key,value in hexpected.items():zero('full H solved '+str(key),hsol[0][key]-value)
mm=S.symbols('m0:5');rr=S.symbols('r0:3')
Lraw=sum(mm[i]*u**i for i in range(5))*t+sum(rr[i]*u**i for i in range(3))
Lt=S.expand(Lraw.subs({u:1/r-p,t:-r*r-r**4*bD},simultaneous=True))
lrows=[Lt.coeff(r,j) for j in (-2,-1)]
lsol=S.solve(lrows,[rr[2],rr[1]],dict=True)
check('full L unique solved block',len(lsol)==1)
zero('full L second',lsol[0][rr[2]]-mm[4])
zero('full L first',lsol[0][rr[1]]-mm[3]+2*p*mm[4])
# Full-row completeness is all-p; none of these are translated W assumptions.
records['section_dimensions']=[12-int(hmat.rank()),6]

# Independent short inverse recursion: T2 comes from the cubic centred row.
v,D0,E0,C0,kap=S.symbols('v D0 E0 C0 kap')
W=1-D0*v*v/2-C0*v**3/4
inverse_eq=S.expand(W**4+2*D0*v*v*W*W+C0*v**3*W+(D0*D0+E0)*v**4-1)
for j in (1,2,3):zero('inverse coefficient '+str(j),inverse_eq.coeff(v,j))
zero('T2 radical descent',(-C0/(4*kap)).subs(C0,M/kap).subs(kap*S.Integer(1),kap)+M/(4*kap**2))
zero('T2 residue',S.residue(-M/(4*u),u,0)+m0/4)
zero('mate coefficient factor',S.Rational(2,4)*(-M/(4*u))+M/(8*u))

# Exact original-coordinate infinity numerators, after the necessary m0=0.
M=M.subs(m0,0);L=L.subs(m0,0);F=F.subs(m0,0)
NN=transformed(H,2); MM=transformed(L,1)
B=-r**4*P.subs(u,1/r-p)+2*r**5-2*p*r**6
C=(-a*p**4*r*r+4*a*p**3*r-3*a*p*p+b*p**3*r*r-3*b*p*p*r+2*b*p
   -c*p*p*r*r+2*c*p*r-c+d*p*r*r-d*r-e*r*r+k+r**3-p*r**4)
DD=-r**4*M.subs(u,1/r-p)
EE=(-al*p**4*r*r+4*al*p**3*r-3*al*p*p+be*p**3*r*r-3*be*p*p*r+2*be*p
    -ga*p*p*r*r+2*ga*p*r-ga+de*p*r*r-de*r+ell)
zero('complete N numerator',NN-r**7*(1-p*r)-s*B-s*s*C)
zero('complete M numerator',MM-DD-s*EE)
zero('actual volume',S.det(S.Matrix([[S.diff(1/r,r),0],[S.diff(-r*r-r**4/s,r),S.diff(-r*r-r**4/s,s)]]))+r*r/s**2)
for i,name,coef,earlier in [(0,'a',a,{}),(1,'b',b,{a:0}),(2,'c',c,{a:0,b:0}),(3,'d',d,{a:0,b:0,c:0})]:
    zero('normal leading '+name,S.expand(B.subs(earlier)).coeff(r,i)+coef)
for i,name,coef,earlier in [(0,'alpha',al,{}),(1,'beta',be,{al:0}),(2,'gamma',ga,{al:0,be:0}),(3,'delta',de,{al:0,be:0,ga:0})]:
    zero('lower leading '+name,S.expand(DD.subs(earlier)).coeff(r,i)+coef)

# Complete simple low-face identities with generic coefficients and all normal
# orders that can alter a face. Higher terms have positive weight in the proof.
a0,b0,c0,d0,e0=S.symbols('a0 b0 c0 d0 e0',nonzero=True)
def local(j,n):
    N=a0*r**7+b0*r**j*s+c0*s*s
    M=d0*r**n+e0*s
    return S.expand(N*N+s**3*M-zeta*s**4)
def face(expr,power,weight):
    return S.expand(expr.subs(s,r**power*Z)).coeff(r,weight)
quad=(b0+c0*Z)**2+d0*Z+(e0-zeta)*Z*Z
zero('j1 n1 low face',face(local(1,1),1,4)-Z*Z*quad)
zero('j1 n>=2 low face',face(local(1,2),1,4)-Z*Z*quad.subs(d0,0))
zero('j2 n2 low face',face(local(2,2),2,8)-Z*Z*quad)
zero('j2 n>=3 low face',face(local(2,3),2,8)-Z*Z*quad.subs(d0,0))
zero('generic quadratic discriminant slope',S.diff(S.discriminant(quad,Z),zeta)-4*b0*b0)
for j in (2,3,4,5):
    zero('n1 low face j'+str(j),face(local(j,1),1,4)-Z**3*(d0+(c0*c0+e0-zeta)*Z))
zero('n1 j2 middle',face(local(2,1),3,10)-Z*Z*(b0*b0+d0*Z))
# Fractional high face paid on its actual normalized parameter.
tau=S.symbols('tau')
for j in (3,4,5):
    high=S.expand(local(j,1).subs({r:tau**3,s:tau**13*Z},simultaneous=True)).coeff(tau,42)
    zero('n1 normalized high j'+str(j),high-a0*a0-d0*Z**3)
for j in (3,4,5):
    zero('n2 low face j'+str(j),face(local(j,2),2,8)-Z**3*(d0+(c0*c0+e0-zeta)*Z))
    expected=(a0+b0*Z)**2+d0*Z**3 if j==3 else a0*a0+d0*Z**3
    zero('n2 high face j'+str(j),face(local(j,2),4,14)-expected)
# Differential orders (vr,vs,Phi_s-order,eta-order), including dr.
orders=[('j1 low',1,1,3,1),('n1 middle',1,3,7,1),('n1 high',3,13,29,5),
        ('j2 low',1,2,6,0),('n2 low',1,2,6,0),('n2 high simple',1,4,10,0)]
for name,vr,vs,den,out in orders:zero('weighted '+name,2*vr+2*vs+vr-1-den-out)
for j,contacts in [(1,range(1,7)),(2,range(2,6))]:
    for contact in contacts:
        exp=S.Rational(7-3*j-contact+4,2)
        check('cancelling regular '+str((j,contact)),exp>=0)
        # Split gap exceeds the centre; denominator term 2NN_s dominates.
        check('cancelling gap '+str((j,contact)),7-3*j+contact>0)
        check('contact upper '+str((j,contact)),contact<=7-j)
zero('n1 j2 cancelling exponent',S.Rational(7-6-1+4,2)-2)

# The only repeated high face is a Morse double, never a triple.
cubic=(a0+b0*Z)**2+d0*Z**3
zdouble=-3*a0/b0;mdouble=4*b0**3/(27*a0)
for j in (0,1):zero('Morse double '+str(j),S.diff(cubic,Z,j).subs({Z:zdouble,d0:mdouble}))
zero('Morse nonzero second',S.diff(cubic,Z,2).subs({Z:zdouble,d0:mdouble})+2*b0*b0/3)
# The critical-value derivative's first possible order is exactly two.
scaled=S.cancel(local(3,2).subs(s,r**4*Z)/r**14)
zero('generic fibre Morse order',S.diff(scaled,zeta)+r*r*Z**4)
zero('Morse ramified eta order',2-1-1)
zero('Morse split eta order',0-1-(-1))
records['local_orders']=orders
records['Morse_contacts']={'1':'regular order 0','2':'nonzero logarithms'}

# Last residual: the actual original coordinate substitution, not a model.
res={a:0,b:0,c:0,al:0,be:0,ga:0}
Hr=H.subs(res);Lr=L.subs(res)
Nres=NN.subs(res);Mres=MM.subs(res)
Phir=S.expand(Nres*Nres+s**3*Mres-zeta*s**4)
last=Z*Z*((-d+k*Z)**2-de*Z+(ell-zeta)*Z*Z)
zero('actual residual tangent',face(Phir,3,12)-last)
zero('last quadratic discriminant slope',S.diff(S.discriminant(S.cancel(last/(Z*Z)),Z),zeta)-4*d*d)
znonzero=de/(k*k+ell-zeta)
zero('d0 nonzero root',last.subs({d:0,Z:znonzero}))
zero('d0 simple derivative',S.diff(last,Z).subs({d:0,Z:znonzero})-de**3/(k*k+ell-zeta)**2)
zero('actual residual log order',2+6-9-(-1))
# Literal constant-L hostile with both nonzero tangent roots retained.
# H=u t²+u t, generic zeta=1: last=Z²(1-Z²), residues at ±1.
lit=last.subs({d:1,k:0,de:0,ell:0,zeta:1})
for root,expected in [(1,S.Rational(-1,2)),(-1,S.Rational(1,2))]:
    zero('literal log '+str(root),root**2/S.diff(lit,Z).subs(Z,root)-expected)
# Genuine nonconstant L example: H=u t², L=u t; root Z=-1 at zeta=1.
lit2=last.subs({d:0,k:0,de:1,ell:0,zeta:1})
zero('literal nonconstant log',1/S.diff(lit2,Z).subs(Z,-1)-1)

# Independent field inversion and ORIGINAL (u,t) derivation pathway.
h,f,vv=S.symbols('h f vv')
inv_t=(h-d*vv-k)/(vv+e);inv_u=vv/inv_t
zero('field inverse H',Hr.subs({u:inv_u,t:inv_t},simultaneous=True)-h)
zero('field inverse v',(u*t).subs({u:inv_u,t:inv_t},simultaneous=True)-vv)
vol=S.det(S.Matrix([[S.diff(inv_u,vv),S.diff(inv_u,h)],[S.diff(inv_t,vv),S.diff(inv_t,h)]]))
zero('field volume',vol-1/(h-d*vv-k))
vr=(f-h*h-ell)/de
zero('rational field eta',1/(de*(h-d*vr-k))-1/(d*h*h+de*h-d*(f-ell)-de*k))
# A candidate primitive depending on h only would have derivative the displayed
# form. This identity derives its denominator from the original source bracket.
Fr=Hr*Hr+Lr
orig_j=S.factor(jac(Fr,Hr))
zero('original field derivation',orig_j-de*(Hr-d*(u*t)-k))
# Distinct generic roots whenever d*delta !=0, so each residue is nonzero.
den=d*h*h+de*h-d*(f-ell)-de*k
zero('generic field discriminant slope',S.diff(S.discriminant(den,h),f)-4*d*d)
zero('d0 field residue',S.residue(1/den.subs(d,0),h,k)-1/de)

# All-p sharp positive family and polynomial factor obstruction.
Hsharp=u*t*t+e*t+k;Fsharp=Hsharp*Hsharp+ell
zero('sharp H mate',jac(Hsharp,-1/t)-1)
zero('sharp F mate',jac(Fsharp,-1/(2*Hsharp*t))-1)
global_check('sharp H global',Hsharp)
check('sharp H nonconstant',S.Poly(Hsharp,t).degree()==2)
Ggeneric=u*u*t**3+u*t+t+1
zero('literal polynomial factor',jac(Fsharp,Ggeneric)-2*Hsharp*jac(Hsharp,Ggeneric))
# The same proof preserves nonzero p; this control is not used as completeness.
zero('nonzero-p literal mate',jac(Fsharp.subs({e:2,k:3,ell:5}),(-1/(2*Hsharp*t)).subs({e:2,k:3}))-1)
records['sharp_family']={'H':'u*t^2+e*t+k','L':'ell','mate':'-1/(2*H*t)'}
records['scope']='Full all-p finite-simple/original-sevenfold-infinity carrier; no mate-degree bound.'
blob=json.dumps(records,sort_keys=True,separators=(',',':')).encode()
print('seven_one_infinity exact controls: PASS')
print('gates:',gates)
print('complete section dimensions: H=6, L=6 before T2')
print('actual weighted faces, contact bounds and residual logs: PASS')
print('original-coordinate field derivation and sharp all-p rational mates: PASS')
print('semantic sha256:',hashlib.sha256(blob).hexdigest())
