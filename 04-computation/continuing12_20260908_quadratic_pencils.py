"""Exact controls for three complete translated global W2 quadratic pencils.

Analytic exhaustion and all-parameter proofs are in the paired report.
No mathematical producer is imported; all gates survive optimization.
"""
from pathlib import Path
from hashlib import sha256
import json,sys
import sympy as S
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve()
DEST=HERE.parent.parent/'05-knowledge/results' if HERE.parent.name=='04-computation' else HERE.parent
u,t,z,w,x,a,d,h,k,s,c,p,q,R,r,b,Y=S.symbols('u t z w x a d h k s c p q R r b Y')
gates=0
def need(ok,msg):
 global gates
 gates+=1
 if not ok:raise ArithmeticError(msg)
def zero(E,msg):need(S.cancel(E)==0,msg)
def jac(F,G,X,T):return S.diff(F,X)*S.diff(G,T)-S.diff(F,T)*S.diff(G,X)

# Whole u^3 span{1,u}: complete global coefficient equations force a critical line.
N3=u**3*(a*u+d)
P3trial=u*u*(p*u*u+q*u+R)
Q3trial=p*(u+h)**2+(q-4*p*h)*(u+h)+s
D3trial=S.expand(P3trial**2-4*N3*Q3trial)
zero(D3trial.coeff(u,8)-p*p,'third-order first forced top coefficient')
zero(D3trial.subs(p,0).coeff(u,6)-q*q,'third-order second forced top coefficient')
F3=N3*t*t+R*u*u*t+s
zero(S.diff(F3,u).subs(u,0),'third-order entire source critical line normal derivative')
zero(S.diff(F3,t).subs(u,0),'third-order entire source critical line tangential derivative')
zero((R*u*u)**2-4*N3*s-u**3*((R*R-4*a*s)*u-4*d*s),'third-order full pencil')

# Whole u^4 span{1,u^2}: its two exhaustive branches fail in different charts.
N4=u**4*(a*u*u+d)
P4trial=u*u*(p*u*u+q*u+R)
Q4trial=(p-a)*u*u+(q+(4*a-2*p)*h)*u+s
zero(S.expand(P4trial**2-4*N4*Q4trial).coeff(u,8)-(p-2*a)**2,'fourth-order top match')
P4=u*u*(2*a*u*u+q*u+R);Q4=a*u*u+q*u+s
D4=S.expand((P4*P4-4*N4*Q4)/u**4)
zero(D4.coeff(u,1)-2*q*(R-2*d),'fourth-order exact branch equation')
zero(a*D4.coeff(u,0)-d*D4.coeff(u,2)-(a*(R-2*d)**2-d*q*q),'fourth-order full pencil rank')
F4a=s-d-k+(a*u*u+d)*w*w+k*w
F4b=s-d+(a*u*u+d)*w*w+q*u*w
zero((N4*t*t+P4*t+Q4).subs({q:0,R:2*d+k})-F4a.subs(w,1+u*u*t),'fourth-order branch A full matching')
zero((N4*t*t+P4*t+Q4).subs(R,2*d)-F4b.subs(w,1+u*u*t),'fourth-order branch B full matching')
zero(S.diff(F4a.subs(w,1+u*u*t),u).subs(u,0),'fourth-order branch A source critical')
zero(S.diff(F4a.subs(w,1+u*u*t),t).subs(u,0),'fourth-order branch A source critical second derivative')
M=2*h-r*(h*h+b*(1-h*r)**2)
Fi4=s-d+a*(1-h*r)**2*M*M+d*r*r*M*M+q*(1-h*r)*M
zero(F4b.subs({u:1/r-h,w:r*M},simultaneous=True)-Fi4,'fourth-order full boundary chart')
zero(S.diff(Fi4,b).subs(r,0),'fourth-order boundary tangential derivative')
zero(S.diff(Fi4,r).subs(r,0)+(4*a*h+q)*(3*h*h+b),'fourth-order exact boundary critical point')
zero(S.diff(Fi4,r).subs({r:0,b:-3*h*h}),'fourth-order actual boundary point always critical')

# Whole u^7 span{1,u}: complete global matching and full field primitive.
N=u**7*(a*u+d)
P=u**4*(2*a*u*u+(2*d-4*a*h)*u+k-4*d*h)
Q=u*(u-2*h)*(a*u*(u-2*h)+d*(u+2*h)+k-4*d*h)+s
Nxx=S.Poly(S.expand(N.subs(u,x-h)),x);Pxx=S.Poly(S.expand(P.subs(u,x-h)),x)
Qxx=S.Poly(S.expand(Q.subs(u,x-h)),x)
zero(Pxx.nth(6)-2*Nxx.nth(8),'seventh-order first mandatory global top coefficient')
zero(Pxx.nth(5)-2*Nxx.nth(7),'seventh-order second mandatory global top coefficient')
for j,E in [(4,Nxx.nth(8)),(3,Nxx.nth(7)),(2,Pxx.nth(4)-Nxx.nth(6)),(1,Pxx.nth(3)-Nxx.nth(5))]:
 zero(Qxx.nth(j)-E,'seventh-order complete global lower matching')
zero(P*P-4*N*(Q-c)-u**7*((k*k+4*a*(c-s))*u+4*d*(c-s)),'seventh-order entire discriminant pencil')
L=(a*u+d)*z+k
f=u*z*L
zs=u-2*h+u**3*t
fs=s+f.subs(z,zs)
zero(fs-(N*t*t+P*t+Q),'seventh-order complete source factorization')
num=2*(8*a*a*u*u-4*a*d*u+3*d*d)*z*z+k*(7*a*u-3*d)*z+k*k
G=num/(30*d**3*u**3*z**3)
Wit=L**3*num/(30*d**3)
zero(f**3*G-Wit,'seventh-order polynomial cubic witness')
zero(jac(fs,G.subs(z,zs),u,t)-1,'seventh-order literal rational source mate')
zero(jac(fs,Wit.subs(z,zs),u,t)-(fs-s)**3,'seventh-order literal polynomial response witness')
Ri=h*h*(3-h*r)+b*(1-h*r)**3
Fi=s+(1-h*r)*Ri*((a*(1-h*r)+d*r)*Ri-k)
zero(fs.subs({u:1/r-h,t:-r*r-r**4*b},simultaneous=True)-Fi,'seventh-order whole second-chart identity')
zero(Fi.subs(r,0)-(s+a*(b+3*h*h)**2-k*(b+3*h*h)),'seventh-order full boundary value')
zero(S.diff(Fi,b).subs(r,0)-(2*a*(b+3*h*h)-k),'seventh-order full tangential derivative')
normal=S.diff(Fi,r).subs(r,0)
zero(normal.subs(k,2*a*(b+3*h*h))-d*(b+3*h*h)**2,'boundary normal derivative on its exact tangential-zero locus')
zero(S.diff(fs,u).subs(u,0)+2*h*(k-2*d*h),'all omitted affine points paid')
Fuy=s+(a+d/u)*Y*Y+k*Y
zero(f.subs(z,Y/u)+s-Fuy,'actual open source coordinates')
zero(S.diff(Fuy,u)+d*Y*Y/u**2,'open critical point forces Y0')
zero(S.diff(Fuy,Y).subs(Y,0)-k,'remaining open derivative nonzero')
zero(jac(u,u*zs,u,t)-u**4,'actual localized coordinate Jacobian')
zero(Fuy.subs(u,d*Y*Y/(c-s-a*Y*Y-k*Y))-c,'rational constant-field inverse')
zero(jac(fs,u*zs,u,t)+d*u*u*(u*zs)**2,'actual derivation in rational coordinate')
Gfy=(c-s)**2/(5*d**3*Y**5)-(c-s)*k/(2*d**3*Y**4)+(k*k-2*a*(c-s))/(3*d**3*Y**3)+a*k/(d**3*Y**2)+a*a/(d**3*Y)
zero(G-Gfy.subs({c:s+f,Y:u*z},simultaneous=True),'single rational primitive in full field')
zero(zs.subs(u,0)+2*h,'Eu Ez disjointness')
zero(L.subs({u:0,z:-2*h})-(k-2*d*h),'Eu EL disjointness')
zero(L.subs(z,0)-k,'Ez EL disjointness')
zero((fs-s-c).subs(u,0)+c,'nonzero source fibres have no lost vertical factor')

# Complete scalar parts via polynomial remainders, with actual Eu jet retained.
A3=(k-2*d*h)**3*(k*k+6*d*h*k+24*d*d*h*h)/(30*d**3)
A2=(k-2*d*h)*(4*a*d*d*h*h+2*a*d*h*k+a*k*k-6*d**3*h)/(3*d**3)
A1=a*a*k/d**3
B3=k**5/(30*d**3);B2=a*k**3/(3*d**3);B1=A1
for variable,ff,pp,coeff in [(u,f.subs(z,zs),Wit.subs(z,zs),(A3,A2,A1)),(z,f,Wit,(B3,B2,B1))]:
 rem=pp-coeff[0]-coeff[1]*ff-coeff[2]*ff**2
 for j in range(3):zero(S.diff(rem,variable,j).subs(variable,0),'complete scalar remainder divisible by local cubic')
zero(A1-B1,'shared simple scalar coefficient')
zero(A2-B2+2*h*(4*a*h*h-6*d*h+3*k)/3,'second coefficient alignment equation')
zero(A3-B3+8*h**3*(12*d*d*h*h-15*d*h*k+5*k*k)/15,'third coefficient alignment equation')
zero(S.discriminant(5*k*k-15*d*h*k+12*d*d*h*h,k)+15*d*d*h*h,'no real cyclicity wall')
need(S.factor((A3*B2-B3*A2).subs(a,0))!=0,'a0 still has two independent coefficient directions')
wallpoly=5*q*q-15*q+12
def reducewall(E):
 num,den=S.fraction(S.cancel(E.subs({d:1,h:1,k:q,a:3*(2-q)/4})))
 return S.rem(num,wallpoly,q),S.rem(den,wallpoly,q)
for E in [A3-B3,A2-B2,A1-B1]:need(reducewall(E)[0]==0,'both exact complex wall branches align')
for E in [a,d,h,k,k-2*d*h,B3]:
 nn,dd=reducewall(E);need(nn!=0 and dd!=0,'complex wall remains in admissible smooth locus')

# Exact order-three PBW compiler, tested as a polynomial operator identity.
T=S.symbols('T');l1,l2,l3=S.symbols('l1 l2 l3')
p0=l1;p1=l2+T*l1;p2=l3+T*l2+T*T*l1/2
zero(p0-l1,'PBW relation first coordinate')
zero(p1-T*p0-l2,'PBW relation second coordinate')
zero(p2-T*p1+T*T*p0/2-l3,'PBW relation third coordinate')
controls=[]
for aa,dd,hh,kk in [(0,1,1,1),(1,1,1,1),(2,1,-1,3),(1,1,1,4)]:
 sub={a:aa,d:dd,h:hh,k:kk}
 C=S.Matrix([[A1.subs(sub),A2.subs(sub),A3.subs(sub)],[B1.subs(sub),B2.subs(sub),B3.subs(sub)]])
 need(dd*hh*kk*(kk-2*dd*hh)!=0,'declared control globally smooth')
 need(C.rank()==2,'declared real controls generate both complete arms')
 need(C[:,2]!=S.zeros(2,1),'declared control retains exact cubic unit order')
 rel=C.nullspace();need(len(rel)==1 and C*rel[0]==S.zeros(2,1),'complete coefficient relation used by PBW compiler')
 controls.append(dict(a=aa,d=dd,h=hh,k=kk,rank=2))

cert=dict(status='FINITE-EXACT controls supporting separate all-parameter proofs',
 closed_pencils=['u^3 span{1,u}: entire source critical line','u^4 span{1,u^2}: source or boundary critical'],
 surviving_pencil='u^7 span{1,u}',global_submersion='d*k*h*(k-2*d*h)!=0; a arbitrary',
 source_components=['u=0','z=u-2h+u^3*t=0','(a*u+d)*z+k=0'],unit_order=3,full_torsion_arms=2,
 rank_one_wall=['4*a*h^2-6*d*h+3*k=0','12*d^2*h^2-15*d*h*k+5*k^2=0'],
 rational_mate=str(G),source_unit_coefficients=[str(A1),str(A2),str(A3)],other_polar_coefficients=[str(B1),str(B2),str(B3)],
 controls=controls,gates=gates,scope='fixed original W2 source ring; other exact pencils and general JC OPEN')
raw=json.dumps(cert,sort_keys=True,separators=(',',':')).encode()+b'\n'
(DEST/(HERE.stem+'_certificate.json')).write_bytes(raw)
print('WHOLE PENCILS: u3 span(1,u) and u4 span(1,u2) have no globally submersive member.')
print('SEVENTH PENCIL: complete global-submersion family; exact source unit order3 and exactly two full torsion arms.')
print('COEFFICIENT RANK: two arms except two admissible complex rank-one branches; no real rank-one wall.')
print('PBW: exact order-three relation compiler retains all coefficient vectors; no Dg3-only claim.')
print('CERTIFICATE_SHA256',sha256(raw).hexdigest())
print('Always-active exact gates:',gates)
