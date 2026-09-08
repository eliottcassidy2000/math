"""Exact controls for the fixed fifth-order quadratic pencil family.

No producer imports. The accompanying proof gives coefficient exhaustion,
all-point global submersion, exact unit order and both Weyl annihilator cases.
"""
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as Q
import json,sys
import sympy as S
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve()
DEST=HERE.parent.parent/'05-knowledge/results' if HERE.parent.name=='04-computation' else HERE.parent
u,t,w,y,c,r,b,a,d,h,k,shift,p,q=S.symbols('u t w y c r b a d h k shift p q')
gates=0
def need(ok,msg):
    global gates
    gates+=1
    if not ok:raise ArithmeticError(msg)
def zero(E,msg):need(S.cancel(E)==0,msg)
def jac(F,G,X,Y):return S.diff(F,X)*S.diff(G,Y)-S.diff(F,Y)*S.diff(G,X)

L=(a*u+d)*w+k
g=u*w*L
F=shift+g
G=(2*(d-2*a*u)*w-k)/(6*d*d*u*u*w*w)
fs=F.subs(w,1+u*u*t);gs=fs-shift;Qs=G.subs(w,1+u*u*t)
N=u**5*(a*u+d)
P=u**3*(2*a*u+2*d+k)
C=a*u*u+(d+k)*u+shift
zero(fs-(N*t*t+P*t+C),'complete source coefficient identity')
zero(P*P-4*N*(C-c)-u**5*((k*k+4*a*(c-shift))*u+4*d*(c-shift)),'entire fixed discriminant pencil')
zero(jac(fs,Qs,u,t)-1,'actual rational source mate')
poly=L*L*(2*(d-2*a*u)*w-k)/(6*d*d)
zero(g*g*G-poly,'actual polynomial upper annihilator')
zero(jac(fs,poly.subs(w,1+u*u*t),u,t)-gs*gs,'polynomial witness bracket')

# Complete coefficient matching after the proved degree and u-adic steps.
Ptrial=u**3*(p*u+q)
Ctrial=(p-a)*u*u+(q-d+(4*a-2*p)*h)*u+shift
boxC=S.expand(N+Ctrial*(u+h)**4-Ptrial*(u+h)**2)
need(S.Poly(boxC,u).degree()<=4,'full global numerator constant box')
need(S.Poly(Ptrial-2*Ctrial*(u+h)**2,u).degree()<=4,'full global numerator linear box')
Dtrial=S.expand(Ptrial**2-4*N*Ctrial)
zero(Dtrial.coeff(u,8)-(p-2*a)**2,'highest discriminant equation forces p=2a')
zero(Ctrial.subs({p:2*a,q:2*d+k})-C,'matched remaining source coefficient')
zero(S.det(S.Matrix([[a,d],[k*k-4*a*shift,-4*d*shift]]))+d*k*k,'full pencil rank criterion')

# The translated root is tested in the actual global chart.
M=2*h-r*(h*h+b*(1-h*r)**2)
Y=(1-h*r)*M
Fi=shift+a*Y*Y+d*r*(1-h*r)*M*M+k*Y
zero(fs.subs({u:1/r-h,t:-r*r-r**4*b},simultaneous=True)-Fi,'complete second-chart identity')
zero(Fi.subs(r,0)-(shift+4*a*h*h+2*k*h),'entire boundary value')
zero(S.diff(Fi,r).subs(r,0)+(4*a*h+k)*(3*h*h+b)-4*d*h*h,'all boundary normal derivatives')
zero(S.diff(Fi,b).subs(r,0),'boundary tangential derivative')
zero(S.diff(Fi,r).subs({r:0,k:-4*a*h})-4*d*h*h,'global submersion under exact matching')
zero(S.diff(fs,u).subs(u,0)-(d+k),'omitted source line derivative')
zero((S.diff(g,w)*u*u).subs(w,0)-k*u**3,'w-zero source derivative')
Fuy=shift+(a+d/u)*y*y+k*y
zero(F.subs(w,y/u)-Fuy,'rational source coordinates off u-zero')
zero(S.diff(Fuy,u)+d*y*y/u**2,'all open-chart critical points force y-zero')
zero(S.diff(Fuy,y).subs(y,0)-k,'remaining source derivative at y-zero')
ui=d*y*y/(c-shift-a*y*y-k*y)
zero(Fuy.subs(u,ui)-c,'full rational constant-field inverse')
zero(jac(fs,(u*w).subs(w,1+u*u*t),u,t)+d*u*(u*(1+u*u*t))**2,'actual derivation on rational fibre coordinate')
Gfy=(c-shift)/(3*d*d*y**3)-a/(d*d*y)-k/(2*d*d*y*y)
zero(G-Gfy.subs({c:F,y:u*w},simultaneous=True),'primitive agrees in full source field')

# Every source component and the entire nonzero-fibre discriminant.
zero(L.subs({u:0,w:1})-(d+k),'Eu and EL separation value')
zero(L.subs(w,0)-k,'Ew and EL separation value')
zero(S.discriminant(g-c,w)-u*((k*k+4*a*c)*u+4*d*c),'whole nonzero-fibre discriminant')
zero((fs-shift-c).subs(u,0)+c,'no extra u component at nonzero fibre')

Aeu=(2*d-k)*(d+k)**2/(6*d*d)
Aew=-k**3/(6*d*d)
Bc=-a*k/(d*d)
Ru=S.cancel(Qs-Aeu/gs**2-Bc/gs)
Rw=S.cancel(G-Aew/g**2-Bc/g)
need(S.factor(S.denom(Ru).subs(u,0))!=0,'complete Eu remainder regular under d(d+k) nonzero')
need(S.factor(S.denom(Rw).subs(w,0))!=0,'complete Ew remainder regular under dku nonzero')
zero(S.limit(gs*gs*Qs,u,0)-Aeu,'independent Eu highest coefficient')
zero(S.limit(gs*(Qs-Aeu/gs**2),u,0)-Bc,'independent Eu simple coefficient')
zero(S.limit(g*g*G,w,0)-Aew,'independent Ew highest coefficient')
zero(S.limit(g*(G-Aew/g**2),w,0)-Bc,'independent Ew simple coefficient')
zero((Aeu-Aew)*Bc+a*k*(2*d+3*k)/(6*d*d),'exact two-direction determinant')
zero((2*d+3*k).subs(k,-4*a*h)-2*(d-6*a*h),'actual global cyclicity wall')
zero((Aeu-Aew).subs(d,-3*k/2),'wall leading directions coincide')
zero((Bc/Aew).subs(d,-3*k/2)-6*a/(k*k),'wall annihilator coefficient')
zero(Aeu.subs(d,k/2),'admitted Eu pole-order drop boundary')

# Independent exact Laurent-module controls at off-wall, wall and pole-drop.
def clean(V):return {j:v for j,v in V.items() if j<0 and any(v)}
def add(U,V):return clean({j:tuple(U.get(j,(Q(0),Q(0)))[i]+V.get(j,(Q(0),Q(0)))[i] for i in range(2)) for j in set(U)|set(V)})
def scale(c,V):return clean({j:tuple(c*x for x in v) for j,v in V.items()})
def T(V):return clean({j+1:v for j,v in V.items()})
def D(V):return clean({j-1:tuple(j*x for x in v) for j,v in V.items()})
def it(op,V,n):
    for _ in range(n):V=op(V)
    return V

# General coefficient-recovery lemma: Euler interpolation retains each arm.
for order in range(1,7):
    packet={-j:(Q(j),Q(j*j)) for j in range(1,order+1)}
    for chosen in range(1,order+1):
        recovered=packet
        for j in range(1,order+1):
            if j!=chosen:
                recovered=scale(Q(1,j-chosen),add(T(D(recovered)),scale(j,recovered)))
        need(recovered=={-chosen:packet[-chosen]},'Euler interpolation separates complete finite coefficient packet')
        need(it(T,recovered,chosen-1)=={-1:packet[-chosen]},'recovered coefficient reaches its first arm level')
records=[]
for av,hv,dv in [(1,1,1),(1,1,6),(1,1,-2),(1,1,2),(2,3,1),(2,3,36),(2,3,-12)]:
    kv=-4*av*hv;sub={a:av,h:hv,d:dv,k:kv}
    Avec=tuple(Q(S.Rational(E.subs(sub))) for E in [Aeu,Aew])
    BB=Q(S.Rational(Bc.subs(sub)));Bvec=(BB,BB)
    theta={-2:Avec,-1:Bvec};wall=dv==6*av*hv
    need(T(theta)!={} and T(T(theta))=={},'unit order exactly two at every declared control')
    need(T(theta)=={-1:Avec},'first leading arm extraction')
    need(add(T(D(theta)),scale(2,theta))=={-1:Bvec},'first lower arm extraction')
    det=Avec[0]*BB-Avec[1]*BB
    need((det==0)==wall,'complete declared wall/off-wall control')
    need(dv*(dv+kv)*kv!=0,'all declared controls pay source submersion')
    for deg in range(7):
        vectors=[it(D,it(T,theta,j),i) for j in [0,1] for i in range(deg+1)]
        Mtx=S.Matrix([[S.Rational(V.get(-power,(Q(0),Q(0)))[axis]) for V in vectors] for power in range(1,deg+3) for axis in range(2)])
        need(Mtx.rank()==(deg+2 if wall else 2*(deg+1)),'whole bounded PBW remainder rank distinguishes wall')
    if wall:
        lam=Q(6*av,kv*kv)
        need(add(add(T(D(theta)),scale(2,theta)),scale(-lam,T(theta)))=={},'exact additional wall annihilator')
        need(Avec[0]==Avec[1] and Avec[0]!=0 and BB!=0,'wall retains exact order-two nonzero unit')
    if dv==-2*av*hv:
        need(Avec[0]==0 and Avec[1]!=0 and BB!=0,'pole-drop remains off-wall and order two')
    if dv==2*av*hv:
        need(S.gcd(N.subs(sub),P.subs(sub))==u**4,'allowed gcd multiplicity jump does not change fibre primitivity')
    records.append(dict(a=av,h=hv,d=dv,k=kv,wall=wall,leading=[str(x) for x in Avec],simple=str(BB)))

# Cheap hostile: arbitrary pencil member is global but need not be submersive.
bad={a:1,d:1,h:1,k:1}
zero(S.diff(Fi,r).subs(bad).subs({r:0,b:-S.Rational(11,5)}),'within-pencil off-matching boundary critical point is real')
zero(S.diff(fs,u).subs({a:1,d:4,h:1,k:-4,u:0}),'forbidden d=4ah is actually source critical')

cert=dict(status='FINITE-EXACT controls; all-parameter theorem in companion report',
    scope='Fixed W2; full discriminant pencil (x-h)^5 span{1,x-h}; a*d*h*(d-4*a*h) nonzero and k=-4*a*h; original source response ring C[x,t]',
    rational_mate=str(G),leading_Eu=str(Aeu),leading_Ew=str(Aew),simple_shared=str(Bc),
    cyclicity_determinant='-a*k*(2*d+3*k)/(6*d^2)',wall='d=6*a*h',
    offwall_left_annihilator='A1*g^2',wall_left_annihilator='A1*g^2 + A1*(g*nabla+2-(6*a/k^2)*g)',
    controls=records,gates=gates)
raw=json.dumps(cert,sort_keys=True,separators=(',',':')).encode()+b'\n'
(DEST/(HERE.stem+'_certificate.json')).write_bytes(raw)
print('FIXED PENCIL: complete global coefficient matching and exact all-point submersion criterion.')
print('UNIT: order two throughout a*d*h*(d-4*a*h)!=0, k=-4*a*h; all three source components retained.')
print('WEYL: cyclic for both complete arms iff d!=6*a*h; wall retains one full arm with the stated two-generator annihilator.')
print('HOSTILE: d=-2*a*h drops only the Eu leading pole; unit order and off-wall cyclicity survive.')
print('CERTIFICATE_SHA256',sha256(raw).hexdigest())
print('Always-active exact gates:',gates)
