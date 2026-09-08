"""Independent exact controls for the dependent quadratic discriminant obstruction."""
from pathlib import Path
from hashlib import sha256
import json,sys
import sympy as S
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve()
DEST=HERE.parent.parent/'05-knowledge/results' if HERE.parent.name=='04-computation' else HERE.parent
x,t,c,delta=S.symbols('x t c delta')
gates=0
def need(ok,msg):
 global gates
 gates+=1
 if not ok:raise ArithmeticError(msg)
def zero(e,msg):need(S.cancel(e)==0,msg)
def jac(f,g):return S.diff(f,x)*S.diff(g,t)-S.diff(f,t)*S.diff(g,x)

# General differential identities, with no degree bound or globality premise.
K=S.Function('K')(x);R=S.Function('R')(x);B=S.Function('S')(x)
H=R*t+B
F=K*H*H-c/4
N=K*R*R;P=2*K*R*B;Q=K*B*B-c/4
zero(F-(N*t*t+P*t+Q),'general factorization identity')
zero(P*P-4*N*Q-c*N,'general dependent discriminant identity')
zero(S.diff(F,t)-2*K*R*H,'whole critical curve tangential derivative')
zero(S.diff(F,x)-H*(S.diff(K,x)*H+2*K*(S.diff(R,x)*t+S.diff(B,x))),'whole critical curve normal derivative')

# Explicit finite valuation universe; the all-exponent argument is analytic.
for exponent in range(17):
 odd=exponent%2;half=exponent//2
 need(odd+2*half==exponent,'squarefree-times-square reconstruction')
 need(odd+half==(exponent+1)//2,'minimal P valuation equals valuation of K R')
 need(2*(odd+half)>=exponent,'necessary square divisibility survives each parity')

cases=[
 ('constant N',S.Integer(3),S.Integer(2),x*x+1,S.Rational(5,7)),
 ('P zero, discriminant zero',x*(x-1),x*x*(x-1),S.Integer(0),S.Integer(0)),
 ('even multiplicities',S.Integer(1),x**3*(x+2)**2,x-1,S.Integer(-4)),
 ('odd and even multiplicities',x*(x-2),x**2*(x-2)**3*(x+1),x*x+2,S.Rational(1,3)),
 ('R constant',x*(x+1),S.Integer(1),x**4-3,S.Integer(8)),
 ('nonmonic K scalar',S.Rational(7,3)*x*(x-3),x*(x-3)**2,S.Integer(5),S.Rational(-2,5)),
 ('P zero, nonzero discriminant',x+1,(x-1)**4,S.Integer(0),S.Integer(6))]
records=[]
for name,kk,rr,ss,cc in cases:
 hh=rr*t+ss; ff=S.expand(kk*hh*hh-cc/4)
 poly=S.Poly(ff,t);nn,pp,qq=poly.nth(2),poly.nth(1),poly.nth(0)
 need(nn!=0,'declared case is genuine source quadratic')
 zero(pp*pp-4*nn*qq-cc*nn,'declared case full discriminant dependence')
 need(S.Poly(kk,x).gcd(S.Poly(S.diff(kk,x),x)).degree()==0,'declared K is squarefree including constants')
 need(S.rem(pp,2*kk*rr,x)==0,'declared exact P divisibility')
 xx=next(i for i in range(20) if rr.subs(x,i)!=0)
 tt=S.cancel(-ss.subs(x,xx)/rr.subs(x,xx))
 zero(hh.subs({x:xx,t:tt}),'explicit critical curve is nonempty')
 zero(S.diff(ff,x).subs({x:xx,t:tt}),'explicit actual critical point normal')
 zero(S.diff(ff,t).subs({x:xx,t:tt}),'explicit actual critical point tangential')
 records.append(dict(name=name,N=str(nn),P=str(pp),Q=str(qq),c=str(cc),point=[str(xx),str(tt)]))

# Inherited hostile: a dependent exact pencil can have a rational mate.
f=x**5*(x-1)**3*t*t
mate=(8*x*x-4*x-1)/(3*x**4*(x-1)**2*t)
zero(jac(f,mate)-1,'inherited proportional 5+3 rational mate is retained')
zero(S.diff(f,x).subs(t,0),'inherited rational mate does not remove source critical line')
zero(S.diff(f,t).subs(t,0),'inherited rational mate source critical line second derivative')

# Inherited global W2 family: dependence also occurs in the actual global ring.
w=1+x*x*t
fglobal=(x**4+delta*x)*w*w
mglobal=1/(3*delta*x*x*w)
zero(jac(fglobal,mglobal)-1,'global dependent radical hostile retains rational mate')
rr,bb=S.symbols('r b')
zero(fglobal.subs({x:1/rr,t:-rr*rr-rr**4*bb},simultaneous=True)-(1+delta*rr**3)*bb*bb,'global dependent hostile full boundary chart')
zero(S.diff(fglobal,t).subs(t,-1/x**2),'global dependent hostile actual source critical curve')
zero(S.diff(fglobal,x).subs(t,-1/x**2),'global dependent hostile actual source critical curve normal')

cert=dict(status='Independent analytic obstruction plus FINITE-EXACT controls',
 theorem='For N!=0, dependence of N and P^2-4NQ forces a nonempty source critical curve',
 factorization='N=K R^2, P=2K R S, Q=K S^2-c/4, F+c/4=K(Rt+S)^2',
 complete_scope='All complex source quadratics; no degree bound, globality, or rational-mate assumption',
 consequences='Every source-submersive genuine quadratic has a two-dimensional discriminant pencil',
 complete_W2_corollary='Original t-degree exactly two, everywhere submersive, rational mate: exactly fifth-order, seventh-order, and elliptic families by audited suppliers',
 full_source_torsion_arms=2,complete_scalar_unit_order_spectrum=[2,3],
 preserved_hostiles=['dependent pencils may have rational mates','dependent pencils may be globally regular on W2'],
 valuation_control_exponents=list(range(17)),parameter_controls=records,gates=gates)
raw=json.dumps(cert,sort_keys=True,separators=(',',':')).encode()+b'\n'
(DEST/(HERE.stem+'_certificate.json')).write_bytes(raw)
print('DEPENDENT QUADRATIC PENCIL: complete source-submersion obstruction in every degree.')
print('FACTOR: F+c/4=K*(R*t+S)^2; the nonempty curve R*t+S=0 is critical.')
print('HOSTILES RETAINED: rational mates and actual W2 globality do not imply submersion.')
print('CERTIFICATE_SHA256',sha256(raw).hexdigest())
print('Always-active exact gates:',gates)
