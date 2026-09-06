"""Original shared-root charts and a genuinely admissible singular wall.
Exact symbolic identities and rational isolations; no producer imports.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb, floor, ceil
import sys, json, hashlib
import sympy as S
sys.stdout.reconfigure(newline='\n')
G=0

def need(ok,label):
    global G
    G+=1
    if not ok:raise ArithmeticError(label)

def I(a,b=None):return (F(a),F(a if b is None else b))
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[1],-a[0])
def mul(a,b):
    v=[x*y for x in a for y in b]
    return(min(v),max(v))
def inverse(a):
    need(a[0]*a[1]>0,'interval denominator excludes zero')
    return (1/a[1],1/a[0])
def power(a,k):
    if k<0:return power(inverse(a),-k)
    ans=I(1)
    for _ in range(k):ans=mul(ans,a)
    return ans
def rational(q):return F(int(S.numer(q)),int(S.denom(q)))
def poly_interval(expr,var,box):
    ans=I(0)
    for c in S.Poly(expr,var).all_coeffs():ans=add(mul(ans,box),I(rational(c)))
    return ans
def fraction_interval(expr,var,box):
    a,b=S.fraction(S.cancel(expr))
    return mul(poly_interval(a,var,box),inverse(poly_interval(b,var,box)))

r,u,x,y,z,v,s,t=S.symbols('r u x y z v s t')
B=v**5-13*v**4+55*v**3-x*v**2+y*v-z
C=v**4-12*v**3+45*v**2-S.Rational(2,3)*x*v+S.Rational(3,7)*y
D=v**4-11*v**3+36*v**2-S.Rational(5,12)*x*v+S.Rational(1,7)*y
phase=z*s**4-S.Rational(12,7)*y*s**3+x*s**2-10*s+S.Rational(1,11)
YC=S.Rational(14,9)*x*r-S.Rational(7,3)*r**2*(r**2-12*r+45)
ZC=r**2*(S.Rational(5,9)*x-r*(4*r**2-45*r+150)/3)
YD=S.Rational(35,12)*x*r-7*r**2*(r**2-11*r+36)
ZD=r**2*(S.Rational(23,12)*x-r*(6*r**2-64*r+197))
AC=5*u**2-24*u+9
AD=23*u**2-60*u+12
HC=44*r**2*u**4-132*r**2*u**3-495*r*u**4+1584*r*u**3-3*r+1650*u**4-5940*u**3+330*u
HD=66*r**2*u**4-132*r**2*u**3-704*r*u**4+1452*r*u**3-r+2167*u**4-4752*u**3+110*u
PC=705672*r**4-15079284*r**3+115450055*r**2-375696750*r+441317250
PD=120434688*r**4-2310536448*r**3+15142183583*r**2-39453584808*r+35022385200
UC=(11286*r**2-109085*r+245025)/(27126*r**2-261360*r+586025)
UD=(446688*r**2-4023865*r+7744176)/(1978416*r**2-17739216*r+34339426)

for name,J,Y,Z,A,H,k,den,P,U in [('C',C,YC,ZC,AC,HC,9,33,PC,UC),('D',D,YD,ZD,AD,HD,12,11,PD,UD)]:
    need(S.expand(B.subs({v:r,y:Y,z:Z}))==0 and S.expand(J.subs({v:r,y:Y}))==0,'exact shared-root chart '+name)
    g=S.expand(phase.subs({y:Y,z:Z,s:u/r}))
    need(S.factor(g-x*u**2*A/(k*r**2)+H/(den*r))==0,'original phase affine coefficient and constant '+name)
    need(S.rem(S.together(A.subs(u,U)).as_numer_denom()[0],P,r)==0,'singular resultant root implies coefficient zero '+name)
    need(S.rem(S.together(H.subs(u,U)).as_numer_denom()[0],P,r)==0,'singular resultant root implies phase constant zero '+name)
    need(S.degree(S.gcd(S.denom(U),P),r)==0,'no lost pole in the singular root reconstruction '+name)
    res=S.factor(S.resultant(A,H,u))
    need(S.degree(S.cancel(res/P),r)==0,'complete degree-four elimination '+name)

isolations={'C':[3528453826,3550824804,6081020695,8208387541],
            'D':[2049574992,2831681376,6130326664,8173391713]}
singular=[]
for name,P,U in [('C',PC,UC),('D',PD,UD)]:
    for a in isolations[name]:
        box=I(F(a,10**9),F(a+1,10**9))
        need(S.Poly(P,r).count_roots(S.Rational(a,10**9),S.Rational(a+1,10**9))==1,'unique exact singular root in rational interval')
        ub=fraction_interval(U,r,box)
        need(ub[0]>0,'reconstructed original phase has u>0')
        singular.append([name,str(box[0]),str(box[1]),str(ub[0]),str(ub[1])])
    need(S.Poly(P,r).count_roots(0,S.oo)==4,'all positive singular roots accounted for')

# One admissible wall; the remaining singular points are not declared admissible.
R=I(F(3528453826758,10**12),F(3528453826759,10**12))
need(S.Poly(PC,r).count_roots(S.Rational(R[0].numerator,R[0].denominator),S.Rational(R[1].numerator,R[1].denominator))==1,'refined unique algebraic witness')
X=S.Rational(3,5)*r*(4*r*r-45*r+150)+S.Rational(1,10)
Y=S.expand(YC.subs(x,X));Z=r*r/18
need(S.expand(ZC.subs(x,X)-Z)==0,'positive exact z=r²/18')
Bw=S.expand(B.subs({x:X,y:Y,z:Z}));Cw=S.expand(C.subs({x:X,y:Y}));Dw=S.expand(D.subs({x:X,y:Y}))
QB,remB=S.div(Bw,v-r,v);QC,remC=S.div(Cw,v-r,v)
need(remB==remC==0,'shared root removed without deleting any other factor')
root_boxes={
 'B':[(18,20),(669,671),(2450,2452),(6331,6333)],
 'C':[(388,390),(1951,1953),(6130,6132)],
 'D':[(181,183),(1493,1495),(3383,3385),(5939,5941)]}
bracket_records=[]
for name,poly in [('B',QB),('C',QC),('D',Dw)]:
    for a,b in root_boxes[name]:
        lo=poly_interval(poly.subs(v,S.Rational(a,1000)),r,R)
        hi=poly_interval(poly.subs(v,S.Rational(b,1000)),r,R)
        need((lo[1]<0<hi[0]) or(hi[1]<0<lo[0]),'uniform opposite endpoint signs for '+name)
        bracket_records.append([name,a,b,-1 if lo[1]<0 else 1])
need(R[0]>F(3385,1000) and R[1]<F(5939,1000),'shared root has the displayed middle ordering')
need(F(20,1000)<F(181,1000)<F(183,1000)<F(388,1000)<F(390,1000)<F(669,1000), 'first B gap contains both distinct interlacer roots')
need(F(671,1000)<F(1493,1000)<F(1495,1000)<F(1951,1000)<F(1953,1000)<F(2450,1000),'second B gap contains both interlacer roots')
need(F(2452,1000)<F(3383,1000)<F(3385,1000)<R[0], 'third B gap contains D root and ends at shared C root')
need(R[1]<F(5939,1000)<F(5941,1000)<F(6130,1000)<F(6132,1000)<F(6331,1000),'fourth B gap contains both final interlacer roots')
xb,yb,zb=poly_interval(X,r,R),poly_interval(Y,r,R),poly_interval(Z,r,R)
ub=fraction_interval(UC,r,R);sb=mul(ub,inverse(R))
need(F(86)<xb[0]<xb[1]<87 and 0<yb[0]<yb[1]<40 and 0<zb[0]<zb[1]<1,
     'necessary x/y/z caps do not exclude the admissible singular wall')

# Reconstruct the ORIGINAL Laurent response, including its inverse carry.
beta={-1:S.Integer(1),0:S.Integer(13),1:S.Integer(55),2:x,3:y,4:z}
cc={-1:S.Integer(1),0:S.Integer(12),1:S.Integer(45),2:S.Rational(2,3)*x,3:S.Rational(3,7)*y}
dd={-1:S.Integer(1),0:S.Integer(11),1:S.Integer(36),2:S.Rational(5,12)*x,3:S.Rational(1,7)*y}
raw={}
for a,c in beta.items():
    for b,d in beta.items():raw[a+b]=raw.get(a+b,0)+c*d
for a,c in cc.items():
    for b,d in dd.items():raw[a+b+1]=raw.get(a+b+1,0)+2*c*d
O={j:comb(14,2*j+1) for j in range(7)};E={j:comb(14,2*j) for j in range(8)}
weights={}
for a,c in O.items():
    for b,d in O.items():weights[a+b]=weights.get(a+b,0)+c*d
for a,c in E.items():
    for b,d in E.items():weights[a+b-1]=weights.get(a+b-1,0)+c*d
q={j:S.expand(c*weights.get(j,0)) for j,c in raw.items() if weights.get(j,0)}
need(sorted(q)==list(range(-1,9)) and q[-1]==28,'entire original response and nonzero lower carry')
need(S.expand(sum(beta[j]*O.get(j,0)*(-s)**j for j in beta)-2002*phase)==0,'original first zero, not a transformed first row')
total=I(0)
for j,c in q.items():
    coeff=I(0)
    for (a,b,d),v0 in S.Poly(c,x,y,z).terms():
        coeff=add(coeff,mul(I(rational(v0)),mul(power(xb,a),mul(power(yb,b),power(zb,d)))))
    term=mul(coeff,power(sb,j))
    total=add(total,neg(term) if j%2 else term)
need(total[1]<0,'full original response remains negative on the admissible singular witness')
coarse=[floor(total[0]),ceil(total[1])]
need(coarse[1]<0,'simple integer response enclosure')
record={'scope':'exact shared-root charts and one admissible singular C wall; general shared-root sign remains OPEN',
        'r_polynomial':str(PC),'r_interval':list(map(str,R)),'u_rational_function':str(UC),
        'x':str(X),'y':str(Y),'z':str(Z),'root_sign_brackets':bracket_records,
        'singular_pairs':singular,'coefficient_bounds':{k:list(map(str,b)) for k,b in [('x',xb),('y',yb),('z',zb),('u',ub),('s',sb)]},
        'original_response_enclosure':list(map(str,total)),'integer_response_enclosure':coarse}
here=Path(__file__).resolve().parent
dest=here.parent/'05-knowledge/results' if here.name=='04-computation' else here
out=dest/(Path(__file__).stem+'_certificate.json')
out.write_bytes((json.dumps(record,indent=2,sort_keys=True)+'\n').encode())
print('CHART_IDENTITIES C,D exact; each has four positive simultaneous coefficient/constant zeros')
print('ADMISSIBLE C singular wall r in',list(map(str,R)),'z=r^2/18; B simple positive; C weak and D strict interlacing')
print('ORIGINAL_PHASE s range',float(sb[0]),float(sb[1]),'original Q integer enclosure',coarse)
print('SCOPE denominator division loses a genuine model wall; this witness is negative, general shared-root sign OPEN')
print('CERTIFICATE_SHA256',hashlib.sha256(out.read_bytes()).hexdigest())
print('PASS',G,'always-active exact gates; actual LF')
