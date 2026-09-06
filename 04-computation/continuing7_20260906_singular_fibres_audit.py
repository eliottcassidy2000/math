"""Independent singular-fibre audit: Laurent coefficients and exact Sturm signs."""
from pathlib import Path
from math import comb
import json,sys,hashlib
import sympy as S
sys.stdout.reconfigure(newline='\n')
gates=0
def need(ok,label):
    global gates
    gates+=1
    if not ok:raise ArithmeticError(label)
here=Path(__file__).resolve().parent
root=here.parent if here.name=='04-computation' else Path('C:/w/s0905')
base=root/'05-knowledge/results' if here.name=='04-computation' else Path('C:/w/continuing7_20260906_moments')
cp=base/'continuing7_20260906_singular_fibres_certificate.json'
raw=cp.read_bytes();need(hashlib.sha256(raw).hexdigest()=='34f178da766757e31968103e47543d71ab8e975100a8ac0d1f4dd1d69c027ee9','pinned producer certificate')
cert=json.loads(raw)
r,u,x,y,z,s,v,t,w=S.symbols('r u x y z s v t w')
B=v**5-13*v**4+55*v**3-x*v**2+y*v-z
C=v**4-12*v**3+45*v**2-S.Rational(2,3)*x*v+S.Rational(3,7)*y
D=v**4-11*v**3+36*v**2-S.Rational(5,12)*x*v+S.Rational(1,7)*y
beta=1/t+13+55*t+x*t**2+y*t**3+z*t**4
cc=1/t+12+45*t+S.Rational(2,3)*x*t**2+S.Rational(3,7)*y*t**3
dd=1/t+11+36*t+S.Rational(5,12)*x*t**2+S.Rational(1,7)*y*t**3
odd=sum(comb(14,2*j+1)*t**j for j in range(7))
even=sum(comb(14,2*j)*t**j for j in range(8))
conv=S.Poly(S.expand(t**2*(beta**2+2*t*cc*dd)),t)
weights=S.Poly(S.expand(t*(odd**2+even**2/t)),t)
Q=0
for j in range(-1,9):
    coefficient=S.expand(conv.nth(j+2)*weights.nth(j+1))
    need(S.expand(coefficient-S.sympify(cert['original_q'][str(j)]))==0,'original Laurent response coefficient')
    Q+=coefficient*(-s)**j
phase=sum(S.expand(beta).coeff(t,j)*comb(14,2*j+1)*(-s)**j for j in range(5))/2002
need(S.expand(phase-(z*s**4-S.Rational(12,7)*y*s**3+x*s**2-10*s+S.Rational(1,11)))==0,'original phase')

def reduce(expr,P):
    n,d=S.fraction(S.cancel(expr));need(S.gcd(d,P)==1,'coefficient denominator unit')
    out=S.rem(n*S.invert(d,P,r),P,r)
    need(S.rem(n-d*out,P,r)==0,'quotient identity')
    return S.expand(out)
def sign(poly,lo,hi,desired,label):
    poly=S.Poly(poly,r)
    need(poly.count_roots(lo,hi)==0,label+' no roots on complete isolator')
    need(S.sign(poly.eval((lo+hi)/2))==desired,label+' rational midpoint sign')
def coeffwire(a):return [str(v) for v in S.Poly(a,r).all_coeffs()]
mom={}
for label,top in [('C',C),('D',D)]:
    # Formal rational expansion at infinity, independent of producer recurrence.
    series=S.series((top/B).subs(v,1/w)/w,w,0,7).removeO().expand()
    mm=[S.expand(series.coeff(w,k)) for k in range(7)]
    need(all(S.expand(a-S.sympify(b))==0 for a,b in zip(mm,cert[label+'_moments'])),'formal Stieltjes moment series')
    mom[label]=mm
def det(mm,n,shift=0):return S.expand(S.Matrix(n,n,lambda i,j:mm[i+j+shift]).det())
H4=det(mom['C'],4)
need(S.expand(H4-S.sympify(cert['C_ordinary_four']))==0,'degree-six Hankel consumer')
cap=S.Rational(7,72)*(x-75)*(135-x)
need(S.expand(det(mom['C'],3)-S.Rational(8,7)*(cap-y))==0,'necessary curved packet')
need(det(mom['C'],2,1)==(x-75)/3,'necessary x lower bound')
data={}
for label,top in [('C',C),('D',D)]:
    ca=cert['charts'][label];P=S.sympify(ca['P']);U=S.sympify(ca['U'])
    solutions=S.solve([B.subs(v,r),top.subs(v,r)],(y,z),dict=True)
    need(len(solutions)==1,'unique shared-root coefficient chart')
    Y,Z=solutions[0][y],solutions[0][z]
    need(S.expand(Y-S.sympify(ca['Y']))==S.expand(Z-S.sympify(ca['Z']))==0,'independent chart solve')
    phasechart=S.Poly(S.factor(phase.subs({y:Y,z:Z,s:u/r})),x)
    coeff=S.factor(phasechart.coeff_monomial(x));constant=S.factor(phasechart.coeff_monomial(1))
    aa=S.Poly(S.fraction(S.cancel(coeff/u**2))[0],u).monic().as_expr()
    hh=S.fraction(constant)[0]
    linear=S.rem(hh,aa,u);need(S.degree(linear,u)==1,'phase wall has linear residual')
    L=S.expand(linear).coeff(u,1);K=S.expand(linear).coeff(u,0)
    need(S.gcd(L,P)==1 and reduce(-K/L-U,P)==0,'unique wall u over each quartic root')
    need(reduce(aa.subs(u,U),P)==reduce(hh.subs(u,U),P)==0,'all reconstructed pairs solve exact wall')
    resultant=S.factor(S.resultant(aa,hh,u))
    need(S.degree(S.cancel(resultant/P),r)==0,'resultant has no omitted positive roots')
    need(S.Poly(P,r).count_roots(0,S.oo)==4,'four positive wall roots exhaust resultant')
    abc=[reduce(a.subs(u,U),P) for a in S.Poly(S.expand(Q.subs({y:Y,z:Z,s:u/r})),x).all_coeffs()]
    need(len(abc)==3 and [coeffwire(a) for a in abc]==ca['response_coefficients_descending_x'],'independent original response reduction')
    data[label]=(P,U,Y,Z,abc)
    for j in range(1,5):
        lo,hi=map(S.Rational,cert['fibres'][label+str(j)]['r_interval'])
        need(S.Poly(P,r).count_roots(lo,hi)==1,'one exact root per producer isolator')
        un,ud=S.fraction(U);sgn=S.sign(ud.subs(r,(lo+hi)/2))
        sign(ud,lo,hi,sgn,'wall denominator');sign(un,lo,hi,sgn,'positive u')
        if j==4:need(lo>S.Rational(71,10),'largest wall violates beta root bound')
need(S.Rational(71,10)**2+(13-S.Rational(71,10))**2/4>59,'Cauchy root cutoff')
for key,limit in {'C1':94,'C3':89,'D1':102,'D2':98,'D3':95}.items():
    P,U,Y,Z,abc=data[key[0]];lo,hi=map(S.Rational,cert['fibres'][key]['r_interval'])
    gap=S.expand(cap-Y)
    sign(gap.subs(x,limit),lo,hi,-1,key+' cap at cutoff')
    sign(S.diff(gap,x).subs(x,limit),lo,hi,-1,key+' cap slope')
    need(S.diff(gap,x,2)==-S.Rational(7,36),'whole excluded x ray by concavity')
    if key=='D3':
        sign(Z.subs(x,limit),lo,hi,-1,'D3 z at cutoff');sign(S.diff(Z,x),lo,hi,1,'D3 z slope')
P,U,Y,Z,abc=data['C'];lo,hi=map(S.Rational,cert['fibres']['C2']['r_interval'])
co=[S.rem(c,P,r) for c in S.Poly(S.expand(H4.subs({y:Y,z:Z}).subs(x,x+88)),x).all_coeffs()]
need(len(co)==5,'all five shifted determinant coefficients')
for c in co:sign(c,lo,hi,-1,'degree-six exclusion beyond88')
for key,limit in {'C1':94,'C2':88,'C3':89,'D1':102,'D2':98}.items():
    P,U,Y,Z,(a,b,c)=data[key[0]];lo,hi=map(S.Rational,cert['fibres'][key]['r_interval'])
    sign(a,lo,hi,-1,key+' curvature')
    sign(2*a*limit+b,lo,hi,1,key+' derivative at cutoff')
    sign(a*limit**2+b*limit+c,lo,hi,-1,key+' response at cutoff')
P,U,Y,Z,(a,b,c)=data['C'];lo,hi=map(S.Rational,cert['fibres']['C2']['r_interval'])
subs={x:92,y:Y.subs(x,92),z:Z.subs(x,92)}
for label,mm in mom.items():
    for shift in [0,1]:
        for n in [1,2,3]:sign(S.rem(det(mm,n,shift).subs(subs),P,r),lo,hi,1,'hostile '+label+' leading minor')
sign(a*92**2+b*92+c,lo,hi,1,'hostile positive original response')
sign(S.rem(H4.subs(subs),P,r),lo,hi,-1,'hostile negative degree-six determinant')

# Independent two-variable Bernstein enclosure for the actual interval.
alpha,beta2=S.symbols('alpha beta2')
def bernstein_sign(expr,rlo,rhi,xlo,xhi):
    poly=S.Poly(S.expand(expr.subs({r:rlo+(rhi-rlo)*alpha,x:xlo+(xhi-xlo)*beta2})),alpha,beta2)
    nr,nx=poly.degree(alpha),poly.degree(beta2)
    values=[]
    for i in range(nr+1):
        for j in range(nx+1):
            values.append(sum(poly.coeff_monomial(alpha**k*beta2**l)*S.Rational(comb(i,k),comb(nr,k))*S.Rational(comb(j,l),comb(nx,l)) for k in range(i+1) for l in range(j+1)))
    need(all(v>0 for v in values) or all(v<0 for v in values),'uniform Bernstein rectangle sign')
    return S.sign(values[0])
control=cert['controls']['C2_actual_full_model_interval'];xlo,xhi=map(S.Rational,control['x_interval'])
BB=S.expand(B.subs({y:Y,z:Z}));CC=S.expand(C.subs(y,Y));DD=S.expand(D.subs(y,Y))
BQ,br=S.div(BB,v-r,v);CQ,cr=S.div(CC,v-r,v)
need(br==cr==0,'known shared root factored exactly')
polys={'B':BQ,'C':CQ,'D':DD}
for label,aa,bb,*_ in control['brackets']:
    sa=bernstein_sign(polys[label].subs(v,S.Rational(aa,1000)),lo,hi,xlo,xhi)
    sb=bernstein_sign(polys[label].subs(v,S.Rational(bb,1000)),lo,hi,xlo,xhi)
    need(sa==-sb,'actual roots at opposite-sign endpoints')
need(bernstein_sign(Y,lo,hi,xlo,xhi)==bernstein_sign(Z,lo,hi,xlo,xhi)==1,'actual interval positive y,z')
# Degree exhaustion plus disjoint bracket orders independently gives interlacing.
brackets={label:sorted((S.Rational(aa,1000),S.Rational(bb,1000)) for lab,aa,bb,*_ in control['brackets'] if lab==label) for label in polys}
betas=brackets['B'][:3]+[(lo,hi)]+brackets['B'][3:]
croots=brackets['C'][:2]+[(lo,hi)]+brackets['C'][2:]
for i,inter in enumerate(croots):
    need(betas[i][1]<=inter[0] or inter==(lo,hi),'C lower interlace boundary')
    need(inter[1]<=betas[i+1][0] or inter==(lo,hi),'C upper interlace boundary')
for i,inter in enumerate(brackets['D']):need(betas[i][1]<inter[0]<inter[1]<betas[i+1][0],'strict D interlace')
print('Independent reconstruction: original phase, all10 response coefficients, formal moments through6')
print('Complete wall classification: eight positive pairs, five entire negative fibres, three excluded')
print('Exact Sturm signs and Bernstein actual-interval controls passed')
print('PASS',gates,'always-active audit gates; actual LF')
