"""Exact algebra/controls for the effective two-anchor tail theorem.

The universal inequalities are proved in the companion note. This program
checks their formal identities, every displayed rational bound, and one
strict positive control plus the rejected weak boundary. No producer imports.
"""
from fractions import Fraction as F
from math import comb
import sympy as S

checks=0
def gate(ok, label):
    global checks
    checks+=1
    if not ok: raise RuntimeError(label)

v,s,e3,e4,e5=S.symbols('v s e3 e4 e5')
B=v**5-13*v**4+55*v**3-e3*v**2+e4*v-e5
C=v**4-12*v**3+45*v**2-S.Rational(2,3)*e3*v+S.Rational(3,7)*e4
D=v**4-11*v**3+36*v**2-S.Rational(5,12)*e3*v+S.Rational(1,7)*e4
gate(S.expand(v*C-S.Rational(2,3)*B-v**3*(v-5)**2/3+S.Rational(5,21)*e4*v-S.Rational(2,3)*e5)==0,'C root identity')
gate(S.expand(v*D-S.Rational(5,12)*B-v**3*(7*v*v-67*v+157)/12+S.Rational(23,84)*e4*v-S.Rational(5,12)*e5)==0,'D root identity')
t,delta,ab=S.symbols('t delta ab')
d=5+delta
ce=55-13*t+t*t-ab-d*(13-t-d)
gate(S.expand(25-5*(13-t-d)+ce-(-3*t+2*delta+t*t+delta*t+delta*delta-ab))==0,'remaining pair identity')
eta=S.symbols('eta')
gate(S.expand(7*(5+eta)**2-67*(5+eta)+157-(-3+3*eta+7*eta**2))==0,'D near five')

epsilon=F(1,100)
gate(118<11**2,'third root lower bound')
gate(9**2>59,'fourth root lower bound')
gate(F(5,7)*epsilon<F(1,10)**2,'fourth root close to five')
f5_bound=3*F(3,100)+2*F(1,10)+F(3,100)**2+F(1,10)*F(3,100)+F(1,10)**2+F(3,100)**2/4
gate(f5_bound==F(2433,8000),'pair bound exact')
gate(f5_bound<F(9,25),'pair bound')
gate(F(9,25)/F(9,5)==F(1,5),'largest root close to five')
quadratic_upper=-3+3*F(1,5)+7*F(1,5)**2
gate(quadratic_upper==-F(53,25),'D quadratic bound')
last_upper=F(5,12)*F(3,2)*epsilon**2/F(24,5)
gate(last_upper==F(1,76800),'small fifth coefficient')
D_upper=F(24,5)**2*quadratic_upper/12+last_upper
gate(D_upper==-F(2544,625)+F(1,76800)<0,'D largest root contradiction')

def conv(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():out[i+j]=out.get(i+j,0)+x*y
    return {i:S.expand(x) for i,x in out.items()}
GB=[1,13,55,e3,e4,e5]
GC=[1,12,45,S.Rational(2,3)*e3,S.Rational(3,7)*e4]
GD=[1,11,36,S.Rational(5,12)*e3,S.Rational(1,7)*e4]
beta={i-1:x for i,x in enumerate(GB)}
cr={i-1:x for i,x in enumerate(GC)}
dr={i-1:x for i,x in enumerate(GD)}
mix=conv(beta,beta)
for j,x in conv(cr,dr).items():mix[j+1]=S.expand(mix.get(j+1,0)+2*x)
q={j:S.expand(comb(28,2*j+2)*mix[j]) for j in range(-1,9)}
gate(q[8]==comb(28,18)*e5**2,'top coefficient')
gate(S.expand(q[7]-comb(28,16)*(2*e4*e5+S.Rational(6,49)*e4**2))==0,'next coefficient')
for j,x in q.items():gate(all(c>=0 for c in S.Poly(x,e3,e4,e5).coeffs()),f'nonnegative raw coefficient {j}')
for j in range(29):gate(comb(28,j)<=comb(28,14),f'binomial maximum {j}')
gate(S.expand(sum(mix.values())-(sum(GB)**2+2*sum(GC)*sum(GD)))==0,'positive coefficient mass')
for cc in (GC,GD):
    gate(all(c>=0 for c in S.Poly(sum(GB)-sum(cc),e3,e4,e5).coeffs()),'contiguous mass bound')
P=182-20020*s+2002*e3*s*s-3432*e4*s**3+2002*e5*s**4
gate(S.expand(P/(S.Integer(2002)*s**3)-(e5*s-S.Rational(12,7)*e4+e3/s-10/s**2+S.Rational(1,11)/s**3))==0,'original root relation')
gate(F(12,7)*comb(28,18)-2*comb(28,16)==-38346750,'top cancellation coefficient')
beta_mass=F(18,5)**5
fifth_bound=F(13,5)**5
delta0=F(6,49)*comb(28,16)*epsilon**2
H=3*comb(28,14)*beta_mass**2+10*comb(28,18)*fifth_bound
threshold=118163898523
gate(delta0==F(2607579,7000),'delta0')
gate(H==F(17194291313831660508,390625),'H')
gate(H/delta0<threshold,'effective tail threshold')
gate(threshold>=1,'tail domain')

# Literal strict model control and exact weak boundary.
roots=[S.Rational(x,5) for x in (1,3,9,22,30)]
poly=S.Poly(S.prod(v-x for x in roots),v)
gate(poly.nth(4)==-13 and poly.nth(3)==55,'strict control anchors')
params={e3:-poly.nth(2),e4:poly.nth(1),e5:-poly.nth(0)}
gate(params[e4]>S.Rational(1,100),'strict control fourth coefficient')
for pp in (C,D):
    for j,x in enumerate(roots):gate((-1)**j*pp.subs(params).subs(v,x)>0,'strict control interlacing signs')
boundary={e3:75,e4:0,e5:0}
gate(S.expand(C.subs(boundary)-v*(v-2)*(v-5)**2)==0,'boundary C still interlaces')
gate(D.subs(boundary).subs(v,5)==-S.Rational(25,4),'boundary D rejects')
gate(S.expand(B.subs(boundary)-v**2*(v-3)*(v-5)**2)==0,'boundary B')

print('PROVED MODEL: both exact anchors and both weak interlacers force e4>1/100.')
print('Root identities, coefficient mass and top cancellation are formal polynomial identities.')
print('D endpoint upper bound:',D_upper)
print('delta0 =',delta0,'; H =',H)
print(f'Every original positive phase s >= {threshold} has Q(-s)<0.')
print('Finite-phase negativity and general actual noncancellation remain OPEN.')
print('PASS:',checks,'always-active exact gates; one strict control and one rejected weak boundary.')
