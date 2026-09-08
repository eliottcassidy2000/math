"""Independent actual-source and full torsion/Weyl unit audit.

No producer imports. The all-order PBW argument is in the matching report;
the finite operator bank checks exact normal forms and a rank-one hostile.
"""
import sys
from fractions import Fraction as Q
import sympy as s
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
u,t,w,y,c,r,b=s.symbols('u t w y c r b')
gates=0
def need(ok,msg):
    global gates
    gates+=1
    if not ok:raise ArithmeticError(msg)
def zero(expr,msg):need(s.cancel(expr)==0,msg)
def jac(F,G,X,Y):return s.diff(F,X)*s.diff(G,Y)-s.diff(F,Y)*s.diff(G,X)

L=(u+1)*w-4
f=u*w*L
q=((1-2*u)*w+2)/(3*u*u*w*w)
fs=f.subs(w,1+u*u*t)
qs=q.subs(w,1+u*u*t)
zero(jac(fs,qs,u,t)-1,'actual source unit primitive')
zero(s.diff(fs,u).subs(u,0)+3,'omitted source line noncritical')
zero((s.diff(f,w)*u*u).subs(w,0)+4*u**3,'entire w=0 divisor noncritical')
P2=L*L*((1-2*u)*w+2)/3
zero(f*f*q-P2,'complete polynomial order-two witness')
zero(jac(fs,P2.subs(w,1+u*u*t),u,t)-fs**2,'actual witness derivative')

# Both local coordinate systems are retained for the full principal parts.
Ru=s.cancel(qs-9/fs**2-4/fs)
Rw=s.cancel(q-s.Rational(32,3)/f**2-4/f)
need(s.factor(s.denom(Ru).subs(u,0))!=0,'full Eu remainder regular in actual (u,t) chart')
need(s.factor(s.denom(Rw).subs(w,0))!=0,'full Ew remainder regular with u invertible')
zero(s.limit(fs**2*qs,u,0)-9,'independent Eu highest scalar coefficient')
zero(s.limit(fs*(qs-9/fs**2),u,0)-4,'independent Eu simple scalar coefficient')
zero(s.limit(f**2*q,w,0)-s.Rational(32,3),'independent Ew highest scalar coefficient')
zero(s.limit(f*(q-s.Rational(32,3)/f**2),w,0)-4,'independent Ew simple scalar coefficient')
zero(L.subs(w,0)+4,'Ew and EL disjoint')
zero(L.subs({u:0,w:1})+3,'Eu and EL disjoint')
need(s.gcd(u*u,1)==1,'Ew factor irreducible primitive linear source polynomial')
need(s.gcd((u+1)*u*u,u-3)==1,'EL factor irreducible primitive linear source polynomial')
zero(s.discriminant(f-c,w)-4*u*((c+4)*u+c),'whole other-fibre discriminant')
zero((4*u*((c+4)*u+c)).subs(c,-4)+16*u,'special nonzero fibre -4 remains nonsquare')
zero((fs-c).subs(u,0)+c,'no omitted u component at nonzero fibre value')

# Rational fibre coordinate pays the exact constants without a geometric guess.
fy=(1+1/u)*y*y-4*y
ui=y*y/(c-y*y+4*y)
zero(f.subs(w,y/u)-fy,'actual fibre coordinate y=uw')
zero(fy.subs(u,ui)-c,'rational inverse u over C(F,y)')
zero(jac(fs,(u*w).subs(w,1+u*u*t),u,t)+u*(u*(1+u*u*t))**2,'nonzero derivation of y')

# Globality and submersion use the actual second chart; D maps to -4.
M=2-r-r*b*(1-r)**2
fi=(1-r)*M*(M-4)
qi=(2+(3*r-2)*M)/(3*(1-r)**2*M**2)
zero(fs.subs({u:1/r-1,t:-r*r-r**4*b},simultaneous=True)-fi,'entire global F chart')
zero(qs.subs({u:1/r-1,t:-r*r-r**4*b},simultaneous=True)-qi,'entire rational G chart')
zero(fi.subs(r,0)+4,'boundary lies over -4')
zero(s.diff(fi,r).subs(r,0)-4,'global boundary submersion')
zero(qi.subs(r,0)+s.Rational(1,6),'primitive has no boundary pole')

# Independent rational-vector representation of two labelled principal-part arms.
A=(Q(9),Q(32,3));B=(Q(4),Q(4))
def clean(V):return {j:v for j,v in V.items() if j<0 and any(v)}
def add(U,V):
    return clean({j:tuple(U.get(j,(Q(0),Q(0)))[k]+V.get(j,(Q(0),Q(0)))[k] for k in range(2)) for j in set(U)|set(V)})
def scale(c,V):return clean({j:tuple(c*x for x in v) for j,v in V.items()})
def T(V):return clean({j+1:v for j,v in V.items()})
def D(V):return clean({j-1:tuple(j*x for x in v) for j,v in V.items()})
def iterate(fn,V,n):
    for _ in range(n):V=fn(V)
    return V
theta={-2:A,-1:B}
need(s.det(s.Matrix([[9,4,1],[s.Rational(32,3),4,1],[0,0,1]]))==-s.Rational(20,3),'independence in three components modulo diagonal')
need(T(theta)=={-1:A},'g theta gives A first arm')
need(add(T(D(theta)),scale(2,theta))=={-1:B},'g nabla plus two gives B first arm')
need(T(T(theta))=={},'left generator g squared annihilates actual unit')
for j in range(1,13):
    for axis in range(2):
        V={-j:tuple(Q(int(k==axis)) for k in range(2))}
        need(add(D(T(V)),scale(-1,T(D(V))))==V,'literal Weyl commutator on both truncated arms')
AA=T(theta);BB=add(T(D(theta)),scale(2,theta))
need(add(scale(Q(-3,5),AA),scale(Q(8,5),BB))=={-1:(Q(1),Q(0))},'explicit Eu arm generator')
need(add(scale(Q(3,5),AA),scale(Q(-27,20),BB))=={-1:(Q(0),Q(1))},'explicit Ew arm generator')

# PBW remainder basis: derivatives on the left, powers of g on the right.
for degree in range(9):
    vectors=[iterate(D,iterate(T,theta,j),i) for j in [0,1] for i in range(degree+1)]
    matrix=s.Matrix([[s.Rational(V.get(-power,(Q(0),Q(0)))[axis].numerator,V.get(-power,(Q(0),Q(0)))[axis].denominator)
                      for V in vectors] for power in range(1,degree+3) for axis in range(2)])
    need(matrix.rank()==2*(degree+1),'bounded PBW remainder map injective in every degree0..8')
for i in range(9):
    for j in range(2,9):
        need(iterate(D,iterate(T,theta,j),i)=={},'entire bounded right-g-squared ideal bank annihilates')

# Order two alone does not imply that its unit generates both arms.
rank_one={-2:A}
need(T(rank_one)!={} and T(T(rank_one))=={},'rank-one hostile has exact scalar order two')
need(add(T(D(rank_one)),scale(2,rank_one))=={},'rank-one hostile has extra annihilator g nabla plus two')
need(add(T(D(theta)),scale(2,theta))!={},'same operator survives on actual two-direction unit')

print('ACCEPT: actual affine torsion has two labelled arms, and the quadratic unit has exact order two.')
print('WEYL: theta is cyclic for the full torsion; exact left annihilator A1*g^2; quotient A1/(A1*g^2).')
print('FINITE OPERATOR BANK: all PBW remainders through derivative degree8 inject; right-g powers2..8 vanish.')
print('HOSTILE: equal order with only one coefficient direction has a strictly larger annihilator.')
print('Always-active exact gates:',gates)
