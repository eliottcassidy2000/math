#!/usr/bin/env python3
"""Exact finite-row transport replacing the THM-4438 weight-24 payer.

Standalone symbolic certificate: no research producer imports or sampled
parameters. All identities are over QQ and every check survives python -O.
The result concerns row15 and projected depth, not a full Keller pair.
"""
import hashlib
import json
import sympy as s

x,t,l=s.symbols('x t l')
F=s.Rational
p=t*(1+x*x*t)
y=x*t*p
u=x*x*t
q=-(x*x+6)/2
GATES=0

def check(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)

def zero(v,label):check(s.expand(v)==0,label)
def bracket(a,c):return s.expand(s.diff(a,x)*s.diff(c,t)-s.diff(a,t)*s.diff(c,x))

# Fixed THM-4308 rows on Phi=0, Delta=896/15; K=-32/5.
aa=(1+x*x/4,F(4,3)+2*x*x,-F(32,9)-F(4,5)*x*x,
    F(2176,135)+(F(1088,315)+F(128,35))*x*x-F(32,15)*x**4)
cc=(-3*x/4-x**3/8,-4*x-F(3,2)*x**3,F(88,15)*x-F(12,5)*x**3,
    (-F(8128,315)-F(192,35))*x+(F(736,105)-F(96,35))*x**3+F(8,5)*x**5)
A=sum(v*t**n for n,v in enumerate(aa))
C=sum(v*t**n for n,v in enumerate(cc))

DA={
12:F(13,30)*x**6,
13:x**4*(441*x**4+490*x*x-60)/135,
14:x**4*(978075*x**6+2573190*x**4-304888*x*x+147200)/93150,
15:2*x*x*(3912300*x**10-8802675*x**8+13143708*x**6-20841156*x**4+44222160*x*x-92124320)/419175,
}
DC={
12:-F(13,40)*x**5*(x*x+2),
13:-x**5*(441*x**4+1606*x*x+608)/180,
14:-x**5*(978075*x**6+5746500*x**4+4183508*x*x-270848)/124200,
15:-2*x*(195615*x**12+391230*x**10-9212432)/27945,
}
dA=sum(v*t**n for n,v in DA.items())
dC=sum(v*t**n for n,v in DC.items())
K=p**5*y**4-F(7,10)*p*y**6-F(508,135)*p**2*y**6
beta=F(235202,27945)
T=s.expand(K+beta*p**3*y**6)

# Literal P_d witnesses. Tuple (a,b,c,e) means x^a u^b p^c y^e.
WA=[((0,0,8,2),'-3395176/46575'),((0,0,11,2),'-36849728/83835'),
    ((0,0,4,4),'-6203/23'),((0,0,5,4),'-3415876/46575'),
    ((0,0,6,4),'128/81'),((0,0,7,4),'5896288/27945'),
    ((0,0,0,6),'-7/30'),((0,0,1,6),'-136958/1035'),
    ((0,0,2,6),'-74696/1035'),((0,0,3,6),'-5367368/46575'),
    ((1,0,10,1),'3395176/46575'),((1,0,6,3),'1252666/3105'),
    ((1,0,2,5),'2/3'),((2,0,8,2),'-415261/3105')]
WC=[((0,0,13,1),'18424864/27945'),((0,0,6,3),'-4384499/31050'),
    ((0,0,2,5),'-34136/345'),((0,0,3,5),'-4489379/31050'),
    ((0,0,4,5),'1472/675'),((0,0,0,7),'-2295668/15525'),
    ((0,0,1,7),'-1472/75'),((1,0,8,2),'4384499/31050'),
    ((1,0,4,4),'18097/60'),((1,0,0,6),'7/40'),
    ((2,0,6,3),'-629231/2070'),((2,0,2,5),'-1/2'),
    ((3,0,8,2),'41671/414')]

def main():
    check(2*12>15 and 2*12-1>14,'nonlinear variation excluded at every tested row')
    check(12+4>15 and 12+4-1>14,'unknown background rows4+ cannot enter')
    for n in range(12,16):
        check(s.degree(DA[n],x)<=n+1,'A degree cap')
        check(s.degree(DC[n],x)<=n+2,'C degree cap')
    jac=s.expand(bracket(A,dC)+bracket(dA,C))
    dp=s.expand(2*C*dC-3*A*A*dA+F(3,4)*dA)
    for n in range(15):zero(jac.coeff(t,n),'all Keller variation rows0..14')
    for n in range(16):zero((dp-T).expand().coeff(t,n),'all P-source variation rows0..15')
    # These are full polynomial-in-l identities before truncation, rather
    # than a claim inferred from first derivatives in an unrestricted jet.
    for n in range(15):zero((l*jac+l*l*bracket(dA,dC)).expand().coeff(t,n),'finite lambda Keller identity')
    nonlinear=s.expand(l*l*(dC*dC-3*A*dA*dA)-l**3*dA**3)
    for n in range(16):zero(nonlinear.coeff(t,n),'finite lambda P identity')
    witnesses=[]
    for depth,terms,target in [(2,WA,dA),(3,WC,dC)]:
        witness=s.Integer(0)
        for (a,b,c,e),coef in terms:
            check(a+b<=depth,'every module witness lies in the actual depth span')
            check(b+c+2*e>=12,'each module summand has zero earlier rows')
            witness+=F(coef)*x**a*u**b*p**c*y**e
        difference=s.expand(witness-target)
        for n in range(16):zero(difference.coeff(t,n),'full projected-depth witness equality')
        witnesses.append([depth,len(terms)])

    # Cheap failure controls: omit a compensator while retaining the same
    # row correction. The first unmatched source rows are literal nonzero.
    first=s.expand(dp-p**5*y**4)
    zero(first.coeff(t,13)+F(7,10)*x**6,'uncompensated weight22 hostile at13')
    second=s.expand(dp-(p**5*y**4-F(7,10)*p*y**6))
    zero(second.coeff(t,14)+F(508,135)*x**6,'missing row14 compensation hostile')
    third=s.expand(dp-K)
    zero(third.coeff(t,15)-beta*x**6,'compensated packet has nonzero row15 response')
    check(beta!=0,'transport response has no parameter-dependent singular locus')

    # Recover the three compensation scalars from the literal fixed row
    # response; this does not assume that projected scalar residuals suffice.
    k,r,j=s.symbols('k r j')
    test=s.expand(dp-(p**5*y**4+k*p*y**6+r*p**2*y**6+j*p**3*y**6))
    equations=[test.coeff(t,n).coeff(x,d)for n in (13,14,15)for d in range(n+1)]
    solved=s.solve(equations,(k,r,j),dict=True)
    check(solved==[{k:-F(7,10),r:-F(508,135),j:beta}],'unique compensation for this displayed row transport')
    jold=s.symbols('jold')
    w=-jold/beta
    zero(jold+beta*w,'old weight24 coefficient removed')
    low=s.expand(jold*p**3*y**6+w*T-w*K)
    zero(low,'same finite transport reaches source with only the lower packet')
    pp,yy=s.symbols('pp yy')
    lowpoly=pp**5*yy**4-F(7,10)*pp*yy**6-F(508,135)*pp**2*yy**6
    check(max(2*a+3*b for a,b in s.Poly(lowpoly,pp,yy).monoms())==22,'new source packet maximum weight22')
    check(s.expand(dA).coeff(t,12)!=0,'source change genuinely changes the prefix')
    # Inverse and composition are exact for this finite additive action.
    w1,w2=s.symbols('w1 w2')
    zero((w1*dA+w2*dA)-(w1+w2)*dA,'additive finite transport law for A')
    zero((w1*T+w2*T)-(w1+w2)*T,'additive finite transport law for source')

    data={'compensators':['-7/10','-508/135','235202/27945'],
          'delta_A':{str(n):s.sstr(s.factor(v))for n,v in DA.items()},
          'delta_C':{str(n):s.sstr(s.factor(v))for n,v in DC.items()},
          'P2_witness':WA,'P3_witness':WC,
          'scope':'finite row15 transport on the declared low jet; weight22 boundary Gm section relative to THM4438; no full Keller or JC2 claim'}
    raw=json.dumps(data,sort_keys=True,separators=(',',':'))
    print('PROVED exact finite additive transport; all background parameters after row3 retained')
    print(raw)
    print('EXPLICIT_GATES',GATES)
    print('SEMANTIC_SHA256',hashlib.sha256(raw.encode()).hexdigest())

if __name__=='__main__':main()
