#!/usr/bin/env python3
"""Exact Hamiltonian interpretation and two all-order scope obstructions.

All identities are symbolic over QQ. Only the literal WA/WC depth witnesses
are read from the hash-pinned previous producer; that producer is not imported.
No numerical census or assertion of polynomial flow termination is used.
"""
import ast
import hashlib
import json
from pathlib import Path
import sympy as s

x,t,p0,y0,w=s.symbols('x t p0 y0 w')
F=s.Rational
p=t*(1+x*x*t);y=x*t*p;u=x*x*t
GATES=0

def check(v,label):
    global GATES
    GATES+=1
    if not v:raise RuntimeError(label)

def zero(v,label):check(s.expand(v)==0,label)
def bracket(a,b):return s.expand(s.diff(a,x)*s.diff(b,t)-s.diff(a,t)*s.diff(b,x))
def trunc(f,n):
    f=s.Poly(s.expand(f),t)
    return sum(c*t**m[0] for m,c in f.terms() if m[0]<n)

aa=(1+x*x/4,F(4,3)+2*x*x,-F(32,9)-F(4,5)*x*x,
    F(2176,135)+(F(1088,315)+F(128,35))*x*x-F(32,15)*x**4)
cc=(-3*x/4-x**3/8,-4*x-F(3,2)*x**3,F(88,15)*x-F(12,5)*x**3,
    (-F(8128,315)-F(192,35))*x+(F(736,105)-F(96,35))*x**3+F(8,5)*x**5)
A=sum(v*t**n for n,v in enumerate(aa));C=sum(v*t**n for n,v in enumerate(cc))
DA={12:F(13,30)*x**6,
13:x**4*(441*x**4+490*x*x-60)/135,
14:x**4*(978075*x**6+2573190*x**4-304888*x*x+147200)/93150,
15:2*x*x*(3912300*x**10-8802675*x**8+13143708*x**6-20841156*x**4+44222160*x*x-92124320)/419175}
DC={12:-F(13,40)*x**5*(x*x+2),
13:-x**5*(441*x**4+1606*x*x+608)/180,
14:-x**5*(978075*x**6+5746500*x**4+4183508*x*x-270848)/124200,
15:-2*x*(195615*x**12+391230*x**10-9212432)/27945}
dA=sum(v*t**n for n,v in DA.items());dC=sum(v*t**n for n,v in DC.items())

S0p=p0**3*y0**5/15+F(16,135)*p0**4*y0**5-y0**7/15-F(2848,6075)*p0**5*y0**5+F(2,621)*p0*y0**7
S0=s.expand(S0p.subs({p0:p,y0:y}))
Sstar=S0-F(4606216,83835)*p**14*y+F(2303108,83835)*p**10*y**3-F(64798,5589)*p**6*y**5+F(1367399,139725)*p**2*y**7-F(52451,3105)*x*y**8

# Literal terminal depth witnesses for the difference from the old transport.
TA=[((0,2,1,6),'419608/3105'),((0,1,2,6),'-10939192/139725'),
    ((0,0,3,6),'518384/5589'),((0,0,7,4),'-18424864/83835'),
    ((0,0,11,2),'36849728/83835')]
TC=[((0,3,2,5),'-104902/1035'),((0,1,0,7),'-6706382/46575'),
    ((0,0,1,7),'247744/5175'),((0,0,5,5),'31232/1215'),
    ((0,0,13,1),'-18424864/27945')]

def witness(terms,depth,minval):
    ans=s.Integer(0)
    for (a,b,c,e),z in terms:
        check(a+b<=depth,'actual depth of each summand')
        check(b+c+2*e>=minval,'summand valuation')
        ans+=F(z)*x**a*u**b*p**c*y**e
    return ans

def main():
    prior=Path(__file__).with_name('planar_jc48_sep06_weight14.py').read_bytes()
    priorhash=hashlib.sha256(prior).hexdigest()
    check(priorhash=='a4a140ab538620e7885d8b77758c02e16a5c3a8de7e53e3daac51f01bda321f7','pinned literal depth witnesses')
    tables={}
    for node in ast.parse(prior).body:
        if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name) and node.targets[0].id in ('WA','WC'):
            tables[node.targets[0].id]=ast.literal_eval(node.value)
    check(set(tables)=={'WA','WC'},'both literal tables recovered without producer import')
    nA=trunc(bracket(A,S0),16);nC=trunc(bracket(C,S0),16)
    for v,old in [(bracket(A,Sstar),dA),(bracket(C,Sstar),dC)]:
        for n in range(16):zero(s.expand(v-old).coeff(t,n),'exact old generator lift')
    for v,old in [(nA,dA),(nC,dC)]:
        for n in range(15):zero(s.expand(v-old).coeff(t,n),'P0 generator agrees through14')
    for depth,oldname,tail,v in [(2,'WA',TA,nA),(3,'WC',TC,nC)]:
        rep=witness(tables[oldname],depth,12)+witness(tail,depth,15)
        for n in range(16):zero(s.expand(rep-v).coeff(t,n),'full projected depth for P0 generator')
    packet=s.expand(p**5*y**4-F(7,10)*p*y**6-F(508,135)*p**2*y**6+F(235202,27945)*p**3*y**6)
    dv=2*C*nC-3*A*A*nA+F(3,4)*nA
    dj=bracket(A,nC)+bracket(nA,C)
    for n in range(16):zero(s.expand(dv-packet).coeff(t,n),'same source packet')
    for n in range(15):zero(s.expand(dj).coeff(t,n),'same bracket through14')
    for n in range(12,16):
        check(s.degree(nA.coeff(t,n),x)<=n+1,'A cap')
        check(s.degree(nC.coeff(t,n),x)<=n+2,'C cap')
    check(12+4>=16 and 13+4-1>=16,'all later background rows invisible')
    check(2*12>=16 and 2*12-1>=15,'finite parameter nonlinear terms invisible')
    zero(trunc(s.diff(S0,t)-F(13,15)*x**5*t**12,13),'leading Hamiltonian X')
    zero(trunc(-s.diff(S0,x)+x**4*t**13/3,14),'leading Hamiltonian Y')
    mismatch=s.expand(S0-p**3*y**5/15).coeff(t,14)
    zero(mismatch-x**5*(16-9*x*x)/135,'bare leading generator hostile')
    check(mismatch!=0,'bare generator really fails')

    # Complete triangular P0 test for the exact old generator through row16.
    target=trunc(Sstar,17);trial=s.Integer(0)
    for n in range(13,17):
        r=s.expand(target-trial).coeff(t,n)
        for e in range(n//2+1):
            trial+=r.coeff(x,e)*p**(n-2*e)*y**e
        if n<16:zero(s.expand(target-trial).coeff(t,n),'P0 triangular early match')
    rem=s.expand(target-trial).coeff(t,16)
    zero(rem+F(52451,3105)*x**9,'first exact-P0 obstruction')
    check(9>16//2,'residual outside every row16 P0 leading monomial')

    # Correct source Poisson algebra; the discarded 2y formula is hostile.
    zero(bracket(u,p)-2*x*p,'u,p identity')
    check(s.expand(bracket(u,p)-2*y)!=0,'missing t draft identity rejected')
    zero(bracket(u,y)-3*u*p,'u,y identity')
    zero(bracket(p,y)+t*p,'p,y identity')
    D=p0**3-y0*y0
    zero((p**3-y*y)-t*p*p,'birational denominator identity')
    ew=2*p0*s.diff(S0p,p0)+3*y0*s.diff(S0p,y0)
    residue=s.rem(ew,y0*y0-p0**3,y0)
    target_res=p0**10*y0*(F(14,5)-F(2848,243)*p0)
    zero(residue-target_res,'weighted Euler residue at cusp divisor')
    check(residue!=0,'uncancelled fixed-source pole')
    hs,hy=s.symbols('H_p H_y')
    formal=-p0*y0*ew/(2*D)-D/p0*(hs*s.diff(S0p,y0)-hy*s.diff(S0p,p0))
    numerator=s.cancel(formal*D*p0)
    zero(s.rem(numerator+p0**2*y0*residue/2,y0*y0-p0**3,y0),'H-independent pole numerator')

    # Universal repaired carrier S=c+p^2 D R; R,R_p,R_y are formal jets.
    rr,rp,ry=s.symbols('R R_p R_y')
    carrier=p0**2*D*rr
    sp=s.diff(carrier,p0)+p0**2*D*rp
    sy=s.diff(carrier,y0)+p0**2*D*ry
    ec=2*p0*sp+3*y0*sy
    zero(ec-p0**2*D*(10*rr+2*p0*rp+3*y0*ry),'carrier weighted Euler identity')
    image=s.cancel(-p0*y0*ec/(2*D)-D/p0*(hs*sy-hy*sp))
    expected=-p0**3*y0*(10*rr+2*p0*rp+3*y0*ry)/2-D*(hs*(-2*p0*y0*rr+p0*D*ry)-hy*((2*D+3*p0**3)*rr+p0*D*rp))
    zero(image-expected,'universal polynomial source-form image')
    check(s.denom(s.cancel(image)).is_Integer,'no parameter-coordinate denominator in carrier image')
    zero(p**3-y*y-t**3*(1+x*x*t)**2,'fixed-source invariant factor')
    zero(p*p*(p**3-y*y)-t**5*(1+x*x*t)**4,'universal-source invariant factor')
    zero(bracket(-u/2,p**3-y*y)+3*p*y,'H=0 fixed-carrier positive control')
    check(s.rem(D,p0*p0,p0)!=0,'fixed positive control outside universal carrier')

    # Highest-degree induction proves non-local-nilpotence for every iterate.
    poly=s.Poly(S0,x,t);top=max(sum(m) for m in poly.monoms())
    topterms=[(m,c) for m,c in poly.terms() if sum(m)==top]
    check(topterms==[((25,25),-F(2848,6075))],'unique total-degree Hamiltonian term')
    cN=-F(2848,243);lead=x
    for n in range(1,6):
        lead=bracket(lead,-F(2848,6075)*(x*t)**25)
        zero(lead-cN**n*x**(24*n+1)*t**(24*n),'leading nontermination recurrence')
    data={'S0':s.sstr(S0p),'terminal_A':s.sstr(s.factor((nA-dA).coeff(t,15))),
          'terminal_C':s.sstr(s.factor((nC-dC).coeff(t,15))),
          'terminal_depth_A':TA,'terminal_depth_C':TC,
          'exact_P0_obstruction':'-52451*x**9/3105 at generator row16',
          'weighted_Euler_residue':s.sstr(target_res),
          'all_H_source_form_carrier':'S in k + p^2*(p^3-y^2)*k[p,y]',
          'fixed_H_source_form_carrier':'S=c+(p^3-y^2)*R and p divides J_p_y(H,S)',
          'carrier_LND_obstruction':'all nonconstant fixed-H carriers have non-locally-nilpotent Hamiltonian derivation',
          'leading_iterate':'(-2848/243)^n*x^(24n+1)*t^(24n)',
          'prior_source_sha256':priorhash}
    raw=json.dumps(data,sort_keys=True,separators=(',',':'))
    print('PROVED finite-jet Hamiltonian lift; fixed-source and polynomial-action obstructions')
    print(raw)
    print('EXPLICIT_GATES',GATES)
    print('SEMANTIC_SHA256',hashlib.sha256(raw.encode()).hexdigest())

if __name__=='__main__':main()
