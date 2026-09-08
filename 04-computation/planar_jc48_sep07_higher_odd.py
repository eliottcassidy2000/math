#!/usr/bin/env python3
"""Exact local exhaustion and exceptional-chart gluing for (4,6) odd cusps.

Topology and connected good-locus statements have a separate analytic proof.
No numerical braid or unproved representative certificate is used here.
"""
import hashlib
import json
import sympy as s

t,z,a,b,c,d,e,g = s.symbols('t z a b c d e g')
F=s.Rational
GATES=0


def check(value,label):
    global GATES
    GATES+=1
    if not bool(value):
        raise RuntimeError(label)


def eq(left,right,label):
    check(s.cancel(left-right)==0,label)


def coeff(poly,n):
    return s.factor(s.expand(poly).coeff(t,n))


def main():
    U=t**4+a*t**3+b*t*t
    V=c*t**6+e*t**5+d*a*t**3+d*b*t*t
    W=s.expand(V-d*U+d*U**2/b**2)
    for n in range(5):
        eq(coeff(W,n),0,'finite low cancellations')
    eq(coeff(W,5),e+2*a*d/b,'finite fifth criterion')
    W=s.expand(W.subs({c:2,e:-2*a*d/b}))
    k3=coeff(W,6)/b**3
    W=s.expand(W-k3*U**3)
    seventh=-a*(6*b*b+(3*a*a+4*b)*d)/b**3
    eq(coeff(W,7),seventh,'general finite seventh')
    eq(coeff(W,6),0,'sixth even term removed')
    un=U.subs(a,1)
    wn=s.expand(W.subs(a,1))
    d9=-6*b*b/(4*b+3)
    ninth_stage=s.expand(wn.subs(d,d9))
    k4=coeff(ninth_stage,8)/b**4
    eq(k4,6*(b+2)/(b**5*(4*b+3)),'finite eighth even coefficient')
    ninth_stage=s.expand(ninth_stage-k4*un**4)
    for n in range(9):
        eq(coeff(ninth_stage,n),0,'finite ninth preceding terms')
    eq(coeff(ninth_stage,9),-44/(b*b*(4*b+3)),'finite ninth cannot vanish')
    eq(seventh.subs({a:1,b:-F(3,4)}),8,'exceptional denominator is finite seven')

    vn=V.subs({a:1,c:2,e:-2*d/b})
    X=s.cancel(z*z*(1+z+b*z*z)/(2-2*d*z/b+d*z**3+d*b*z**4))
    Z=s.cancel(z**6/(2-2*d*z/b+d*z**3+d*b*z**4))
    inf7=s.limit((Z-4*X**3)/z**7,z,0)
    eq(inf7,-(2*d+3*b)/(2*b),'general infinity seventh')
    eq(inf7.subs(d,d9),-9/(2*(4*b+3)),'finite nine forces infinity seven')
    d79=-3*b/2
    X9=X.subs(d,d79);Z9=Z.subs(d,d79)
    rho=s.cancel(Z9/X9**3-4)
    k_inf4=-6-24*b
    eq(s.limit((rho/X9-k_inf4)/z,z,0),7,'infinity ninth supplier constant')
    eq(s.limit((Z9-4*X9**3-k_inf4*X9**4)/z**9,z,0),F(7,16),'actual infinity ninth')
    eq(seventh.subs({a:1,d:d79}),9/(2*b*b),'finite seven on infinity nine')

    # The exceptional finite chart is not a b-division limit assumption.
    ue=t**4;ve=2*t**6+e*t**5+g*t*t
    we=s.expand(ue-ve**2/g**2)
    for n in range(7):
        eq(coeff(we,n),0,'exceptional cusp low rows')
    eq(coeff(we,7),-2*e/g,'exceptional finite seventh')
    xe=z*z/(2+e*z+g*z**4);ze=z**6/(2+e*z+g*z**4)
    eq(s.limit((ze-4*xe**3)/z**7,z,0),e/2,'exceptional infinity seventh')
    check(all(coeff(ve.subs(e,0),n)==0 for n in [1,3,5]),'exceptional even-map hostile')

    # One irreducible coefficient family contains regular and exceptional charts.
    relation=2*a*g+e*b*b
    familyV=2*t**6+e*t**5-e*b*t**3/2+g*t*t
    eq(familyV.subs(e,-2*a*g/b**2),V.subs({c:2,d:g/b,e:-2*a*g/b**2}),
       'b-chart is exact general normal form')
    factors=s.factor_list(relation)[1]
    check(len(factors)==1 and factors[0][1]==1,'irreducible defining polynomial control')
    eq(s.diff(relation,a),2*g,'smooth g chart')
    eq(s.diff(relation,e),b*b,'smooth b chart')
    ag=-e*b*b/(2*g)
    Ug=U.subs(a,ag);Vg=familyV
    Wg=s.expand(Ug-b*Vg/g-Vg**2/g**2)
    kg3=coeff(Wg,6)/g**3
    eq(kg3,-b*(b*e*e+8*g)/(4*g**5),'g-chart even coefficient')
    Wg=s.expand(Wg-kg3*Vg**3)
    for n in range(7):
        eq(coeff(Wg,n),0,'g-chart preceding terms')
    g7=-e*(3*b**3*e**2+24*b*b*g+16*g*g)/(8*g**3)
    eq(coeff(Wg,7),g7,'g-chart finite seventh')
    eq(g7,-b/g*seventh.subs({a:ag,d:g/b}),'overlap nonvanishing agrees')
    eq(g7.subs(b,0),-2*e/g,'exceptional section is interior to cusp-seven locus')
    eq((e-3*ag)/2,e*(3*b*b+2*g)/(4*g),'g-chart infinity seventh')
    eq(Ug.subs(b,0),ue,'actual exceptional U in family')
    eq(Vg.subs(b,0),ve,'actual exceptional V in family')
    eq(seventh.subs({a:1,b:1,d:0}),-6,'b-only chart positive local control')
    eq(inf7.subs({b:1,d:0}),-F(3,2),'b-only infinity-seven control')
    eq(g7.subs({b:0,e:1,g:1}),-2,'exceptional local control')
    eq(((e-3*ag)/2).subs({b:0,e:1,g:1}),F(1,2),'exceptional local infinity control')
    # a=0 in the b-chart gives an even map and cannot be birational.
    evenU=U.subs(a,0);evenV=V.subs({a:0,c:2,e:0})
    eq(evenU.subs(t,-t),evenU,'excluded a-zero U is even')
    eq(evenV.subs(t,-t),evenV,'excluded a-zero V is even')

    pairs=[(7,7,4),(7,9,3),(9,7,3)]
    for finite,infinity,nodes in pairs:
        check((finite-1)//2+(infinity-1)//2+nodes==10,'exact rational sextic genus')
        check(nodes>=2,'declared node range')
    data={'finite_infinity_nodes':pairs,'finite_seventh':str(seventh),
          'finite_ninth':'-44/(b^2*(4b+3))','infinity_ninth':'7/16',
          'family_equation':str(relation),'g_chart_finite_seventh':str(g7),
          'g_chart_infinity_seventh':str(e*(3*b*b+2*g)/(4*g)),
          'scope':'local exhaustion and connected-family mechanism; actual certified representatives are separate dependencies'}
    print('STATUS: exact local classification; independent audit and representative certificates pending')
    print(json.dumps(data,sort_keys=True,indent=2))
    print('NO_FINITE_ELEVEN: finite ninth numerator is the nonzero constant -44')
    print('EXCEPTIONAL_RETAINED: b=0,g*e!=0 is interior to one smooth irreducible finite7/infinity7 family')
    print('SEMANTIC_SHA256',hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest())
    print('PASS gates='+str(GATES))


if __name__=='__main__':
    main()
