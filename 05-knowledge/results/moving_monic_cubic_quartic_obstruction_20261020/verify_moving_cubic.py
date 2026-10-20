#!/usr/bin/env python3
"""Exact controls for the entire moving cubic/quartic proof.

Symbolic identities are universal in the named variables. The finite controls
are explanatory witnesses only, not substitutes for the proof of existence
or nonexistence of arbitrary generators. Every check remains live under -O.
"""
import json
import sympy as S

x,y,t,s,U,V,a,b,c,ap,bp,cp,j,lam,B = S.symbols('x y t s U V a b c ap bp cp j lam B')
checks = []
def eq(name, left, right=0):
    value = S.factor(S.together(left-right))
    if value != 0:
        raise RuntimeError(f'{name}: {value}')
    checks.append(name)
def truth(name, value):
    if not bool(value):
        raise RuntimeError(name)
    checks.append(name)
def rem(poly, basis):
    out = S.expand(poly)
    coeffs=[]
    for g in sorted(basis, key=lambda z:S.degree(z,y), reverse=True):
        d=S.degree(g,y)
        co=S.Poly(out,y).coeff_monomial(y**d)/S.Poly(g,y).LC()
        out=S.expand(out-co*g)
        coeffs.append((str(g),str(S.factor(co))))
    return S.factor(out), coeffs

def formal_dx(poly):
    return S.diff(poly,a)*ap+S.diff(poly,b)*bp+S.diff(poly,c)*cp

h=y**3+a*y
k=y**4+b*y**2+c*y
F=-S.resultant(h-U,k-V,y)
I=3*c**2+4*a*(b-S.Rational(2,3)*a)**2
C=y**6+(a+b)*y**4+2*c*y**3+(a*a-a*b+b*b)*y*y+(a+b)*c*y+a*(a-b)**2+c*c
H=s**3-(b-2*a)*s-c
pair=y*y-s*y+s*s+a

eq('raw immersion resultant is 9I',S.resultant(S.diff(h,y),S.diff(k,y),y),9*I)
eq('plane relation evaluates to zero',F.subs({U:h,V:k}))
truth('plane relation is monic degree three in V',S.degree(F,V)==3 and S.Poly(F,V).LC()==1)
eq('conductor quotient identity without division',S.diff(F,V).subs({U:h,V:k}),S.diff(h,y)*C)
eq('global conductor membership C',C,h*h+(b-a)*k+2*c*h+a*(a-b)**2+c*c)
eq('global conductor membership Cy',C*y,h*k+c*k+(a-b)**2*h)
eq('global conductor membership Cy²',C*y*y,k*k+(a-b)*h*h-a*(a-b)*k+c*(a-b)*h)
truth('conductor remains monic degree six',S.degree(C,y)==6 and S.Poly(C,y).LC()==1)
eq('collision-product conductor',S.resultant(H,pair,s),C)
eq('separation resultant',S.resultant(H,-3*s*s-4*a,s),-9*I)
eq('H discriminant',S.discriminant(H,s),4*(b-2*a)**3-27*c*c)
# Pairwise coordinate differences in elementary symmetric variables.
ht=t**3+a*t
kt=t**4+b*t*t+c*t
hd=S.cancel((h-ht)/(y-t))
kd=S.cancel((k-kt)/(y-t))
eq('cubic divided difference',hd,(y+t)**2-y*t+a)
eq('quartic divided difference',kd,(y+t)*((y+t)**2-2*y*t+b)+c)
eq('cubic mod pair gives common h',S.rem(h,pair,y),-s*(s*s+a))
eq('common h after H reduction',S.rem(-s*(s*s+a)+c+(b-a)*s,H,s))
eq('quartic mod pair independent of y modulo H',S.rem(S.Poly(S.rem(k,pair,y),y).coeff_monomial(y),H,s))
eq('pair root separation',s*s-4*(s*s+a),-3*s*s-4*a)
tangent=(3*y*y+a)*(4*(s-y)**3+2*b*(s-y)+c)-(3*(s-y)**2+a)*(4*y**3+2*b*y+c)
tangent_quotient=S.cancel(tangent/(s-2*y))
eq('direct node transversality',S.rem(S.rem(tangent_quotient-(3*s*s+4*a)*S.diff(H,s),pair,y),H,s))
eq('b=a unit identity',9*I.subs(b,a),4*a**3+27*c*c)
eq('repeated H factors',H.subs({b:2*a+3*t*t,c:-2*t**3}),(s-t)**2*(s+2*t))
eq('repeated H immersion factor',I.subs({b:2*a+3*t*t,c:-2*t**3}),S.Rational(4,9)*(a+3*t*t)*(4*a+3*t*t)**2)
# Relative Jacobian and gap equations.
J1=formal_dx(h)*S.diff(k,y)-S.diff(h,y)*formal_dx(k)
J_expected=(4*ap-3*bp)*y**4-3*cp*y**3+(2*b*ap-a*bp)*y*y+(c*ap-a*cp)*y
eq('standard relative Jacobian',J1,J_expected)
T=j*y-lam*S.integrate(J1,y)
eq('gap five coefficient',S.Poly(T,y).coeff_monomial(y**5),-lam*(4*ap-3*bp)/5)
Tlow=S.expand(T.subs(bp,S.Rational(4,3)*ap))
reduced,_=rem(Tlow,[k,h,1])
eq('gap two coefficient',S.Poly(reduced,y).coeff_monomial(y**2),-lam*((3*b-2*a)*cp+2*c*ap)/4)
eq('standard y coefficient',S.Poly(reduced,y).coeff_monomial(y),j+lam*a*(2*b*ap-S.Rational(4,3)*a*ap)/3-S.Rational(3,4)*lam*c*cp)
eq('affine b substitution product derivative',((3*b-2*a)*cp+2*c*ap).subs(b,S.Rational(4,3)*a+B),(2*a+3*B)*cp+2*c*ap)
eq('unit I on standard differential locus',I.subs(b,S.Rational(4,3)*a+B),3*c*c+S.Rational(4,9)*a*(2*a+3*B)**2)
# Second reference exists globally once c is a nonzero constant.
k2=k.subs(b,S.Rational(2,3)*a)
q=h*h+S.Rational(3,2)*c*h
recovery=q*q-k2**3-2*a*k2*q+S.Rational(2,3)*a*a*k2*k2+(S.Rational(2,27)*a**3+S.Rational(3,4)*c*c)*q-(a**4/S.Integer(9)+a*c*c)*k2
eq('second-reference inverse identity',recovery,c**3*h/8)
eq('second-locus I is 3c²',I.subs(b,S.Rational(2,3)*a),3*c*c)
J2=S.diff(k2,a)*ap*S.diff(q,y)-S.diff(k2,y)*S.diff(q,a)*ap
J1_second=ap*(2*y**4+S.Rational(2,3)*a*y*y+c*y)
eq('second reference Jacobian',J2,-(2*h+S.Rational(3,2)*c)*J1_second)
V2=S.integrate(S.cancel(J2/ap),y)
V2_expected=-k2*k2/2-S.Rational(2,9)*a*h*h+a*a*k2/3-a*c*h/3-c*c*y*y/4
eq('second primitive full identity',V2,V2_expected)
r2,coeffs2=rem(V2,[k2*k2,h*k2,h*h,k2,h,1])
eq('second primitive gap remainder',r2,-c*c*y*y/4)
rT,_=rem(j*y-lam*ap*V2,[k2*k2,h*k2,h*h,k2,h,1])
eq('second required primitive gap remainder',rT,j*y+lam*ap*c*c*y*y/4)
# Hostile and positive finite controls, using rational exact values.
models={
 'three_ordinary_nodes':{a:-1,b:0,c:0},
 'ordinary_triple_point':{a:0,b:0,c:1},
 'repeated_H_immersed':{a:-S.Rational(9,4),b:-S.Rational(3,2),c:-2},
 'ramified_cusp':{a:0,b:0,c:0},
}
for name,sub in models.items():
    C0=S.Poly(C.subs(sub),y)
    eq(name+' conductor quotient',S.diff(F,V).subs({U:h,V:k}).subs(sub),S.diff(h,y).subs(sub)*C0.as_expr())
    truth(name+' degree six',C0.degree()==6 and C0.LC()==1)
    truth(name+' relation monic',S.Poly(F.subs(sub),V).LC()==1)
truth('three-node model squarefree conductor',S.degree(S.gcd(C.subs(models['three_ordinary_nodes']),S.diff(C,y).subs(models['three_ordinary_nodes'])),y)==0)
eq('triple-point conductor',(C.subs(models['ordinary_triple_point'])),(y**3+1)**2)
truth('triple-point H remains squarefree',S.discriminant(H.subs(models['ordinary_triple_point']),s)!=0)
eq('repeated-H control discr zero',S.discriminant(H.subs(models['repeated_H_immersed']),s))
truth('repeated-H control immersed',I.subs(models['repeated_H_immersed'])!=0)
eq('ramified model I zero',I.subs(models['ramified_cusp']))
# The moving family is immersed and two-generated in both classes,
# but its displayed Jacobians have unavoidable y dependence.
moving={a:x,b:S.Rational(2,3)*x,c:1}
eq('nonconstant moving unit-I hostile',I.subs(moving),3)
standard_moving={a:S.Rational(3,4)*(1-x*x),b:(3-x*x)/2,c:x}
eq('standard-only moving unit-I hostile',I.subs(standard_moving),3)
eq('standard-only hostile excludes second generator',(2*a-3*b).subs(standard_moving),-3)
eq('moving second reference inverse identity',recovery.subs({a:x,c:1}),h.subs(a,x)/8)
truth('moving standard Jacobian nonconstant in y',S.degree(J1.subs({**moving,ap:1,bp:S.Rational(2,3),cp:0}),y)==4)
truth('moving second Jacobian nonconstant in y',S.degree(J2.subs({a:x,c:1,ap:1}),y)==7)
# An entire ambient target automorphism provides scalar determinant 1.
G=U+x*V*V
L=V+x
eq('entire ambient-generator determinant',S.diff(G,U)*S.diff(L,V)-S.diff(G,V)*S.diff(L,U),1)
# Generic but not whole generation allows a nonunit scalar x.
Gbad=x*U
Lbad=V
eq('generic-only determinant hostile',S.diff(Gbad,U)*S.diff(Lbad,V)-S.diff(Gbad,V)*S.diff(Lbad,U),x)
eq('generic-only first generator vanishes at special fiber',Gbad.subs(x,0))
# Free-basis residues and exact bounded low-degree filtration.
truth('standard residue-separated basis',[d%3 for d in [0,4,8]]==[0,1,2])
semigroup={3*i+4*j for i in range(10) for j in range(10)}
truth('exact gaps through conductor',[i for i in range(12) if i not in semigroup]==[1,2,5])
truth('low degree four span',[i for i in range(5) if i in semigroup]==[0,3,4])
truth('low degree eight span',[i for i in range(9) if i in semigroup]==[0,3,4,6,7,8])
output={
 'status':'PASS','exact_checks':len(checks),'checks':checks,
 'scope':'Universal symbolic identities plus named exact controls; the all-generating-pair theorem also uses the written proof and pinned classical degree bound.',
 'second_primitive_reduction':coeffs2,
 'conductor':str(C),
 'hostile':'a=x,b=2x/3,c=1 has I=3 and both global reference generating pairs, but neither has constant nonzero Jacobian. Generic-only pair (xh,k) has ambient determinant x and fails whole generation at x=0.'
}
print(json.dumps(output,indent=2,sort_keys=True))
