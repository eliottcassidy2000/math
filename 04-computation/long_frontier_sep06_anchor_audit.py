#!/usr/bin/env python3
"""Independent ordinary-carrier audit of the C-only first-phase theorem."""
import hashlib
import sympy as S

a,b,f,s,u,w = S.symbols('a b f s u w')
GATES=0
TRACE=hashlib.sha256()


def gate(ok,label):
    global GATES
    if ok != True:
        raise RuntimeError(label+': '+str(ok))
    GATES+=1
    TRACE.update((label+'\n').encode())


def zero(expr,label):
    gate(S.cancel(expr)==0,label)


def main():
    B=[-f,b,-a,55,-13,1]
    C=[3*b/7,-2*a/3,45,-12,1]
    D=[b/7,-5*a/12,36,-11,1]
    def carrier(coeff):
        n=len(coeff)-1
        return S.Poly(S.expand((1+u)**14*sum(c*s**(n-j)*u**(2*j)
                                  for j,c in enumerate(coeff))),u)
    HB,HC,HD=map(carrier,(B,C,D))
    original=182-20020*s+2002*a*s**2-3432*b*s**3+2002*f*s**4
    zero(HB.nth(9)+s*original,'ordinary original P extraction')
    response=S.cancel(((HB*HB).nth(18)-2*s*(HC*HD).nth(16))/s)
    # Independent expansion at infinity of the rational interlacer quotient.
    denominator=sum(B[5-j]*w**j for j in range(6))
    numerator=sum(C[4-j]*w**j for j in range(5))
    moments=S.Poly(S.series(numerator/denominator,w,0,5).removeO(),w)
    mu=[moments.nth(j) for j in range(5)]
    expected=[1,1,3,a/3-16,16*a/3-373-4*b/7]
    for j in range(5):
        zero(mu[j]-expected[j],'independent residue moment '+str(j))
    zero(S.det(S.Matrix(3,3,lambda i,j:mu[i+j]))
         -((a-75)*(135-a)/9-8*b/7),'independent full Hankel determinant')
    replacement=10/s-1/(11*s**2)+12*b*s/7-f*s**2
    zero(original.subs(a,replacement),'original zero retained')
    reduced=S.Poly(S.expand(response.subs(a,replacement)),b,f)
    H=9335990950*s**4+19125419885*s**3-1420269695*s**2+19755450*s-63866
    terms={
      (2,0):-39330*s**7*(19601*s+31920)/49,
      (1,1):26220*s**8*(9605*s+24708)/7,
      (1,0):30*s**3*(600168371*s**3+1898076564*s**2-166142340*s+1753752)/77,
      (0,2):-2185*s**9*(7735*s+24708),
      (0,1):-5*s**4*(1724814091*s**3+5260283632*s**2-716499174*s+13585572)/11,
      (0,0):-7*H/121}
    gate(set(reduced.monoms())==set(terms),'all six response monomials retained')
    for powers,value in terms.items():
        zero(reduced.coeff_monomial(b**powers[0]*f**powers[1])-value,
             'ordinary-carrier eliminated coefficient '+str(powers))
    lo,hi=S.Rational(1,110),S.Rational(1,90)
    gate(S.Rational(13,5)**5<119,'AM-GM endpoint')
    gate(75-S.Rational(12,7)*S.Rational(175,2)*lo>0,'first phase left margin')
    gate(182-20020*hi+2002*135*hi**2+2002*119*hi**4<0,'first phase right margin')
    gate(-20020+4004*135*hi+8008*119*hi**3<0,'first phase strictly decreasing')
    gate(12*9335990950*hi**2+6*19125419885*hi-2*1420269695<0,'H second derivative upper bound')
    gate(S.diff(H,s).subs(s,lo)<0,'H first derivative left endpoint')
    gate(H.subs(s,hi)>6000,'H full interval minimum')
    gate(13585572-716499174*hi>0,'negative f coefficient on full interval')
    gate(S.Rational(30,77)*hi**3*(600168371*hi**3+1898076564*hi**2+1753752)<2,'positive b coefficient upper bound')
    gate(S.Rational(26220,7)*hi**8*(9605*hi+24708)<S.Rational(1,1000),'positive bf coefficient upper bound')
    gate(-S.Rational(42000,121)+175+S.Rational(175,2)*119/1000<-160,'full uniform response margin')
    alternative=12*b/(7*s)-a/s**2+10/s**3-1/(11*s**4)
    alt=S.expand(response.subs(f,alternative))
    det=S.det(S.hessian(alt,(a,b)))
    zero(det+S.Rational(1031232600,49)*s**12*(10966105*s**2+72692884*s+144097056),
         'alternate coefficient elimination remains indefinite')
    v=S.Symbol('v')
    boundary={a:75,b:0,f:0}
    bb=sum(c*v**j for j,c in enumerate(B))
    cc=sum(c*v**j for j,c in enumerate(C))
    zero((cc/bb).subs(boundary)-(S.Rational(2,3)/v+S.Rational(1,3)/(v-3)),
         'weak repeated boundary positive residue measure')
    print('UNIVERSE formal a,b,f,s coefficient identities; exact interval bounds; weak repeated boundary')
    print('INDEPENDENT_PATH ordinary (1+u)^14 carriers and coefficient extractions 9,18,16; rational series at infinity')
    print('CLAIM C-only interlacing and e1=13,e2=55: unique smallest phase in (1/110,1/90); sQ(-s)<-160')
    print('SCOPE other phases OPEN; no D interlacing needed; no finite parameter extrapolation')
    print('PASS explicit_gates='+str(GATES)+' semantic_sha256='+TRACE.hexdigest())


if __name__=='__main__':
    main()
