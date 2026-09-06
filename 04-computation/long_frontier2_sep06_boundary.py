#!/usr/bin/env python3
"""Exact indefinite-curvature direction for fixed original-phase fibres.

The compactness, root-stability and boundary argument is analytic; this
source checks the complete carrier identity and its exact curvature.
"""
import hashlib
import json
import sympy as S

a,b,f,s,u=S.symbols('a b f s u')
GATES=0


def gate(ok,label):
    global GATES
    if ok != True:
        raise RuntimeError(label+': '+str(ok))
    GATES+=1


def zero(value,label):
    gate(S.cancel(value)==0,label)


def main():
    beta=[-f,b,-a,55,-13,1]
    c=[3*b/7,-2*a/3,45,-12,1]
    d=[b/7,-5*a/12,36,-11,1]
    def carrier(cs):
        n=len(cs)-1
        return S.Poly(S.expand((1+u)**14*sum(co*s**(n-j)*u**(2*j)
                                          for j,co in enumerate(cs))),u)
    hb,hc,hd=map(carrier,(beta,c,d))
    p=182-20020*s+2002*a*s**2-3432*b*s**3+2002*f*s**4
    zero(hb.nth(9)+s*p,'complete original ordinary carrier')
    qbar=S.cancel(((hb*hb).nth(18)-2*s*(hc*hd).nth(16))/s)
    fzero=12*b/(7*s)-a/s**2+10/s**3-1/(11*s**4)
    zero(p.subs(f,fzero),'same original phase retained')
    h=S.Poly(S.cancel(-S.Rational(11,14)*qbar.subs(f,fzero)),a,b,s)
    gate(len(h.terms())==19 and h.degree_list()==(2,2,8),'complete eliminated support')
    matrix=S.hessian(h.as_expr(),(a,b))
    aa=26558675*s**6+S.Rational(593856780,7)*s**5
    ab=-S.Rational(845791650,49)*s**7-S.Rational(3563140680,49)*s**6
    bb=S.Rational(286833690,49)*s**8+S.Rational(1972792800,49)*s**7
    for actual,expected in zip(matrix,[aa,ab,ab,bb]):
        zero(actual-expected,'exact Hessian entry')
    determinant=-S.Rational(31194786150,2401)*s**12*(10966105*s**2+72692884*s+144097056)
    zero(matrix.det()-determinant,'strictly negative determinant for s positive')
    direction=S.Matrix([-ab,aa])
    zero((direction.T*matrix*direction)[0]-aa*determinant,'explicit negative H curvature')
    df=12*direction[1]/(7*s)-direction[0]/s**2
    zero(S.diff(p,a)*direction[0]+S.diff(p,b)*direction[1]+S.diff(p,f)*df,
         'lifted direction preserves original phase identically')
    gate(all(x>0 for x in S.Poly(aa,s).coeffs()),'positive first diagonal')
    gate(all(x>0 for x in S.Poly(-ab,s).coeffs()),'positive first direction coordinate')
    gate(all(x>0 for x in S.Poly(-determinant,s).coeffs()),'strict negative determinant signs')
    # Ordinary polynomial response has negative Hessian determinant too;
    # the chosen direction has strictly positive Qbar second derivative.
    zero(S.det(S.hessian(qbar.subs(f,fzero),(a,b)))-S.Rational(196,121)*determinant,
         'response Hessian scaling')
    curvature=S.factor(-S.Rational(14,11)*aa*determinant)
    gate(all(x>0 for x in S.Poly(curvature,s).coeffs()),'strict positive response curvature')
    for phase in [S.Rational(1,100),S.Rational(1,8),S.Rational(3,2),S.Integer(10)]:
        gate(determinant.subs(s,phase)<0 and curvature.subs(s,phase)>0,
             'named exact positive phase control')
    # A weak, repeated/zero B boundary must stay in the boundary family.
    v=S.Symbol('v')
    zero((v**5-13*v**4+55*v**3-a*v*v+b*v-f).subs({a:75,b:0,f:0})
         -v*v*(v-3)*(v-5)**2,'actual nonnegative repeated and zero boundary')
    bank={'status':'PROVED fixed-phase boundary reduction; boundary sign remains OPEN',
          'universe':'formal complete carriers and exact Hessian; four rational phase controls; one actual B boundary',
          'Hessian_determinant':str(determinant),
          'direction_a':str(-ab),'direction_b':str(aa),'direction_f':str(S.factor(df)),
          'Qbar_second_derivative':str(curvature),
          'boundary_B_only':'B has a zero root or a repeated root',
          'boundary_two_interlacer':'B boundary, or gcd(B,C) nonconstant, or gcd(B,D) nonconstant'}
    data=json.dumps(bank,sort_keys=True,separators=(',',':'))
    print('EXACT FIXED-PHASE EXTREMA OCCUR ON NAMED COEFFICIENT BOUNDARIES')
    print(data)
    print('PASS explicit_gates='+str(GATES)+' semantic_sha256='+hashlib.sha256(data.encode()).hexdigest())


if __name__=='__main__':
    main()
