#!/usr/bin/env python3
"""All-h beta-step identities and an exact obstruction to sign propagation.

The obstruction is a model with genuine first rows and carry/top anchors,
not an actual Laurent counterexample. The eight genuine h7/h8 controls
retain every count and negative auxiliary index. No endpoint table compiler.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import factorial, comb, gcd
import json
import sympy as S

GATES=0
TRACE=[]


def gate(ok,label,data=None):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(f'{label}: {data}')
    TRACE.append((label,data))


def rows(h,x):
    p=[F(factorial(2*h+1)*factorial(x+h),
         factorial(x+j)*factorial(3*h-3*j)*factorial(1+2*j))
       for j in range(h+1)]
    q=[F(factorial(2*x+2*h),
         factorial(2*x+e)*factorial(6*h-3*e)*factorial(2+2*e))
       for e in range(-1,2*h+1)]
    return p,q


def remainder(a,p):
    degree=len(p)-1
    a=list(a)
    while len(a)>degree:
        c=a.pop()/p[-1]
        offset=len(a)-degree
        for j in range(degree):
            a[offset+j]-=c*p[j]
    return a+[F(0)]*(degree-len(a))


def response_characteristic(a,p):
    degree=len(p)-1
    rr=remainder(a[1:],p)
    rr=[rr[j]-a[0]*p[j+1]/p[0] for j in range(degree)]
    cols=[remainder([F(0)]*j+rr,p) for j in range(degree)]
    mat=S.Matrix(degree,degree,lambda i,j:cols[j][i])
    return [F(v) for v in mat.charpoly().all_coeffs()]


def literal(h,x,mass):
    g=x+3*h+1
    a=6*h+3
    out={}
    for nc in range(mass+1):
        for nb in range(mass-nc+1):
            na=mass-nb-nc
            if -a*na+(2*g-a)*nb+(3*g-a)*nc==0:
                out[(na,nb,nc)]=F(factorial(mass),
                    factorial(na)*factorial(nb)*factorial(nc))
    return out


def symbolic_step_checks():
    x=S.Symbol('x')
    def fall(a,n):
        return S.prod(a-j for j in range(n))
    for h in (1,2,7,8):
        for j in range(h+1):
            p=fall(x+h,h-j)*S.Rational(factorial(2*h+1),
                                      factorial(3*h-3*j)*factorial(1+2*j))
            expr=(x+j+1)*p.subs(x,x+1)-(x+h+1)*p
            gate(S.Poly(expr,x).is_zero,'symbolic first beta step',(h,j))
        for e in range(-1,2*h+1):
            q=fall(2*x+2*h,2*h-e)*S.Rational(1,
                       factorial(6*h-3*e)*factorial(2+2*e))
            expr=(2*x+e+1)*(2*x+e+2)*q.subs(x,x+1)
            expr-=(2*x+2*h+1)*(2*x+2*h+2)*q
            gate(S.Poly(expr,x).is_zero,'symbolic carried second beta step',(h,e))
    print('IDENTITIES all-h proof by coefficient ratios; exact symbolic controls h=1,2,7,8')


def model_obstruction():
    t=S.Symbol('t')
    old=S.Poly((t+S.Rational(2,43))*(t+S.Rational(1,2))
               *(t+S.Rational(43,21)),t)
    new=S.Poly(sum(30*old.nth(k)*t**k/((k+2)*(k+3))
                   for k in range(4)),t)
    expected_old=[S.Rational(1,21),S.Rational(2063,1806),
                  S.Rational(4685,1806),S.Integer(1)]
    expected_new=[S.Rational(5,21),S.Rational(10315,3612),
                  S.Rational(4685,1204),S.Integer(1)]
    gate([old.nth(k) for k in range(4)]==expected_old,'exact model starting cubic')
    gate([new.nth(k) for k in range(4)]==expected_new,'exact model transformed cubic')
    p1,q1=rows(1,1)
    p2,q2=rows(1,2)
    gate(p1==[F(2),F(1)] and p2==[F(3),F(1)],'both first rows are genuine')
    gate([F(expected_old[k]) for k in (0,3)]==[720*q1[k] for k in (0,3)],
         'original model carry and leading coefficients are genuine')
    gate([F(expected_new[k]) for k in (0,3)]==[720*q2[k] for k in (0,3)],
         'transformed model carry and leading coefficients are genuine')
    gate(any(F(expected_old[k])!=720*q1[k] for k in (1,2)),
         'model interior coefficients explicitly not genuine')
    gate(old.eval(-2)/(-2)==-S.Rational(3,43),'initial same-root response negative')
    gate(new.eval(-3)/(-3)==S.Rational(557,5418),'next same-root response strictly positive')
    gate(all(c>0 for c in expected_old+expected_new),'all model coefficients strictly positive')
    brackets=((-3,-2),(-1,S.Rational(-1,2)),(S.Rational(-1,2),0))
    endpoint_values=[]
    for lo,hi in brackets:
        lv,hv=new.eval(lo),new.eval(hi)
        gate(lv*hv<0,'each transformed negative root has its own sign bracket',(str(lo),str(hi)))
        endpoint_values.append([str(lo),str(hi),str(lv),str(hv)])
    gate(S.discriminant(new,t)==S.Rational(1152093393330275,56737436781312)>0,
         'independent positive transformed cubic discriminant')
    # Literal reconstruction of both actual h1 rows verifies all anchors.
    for xx,pp,qq in ((1,p1,q1),(2,p2,q2)):
        g=xx+4
        first=literal(1,xx,g);second=literal(1,xx,2*g)
        gate(set(first)=={(xx+j,3-3*j,1+2*j) for j in range(2)},
             'full actual h1 first fibre',xx)
        gate(set(second)=={(2*xx+e,6-3*e,2+2*e) for e in range(-1,3)},
             'full actual h1 carried fibre',xx)
        gate([first[(xx+j,3-3*j,1+2*j)]/comb(g,3) for j in range(2)]==pp,
             'literal genuine first coefficient laws',xx)
        scale=F(factorial(2*g),factorial(2*xx+2))
        gate([second[(2*xx+e,6-3*e,2+2*e)]/scale for e in range(-1,3)]==qq,
             'literal genuine second coefficient laws',xx)
    data=dict(h=1,x=1,first_rows=['t+2','t+3'],
              model_tQ_before=[str(v) for v in expected_old],
              model_tQ_after=[str(v) for v in expected_new],
              response_before='-3/43',response_after='557/5418',
              genuine_scaled_tQ_before=[str(720*v) for v in q1],
              genuine_scaled_tQ_after=[str(720*v) for v in q2],
              transformed_root_brackets=endpoint_values)
    print('MODEL_OBSTRUCTION '+json.dumps(data,separators=(',',':')))
    # A second obstruction has both adjacent designated masses primitive.
    old3=S.Poly((t+S.Rational(1,6))*(t+S.Rational(100,101))
                 *(t+S.Rational(101,25)),t)
    new3=S.Poly(sum(90*old3.nth(k)*t**k/((k+6)*(k+7))
                    for k in range(4)),t)
    p3,q3=rows(1,3)
    p4,q4=rows(1,4)
    gate(p3==[F(4),F(1)] and p4==[F(5),F(1)],'primitive-mass model first rows')
    gate(gcd(7,9)==gcd(8,9)==1,'both adjacent designated masses are primitive')
    gate([F(old3.nth(k)) for k in (0,3)]==[720*q3[k] for k in (0,3)],
         'primitive-mass original carry and top anchors')
    gate([F(new3.nth(k)) for k in (0,3)]==[720*q4[k] for k in (0,3)],
         'primitive-mass transformed carry and top anchors')
    gate(old3.eval(-4)/(-4)==-S.Rational(874,7575), 'primitive-mass initial response')
    gate(new3.eval(-5)/(-5)==S.Rational(221,21210), 'primitive-mass reversed response')
    gate(all(new3.nth(k)>0 for k in range(4)) and S.discriminant(new3,t)>0,
         'primitive-mass transformed cubic has distinct negative roots')
    print('PRIMITIVE_MASS_MODEL '+json.dumps(dict(h=1,x=3,g_before=7,g_after=8,
        old_tQ=[str(old3.nth(k)) for k in range(4)],
        new_tQ=[str(new3.nth(k)) for k in range(4)],
        response_before='-874/7575',response_after='221/21210',
        new_discriminant=str(S.discriminant(new3,t))),separators=(',',':')))
    print('SCOPE model route refuted; actual interior factorial amplitudes differ; no actual Laurent counterexample')


def genuine_named_controls():
    controls=[]
    for h in (7,8):
        for xx in (1,2,10,100):
            g=xx+3*h+1
            p,q=rows(h,xx)
            B={j:comb(xx+3*h-2*j,xx+j) for j in range(-xx,h+1)}
            # Complete negative support is retained before convolution.
            K=F(factorial(2*g),factorial(2*xx+2*h))
            W=[F(comb(2*g,2+2*e)*sum(b*B.get(e-j,0) for j,b in B.items()),K)
               for e in range(-1,2*h+1)]
            skip=[a-b for a,b in zip(q,W)]
            gate(all(a>=0 for a in skip),'actual skip coefficient dominance',(h,xx))
            cq=response_characteristic(q,p)
            cs=response_characteristic(skip,p)
            gate(all(a>0 for a in cq),'actual carried response negative at all first roots',(h,xx))
            gate(all(a>0 for a in cs),'actual skip response negative at all first roots',(h,xx))
            first=literal(h,xx,g);second=literal(h,xx,2*g)
            fcounts=[(xx+j,3*h-3*j,1+2*j) for j in range(h+1)]
            qcounts=[(2*xx+e,6*h-3*e,2+2*e) for e in range(-1,2*h+1)]
            gate(set(first)==set(fcounts) and set(second)==set(qcounts),
                 'named complete literal fibres',(h,xx))
            gate([first[c]/comb(g,2*h+1) for c in fcounts]==p,
                 'named literal first coefficient row',(h,xx))
            gate([second[c]/K for c in qcounts]==q,
                 'named literal full carried coefficient row',(h,xx))
            controls.append(dict(h=h,x=xx,g=g,first_channels=h+1,
                doubled_channels=2*h+2,negative_first_root_values=h,
                full_characteristic=[str(c) for c in cq],
                skip_characteristic=[str(c) for c in cs]))
    for row in controls:
        print('GENUINE_CONTROL '+json.dumps(row,separators=(',',':')))
    print('CONTROL_SCOPE eight named actual coefficient rows, including nonprimitive rows; no all-h sign claim')


def main():
    symbolic_step_checks()
    model_obstruction()
    genuine_named_controls()
    print('PASS explicit_gates='+str(GATES))
    print('SEMANTIC_SHA256 '+sha256(json.dumps(TRACE,separators=(',',':')).encode()).hexdigest())


if __name__=='__main__':
    main()
