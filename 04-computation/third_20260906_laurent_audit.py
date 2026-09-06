#!/usr/bin/env python3
"""Independent endpoint-39 audit: literal fibers, resultants, 58-node interpolation.

No producer or repository implementation imports. SymPy computes only exact
rational resultants; finite differences independently recover power coefficients.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import comb, factorial, gcd
from pathlib import Path
import json
import sympy as S

T,Z=S.symbols('T Z')
GATES=0
TRACE=sha256()
ROOT=Path(__file__).resolve().parents[1]

def gate(ok,label,data=None):
    global GATES
    if not ok:raise RuntimeError(f'{label}: {data!r}')
    GATES+=1
    TRACE.update((label+':'+repr(data)+'\n').encode())

def ff(a,n):
    v=1
    for i in range(n):v*=a-i
    return v

def derivative_ff(a,n,multiplier=1):
    return multiplier*sum(ff(a,i)*ff(a-i-1,n-i-1) for i in range(n))

def boundary_controls():
    count=0
    for h in range(1,10):
        for r in range(1,h+1):
            a=h-r;b=2*h-2*r
            gate(all(ff(a,h-j)==0 for j in range(r)) and ff(a,h-r)!=0,'small-zero-root-cluster',(h,r))
            gate(derivative_ff(a,h)!=0 and derivative_ff(b,2*h+1,2)!=0,'simple-parameter-zeros',(h,r))
            gate(all(ff(b,2*h-e)==0 for e in range(2*r)) and ff(b,2*h-2*r)!=0,'regular-response-orders',(h,r))
            for j in range(1,r+1):
                gate((j*(r-1)+r-1)//r==min(j,r-1),'integral-order-rounding',(h,r,j))
            count+=1
        for k in range(1,h+1):
            exponents=[max(0,min(k-h+r,r-1)) for r in range(2,h+1)]
            expected=k*(k+1)//2 if k<h else h*(h-1)//2
            gate(sum(exponents)==expected,'factor-degree',(h,k))
    print('BOUNDARY_AUDIT h=1..9 clusters='+str(count)+' includes_r1=True')

def fibers(g,m):
    rows=[]
    for c in range(m+1):
        numerator=39*m-3*g*c
        if numerator%(2*g):continue
        b=numerator//(2*g);a=m-b-c
        if min(a,b,c)<0:continue
        rows.append(((a,b,c),F(factorial(m),factorial(a)*factorial(b)*factorial(c))))
    return rows

def node(x):
    g=x+19;p={};q={}
    for counts,v in fibers(g,g):
        a,b,c=counts
        gate((c-1)%2==0,'first-phase-integrality',(x,counts))
        p[(c-1)//2]=v/comb(g,13)
    for counts,v in fibers(g,2*g):
        a,b,c=counts
        gate((c-2)%2==0,'second-phase-integrality',(x,counts))
        q[(c-2)//2]=v/F(factorial(2*g),factorial(2*g-26))
    gate(set(p)==set(range(7)) and p[6]==1,'complete-monic-first-fiber',x)
    gate(set(q)==set(range(-1,13)) and q[-1]>0,'complete-carried-second-fiber',x)
    P=S.Poly.from_dict({(j,):S.Rational(v.numerator,v.denominator) for j,v in p.items()},T)
    # Multiplication by T clears the actual Laurent carry. The resultant
    # is product_i[T_i*(Z-q(T_i))]; degree six gives product_i T_i=P(0).
    Q=S.Poly.from_dict({(j+1,):S.Rational(v.numerator,v.denominator) for j,v in q.items()},T)
    quotient_response=Q.rem(P)
    result=S.Poly(S.resultant(P.as_expr(),Z*T-quotient_response.as_expr(),T)/P.TC(),Z)
    coeffs=[F(a) for a in result.all_coeffs()]
    gate(result.degree()==6 and coeffs[0]==1,'normalized-resultant',x)
    return coeffs[1:]

def factor(k,x):
    answer=1
    for r in range(2,7):answer*=(x+r)**max(0,min(k-6+r,r-1))
    return answer

def interpolate(values):
    differences=list(values);basis=[F(1)];answer=[F(0)]*len(values)
    for j in range(len(values)):
        for i,v in enumerate(basis):answer[i]+=differences[0]*v
        differences=[b-a for a,b in zip(differences,differences[1:])]
        if j+1<len(values):
            updated=[F(0)]*(len(basis)+1)
            for i,v in enumerate(basis):
                updated[i]-=v
                updated[i+1]+=v/F(j+1)
            # Current basis is binomial(x-1,j), next is basis*(x-j-1)/(j+1).
            basis=updated
    while answer and answer[-1]==0:answer.pop()
    return answer

def evaluate(poly,x):
    value=F(0)
    for c in reversed(poly):value=value*x+c
    return value

def certificate_audit():
    values=[node(x) for x in range(1,59)]
    polynomials=[]
    for k in range(1,7):
        poly=interpolate([row[k-1]/factor(k,x) for x,row in enumerate(values,1)])
        degree=12*k-(k*(k+1)//2 if k<6 else 15)
        gate(len(poly)==degree+1,'deflated-degree',(k,degree))
        for i,c in enumerate(poly):gate(c>0,'positive-power-coefficient',(k,i))
        for x,row in enumerate(values,1):gate(evaluate(poly,x)*factor(k,x)==row[k-1],'interpolation-replay',(k,x))
        polynomials.append(poly)
    for x in (59,71,83):
        row=node(x)
        for k,poly in enumerate(polynomials,1):gate(evaluate(poly,x)*factor(k,x)==row[k-1],'independent-heldout-parameter',(k,x))
    frozen=(ROOT/'05-knowledge/results/third_20260906_laurent.out').read_bytes()
    gate(sha256(frozen).hexdigest()=='4020a649ea09986521314e5a7d36c1b8df6f9e8883e881ebb23f89487c635221','frozen-certificate-hash')
    records=[json.loads(line.split(' ',1)[1]) for line in frozen.decode().splitlines() if line.startswith('DEFLATED_CERTIFICATE ')]
    gate(len(records)==6,'six-certificate-records')
    for k,(poly,record) in enumerate(zip(polynomials,records),1):
        gate(record['h']==6 and record['k']==k,'certificate-identity',k)
        expected=[F(record['content'])*v for v in record['coefficients_ascending']]
        gate(poly==expected,'independent-all-coefficients-match',k)
        exponents=[[r,max(0,min(k-6+r,r-1))] for r in range(2,7) if max(0,min(k-6+r,r-1))]
        gate(record['factor_exponents']==exponents,'independent-factor-exponents',k)
        payload=json.dumps([str(v) for v in poly],separators=(',',':')).encode()
        print('RECONSTRUCTED k='+str(k)+' degree='+str(len(poly)-1)+' positive_coefficients='+str(len(poly))+' rational_coefficient_sha256='+sha256(payload).hexdigest())
    gate(sum(map(len,polynomials))==208,'all208-coefficients')
    print('INTERPOLATION_AUDIT nodes=1..58 heldout=59,71,83 route=literal_charge_fibers_then_resultant_then_finite_differences all208_match=True')

def admissibility_controls():
    for g in (20,22,23,25,38,40,41,58,61,77):
        gate(gcd(g,39)==1,'primitive-family-control',g)
        first=next(m for m in range(1,g+1) if fibers(g,m))
        gate(first==g,'literal-first-admissible-mass',g)
    first=next(m for m in range(1,22) if fibers(21,m))
    gate(first==7,'g21-first-admissible-mass-hostile',first)
    print('ADMISSIBILITY_AUDIT ten_primitive_parameters=True nonprimitive_g21_first_mass=7')

def quadratic_controls():
    x=S.Symbol('x');h=2
    p=sum(S.Rational(factorial(5),factorial(6-3*j)*factorial(1+2*j))*S.prod(x+2-i for i in range(2-j))*T**j for j in range(3))
    regular=sum(S.Rational(1,factorial(12-3*e)*factorial(2+2*e))*S.prod(2*x+4-i for i in range(4-e))*T**e for e in range(5))
    carry=S.prod(2*x+4-i for i in range(5))/factorial(15)
    inv=S.invert(T,p,T)
    response=S.rem(regular+carry*inv,p,T)
    actual=S.factor(S.resultant(p,Z-response,T))
    removed=S.factor(S.resultant(p,Z-regular,T))
    determinant=S.Poly(actual,Z).TC();no_carry=S.Poly(removed,Z).TC()
    gate(determinant.subs(x,-1)!=0,'retained-carry-xminus1')
    gate(S.cancel(determinant/(x+2)).subs(x,-2)!=0,'retained-carry-xminus2-simple')
    gate(S.rem(no_carry,(x+1)*(x+2)**2,x)==0,'deleted-carry-extra-forced-factor')
    y=S.Symbol('y')
    partner=S.rem(-S.diff(p,T)*response,p,T)
    bezout=S.cancel((p*partner.subs(T,y)-p.subs(T,y)*partner)/(T-y))
    bezout=S.Poly(bezout,T,y)
    matrix=S.Matrix([[S.expand(bezout.coeff_monomial(T**i*y**j)).coeff(x,2-i-j) for j in range(2)] for i in range(2)])
    scale=matrix[0,0]/49138290
    primitive=matrix/scale
    gate(scale>0 and primitive==S.Matrix([[49138290,357162621],[357162621,2402157940]]),'Bezout-coefficient-matrix')
    gate(primitive.det()==-9527204358067041,'Bezout-indefinite-coefficient')
    print('QUADRATIC_AUDIT actual_carry_orders=(0,1) removed_carry_forces=(1,2) coefficient_matrix_determinant=-9527204358067041')

if __name__=='__main__':
    boundary_controls();certificate_audit();admissibility_controls();quadratic_controls()
    print('STATUS: independent analytic acceptance plus exact 58-node reconstruction; all-height residual positivity remains OPEN')
    print('PASS explicit_gates='+str(GATES)+' semantic_sha256='+TRACE.hexdigest())
