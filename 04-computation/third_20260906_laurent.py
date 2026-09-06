#!/usr/bin/env python3
"""All-height forced-factor gates, endpoint39 certificate, and PSD hostile.

No repository producer imports. All load-bearing arithmetic is rational.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import comb, factorial, gcd
import json
import sympy as S

x=S.Symbol('x')
ZERO=S.Poly(0,x,domain=S.QQ)
ONE=S.Poly(1,x,domain=S.QQ)
GATES=0
TRACE=sha256()

def gate(ok,label,data=None):
    global GATES
    if not ok:raise RuntimeError(f'{label}: {data}')
    GATES+=1
    TRACE.update((label+':'+repr(data)+'\n').encode())

def ff(a,n):
    p=ONE
    for j in range(n):p*=S.Poly(a-j,x,domain=S.QQ)
    return p

def rem(a,p):
    a=list(a)
    while len(a)>=len(p):
        v=a.pop()
        for j,pj in enumerate(p[:-1]):a[len(a)-len(p)+1+j]-=v*pj
    return a+[ZERO]*(len(p)-1-len(a))

def rows(h):
    p=[ff(x+h,h-j).mul_ground(S.Rational(factorial(2*h+1),factorial(3*h-3*j)*factorial(1+2*j))) for j in range(h+1)]
    q=[ff(2*x+2*h,2*h-e).mul_ground(S.Rational(1,factorial(6*h-3*e)*factorial(2+2*e))) for e in range(-1,2*h+1)]
    carry=S.Poly(x,x,domain=S.QQ).mul_ground(S.Rational(2**(h+1)*factorial(3*h),factorial(2*h+1)*factorial(6*h+3)))
    for j in range(h):carry*=S.Poly(2*x+2*j+1,x,domain=S.QQ)
    gate(q[0]==carry*p[0],'generic-carry-cancellation',h)
    r=rem(q[1:],p)
    r=[r[j]-carry*p[j+1] for j in range(h)]
    for j,a in enumerate(r):gate(a.degree()<=2*h-j,'generic-response-degree',(h,j))
    return p,q,r

def mm(a,b):
    n=len(a)
    return [[sum((a[i][j]*b[j][k] for j in range(n)),ZERO) for k in range(n)] for i in range(n)]

def charpoly(p,r):
    h=len(r)
    cols=[rem([ZERO]*j+r,p) for j in range(h)]
    mat=[[cols[j][i] for j in range(h)] for i in range(h)]
    power=[[ONE if i==j else ZERO for j in range(h)] for i in range(h)]
    traces=[]
    for _ in range(h):
        power=mm(power,mat);traces.append(sum((power[i][i] for i in range(h)),ZERO))
    c=[ONE]
    for k in range(1,h+1):
        c.append(sum((c[k-i]*traces[i-1] for i in range(1,k+1)),ZERO).mul_ground(S.Rational(-1,k)))
        gate(c[k].degree()<=2*h*k,'generic-characteristic-degree',(h,k))
    total=[[ONE if i==j else ZERO for j in range(h)] for i in range(h)]
    for k in range(1,h+1):
        total=mm(total,mat)
        for i in range(h):total[i][i]+=c[k]
    gate(all(v.is_zero for row in total for v in row),'Cayley-Hamilton',h)
    return c

def exponent(h,k,r):return max(0,min(k-h+r,r-1))

def forced_factor(h,k):
    d=ONE
    for r in range(2,h+1):d*=S.Poly((x+r)**exponent(h,k,r),x,domain=S.QQ)
    return d

def boundary_gates(h,p,q,c):
    for r in range(1,h+1):
        gate(all(p[j].eval(-r)==0 for j in range(r)) and p[r].eval(-r)!=0,'zero-root-cluster',(h,r))
        gate(p[0].diff().eval(-r)!=0 and q[0].diff().eval(-r)!=0,'simple-parameter-defects',(h,r))
        gate(all(q[e+1].eval(-r)==0 for e in range(2*r)),'regular-response-boundary',(h,r))
        gate(q[2*r+1].eval(-r)!=0,'regular-response-first-survivor',(h,r))
    for k,a in enumerate(c[1:],1):
        d=forced_factor(h,k)
        b=a.exquo(d)
        for r in range(2,h+1):
            # Exact orders are FINITE-EXACT controls, not the all-h theorem.
            gate(b.eval(-r)!=0,'finite-exact-forced-orders',(h,k,r))
    print('BOUNDARY_CONTROL h='+str(h)+' exact_orders_match=True')

def literal_controls(p,q):
    h=6;primitive=0
    for xx in range(1,25):
        g=xx+19
        first=[F(factorial(g),factorial(xx+j)*factorial(18-3*j)*factorial(1+2*j)) for j in range(7)]
        gate(first==[comb(g,13)*F(a.eval(xx)) for a in p],'literal-first-row',xx)
        K=factorial(2*g)//factorial(2*g-26)
        second=[F(factorial(2*g),factorial(2*xx+e)*factorial(36-3*e)*factorial(2+2*e)) for e in range(-1,13)]
        gate(second==[K*F(a.eval(xx)) for a in q],'literal-full-doubled-row',xx)
        for j in range(7):
            counts=(xx+j,18-3*j,1+2*j)
            gate(sum(counts)==g and -39*counts[0]+(2*g-39)*counts[1]+(3*g-39)*counts[2]==0,'literal-first-counts',(xx,j))
        gate(second[0]>0 and len(second)==14,'actual-lower-carry',xx)
        primitive+=gcd(g,39)==1
    source=(-39,3,24)
    first=next(m for m in range(1,22) if any(source[0]*a+source[1]*b+source[2]*(m-a-b)==0 for a in range(m+1) for b in range(m-a+1)))
    gate(first==7,'g21-gcd-hostile',first)
    print('LITERAL_CONTROLS x=1..24 primitive_rows='+str(primitive)+' gcd_hostile_g21_first_admissible_mass=7')

def psd_hostile(p,r):
    h=2
    derivative=[p[j+1].mul_ground(j+1) for j in range(h)]
    product=[ZERO]*(2*h-1)
    for i,a in enumerate(derivative):
        for j,b in enumerate(r):product[i+j]-=a*b
    q=rem(product,p)
    be=[[ZERO for _ in range(h)] for _ in range(h)]
    for i in range(h+1):
        for j in range(h):
            if i>j:
                for a in range(i-j):be[i-1-a][j+a]+=p[i]*q[j]
            elif i<j:
                for a in range(j-i):be[j-1-a][i+a]-=p[i]*q[j]
    # Normalized matrix entry w^(6-i-j) B_ij(1/w), at w^4.
    matrix=S.Matrix([[be[i][j].nth(2-i-j) for j in range(h)] for i in range(h)])
    denominator=S.ilcm(*[S.denom(a) for a in matrix]);integer=matrix*denominator
    content=S.igcd(*[int(a) for a in integer]);primitive=integer/content
    expected=S.Matrix([[49138290,357162621],[357162621,2402157940]])
    gate(primitive==expected and primitive.det()==-9527204358067041,'coefficient-PSD-hostile')
    print('PSD_HOSTILE h=2 power_w=4 primitive_matrix=[[49138290,357162621],[357162621,2402157940]] determinant=-9527204358067041')

def main():
    saved={}
    for h in range(2,7):
        p,q,r=rows(h);c=charpoly(p,r)
        boundary_gates(h,p,q,c);saved[h]=(p,q,r,c)
    p,q,r,c=saved[6]
    count=0
    for k,a in enumerate(c[1:],1):
        d=forced_factor(6,k);b=a.exquo(d)
        content,primitive=b.primitive();coefficients=list(reversed(primitive.all_coeffs()))
        gate(content>0,'positive-deflated-content',k)
        for j,v in enumerate(coefficients):gate(v>0,'positive-deflated-coefficient',(k,j));count+=1
        gate(b.degree()==12*k-d.degree(),'deflated-degree',k)
        cert={'h':6,'k':k,'degree':b.degree(),'factor_exponents':[[r,exponent(6,k,r)] for r in range(2,7) if exponent(6,k,r)],'content':str(content),'coefficients_ascending':[int(v) for v in coefficients]}
        print('DEFLATED_CERTIFICATE '+json.dumps(cert,separators=(',',':')))
    gate(count==208,'certificate-coefficient-count',count)
    literal_controls(p,q)
    pp,qq,rr,cc=saved[2]
    psd_hostile(pp,rr)
    no_carry=charpoly(pp,rem(qq[1:],pp))[-1]
    gate(no_carry.rem(S.Poly((x+1)*(x+2)**2,x)).is_zero,'omitted-carry-extra-forced-factor')
    gate(cc[-1].eval(-1)!=0 and cc[-1].exquo(S.Poly(x+2,x)).eval(-2)!=0,'retained-carry-removes-extra-orders')
    print('CARRY_CONTROL h2 no_carry_norm_divisible_by=(x+1)*(x+2)^2 actual_norm_has_orders=(0,1)')
    print('CLAIM endpoint39 g=x+19>=20 gcd(g,39)=1 first_nonzero_mass=g_or_2g; structural_divisibility_all_h; exact_factor_orders_only_h2..6')
    print('PASS explicit_gates='+str(GATES)+' semantic_sha256='+TRACE.hexdigest())

if __name__=='__main__':main()
