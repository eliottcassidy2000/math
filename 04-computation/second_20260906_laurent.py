#!/usr/bin/env python3
"""Endpoint-33 actual carried return certificate and observer hostiles.

Run from the repository root. Uses exact rational polynomial arithmetic;
no repository mathematical producer or numerical root locator is imported.
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
    if not ok: raise RuntimeError(f'{label}: {data}')
    GATES+=1
    TRACE.update((label+':'+repr(data)+'\n').encode())

def fall(a,n):
    z=ONE
    for j in range(n): z*=S.Poly(a-j,x,domain=S.QQ)
    return z

def rem(a,p):
    a=list(a)
    while len(a)>=len(p):
        v=a.pop()
        for j,pj in enumerate(p[:-1]): a[len(a)-len(p)+1+j]-=v*pj
    return a+[ZERO]*(len(p)-1-len(a))

def symbolic_rows():
    p=[fall(x+5,5-j).mul_ground(S.Rational(factorial(11),factorial(15-3*j)*factorial(1+2*j))) for j in range(6)]
    q=[fall(2*x+10,10-e).mul_ground(S.Rational(1,factorial(30-3*e)*factorial(2+2*e))) for e in range(-1,11)]
    carry=fall(x,1).mul_ground(S.Rational(64*factorial(15),factorial(33)*factorial(11)))
    for j in range(5): carry*=S.Poly(2*x+2*j+1,x,domain=S.QQ)
    gate(q[0]==carry*p[0],'exact-carry-denominator-cancellation')
    r=rem(q[1:],p)
    r=[r[j]-carry*p[j+1] for j in range(5)]
    for j,a in enumerate(r): gate(a.degree()<=10-j,'response-weight-bound',j)
    return p,q,r

def mm(a,b):
    return [[sum((a[i][j]*b[j][k] for j in range(5)),ZERO) for k in range(5)] for i in range(5)]

def characteristic(p,r):
    cols=[rem([ZERO]*j+r,p) for j in range(5)]
    mat=[[cols[j][i] for j in range(5)] for i in range(5)]
    powers=[[ONE if i==j else ZERO for j in range(5)] for i in range(5)]
    traces=[]
    for _ in range(5):
        powers=mm(powers,mat)
        traces.append(sum((powers[i][i] for i in range(5)),ZERO))
    cs=[ONE]
    for k in range(1,6):
        cs.append(sum((cs[k-i]*traces[i-1] for i in range(1,k+1)),ZERO).mul_ground(S.Rational(-1,k)))
        gate(cs[k].degree()==10*k,'characteristic-degree',k)
    # Cayley-Hamilton as a second symbolic consequence-object identity.
    total=[[ONE if i==j else ZERO for j in range(5)] for i in range(5)]
    for k in range(1,6):
        total=mm(total,mat)
        for i in range(5): total[i][i]+=cs[k]
    gate(all(a.is_zero for row in total for a in row),'symbolic-Cayley-Hamilton')
    return cs

def numeric_rem(a,p):
    a=list(map(F,a));p=list(map(F,p))
    while len(a)>=len(p):
        v=a.pop()/p[-1]
        for j,pj in enumerate(p[:-1]):a[len(a)-len(p)+1+j]-=v*pj
    return a+[F(0)]*(len(p)-1-len(a))

def literal_controls(p,q):
    primitive=0
    for xx in range(1,25):
        g=xx+16
        pp=[F(a.eval(xx)) for a in p]
        scale=comb(g,11)
        first=[F(factorial(g),factorial(xx+j)*factorial(15-3*j)*factorial(1+2*j)) for j in range(6)]
        gate(first==[scale*a for a in pp],'literal-first-row',xx)
        K=factorial(2*g)//factorial(2*g-22)
        second=[F(factorial(2*g),factorial(2*xx+e)*factorial(30-3*e)*factorial(2+2*e)) for e in range(-1,11)]
        gate(second==[K*F(a.eval(xx)) for a in q],'literal-full-doubled-row',xx)
        gate(second[0]>0 and len(second)==12,'retained-lower-carry',xx)
        for j in range(6):
            counts=(xx+j,15-3*j,1+2*j)
            gate(sum(counts)==g and -33*counts[0]+(2*g-33)*counts[1]+(3*g-33)*counts[2]==0,'first-count-charge',(xx,j))
        primitive+=gcd(g,33)==1
    # Omission of gcd(g,33)=1 does not preserve the first mass.
    support=(-33,3,21)
    reachable=[]
    for m in range(1,19):
        if any(support[0]*a+support[1]*b+support[2]*(m-a-b)==0 for a in range(m+1) for b in range(m-a+1)):
            reachable.append(m)
    gate(reachable[0]==6,'gcd-hostile-g18-first-mass6',reachable)
    print('LITERAL_CONTROLS x=1..24 indexed_rows=24 primitive_rows='+str(primitive)+' gcd_hostile_g18_first_mass=6')

def cone_hostile():
    # H=(1+u)^8*(u^6+7t*u^4+10t^2*u^2+t^3), k=5.
    p=[F(1),F(10),F(1)]
    hc=[]
    for j in range(11):
        a=[F(0)]*4
        for l,b in enumerate((1,7,10,1)):
            jj=j-6+2*l
            if 0<=jj<=8:a[l]=F(b*comb(8,jj))
        hc.append(a)
    def convolution(a,b):
        out=[F(0)]*(len(a)+len(b)-1)
        for i,u in enumerate(a):
            for j,v in enumerate(b):out[i+j]+=u*v
        return out
    products=[numeric_rem(convolution(hc[j],hc[10-j]),p) for j in range(11)]
    rays=[]
    for r in range(5):
        def ff(n):return factorial(n)//factorial(n-r) if n>=r else 0
        z=[sum(ff(j)*ff(10-j)*products[j][l] for j in range(11)) for l in range(2)]
        gate(all(a.denominator==1 for a in z),'cone-integral-ray',r)
        z=list(map(int,z));gg=gcd(*z);rays.append([a//gg for a in z])
    actual=[F(0)]*7
    for e in range(-1,5):actual[e+2]=F(factorial(16),factorial(2+e)*factorial(12-3*e)*factorial(2+2*e))
    target=numeric_rem(actual,p);gg=gcd(*map(int,target));target=[int(a)//gg for a in target]
    det=lambda a,b:a[0]*b[1]-a[1]*b[0]
    tests=[det(rays[0],a) for a in rays]
    rejected=det(rays[0],target)
    gate(all(a>=0 for a in tests) and rejected<0,'actual-row-outside-complete-derivative-cone',(tests,rejected))
    print('CONE_HOSTILE primitive_rays='+json.dumps(rays,separators=(',',':'))+' target='+json.dumps(target,separators=(',',':')))
    print('CONE_SEPARATOR values='+str(tests)+' actual='+str(rejected))

def phase_pade_hostile(p,q):
    # Exact bounded 5x5 solve at x=10000. All sums have fixed length.
    xx=10000;g=xx+16;qq=xx+5
    pp=[F(a.eval(xx)) for a in p]
    qraw=[F(a.eval(xx)) for a in q]
    qr=numeric_rem(qraw[1:],pp)
    carry=qraw[0]/pp[0]
    qr=[qr[j]-carry*pp[j+1] for j in range(5)]
    beta=[F(comb(qq+2*a,3*a)) for a in range(12)]
    wr=[]
    for e in range(-1,11):
        ss=10-e
        wr.append(F(comb(2*g,2+2*e))*sum(beta[a]*beta[ss-a] for a in range(ss+1)))
    w=numeric_rem(wr[1:],pp);carry=wr[0]/pp[0]
    w=[w[j]-carry*pp[j+1] for j in range(5)]
    wa=[numeric_rem([F(0)]*(2*i)+w,pp) for i in range(3)]
    qa=[numeric_rem([F(0)]*(2*i)+qr,pp) for i in range(1,3)]
    A=S.Matrix([[wa[i][l] for i in range(3)]+[-qa[i][l] for i in range(2)] for l in range(5)])
    b=S.Matrix(qr)
    solution=A.inv()*b
    gate(A*solution==b,'phase-Pade-exact-identity')
    signs=[int(S.sign(a)) for a in solution]
    gate(signs==[1,-1,-1,-1,-1],'square-phase-Pade-scale-hostile',signs)
    print('PADE_HOSTILE x=10000 coefficients_order=B0,B1,B2,A1,A2 signs='+str(signs))

def main():
    p,q,r=symbolic_rows()
    cs=characteristic(p,r)
    certificates=[]
    for k in range(1,6):
        shifted=S.Poly(cs[k].as_expr().subs(x,x+1),x,domain=S.QQ)
        content,primitive=shifted.primitive()
        coefficients=list(reversed(primitive.all_coeffs()))
        gate(content>0,'positive-characteristic-content',k)
        for i,a in enumerate(coefficients):gate(a>0,'positive-shifted-characteristic-coefficient',(k,i))
        certificate={'k':k,'degree':10*k,'content':str(content),'coefficients_ascending':[int(a) for a in coefficients]}
        certificates.append(certificate)
        print('CHAR_CERTIFICATE '+json.dumps(certificate,separators=(',',':')))
    print('CHAR_SCHEMA variable=v=x-1; c_k(x)=content*sum(coefficients_ascending[j]*(x-1)^j)')
    literal_controls(p,q)
    cone_hostile()
    phase_pade_hostile(p,q)
    print('CLAIM endpoint33 g=x+16>=17 gcd(g,33)=1; first_nonzero_mass=g_or_2g; all_5_canceling_phases_detect_at_2g')
    print('PASS explicit_gates='+str(GATES)+' semantic_sha256='+TRACE.hexdigest())

if __name__=='__main__':main()
