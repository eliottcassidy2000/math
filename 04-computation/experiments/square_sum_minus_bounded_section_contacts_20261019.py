#!/usr/bin/env python3
"""Exact finite atlas for a distinct minus-port Pythagorean contact.

The cusp uses 3Q on the negative-z anchor sheet. Outside it, use the
positive-z anchor and the inherited even signed recurrence through 256.
SymPy is needed only for symbolic checks. All checks survive python -O.
"""
from collections import Counter
from fractions import Fraction as F
from math import gcd, isqrt
import json
import sys
import sympy as sp

sys.set_int_max_str_digits(0)
CHECKS = Counter()
WCOEFFS = [1,-8,27,-112,743,-3368,8293,-8000,-11894,44912,-39122,
           -48288,150510,-132112,-7414,100544,-66331,24216,-26297,
           21392,-6469,824,1]


def check(condition, name, kind='exact_rational'):
    CHECKS[kind] += 1
    if not condition:
        raise RuntimeError(name)


def spread(t):
    return 4*t*(t*t-1)/(t*t+1)**2


def k_of(t):
    # Positive on the selected low cone 1<t<4/3.
    return (-t*t+2*t+1)/(t*t+1)


def swap(t):
    return (t+1)/(t-1)


def components(t):
    A=t**4+2*t*t-8*t+1
    B=3*t**6-12*t**5+17*t**4-23*t*t+12*t-5
    C=t**5-2*t**4+2*t**3+4*t*t+t-2
    D=t**6-4*t**5+3*t**4+16*t**3-29*t*t+20*t+1
    return A,B,C,D


def horner(coeffs,t):
    out=0
    for coefficient in coeffs:
        out=out*t+coefficient
    return out


def section(t):
    A,B,C,D=components(t)
    N,E=C*D,A*B
    return -N/E, -horner(WCOEFFS,t)/(E*E)


def on_curve(R,k):
    return R is None or R[1]**2==(R[0]-k*k)*(R[0]**2-1)


def add(R,S,k):
    if R is None: return S
    if S is None: return R
    x,y=R;u,v=S
    if x==u and y==-v: return None
    slope=(3*x*x-2*k*k*x-1)/(2*y) if R==S else (v-y)/(u-x)
    xx=slope*slope+k*k-x-u
    return xx,-y+slope*(x-xx)


def neg(R):
    return None if R is None else (R[0],-R[1])


def anchor(t, sheet=-1):
    k=k_of(t);K=k*k;u=-t;z=sheet*(t*t+1)
    x=(k*z+K+2*u+K*u*u)/(2*u*u)
    return x,(x-K)*(x*u-1)/k


def inverse(R,k):
    if R is None: return None
    x,y=R
    if x==0 or x==k*k: return None
    u=(k*y+x-k*k)/(x*(x-k*k))
    z=(2*u*u*x-k*k-2*u-k*k*u*u)/k
    check(z*z==k*k*(u*u+1)**2-4*u*(u*u-1),'inverse quartic')
    return u,z


def normalize(t):
    check(t>1,'positive parameter')
    n,m=t.numerator,t.denominator
    d=2 if (n+m)%2==0 else 1
    a,b,c=(n*n-m*m)//d,2*n*m//d,(n*n+m*m)//d
    if a%2==0: a,b=b,a
    check(a>0 and b>0 and a*a+b*b==c*c,'positive PPT')
    check(gcd(a,b)==gcd(a,c)==gcd(b,c)==1,'primitive PPT')
    check(a%2==1 and b%2==0,'parity-oriented PPT')
    check(F(2*a*b,c*c)==spread(t),'normalization spread')
    return (a,b,c),d


def certify_output(t,u,z):
    v=-u;f0,f1=spread(t),spread(v)
    check(v>1 and 0<f1<F(1,2),'admissible output')
    check(f1!=f0,'changed spread')
    r=z/(v*v+1)
    check(r*r==1-f0+f1,'minus-plus square contact')
    old,_=normalize(t)
    new,d=normalize(v)
    check(new!=old,'distinct primitive petal')
    return new,d,abs(r)


def construct(t):
    check(t>1 and 0<spread(t)<F(1,2),'admissible input')
    t=min(t,swap(t))
    check(F(1)<t<F(4,3),'low-cone representative')
    if t<=F(17,16):
        u,z=section(t)
        A,B,C,D=components(t)
        check(A<0 and B<0 and C>0 and D>0,'local denominator signs')
        check(t<-u<F(5,4),'local strict interval')
        k=k_of(t);Q=anchor(t,-1)
        R=add(Q,add(Q,Q,k),k)
        check(on_curve(R,k),'3Q section on curve')
        check(inverse(R,k)==(u,z),'section inverse formula')
        new,d,r=certify_output(t,u,z)
        check(new[2] <= 776516*(2*normalize(t)[0][2])**11,
              'explicit cusp hypotenuse upper bound')
        return {'chart':'3Q_negative_sheet','t':str(t),'u':str(u),'triple':list(new),
                'divisor':d,'root':str(r),'q':None}
    f0=spread(t)
    check(f0>F(35904,297025)>F(1,10),'noncusp spread lower bound')
    J=0
    while 2**J*f0<2: J+=1
    budget=2*(100*(4+J)//7)
    check(J<=5 and budget<=256,'uniform exact fallback budget')
    k=k_of(t);P=(F(0),k);Q=anchor(t,+1);A=add(P,P,k);R=None
    for j in range(1,budget//2+1):
        R=add(R,A,k)
        check(on_curve(R,k),'even translate on curve')
        for sigma in (1,-1):
            pair=inverse(add(Q,R if sigma==1 else neg(R),k),k)
            if pair is None: continue
            u,z=pair
            if u>=-1 or z<=0: continue
            f1=spread(-u)
            if not(0<f1<F(1,2)) or f1==f0: continue
            new,d,r=certify_output(t,u,z)
            return {'chart':'Q_positive_sheet+qP','t':str(t),'u':str(u),
                    'triple':list(new),'divisor':d,'root':str(r),'q':2*sigma*j,
                    'J':J,'budget':budget}
    raise RuntimeError('certified fallback budget exhausted')


def taylor_rat(expr,t,h,degree):
    n,d=sp.fraction(sp.cancel(expr))
    nc=sp.Poly(n.subs(t,1+h),h);dc=sp.Poly(d.subs(t,1+h),h)
    coeff=[]
    for j in range(degree+1):
        coeff.append(sp.cancel((nc.nth(j)-sum(dc.nth(i)*coeff[j-i]
                      for i in range(1,j+1)))/dc.nth(0)))
    return coeff


def symbolic_controls():
    t,h=sp.symbols('t h')
    A,B,C,D=components(t);N,E=sp.expand(C*D),sp.expand(A*B)
    k=k_of(t);v=N/E;W=horner(WCOEFFS,t)
    def identity(x,name):
        check(sp.cancel(x)==0,name,'exact_symbolic')
    identity(k*k-1+spread(t),'anchor square identity')
    identity(spread(swap(t))-spread(t),'leg-swap spread')
    identity(k_of(swap(t))+k,'leg-swap signed k')
    # Reciprocal port duality with the companion plus section, checked
    # algebraically here rather than by importing another script.
    Ap=t**4-8*t**3+2*t*t+1
    Bp=5*t**6-12*t**5+23*t**4-17*t*t+12*t-3
    Cp=2*t**5-t**4-4*t**3-2*t*t+2*t-1
    Dp=t**6+20*t**5-29*t**4+16*t**3+3*t*t-4*t+1
    Np=t*Ap*Bp;Ep=Cp*Dp
    Zp=horner(list(reversed(WCOEFFS)),t)
    identity(t**11*Np.subs(t,1/t)+E,'reciprocal numerator reversal')
    identity(t**11*Ep.subs(t,1/t)+N,'reciprocal denominator reversal')
    identity(t**22*Zp.subs(t,1/t)-W,'reciprocal square reversal')
    identity(v*(Np/Ep).subs(t,1/t)-1,'reciprocal section identity')
    identity(spread(1/t)+spread(t),'reciprocal spread sign')
    identity(((t*t+2*t-1)/(t*t+1)).subs(t,1/t)-k,'reciprocal square root')
    xx,yy,kk,uu,zz=sp.symbols('xx yy kk uu zz')
    Fquart=lambda a,b: b*b-kk*kk*(a*a+1)**2+4*a*(a*a-1)
    identity(Fquart(-1/uu,-zz/(uu*uu))-Fquart(uu,zz)/uu**4,
             'quartic reciprocal involution')
    ui=(kk*yy+xx-kk*kk)/(xx*(xx-kk*kk))
    ximage=kk*kk+ui*ui*(kk*kk-xx)
    yimage=(ximage-kk*kk)*(-ximage/ui-1)/kk
    slope=(kk-yy)/xx
    xexpected=slope*slope+kk*kk-xx
    yexpected=yy+slope*(xx-xexpected)
    curve=sp.Poly(yy*yy-(xx-kk*kk)*(xx*xx-1),yy)
    for actual,expected in [(ximage,xexpected),(yimage,yexpected)]:
        numerator=sp.fraction(sp.factor(actual-expected))[0]
        remainder=sp.rem(sp.Poly(numerator,yy),curve).as_expr()
        check(sp.factor(remainder)==0,'cubic involution is -R-P','exact_symbolic')
    identity((t*t+1)**2*W**2-(t*t-2*t-1)**2*(N*N+E*E)**2
             -4*N*E*(N*N-E*E)*(t*t+1)**2,'polynomial contact certificate')
    # Derive the section from actual rational-function group operations.
    Q=anchor(t,-1)
    Qexpected=((t*t-4*t+1)/(t*t+1),-8*t*(t-1)/(t*t+1)**2)
    for actual,expected in zip(Q,Qexpected): identity(actual-expected,'negative-sheet Q formula')
    def sadd(V,W): return tuple(sp.factor(x) for x in add(V,W,k))
    S=sadd(Q,Q)
    identity(S[0]-(t*t+1)*(t*t-4*t+5)/(4*(t-1)**2),'2Q x formula')
    R=sadd(Q,S);x,y=R
    uu=sp.factor((k*y+x-k*k)/(x*(x-k*k)))
    zz=sp.factor((2*uu*uu*x-k*k-2*uu-k*k*uu*uu)/k)
    identity(uu+v,'section parameter formula')
    identity(zz+W/E**2,'section z formula')
    check(taylor_rat(v,t,h,3)==[1,3,6,14],'v Taylor coefficients','exact_symbolic')
    check(taylor_rat(spread(v),t,h,3)==[0,6,3,-35],
          'output spread Taylor coefficients','exact_symbolic')
    check(taylor_rat(spread(t),t,h,3)==[0,2,-1,-1],
          'input spread Taylor coefficients','exact_symbolic')
    check(sp.limit((t-1)**2*S[0],t,1)==1,'2Q specializes to O','exact_symbolic')
    EP=sp.Poly(E.subs(t,h+1),h)
    check([EP.nth(i) for i in range(11)]==[32,32,-80,-128,-16,48,40,64,50,18,3],
          'E shifted coefficients','exact_symbolic')
    etail=sum(F(-int(EP.nth(i)),16**i) for i in range(11) if EP.nth(i)<0)
    check(etail==F(1409,4096)<1,'uniform E positive certificate')
    NP=sp.Poly(N.subs(t,h+1),h)
    check(all(x>0 for x in NP.all_coeffs()),'N shifted positive coefficients','exact_symbolic')
    L=sp.Poly(sp.cancel((N-t*E)/(t-1)).subs(t,h+1),h)
    check([L.nth(i) for i in range(11)]==[64,256,480,384,96,0,-48,-96,-60,-16,-2],
          'v-t shifted coefficients','exact_symbolic')
    ltail=sum(F(-int(L.nth(i)),16**i) for i in range(11) if L.nth(i)<0)
    check(ltail==F(1777281,549755813888)<1,'uniform v-t positive certificate')
    G=sp.Poly((5*E-4*N).subs(t,h+1),h)
    check([G.nth(i) for i in range(12)]==[32,-352,-1232,-1728,-1040,-272,-152,96,178,58,-5,-4],
          'G shifted coefficients','exact_symbolic')
    gtail=sum(F(-int(G.nth(i)),16**i) for i in range(12) if G.nth(i)<0)
    check(gtail==F(119849021866005,4398046511104)<28,'uniform G positive certificate')
    check(sum(abs(x) for x in sp.Poly(N,t).all_coeffs())==446,
          'numerator coefficient norm','exact_symbolic')
    check(sum(abs(x) for x in sp.Poly(E,t).all_coeffs())==760,
          'denominator coefficient norm','exact_symbolic')
    check(spread(F(17,16))==F(35904,297025)>F(1,10),'rational cusp boundary')
    check(spread(F(5,4))==F(720,1681)<F(1,2),'strict admissibility cap')
    bad_u,_=section(F(5,4))
    check(not(0<spread(-bad_u)<F(1,2)),'hostile nonlocal section control')
    check(abs(k_of(F(17,16)))>0,'positive low-cone k')


def petals(bound):
    rows=[]
    for n in range(2,isqrt(bound)+1):
        for m in range(1,n):
            a,b,c=n*n-m*m,2*n*m,n*n+m*m
            if c<=bound and gcd(n,m)==1 and (n+m)%2 and c*c>4*a*b:
                rows.append((a,b,c,n,m))
    return sorted(rows,key=lambda row:(row[2],row[0],row[1]))


def main():
    symbolic_controls()
    universe=petals(10000)
    check(len(universe)==528,'complete c<=10000 census')
    hist=Counter();divs=Counter();qhist=Counter();maxbits=0
    for a,b,c,n,m in universe:
        row=construct(F(n,m))
        hist[row['chart']]+=1;divs[row['divisor']]+=1
        if row['q'] is not None: qhist[row['q']]+=1
        maxbits=max(maxbits,row['triple'][2].bit_length())
        check(normalize(F(n,m))[0]==(a,b,c),'census input identity')
        check(construct(swap(F(n,m)))==row,'census swapped-leg invariance')
    probes=[]
    for t in [F(17,16),F(65,64),F(101,100),F(1000001,1000000),
              F(10**30+1,10**30),F(10**30),F(5,4),F(9,7),F(8)]:
        row=construct(t)
        check(construct(swap(t))==row,'global swapped-leg invariance')
        probes.append({'original_t':str(t),'low_t':row['t'],'chart':row['chart'],
                       'q':row['q'],'divisor':row['divisor'],
                       'input_c_bits':normalize(t)[0][2].bit_length(),
                       'output_c_bits':row['triple'][2].bit_length(),
                       'spread_in':str(spread(t)),
                       'spread_out':str(F(2*row['triple'][0]*row['triple'][1],row['triple'][2]**2))})
    data={'status':'FINITE-EXACT controls; general theorem proved in companion note',
          'universe':'all primitive odd/even oriented admissible PPTs with c<=10000',
          'anchors':len(universe),'chart_counts':dict(hist),'normalization_divisors':dict(divs),
          'fallback_signed_index_histogram':dict(sorted(qhist.items())),
          'largest_output_c_bits_in_census':maxbits,'selected':probes,'checks':dict(CHECKS)}
    print(json.dumps(data,sort_keys=True,indent=2))
    if len(sys.argv)>1:
        with open(sys.argv[1],'w') as out:
            json.dump(data,out,sort_keys=True,indent=2);out.write('\n')


if __name__=='__main__': main()
