#!/usr/bin/env python3
"""Exact finite atlas for a distinct plus-port Pythagorean contact.

The cusp chart is a rational section; the other chart searches even signed
translates through 256. SymPy is needed for symbolic identities. All checks
remain enabled with python -O. No floating point or external data are used.
"""
from collections import Counter
from fractions import Fraction as F
from math import gcd, isqrt
import json
import sys
import sympy as sp

sys.set_int_max_str_digits(0)
CHECKS = Counter()
ZCOEFFS = [1,824,-6469,21392,-26297,24216,-66331,100544,-7414,
           -132112,150510,-48288,-39122,44912,-11894,-8000,8293,
           -3368,743,-112,27,-8,1]


def check(condition, name, kind='exact_rational'):
    CHECKS[kind] += 1
    if not condition:
        raise RuntimeError(name)


def spread(t):
    return 4*t*(t*t-1)/(t*t+1)**2


def k_of(t):
    return (t*t+2*t-1)/(t*t+1)


def swap(t):
    return (t+1)/(t-1)


def components(t):
    A=t**4-8*t**3+2*t*t+1
    B=5*t**6-12*t**5+23*t**4-17*t*t+12*t-3
    C=2*t**5-t**4-4*t**3-2*t*t+2*t-1
    D=t**6+20*t**5-29*t**4+16*t**3+3*t*t-4*t+1
    return A,B,C,D


def horner(coeffs,t):
    out=0
    for coefficient in coeffs:
        out=out*t+coefficient
    return out


def section(t):
    A,B,C,D=components(t)
    N,E=t*A*B,C*D
    return N/E, horner(ZCOEFFS,t)/(E*E)


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


def anchor(t):
    return (t*t+4*t+1)/(t*t+1),8*t*(t+1)/(t*t+1)**2


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
    f0,f1=spread(t),spread(u)
    check(u>1 and 0<f1<F(1,2),'admissible output')
    check(f1!=f0,'changed spread')
    r=z/(u*u+1)
    check(r*r==1+f0-f1,'plus-minus square contact')
    old,_=normalize(t)
    new,d=normalize(u)
    check(new!=old,'distinct primitive petal')
    return new,d,abs(r)


def construct(t):
    check(t>1 and 0<spread(t)<F(1,2),'admissible input')
    t=min(t,swap(t))
    check(F(1)<t<F(4,3),'low-cone representative')
    if t<=F(17,16):
        u,z=section(t)
        A,B,C,D=components(t)
        check(A<0 and B>0 and C<0 and D>0,'local denominator signs')
        check(t<u<F(5,4),'local strict interval')
        k=k_of(t);P=(F(0),k);Q=anchor(t)
        # Two independent group-law expressions for the new section.
        R1=add(add(add(Q,Q,k),Q,k),add(P,P,k),k)
        V=add(Q,P,k);R2=add(Q,add(V,V,k),k)
        check(R1==R2 and on_curve(R1,k),'3Q+2P section group identity')
        check(inverse(R1,k)==(u,z),'section inverse formula')
        new,d,r=certify_output(t,u,z)
        check(new[2] <= 776516*(2*normalize(t)[0][2])**11,
              'explicit cusp hypotenuse upper bound')
        return {'chart':'3Q+2P','t':str(t),'u':str(u),'triple':list(new),
                'divisor':d,'root':str(r),'q':None}
    f0=spread(t)
    check(f0>=F(35904,297025)>F(1,10),'noncusp spread lower bound')
    J=0
    while 2**J*f0<2+f0: J+=1
    budget=2*(100*(4+J)//7)
    check(J<=5 and budget<=256,'uniform exact fallback budget')
    k=k_of(t);P=(F(0),k);Q=anchor(t);A=add(P,P,k);R=None
    for j in range(1,budget//2+1):
        R=add(R,A,k)
        check(on_curve(R,k),'even translate on curve')
        for sigma in (1,-1):
            pair=inverse(add(Q,R if sigma==1 else neg(R),k),k)
            if pair is None: continue
            u,z=pair
            if u<=1 or z<=0: continue
            f1=spread(u)
            if not(0<f1<F(1,2)) or f1==f0: continue
            new,d,r=certify_output(t,u,z)
            return {'chart':'Q+qP','t':str(t),'u':str(u),'triple':list(new),
                    'divisor':d,'root':str(r),'q':2*sigma*j,
                    'J':J,'budget':budget}
    raise RuntimeError('certified fallback budget exhausted')


def taylor_rat(expr,t,h,degree):
    n,d=sp.fraction(sp.cancel(expr))
    nc=sp.Poly(n.subs(t,1+h),h)
    dc=sp.Poly(d.subs(t,1+h),h)
    coeff=[]
    for j in range(degree+1):
        coeff.append(sp.cancel((nc.nth(j)-sum(dc.nth(i)*coeff[j-i]
                      for i in range(1,j+1)))/dc.nth(0)))
    return coeff


def symbolic_controls():
    t,h=sp.symbols('t h')
    A,B,C,D=components(t);N,E=sp.expand(t*A*B),sp.expand(C*D)
    k=k_of(t);u=N/E;Z=horner(ZCOEFFS,t)
    def identity(x,name):
        check(sp.cancel(x)==0,name,'exact_symbolic')
    identity(k*k-1-spread(t),'anchor square identity')
    identity(spread(swap(t))-spread(t),'leg-swap spread')
    identity(k_of(swap(t))-k,'leg-swap k')
    identity((t*t+1)**2*Z**2-(t*t+2*t-1)**2*(N*N+E*E)**2
             +4*N*E*(N*N-E*E)*(t*t+1)**2,'polynomial contact certificate')
    # Derive the section from actual rational-function group operations.
    P=(sp.Integer(0),k);Q=anchor(t)
    def sadd(V,W): return tuple(sp.factor(x) for x in add(V,W,k))
    V=sadd(Q,P)
    Vexpected=((t*t-4*t+1)/(t*t+1),-8*t*t*(t-1)/(t*t+1)**2)
    for actual,expected in zip(V,Vexpected): identity(actual-expected,'Q+P formula')
    S=sadd(V,V)
    SX=(t*t+1)*(5*t*t-4*t+1)/(4*t*t*(t-1)**2)
    identity(S[0]-SX,'2(Q+P) x formula')
    R=sadd(Q,S);x,y=R
    uu=sp.factor((k*y+x-k*k)/(x*(x-k*k)))
    zz=sp.factor((2*uu*uu*x-k*k-2*uu-k*k*uu*uu)/k)
    identity(uu-u,'section parameter formula')
    identity(zz-Z/E**2,'section z formula')
    check(taylor_rat(u,t,h,3)==[1,3,0,2],'u Taylor coefficients','exact_symbolic')
    check(taylor_rat(spread(u),t,h,3)==[0,6,-9,-23],
          'output spread Taylor coefficients','exact_symbolic')
    check(taylor_rat(spread(t),t,h,3)==[0,2,-1,-1],
          'input spread Taylor coefficients','exact_symbolic')
    check(sp.limit((t-1)**2*S[0],t,1)==1,'S specializes to O','exact_symbolic')
    check(sp.limit((-S[0]/S[1])/(t-1),t,1)==-1,
          'nonzero formal leading coefficient','exact_symbolic')
    for name,poly in [('B',B),('D',D)]:
        check(all(x>0 for x in sp.Poly(poly.subs(t,1+h),h).all_coeffs()),
              name+' shifted positive coefficients','exact_symbolic')
    identity(A.subs(t,h+1)-(h**4-4*h**3-16*h*h-16*h-4),'A shift')
    identity(C.subs(t,h+1)-(2*h**5+9*h**4+12*h**3-8*h-4),'C shift')
    check(2*F(1,16)**5+9*F(1,16)**4+12*F(1,16)**3<1,'C positive tail below 1')
    L=sp.Poly(sp.cancel(-(N-t*E)/(t-1)).subs(t,h+1),h)
    check(all(x>0 for x in L.all_coeffs()),'u-t sign polynomial','exact_symbolic')
    G=sp.Poly((4*N-5*E).subs(t,h+1),h)
    check(G.nth(0)==32 and G.nth(11)==10,'G endpoint coefficients','exact_symbolic')
    check(all(G.nth(i)<0 for i in range(1,11)),'G negative tail','exact_symbolic')
    negtail=sum(F(-int(G.nth(i)),16**i) for i in range(1,11))
    check(negtail<20,'uniform G positive certificate')
    check(sum(abs(x) for x in sp.Poly(N,t).all_coeffs())==760,
          'numerator coefficient norm','exact_symbolic')
    check(sum(abs(x) for x in sp.Poly(E,t).all_coeffs())==446,
          'denominator coefficient norm','exact_symbolic')
    check(spread(F(17,16))==F(35904,297025)>F(1,10),'rational cusp boundary')
    check(spread(F(5,4))==F(720,1681)<F(1,2),'strict admissibility cap')
    # The old chart's hostile fixed index q=0 is the anchor, while the new
    # section is outside the admissible region at this larger parameter.
    bad_u,_=section(F(5,4))
    check(not(0<spread(bad_u)<F(1,2)),'hostile nonlocal section control')


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
