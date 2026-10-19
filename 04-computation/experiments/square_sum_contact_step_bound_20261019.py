#!/usr/bin/env python3
"""Exact finite probe of the specified elliptic contact-step search.

No external packages. Explicit checks remain enabled under python -O.
The two acceptance modes are deliberately kept separate: strict u>1 versus
all positive-spread real charts, including -1<u<0 via t=-1/u.
This census is finite; the companion note proves unbounded search indices.
"""
from collections import Counter
from fractions import Fraction as F
from math import gcd, isqrt
import json
import sys

sys.set_int_max_str_digits(0)
CHECKS = 0


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def spread(t):
    return 4*t*(t*t-1)/(t*t+1)**2


def petals(bound):
    out = []
    for n in range(2, isqrt(bound)+1):
        for m in range(1, n):
            a,b,c = n*n-m*m, 2*n*m, n*n+m*m
            if c <= bound and gcd(m,n) == 1 and (m+n)%2 and c*c > 4*a*b:
                out.append((a,b,c,n,m))
    return sorted(out, key=lambda row: (row[2],row[0],row[1]))


def add(P,Q,k):
    if P is None: return Q
    if Q is None: return P
    x,y = P
    u,v = Q
    if x == u and y == -v: return None
    slope = (3*x*x-2*k*k*x-1)/(2*y) if P == Q else (v-y)/(u-x)
    xx = slope*slope+k*k-x-u
    return xx, -y+slope*(x-xx)


def neg(P):
    return None if P is None else (P[0],-P[1])


def on_curve(Q,k):
    return Q is None or Q[1]**2 == (Q[0]-k*k)*(Q[0]**2-1)


def anchor(t,k):
    x = (t*t+4*t+1)/(t*t+1)
    y = 8*t*(t+1)/(t*t+1)**2
    check(on_curve((x,y),k), 'anchor on curve')
    return x,y


def inverse(Q,k):
    if Q is None: return None,'O'
    x,y = Q
    if x == 0: return None,'X=0'
    if x == k*k: return None,'X=k^2'
    u = (k*y+x-k*k)/(x*(x-k*k))
    z = (2*u*u*x-k*k-2*u-k*k*u*u)/k
    check(z*z == k*k*(u*u+1)**2-4*u*(u*u-1), 'inverse quartic equation')
    return (u,z),None


def normalize(t):
    check(t>1, 'normalized parameter positive')
    n,m = t.numerator,t.denominator
    divisor = 2 if (n+m)%2 == 0 else 1
    a,b,c = (n*n-m*m)//divisor,2*n*m//divisor,(n*n+m*m)//divisor
    if a%2 == 0: a,b = b,a
    check(a>0 and b>0 and a*a+b*b == c*c, 'positive PPT')
    check(gcd(a,b)==gcd(b,c)==gcd(a,c)==1, 'primitive PPT')
    check(a%2==1 and b%2==0, 'parity convention')
    check(F(2*a*b,c*c)==spread(t), 'normalization preserves spread')
    return (a,b,c),divisor


def probe(t,cap,all_charts):
    check(t>1 and 0<spread(t)<F(1,2), 'admissible anchor')
    k=(t*t+2*t-1)/(t*t+1)
    f0=spread(t)
    Q=anchor(t,k)
    P=(F(0),k)
    R=None
    exceptions=[]
    for q in range(1,cap+1):
        R=add(R,P,k)
        check(on_curve(R,k), 'multiple on curve')
        found=[]
        for sign in (1,-1):
            S=add(Q,R if sign==1 else neg(R),k)
            check(on_curve(S,k), 'translate on curve')
            pair,exception=inverse(S,k)
            if exception:
                exceptions.append((sign*q,exception))
                continue
            u,z=pair
            r=z/(u*u+1)
            f=spread(u)
            check(r*r==1+f0-f, 'exact contact equation before normalization')
            negative_chart=all_charts and -1<u<0
            t1=-1/u if negative_chart else u
            if t1<=1 or not(0<f<F(1,2)) or f==f0:
                continue
            check(spread(t1)==f, 'negative chart preserves spread')
            abc,divisor=normalize(t1)
            found.append({'q':sign*q,'negative_chart':negative_chart,
                          'side':'up' if f>f0 else 'down',
                          'normalization_divisor':divisor,
                          'neighbor_c_digits':len(str(abc[2]))})
        if found:
            return {'anchor_t':str(t),'first_abs_q':q,'found':found,'exceptions':exceptions}
    return {'anchor_t':str(t),'first_abs_q':None,'tested_through':cap,'found':[],'exceptions':exceptions}


def near_wall(m,lower):
    # f increases on (1,4/3); decreases on (7,8).
    lo,hi=(m,4*m//3) if lower else (7*m,8*m)
    while hi-lo>1:
        mid=(lo+hi)//2
        if (spread(F(mid,m))<F(1,2)) == lower:
            lo=mid
        else:
            hi=mid
    t=F(lo if lower else hi,m)
    check(0<spread(t)<F(1,2), 'wall anchor admissible')
    return t


def run():
    universe=petals(10000)
    check(len(universe)==528, 'complete c<=10000 anchor count')
    results={}
    for all_charts,cap,label in [(False,22,'strict_positive_chart'),(True,11,'all_positive_spread_charts')]:
        rows=[]
        for a,b,c,n,m in universe:
            row=probe(F(n,m),cap,all_charts)
            row['ppt']=[a,b,c]
            rows.append(row)
        check(all(row['first_abs_q'] is not None for row in rows), label+' census resolved')
        check(not any(row['exceptions'] for row in rows), label+' no inverse-chart exceptions')
        histogram=Counter(row['first_abs_q'] for row in rows)
        results[label]=rows
        print(label, 'histogram', json.dumps(dict(sorted(histogram.items()))))
        top=max(histogram)
        print(label, 'hardest', json.dumps([row for row in rows if row['first_abs_q']==top]))
        print(label, 'odd/odd successes',sum(r['normalization_divisor']==2 for row in rows for r in row['found']))
    check(normalize(F(9,7))[0]==normalize(F(8))[0], 'odd/odd swapped-leg anchor regression')
    selected=[]
    for t,cap in [(F(100),20),(F(1000),20),(F(10000),20),
                  (F(10**6),20),(F(10**9),20),(F(10**12),20),
                  (F(10**6+1,10**6),20),
                  (F(10**18),12),(F(10**30),12),
                  (F(10**12+1,10**12),12),(F(10**30+1,10**30),12)]:
        row=probe(t,cap,True)
        row['kind']='parabolic_spine'
        selected.append(row)
    for digits in (4,8,12):
        for lower in (True,False):
            row=probe(near_wall(10**digits,lower),12,True)
            row['kind']='lower_wall' if lower else 'upper_wall'
            selected.append(row)
    check(not any(row['exceptions'] for row in selected), 'selected no inverse-chart exceptions')
    results['selected_larger_anchors']=selected
    print('SELECTED',json.dumps(selected))
    # Exact Pell certificates for both singular-limit families.
    p0,p1=0,1
    for q in range(1,51):
        p2=2*p1+p0
        U=F(p2,p1)
        check(F(2)<=U<=F(5,2), 'large-t positive Pell limit cone')
        check(spread(U)==1-F(1,(p2*p2+p1*p1)**2), 'Pell spread deficit')
        check(spread(U)>=F(24,25), 'large-t Pell cone exclusion')
        if q>=2:
            V=-F(p0,p1)
            check(spread(V)==spread(-1/V), 'negative Pell chart equivalence')
            check(spread(V)>=F(24,25), 'large-t negative Pell cone exclusion')
        p0,p1=p1,p2
    # t -> 1 limits from the node's other component.
    plus,minus=F(1),F(1)
    for q in range(1,51):
        plus,minus=1/(plus-2),2+1/minus
        check(F(2)<minus<=3 and spread(minus)>=F(24,25), 'small-t minus limit excluded')
        if q>=2:
            check(-1<plus<0 and spread(plus)>=F(24,25), 'small-t plus limit excluded')
    print('EXPLICIT CHECKS',CHECKS)
    return results

if __name__=='__main__':
    data=run()
    if len(sys.argv)==2:
        with open(sys.argv[1],'w') as stream:
            json.dump(data,stream,indent=2)
            stream.write('\n')
