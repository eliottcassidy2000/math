#!/usr/bin/env python3
"""Independent supersingular checks via Hasse coefficients and elimination.
No imports from producer; no Fp2 enumeration; no class-number code.
"""
from math import comb, isqrt
import json
import sys
from pathlib import Path

GATES=0

def need(ok, msg):
    global GATES
    GATES+=1
    if not ok: raise RuntimeError(msg)


def trim(a,p):
    a=[x%p for x in a]
    while a and a[-1]==0:a.pop()
    return a

def add(a,b,p):
    c=[0]*max(len(a),len(b))
    for i,x in enumerate(a):c[i]+=x
    for i,x in enumerate(b):c[i]+=x
    return trim(c,p)

def scale(a,k,p):return trim([k*x for x in a],p)

def mul(a,b,p):
    c=[0]*(len(a)+len(b)-1) if a and b else []
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return trim(c,p)

def divmodp(a,b,p):
    need(bool(b),'polynomial division by zero')
    a=trim(a,p);b=trim(b,p);q=[0]*max(0,len(a)-len(b)+1)
    ib=pow(b[-1],-1,p)
    while len(a)>=len(b):
        j=len(a)-len(b); k=a[-1]*ib%p;q[j]=k
        a=add(a,[0]*j+scale(b,-k,p),p)
    return trim(q,p),a

def monic(a,p):return scale(a,pow(a[-1],-1,p),p) if a else []

def gcdp(a,b,p):
    while b:a,b=b,divmodp(a,b,p)[1]
    return monic(a,p)

def exactdiv(a,b,p):
    q,r=divmodp(a,b,p);need(not r,'not an exact polynomial quotient');return q

def deriv(a,p):return trim([i*a[i] for i in range(1,len(a))],p)

def power(a,k,f,p):
    z=[1]
    while k:
        if k&1:z=divmodp(mul(z,a,p),f,p)[1]
        a=divmodp(mul(a,a,p),f,p)[1];k//=2
    return z

def evaluate(a,x,p):
    v=0
    for k in reversed(a):v=(v*x+k)%p
    return v


def supersingular_j(p):
    m=(p-1)//2;out=[]
    for j in range(p):
        if j==0:a,b=0,1
        elif j==1728%p:a,b=1,0
        else:
            k=j*pow((1728-j)%p,-1,p)%p;a,b=3*k%p,2*k%p
        h=0
        for cubics in range(m+1):
            linears=p-1-3*cubics
            constants=m-cubics-linears
            if linears>=0 and constants>=0:
                h+=comb(m,cubics)*comb(m-cubics,linears)*pow(a,linears,p)*pow(b,constants,p)
        if h%p==0:out.append(j)
    return out


def eliminated_polynomial(p,j):
    # Remainder of (u+256)^3-j*u^2 modulo u^2-s*u+4096 is A(s)u+B(s).
    A=[192512,768-j,1]
    B=[4096*(j+3328),-4096]
    # Norm(Au+B)=B^2+s*A*B+4096*A^2.
    return monic(add(add(mul(B,B,p),[0]+mul(A,B,p),p),scale(mul(A,A,p),4096,p),p),p)


def level_two(p,ss):
    f=[1]
    for j in ss:
        h=eliminated_polynomial(p,j)
        need(len(h)>1,'constant elimination polynomial')
        dh=deriv(h,p)
        need(bool(dh),'inseparable cubic needs separate handling')
        h=exactdiv(h,gcdp(h,dh,p),p)
        f=monic(mul(f,exactdiv(h,gcdp(f,h,p),p),p),p)
    rational_poly=gcdp(f,add(power([0,1],p,f,p),[0,-1],p),p)
    all_degree2=gcdp(f,add(power([0,1],p*p,f,p),[0,-1],p),p)
    need(all_degree2==f,'some level-two points require larger than Fp2')
    need(gcdp(f,deriv(f,p),p)==[1],'combined polynomial not squarefree')
    rational=[s for s in range(p) if evaluate(f,s,p)==0]
    need(len(rational)==len(rational_poly)-1,'rational root count')
    remaining=monic(exactdiv(f,rational_poly,p),p)
    need((len(f)-len(rational_poly))%2==0,'quadratic orbit parity')
    return {'p':p,'supersingular_j':ss,'squarefree_polynomial_low_first':f,
            'total_points':len(f)-1,'rational_points':rational,
            'nonrational_factor_low_first':remaining,'minus_dimension':(len(f)-len(rational_poly))//2}


def main():
    if '--test-failing-gate' in sys.argv:need(False,'intentional failure control')
    # The producer receipt is used solely as a comparison target after independent computation.
    baseline=Path(__file__).resolve().parent.parent/'moonshine'/'results.json'
    source=json.loads(baseline.read_text())
    point_results={}
    for p in range(5,200):
        if any(p%d==0 for d in range(2,isqrt(p)+1)):continue
        ss=supersingular_j(p)
        need(ss==source['direct_pointcount_ss_j'][str(p)],f'coefficient vs point-count mismatch p={p}')
        point_results[p]=ss
    rows=[]
    for p in [5,7,11,13,17,19,23,29,31,41,47,59,71]:
        row=level_two(p,point_results[p]);rows.append(row)
        producer=next(x for x in source['level_two'] if x['p']==p)
        need(row['total_points']==producer['total_points'],f'all roots p={p}')
        need(row['rational_points']==producer['rational_points'],f'rational roots p={p}')
        need(row['minus_dimension']==producer['minus_dimension'],f'odd rank p={p}')
        if row['minus_dimension']:
            need(row['nonrational_factor_low_first']==list(reversed(producer['irreducible_quadratics'][0])),f'quadratic p={p}')
    need({r['p'] for r in rows if r['minus_dimension']}=={29,41,59,71},'four-prime odd shell')
    print(json.dumps({'status':'PASS','gates':GATES,'point_count_comparisons':len(point_results),
        'point_scope':'every prime 5 <= p <= 199','level_two_scope':'13 Ogg primes >= 5',
        'method':'Hasse coefficients and squarefree polynomial resultant gcd; no Fp2 enumeration or class numbers',
        'level_two':rows},sort_keys=True,indent=2))

if __name__=='__main__':main()
