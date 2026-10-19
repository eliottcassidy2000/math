#!/usr/bin/env python3
"""Exact controls for the expanding-bank / reset-bank bounded-head theorem."""
from fractions import Fraction as F
from itertools import product
import json

COUNT = 0

def check(ok, message):
    global COUNT
    COUNT += 1
    if not ok:
        raise RuntimeError(message)

def ceil(x):
    return -(-x.numerator // x.denominator)

def coeff(bits):
    a, b, m = 1, 0, 1
    for bit in bits:
        if bit == '1':
            a, b = 3*a, 3*b+m
        elif bit != '0':
            raise ValueError('nonbinary word')
        m *= 2
    return a, b, m

def data(lo, hi):
    a,b,m=coeff(lo)
    check(coeff(hi)==(a,b+m-a,m), 'independent companion carry')
    rho=(-b*pow(a,-1,m))%m
    au,bu,mu=coeff(hi)
    check((-bu*pow(au,-1,mu))%mu==(rho+1)%m,'companion source guard')
    return dict(lo=lo,hi=hi,a=a,b=b,m=m,rho=rho)

def walk(x,bits):
    y=x
    for bit in bits:
        check(y%2==int(bit),'literal parity guard')
        y=(3*y+1)//2 if y%2 else y//2
    return y

KLO='110110110110110101011011101010'
KHI='001111111111111101010110000100'
BANK={
 'K':data(KLO,KHI),
 'V':data(KLO*4+'10',KHI*4+'01'),
 'F':data('10','01'),
 'H':data('110101000','001101001')}
C=max(F(BANK[k]['b'], BANK[k]['a']-BANK[k]['m']) for k in 'KV')
P=max(F(BANK[k]['b'], BANK[k]['m']-BANK[k]['a']) for k in 'FH')
check(BANK['K']['rho']>3*C+4*P,'uniform source floor for optional extra reset')
HEAD=ceil(C+P)
for k in 'KV':
    w=BANK[k]
    check(w['a']>w['m'],'expansion')
    for t in range(1,len(w['lo'])+1):
        aa,bb,mm=coeff(w['lo'][:t])
        check(aa>mm,'every lower prefix has expanding multiplier')
    c=F(w['b'],w['a']-w['m'])
    check(F(coeff(w['hi'])[1],w['a']-w['m'])==c-1,'upper expansion center shift')
for k in 'FH':
    w=BANK[k]
    check(w['a']<w['m'],'contraction')
    p=F(w['b'],w['m']-w['a'])
    check(F(coeff(w['hi'])[1],w['m']-w['a'])==p+1,'upper contraction center shift')
check(F(BANK['V']['b'],BANK['V']['a']-BANK['V']['m'])>F(BANK['K']['b'],BANK['K']['a']-BANK['K']['m']), 'genuinely different centers')

rows=[]
sequences=[''.join(z) for n in range(1,6) for z in product('KV',repeat=n)]
sequences += ['KV'*8,'VK'*8,'K'*32,'V'*16]
for seq in sequences:
    elo=''.join(BANK[k]['lo'] for k in seq)
    ehi=''.join(BANK[k]['hi'] for k in seq)
    ea,eb,em=coeff(elo)
    A,B=F(ea,em),F(eb,em)
    check(B<=C*(A-1),'expanding bank telescope')
    for pattern in ['F','H','FH','HF','FFH','HHF']:
        r,d=F(1),F(0)
        resets=''
        while r*A>=1:
            k=pattern[len(resets)%len(pattern)]
            w=BANK[k]
            u,v=F(w['a'],w['m']),F(w['b'],w['m'])
            r,d=u*r,u*d+v
            resets+=k
        check(d<=P*(1-r),'reset bank telescope')
        check(r*B+d<C+P,'uniform normalized intercept')
        lo=elo+''.join(BANK[k]['lo'] for k in resets)
        hi=ehi+''.join(BANK[k]['hi'] for k in resets)
        w=data(lo,hi)
        a,b,m,rho=(w[k] for k in ('a','b','m','rho'))
        delta=m-a
        check(delta>=1 and 0<rho<m-1,'integer slope gap and paired nonwrapping residues')
        check(F(b,m)==r*B+d,'exact composed intercept')
        tmin=max(0,(b-delta*rho)//(delta*m)+1)
        check(tmin<=HEAD,'uniform excluded head')
        # Exact rational endpoint test around cutoff, independently literal replay.
        for t in sorted(set([0,max(0,tmin-1),tmin,HEAD])):
            x=rho+t*m
            y=walk(x,lo)
            yu=walk(x+1,hi)
            check(y==(a*x+b)//m and yu==y+1,'paired exact endpoint')
            check((y<x)==(t>=tmin),'strict exact lift cutoff')
        for extra in 'FH':
            more=data(lo+BANK[extra]['lo'],hi+BANK[extra]['hi'])
            aa,bb,mm,rr=(more[k] for k in ('a','b','m','rho'))
            check((mm-aa)*rr>bb,'one extra reset removes all positive-lift exceptions')
            yy=walk(rr,more['lo'])
            check(walk(rr+1,more['hi'])==yy+1 and yy<rr,'extra reset literal paired receipt')
        rows.append(dict(expanding=seq,resets=resets,steps=len(lo),t_min=tmin))

# Sharp excluded-head equality from the existing one-letter theorem.
a,b,m=coeff('10');rho=(-b*pow(a,-1,m))%m
check(max(0,(b-(m-a)*rho)//((m-a)*m)+1)==1,'sharp ceiling control')
# Wrong-order control: arbitrarily large normalized intercept while slope<1.
wrong=[]
for h in (4,8,16,32):
    a,b,m=coeff('0'*h+'1'*h)
    check(F(a,m)<1 and F(b,m)==F(3,2)**h-1,'mixed order exact identity')
    check(F(b,m)>1,'mixed order violates C+P=1 intercept bound')
    wrong.append(dict(h=h,slope=str(F(a,m)),intercept=str(F(b,m))))
summary=dict(status='PROVED conditional cylinder theorem; no universal source entry',
 checks=COUNT,bank={k:{z:(str(v) if isinstance(v,F) else v) for z,v in w.items() if z not in ('lo','hi')} for k,w in BANK.items()},
 C=str(C),P=str(P),uniform_head=HEAD,cylinders=len(rows),max_observed_t_min=max(x['t_min'] for x in rows),
 sample_rows=rows[:6]+rows[-6:],wrong_order_controls=wrong)
print(json.dumps(summary,indent=2,sort_keys=True))
