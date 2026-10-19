#!/usr/bin/env python3
"""Exact controls for the guarded compactification note; no network or assertions."""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import gcd

checks=0

def need(test,label):
    global checks
    checks+=1
    if not test:
        raise RuntimeError(label)

def branch(word):
    a,b,M=1,0,1
    for z in word:
        if z=='1':
            a,b=3*a,3*b+M
        elif z!='0':
            raise ValueError(z)
        M*=2
    return a,b,M

def step(x):
    return (3*x+1)//2 if x&1 else x//2

def replay(x,word):
    for bit in word:
        if x%2!=int(bit):
            return False,x
        x=step(x)
    return True,x

def parity(x,n):
    w=''
    for _ in range(n):
        w+=str(x%2)
        x=step(x)
    return w,x

def frac(x):
    return f'{x.numerator}/{x.denominator}'

def v2(x):
    x=abs(x)
    if x==0:
        return None
    return (x & -x).bit_length()-1

def prime(n):
    if n<2:return False
    p=2
    while p*p<=n:
        if n%p==0:return False
        p+=1
    return True

# Every binary word through length 10, its canonical repeated cylinders,
# three positive lifts, and an adjacent nonmember. Direct integer replay is
# separate from the closed repeated-word residue formula.
digest=sha256(); word_count=contracting_count=repeat_lifts=zero_returns=0
for n in range(1,11):
    for code in range(1<<n):
        w=format(code,f'0{n}b');a,b,M=branch(w);h=w.count('1')
        need(a==3**h and M==2**n,'word coefficient reconstruction')
        need(F(b,M)<=F(a-1,4),'universal carry bound')
        word_count+=1
        if a<M:
            contracting_count+=1
            need(b<M*M,'contracting carry below source modulus')
        for k in (1,2,3):
            Mk=M**k
            r=(-b*pow(a-M,-1,Mk))%Mk
            aa,bb,MM=branch(w*k)
            need(MM==Mk and aa==a**k,'repeated coefficient')
            need((a-M)*bb==b*(a**k-M**k),'repeated carry identity')
            need(r==(-bb*pow(aa,-1,MM))%MM,'independent repeated source residue')
            for u in range(3):
                x=r+u*Mk
                if x==0:continue
                ok,y=replay(x,w*k)
                need(ok,'canonical lift direct parity replay')
                need(y==F(aa*x+bb,MM),'repeated endpoint equality')
                ok,z=replay(x,w)
                need(ok and (z-x)%M**(k-1)==0,'displacement divisibility')
                mismatch=(a-M)*x+b
                need(mismatch%Mk==0,'repetition fuel divisibility')
                if mismatch:
                    need(v2(mismatch)//n>=k,'valuation repetition budget')
                if a<M and k>=2:
                    need(z<=x,'square-guard descent or exact return')
                    if z==x:zero_returns+=1
                if a>M and k>=2:
                    need((a-M)*x+b>=Mk,'expanding crossover exclusion')
                    if b>0 and b<Mk:
                        need(F(b,(a-M)*x)<=F(b,Mk-b),'scale-ratio exclusion')
                repeat_lifts+=1
                digest.update(f'{w},{k},{x},{z},{y}\n'.encode())
            need(not replay(r+1,w*k)[0],'adjacent source is not the guarded cylinder')

# Genuine returns and a paradoxical hostile whose parity does not repeat.
need(replay(1,'1010')==(True,1),'root-cycle return guard')
need(replay(2,'0101')==(True,2),'other root-cycle return guard')
wp,yp=parity(7,8)
ap,bp,Mp=branch(wp)
need((wp,yp)==('11101001',8),'positive overshoot hostile')
need(ap<Mp and not replay(7,wp*2)[0],'overshoot lacks the repeated guard')

K='110110110110110101011011101010'
P=3**19;Q=2**30;C0=4912087879;D=P-Q;c=F(C0,D)
need(branch(K)==(P,C0,Q),'K branch')
need(D==193*458651,'phase denominator factorization')
need(prime(193) and prime(458651) and prime(9173),'phase factor primality')
need(gcd(C0,D)==1,'unit initial phase')
order=146768
need(order==2**4*9173,'order factorization')
need(pow(Q,order,D)==1,'claimed modular order is a multiple')
need(pow(Q,order//2,D)!=1 and pow(Q,order//9173,D)!=1,'claimed modular order is minimal')
phases=[]
for L in range(1,13):
    mod=Q**L;r=(-C0*pow(D,-1,mod))%mod
    num=D*r+C0
    need(num%mod==0,'phase scaling identity')
    t=num//mod
    need(t%D==C0*pow(Q,-L,D)%D,'phase modular recurrence value')
    if L>=2:need(1<=t<D,'canonical compact phase interval')
    if L>=3:need(t==phases[-1]['phase']*pow(Q,-1,D)%D,'phase consecutive recurrence')
    need(F(r,mod)==F(t,D)-c/mod,'compact real scaling formula')
    need(D*r>=mod-C0,'uniform minimal source bound')
    need(replay(r,K*L)[0],'K source replay')
    phases.append({'L':L,'source_residue':r,'phase':t})

# Formal near-neutral maps versus exact legal-repeat crossover obstruction.
neutral=[]
for L,j in [(11,3),(40,11),(69,19)]:
    w=K*L+'10'*j;a,b,M=branch(w)
    s=F(a,M);e=F(b,M);mu=F(3,4)
    need(s==F(P,Q)**L*mu**j,'near-neutral exact slope')
    need(e==c*s+1-(c+1)*mu**j,'near-neutral exact carry')
    need(1<s<F(101,100),'declared expanding near-neutral examples')
    need(0<e<64 and M>2*e,'uniform example carry control')
    need(e/(M-e)<F(1,2**(len(w)-7)),'explicit nonlinear-ratio upper bound')
    r=(-b*pow(a-M,-1,M*M))%(M*M)
    need(replay(r,w*2)[0],'large near-neutral square guard')
    need(F((a-M)*r+b,M)>=M,'near-neutral guarded displacement gap')
    for t in [F(1,10),F(1,100),F(1,1000)]:
        need(1/(s/t+e)==t/(s+e*t),'reciprocal conjugacy')
        t0=t/(1+(c+1)*t)
        need(1/((c+1)*t0)==1/((c+1)*t)+1,'exact parabolic translation coordinate')
    neutral.append({'L':L,'F_count':j,'Terras_length':len(w),
                    'slope_minus_one':frac(s-1), 'carry':frac(e),
                    'nonlinear_ratio_less_than_power_of_two':7-len(w),
                    'affine_data_sha256':sha256(f'{a},{b},{M}'.encode()).hexdigest()})

# Rational-anchor tests: exact declared examples only, not a census theorem.
exceptions=[14196,15364,18436,49324,128308,215276,266092,288932,295660,
            310916,390388,394148,480032,522404,529652]
anchor=[]
for d in [4,12852]+exceptions:
    z=D*d-C0;a,b,M=1,0,1;w=''
    for q in range(1,10001):
        bit=z%2;w+=str(bit)
        if bit:
            z=(3*z+D)//2;a,b=3*a,3*b+M
        else:z//=2
        M*=2
        if a<M:break
    else:raise RuntimeError('declared rational control has not stopped')
    need(branch(w)==(a,b,M),'independent rational parity branch')
    R=F(a*d+b,M-a)
    L=(q+29)//30
    while Q**L<=D*R+C0:L+=1
    mod=Q**L;r=(-C0*pow(D,-1,mod))%mod
    need(30*L>=q and mod>D*R+C0,'quantitative guard-depth criterion')
    need(r>R,'minimal lift exceeds original-source threshold')
    for u in (0,1,3):
        x=r+u*mod
        ok,y=replay(x+d,w)
        need(ok and y<x,'rational-anchor exact lifted repair')
    anchor.append({'d':d,'coefficient_stop':q,'sufficient_K_blocks':L,'threshold':frac(R)})
need(anchor[0]['coefficient_stop']==92 and anchor[0]['sufficient_K_blocks']==4,'known d=4 control')
need(anchor[1]['coefficient_stop']==135 and anchor[1]['sufficient_K_blocks']==5,'refined d=12852 control')

# Endpoint filtering identities on all words through length 7 and offsets 1..8.
filter_checks=0
for n in range(1,8):
    for code in range(1<<n):
        w=format(code,f'0{n}b');a,b,M=branch(w)
        if a>=M:continue
        for d in range(1,9):
            num=b-M*d;den=M-a
            candidate=num>0 and num%den==0
            if candidate:
                x=num//den
                ok,y=replay(x,w)
                need(ok and y==x+d,'exact positive endpoint-source filter')
            filter_checks+=1

# No globally uniform repair clock: explicit all-odd shifted cylinders.
for d in (1,4,12852):
    for J in range(1,65):
        q=d+2;x=(1<<J)*q-d-1
        word,y=parity(x+d,J)
        need(word=='1'*J and y>x,'all-odd boundary hostile')

print(json.dumps({'status':'PASS','explicit_checks':checks,
 'scope':'scoped symbolic proofs; finite exact controls do not assert pointwise Collatz coverage',
 'all_words_max_length':10,'word_count':word_count,
 'contracting_word_count':contracting_count,'positive_repeated_lifts':repeat_lifts,
 'zero_return_controls':zero_returns,'repeat_controls_digest':digest.hexdigest(),
 'phase_modulus':D,'phase_period':order,'phase_controls':phases,
 'near_neutral_examples':neutral,'rational_anchor_controls':anchor,
 'endpoint_filter_word_max_length':7,'endpoint_filter_word_offset_pairs':filter_checks,
 'overshoot_hostile':{'source':7,'word':wp,'endpoint':yp,'square_guard':False}},
 indent=2,sort_keys=True))
