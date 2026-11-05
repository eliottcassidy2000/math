#!/usr/bin/env python3
"""All-modulus exterior banks and conditioned size/residue laws.
Standalone exact controls. No gates disappear under python -O.
"""
from fractions import Fraction
from itertools import product
from collections import Counter
from pathlib import Path
import hashlib,json,math
HERE=Path(__file__).resolve().parent
CHECKS={}
def require(ok,msg,key):
    if not ok: raise RuntimeError(msg)
    CHECKS[key]=CHECKS.get(key,0)+1

def bits(a,n): return tuple((a>>i)&1 for i in range(n-1,-1,-1))
def A(w):
    a=0
    for b in w:a=2*a+b
    return a

def C(w):return sum((1-b)<<i for i,b in enumerate(w))
def p(w):return w[-1] if w else 0
def eps(w):return int(bool(w) and w[0]==0)
def count(w):
    h,a,last=1,0,0
    for b in w:
        if not b:h+=a+last
        a=2*a+b;last=b
    return h

def slot(w):
    if not w:return 1
    return w[0]+1-w[-1]+sum(x==0 and y==1 for x,y in zip(w,w[1:]))+sum(1<<(i-j-1) for j in range(len(w)) for i in range(j+1,len(w)) if w[j] and not w[i])

def rows_one(w):
    n=len(w)+1;rows=[0]*n
    for i,b in enumerate(w,1):
        if b:rows[0]|=1<<i
        else:rows[i]|=1
        for j in range(i+1,n):rows[i]|=1<<j
    return rows

def rows_two(w):
    n=len(w)+2;rows=[0]*n;rows[0]|=2
    for i,c in enumerate(w,2):
        for a in (0,1):
            if (c>>a)&1:rows[a]|=1<<i
            else:rows[i]|=1<<a
        for j in range(i+1,n):rows[i]|=1<<j
    return rows

def graph_count(rows):
    n=len(rows);dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i in range(n):
            value=dp[mask][i]
            if value:
                more=rows[i]&~mask
                while more:
                    bit=more&-more;more-=bit
                    dp[mask|bit][bit.bit_length()-1]+=value
    return sum(dp[-1])

def strong(rows):
    n=len(rows)
    for i in range(n):
        seen=1<<i;old=-1
        while old!=seen:
            old=seen
            for j in range(n):
                if (seen>>j)&1:seen|=rows[j]
        if seen!=(1<<n)-1:return False
    return True

def parameters(M):
    odd=M;k=0
    while not odd&1:k+=1;odd//=2
    ell=max(1,k)
    while pow(2,ell,odd)!=1%odd:ell+=1
    g=math.gcd(2,M);T=M//g;Q=1<<ell
    return {'M':M,'k':k,'oddpart':odd,'ell':ell,'g':g,'T':T,'L':2*ell*T,'D':(Q-1)**2+1}

def bank(par,t):
    ell=par['ell'];T=par['T']
    return ((1,)*ell+(0,)*ell)*t+(0,)*(2*ell*(T-t))

def admissible(M):return list(range(M)) if M&1 else list(range(1,M,2))
def choose(par,u,w,r):
    M=par['M'];g=par['g'];T=par['T'];L=par['L'];D=par['D']
    K=count(u)+count(w)-1+(1<<L)*A(u)*C(w)
    d=(r-K)%M
    require(d%g==0,'invalid residue parity','decoder_parity')
    if T==1:t=1
    else:
        t=((d//g)*pow(D//g,-1,T))%T
        if t==0:t=T
    return u+bank(par,t)+w,t

def ends(n,L,h):
    a=(n-L)//2;b=n-L-a
    def root(power):return math.isqrt(h<<power) if power>=0 else math.isqrt(h>>(-power))
    aa=min((1<<a)-1,root(2*a-n));bb=min((1<<b)-1,root(2*b-n))
    return bits(aa,a),tuple(1-x for x in bits(bb,b)[::-1]),(1<<(n-a))+(1<<(n-b))+n//2+1

# Independent physical graph counters and the exact concatenation interface.
for n in range(9):
    for w in product((0,1),repeat=n):
        require(count(w)==slot(w),'slot count mismatch','slot')
        require(count(w)==graph_count(rows_one(w)),'physical one-apex count','graph_one')
        require(count(w)&1,'parity','oddness')
short=[w for n in range(4) for w in product((0,1),repeat=n)]
for u,v,w in product(short,repeat=3):
    require(count(u+v)==count(u)+count(v)-1+A(u)*C(v)+p(u)*eps(v),'two-word law','concat_two')
    if v:
        expected=count(u)+count(v)+count(w)-2+A(u)*C(v)+((1<<len(v))*A(u)+A(v))*C(w)+p(u)*eps(v)+p(v)*eps(w)
        require(count(u+v+w)==expected,'three-word law','concat_three')

bank_rows=[]
exteriors=[((),()),((1,),(0,)),((0,1,1),(1,0,0)),((1,0,1,1),(0,1,1,0))]
for M in range(1,65):
    par=parameters(M);L=par['L'];T=par['T'];D=par['D'];g=par['g']
    require(math.gcd(D,M)==g,'increment gcd','bank_gcd')
    vals=[];words=[]
    for t in range(1,T+1):
        v=bank(par,t);h=count(v);words.append(v);vals.append(h%M)
        require(len(v)==L and v[0]==1 and v[-1]==0,'bank physical type','bank_type')
        require(A(v)%M==C(v)%M==0,'boundary signature vanishing','bank_boundaries')
        require(h%M==(1+t*D)%M,'bank count progression','bank_counter')
        if M<=12:require(h==slot(v),'independent bank slot count','bank_slot')
        if L<=12:require(h==graph_count(rows_one(v)),'graph bank witness','bank_graph')
    require(len(set(words))==T,'bank members not distinct','bank_distinct')
    require(sorted(vals)==admissible(M),'empty-exterior onto','bank_empty_onto')
    for u,w in exteriors:
        K=count(u)+count(w)-1+(1<<L)*A(u)*C(w)
        observed=[]
        for t,v in enumerate(words,1):
            hv=count(u+v+w);observed.append(hv%M)
            require(hv%M==(K+t*D)%M,'universal exterior law','bank_exterior')
        require(sorted(observed)==admissible(M),'conditional uniform interface','bank_uniform')
        for r in admissible(M):
            value,t=choose(par,u,w,r)
            require(count(value)%M==r,'residue decoder','bank_decode')
    bank_rows.append({k:par[k] for k in ('M','k','oddpart','ell','T','L')}|{'residue_count':T,'increment_mod_M':D%M})

# A hostile for the tempting t=0 endpoint in an even modulus.
par=parameters(6);zero=(0,)*par['L']
require(C(zero)%6!=0,'t=0 hostile absent','zero_hostile')

# Exact simultaneous selectors. Big counts are handled as integers, not floats.
joint=[]
for M in list(range(1,25))+[31,32,45,64,65,81,97]:
    par=parameters(M);L=par['L']
    for extra in (2,3,12,25):
        n=L+extra
        for h in (1<<(n-2),3*(1<<(n-2)),(1<<n)-1,1<<n):
            u,w,error=ends(n,L,h)
            residues=admissible(M) if M<=24 else admissible(M)[::max(1,len(admissible(M))//7)]
            for r in residues:
                full,t=choose(par,u,w,r);H=count(full)
                require(len(full)==n and full[0]==1 and full[-1]==0,'strong joint type','joint_type')
                require(H%M==r and abs(H-h)<=error,'joint magnitude/residue','joint')
        joint.append({'M':M,'n':n,'L':L,'error_bits':error.bit_length()})

# Complete small uniform strong-word histograms, with exact TV inequalities.
histograms=[]
pars={M:parameters(M) for M in range(1,17)}
for n in range(2,17):
    hist={M:Counter() for M in pars};total=1<<(n-2)
    for inner in range(total):
        w=(1,)+bits(inner,n-2)+(0,);H=count(w)
        for M in hist:hist[M][H%M]+=1
    records=[]
    for M,par in pars.items():
        T=par['T'];L=par['L'];q=(n-2)//L
        tv=sum(abs(Fraction(hist[M][r],total)-Fraction(1,T)) for r in admissible(M))/2
        eta=(1-Fraction(T,1<<L))**q
        require(tv<=eta,'residue TV bound','mixing_histogram')
        if M in (3,4,5,6,8):records.append({'M':M,'TV':str(tv),'q':q})
    histograms.append({'n':n,'presentations':total,'diagnostics':records})

# Two-apex transfer: independent nonnegative run-slot implementation.
TRANS=[[[] for i in range(36)] for c in range(4)];ACCEPT=[]
for order in range(2):
    for before,middle,after in product(range(3),range(3),range(2)):
        i=18*order+6*before+2*middle+after
        if before!=2 and (middle==1 or (middle==0 and order==0)):ACCEPT.append(i)
        for c in range(4):
            first=(c>>order)&1;last=(c>>(1-order))&1
            TRANS[c][i].append(18*order+6*(2 if first else 1)+2*middle+after)
            if middle or first:TRANS[c][i].append(18*order+6*before+2*(2 if last else 1)+after)
            if after or last:TRANS[c][i].append(18*order+6*before+2*middle+1)
def two_count(w):
    state=[0]*36;state[0]=state[18]=1
    for c in w:
        nxt=[0]*36
        for i,a in enumerate(state):
            for j in TRANS[c][i]:nxt[j]+=a
        state=nxt
    return sum(state[i] for i in ACCEPT)
for n in range(6):
    for w in product(range(4),repeat=n):
        require(two_count(w)==graph_count(rows_two(w)),'physical two-apex count','graph_two')

def boundary_prefix(x,k):
    """Exactly represent x through k inverse branches of ratio 2/3."""
    out=[];tail=x
    for _ in range(k):
        bit=0 if tail<=Fraction(2,3) else 1
        out.append(1+2*bit);tail=(3*tail-bit)/2
        require(0<=tail<=1,'boundary inverse left its domain','boundary_inverse')
    return tuple(out),tail

def G(w):
    # Valid first-slot assignments with all three labels seen.
    state=[0]*8;state[0]=1
    for c in w:
        nxt=[0]*8
        for seen,value in enumerate(state):
            for j in range(3):
                if j and not ((seen>>j)&1) and not ((c>>(j-1))&1):continue
                nxt[seen|(1<<j)]+=value
        state=nxt
    return Fraction(state[7],3**len(w))
net=[]
for k in range(1,17):
    m=2*k+1
    for j in range(33):
        target=Fraction(j,16);x=target/2
        u,tail=boundary_prefix(x,k)
        # The surviving one-slot sum plus its exact tail is x.
        partial=sum(Fraction((c>>1)&1,3)*Fraction(2,3)**i for i,c in enumerate(u))
        require(partial+Fraction(2,3)**k*tail==x,'boundary representation','boundary_representation')
        require(0<=x-G(u)<=3*Fraction(2,3)**k,'coupled boundary G error','boundary_G')
        w=u+(2,)+(0,)*k;H=two_count(w)
        error=12*Fraction(2,3)**k
        require(w[0]!=0 and w[-1]!=3 and 2 in w,'strong two-apex word type','net_type')
        require(abs(Fraction(H,3**m)-target)<=error,'two-apex normalized net','two_net')
        if k<=4:
            require(H==graph_count(rows_two(w)),'graph-level net witness','net_graph')
            require(strong(rows_two(w)),'graph strongness','net_strong')
        if k in (4,8,16) and j in (0,8,16,24,32):
            net.append({'k':k,'m':m,'target':str(target),'word':''.join(map(str,w)),'H':H,'error':str(error)})

# r=1 strong support is not [0,1]; named H13 one-apex finite hostile.
require({count(w) for n in range(2,6) for w in product((0,1),repeat=n) if w[0] and not w[-1]}=={3,5,9,11,15,17,19,23},'small one-apex control','one_apex_control')

receipt={'status':'PASS','scope':'All-order proofs plus finite exact controls; no integer-support completion','checks':CHECKS,'arbitrary_modulus_banks':bank_rows,'simultaneous_tests':joint,'finite_mixing_diagnostics':histograms,'two_apex_net_examples':net,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':'PASS','checks':CHECKS,'receipt_sha256':hashlib.sha256((HERE/'receipt.json').read_bytes()).hexdigest()},sort_keys=True))
