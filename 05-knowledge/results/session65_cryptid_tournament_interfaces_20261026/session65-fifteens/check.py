#!/usr/bin/env python3
"""Exact bounded controls; no unbounded conjecture claimed."""
from collections import Counter
from itertools import combinations
from math import comb,gcd,isqrt
from pathlib import Path
import hashlib,json

def need(c,msg):
    if not c: raise RuntimeError(msg)
def ternary(n):
    if n==0:return '0'
    out=''
    while n:out=str(n%3)+out;n//=3
    return out
def primes(lim):
    return [p for p in range(2,lim+1) if all(p%d for d in range(2,isqrt(p)+1))]
def classno(D):
    forms=[]
    for a in range(1,isqrt((-D)//3)+1):
        for b in range(-a,a+1):
            if (b*b-D)%(4*a):continue
            c=(b*b-D)//(4*a)
            if a>c or gcd(gcd(a,abs(b)),c)>1:continue
            if (abs(b)==a or a==c) and b<0:continue
            forms.append((a,b,c))
    return len(forms)

def defect_distribution(n,b):
    sig=list(range(n))
    for j in range(b):sig[2*j],sig[2*j+1]=2*j+1,2*j
    pairs=list(combinations(range(n),2)); idx={e:i for i,e in enumerate(pairs)}
    trans=[]
    for x,y in pairs:
        u,v=sig[x],sig[y]
        trans.append((idx[min(u,v),max(u,v)],u>v))
    cnt=Counter()
    for mask in range(1<<len(pairs)):
        d=sum(((mask>>i)&1)!=(((mask>>j)&1)^flip) for i,(j,flip) in enumerate(trans))
        cnt[d]+=1
    f=n-2*b;m=b*(f+b-1)
    expected={b+2*j:(1<<(comb(f,2)+b+m))*comb(m,j) for j in range(m+1)}
    need(dict(cnt)==expected,('defect',n,b,cnt,expected))
    return {'n':n,'fixed':f,'pairs':b,'pair_orbits':m,'histogram':dict(sorted(cnt.items()))}

def field37():
    p=37;D=next(x for x in range(2,p) if pow(x,(p-1)//2,p)==p-1)
    def add(x,y):return ((x[0]+y[0])%p,(x[1]+y[1])%p)
    def neg(x):return ((-x[0])%p,(-x[1])%p)
    def mul(x,y):return ((x[0]*y[0]+D*x[1]*y[1])%p,(x[0]*y[1]+x[1]*y[0])%p)
    def inv(x):
        d=pow((x[0]*x[0]-D*x[1]*x[1])%p,-1,p)
        return (x[0]*d%p,-x[1]*d%p)
    def pw(x,n):
        a=(1,0)
        while n:
            if n&1:a=mul(a,x)
            x=mul(x,x);n//=2
        return a
    coeff=[comb((p-1)//2,j)**2%p for j in range((p-1)//2+1)]
    roots=[];js=set()
    for a in range(p):
        for b in range(p):
            x=(a,b);h=(0,0)
            for c in reversed(coeff):h=add(mul(h,x),(c,0))
            if h!=(0,0):continue
            roots.append(x);xm1=add(x,(-1,0));num=add(add(mul(x,x),neg(x)),(1,0))
            j=mul((256%p,0),mul(pw(num,3),inv(mul(pw(x,2),pw(xm1,2)))))
            js.add(j)
    js=sorted(js);sigma=[js.index(pw(j,p)) for j in js]
    rational=[a for a,b in js if b==0]
    direct=[]
    for j in range(p):
        if j==0:a,b=0,1
        elif j==1728%p:a,b=1,0
        else:
            k=j*pow((1728-j)%p,-1,p)%p;a,b=3*k%p,2*k%p
        character_sum=0
        for x in range(p):
            y=(x*x*x+a*x+b)%p;chi=pow(y,(p-1)//2,p)
            character_sum+=-1 if chi==p-1 else chi
        if character_sum==0:direct.append(j)
    need(len(js)==3 and len(rational)==1 and rational==direct,'p37 count')
    return {'p':p,'field':'F_p[t]/(t^2-D)','D':D,'Hasse_roots':roots,'j_invariants':js,'frobenius_permutation':sigma,'rational_j':rational,'independent_trace_zero_j':direct,'b':1}

O=[2,3,5,7,11,13,17,19,23,29,31,41,47,59,71]
rows=[]
for p in primes(211):
    if p<=3:s=r=1
    else:
        s=p//12+{1:0,5:1,7:1,11:2}[p%12]
        r=(classno(-4*p)+(classno(-p) if p%4==3 else 0))//2
    need(s>=r and (s-r)%2==0,'ss arithmetic')
    rows.append({'p':p,'s':s,'r':r,'b':(s-r)//2})
need([x['p'] for x in rows if x['b']==0]==O,'Ogg finite list')

digits={d:[] for d in '012'}
for n in range(1001):
    s=ternary(1<<n)
    for d in digits:
        if d not in s:digits[d].append(n)
need(digits=={'0':[0,1,2,3,4,15],'1':[1,3,9],'2':[0,2,8]},'digit test')
need(ternary(2**15)=='1122221122' and ternary(2**6)=='2101','persistence')

lifts=[]
for k in range(1,9):
    mod=3**k;u=2*3**(k-1);seq=[pow(2,n,mod) for n in range(u)]
    need(len(set(seq))==u,'primitive root')
    need(pow(2,u,3*mod)==1+mod,'principal unit')
    tails={d:[] for d in [0,2]}
    for n,x in enumerate(seq):
        ds=ternary(x).zfill(k)
        for d in tails:
            if str(d) not in ds:tails[d].append(n)
        for a in range(3):
            target=pow(2,n+a*u,mod*3)//mod
            expected=(pow(2,n,mod*3)//mod+a*(x%3))%3
            need(target==expected,'lift digit')
    need(len(tails[0])==2**k and len(tails[2])==2**(k-1),'tail census')
    rep=(mod-1)//2;exp=seq.index(rep)
    lifts.append({'k':k,'period':u,'zero_free_tails':len(tails[0]),'two_free_tails':len(tails[2]),'all_ones_tail_exponent':exp,'arbitrarily_large_exponents':'n = exponent + t*period'})

# Recovered 3-adic tournament: low nonzero ternary digit controls orientation.
def arc(x,y,mod):
    d=(y-x)%mod
    need(d!=0,'diagonal')
    while d%3==0:d//=3
    return d%3==1
signed=[]
for k in range(1,6):
    q=3**k;count=0
    for x in range(q):
        for y in range(x+1,q):
            a=arc(x,y,q)
            need(arc(2*x%q,2*y%q,q)!=a,'times2 anti')
            need(arc(4*x%q,4*y%q,q)==a,'times4 auto')
            count+=1
    signed.append({'k':k,'q':q,'unordered_pairs_checked':count})

# Finite no-zero tails are not invariant under the C3 kernel; k=1->2.
kernel_cycle=[1,4,7]
need([ternary(x).zfill(2) for x in kernel_cycle]==['01','11','21'],'kernel hostile')
need([0 not in [int(c) for c in ternary(x).zfill(2)] for x in kernel_cycle]==[False,True,True],'sieve not invariant')

machine='2RB2RB1RA3RA_0LC1RB3RB2RA_0LE1LD2LC3LC_0RB1LD2LD3LD_1RD2RD3RD0LF_1RE2RZ------'
parsed=[[r[i:i+3] for i in range(0,len(r),3)] for r in machine.split('_')]
need(len(parsed)==6 and all(len(r)==4 for r in parsed),'machine shape')
flat=sum(parsed,[])
need(sum(t=='---' for t in flat)==2 and sum(t.endswith('Z') for t in flat)==1,'instruction shape')
result={'status':'FINITE-EXACT; no conjecture or machine-correctness proof','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'ss_class_number_formula_range':'all primes <=211; classical formulas imported, not proved by this census','ss_rows':rows,'defect_distributions':[defect_distribution(n,b) for n in range(7) for b in range(n//2+1)],'p37':field37(),'digit_avoidance_exponents_n_le_1000':digits,'sloane_persistence_chain_decimal':[3**15-1,2**15,64,0],'sloane_persistence_chain_ternary':['2'*15,'1122221122','2101','0'],'modular_lifts':lifts,'signed_tournaments':signed,'sieve_kernel_hostile':{'cycle':kernel_cycle,'words':['01','11','21']},'sloane_machine_syntax_only':{'machine':machine,'states':6,'symbols':4,'nonhalting_instructions':21,'halting_instructions':1,'undefined_entries':2}}
print(json.dumps(result,sort_keys=True,indent=2))
