#!/usr/bin/env python3
"""Exact verifier; Python standard library only. Optional --generate needs SymPy.
All primality assertions use recursive Lucas certificates, never probable-prime output.
"""
import argparse, functools, hashlib, json, math
from pathlib import Path
ROOT = Path(__file__).resolve().parent
OGG = {2,3,5,7,11,13,17,19,23,29,31,41,47,59,71}
EXPECTED_F = [3,5,7,13,23,17,19,23,37,61,67,61,71,47,107,59,61,109,89]
EXPECTED_L = [None,3,7,11,13,17,29,23,43,41,73,59,47,89,67,73,107,89,101]

def require(ok, msg):
    if not ok: raise ValueError(msg)

def small_prime(n):
    return n >= 2 and all(n % d for d in range(2, math.isqrt(n)+1))

def mr_failure(n, a):
    """True is an exact compositeness witness, not a probable-primality claim."""
    if n < 3 or n % 2 == 0: return False
    d, s = n-1, 0
    while d % 2 == 0: d //= 2; s += 1
    x = pow(a, d, n)
    if x in (1,n-1): return False
    for _ in range(s-1):
        x = x*x % n
        if x == n-1: return False
    return True

def generate():
    from sympy import factorint, isprime
    certs = {}
    def cert(n):
        n = int(n)
        if str(n) in certs: return
        if n == 2: certs['2'] = {'factors': [], 'base': 1}; return
        print('certifying',n,flush=True)
        factors = sorted((int(q),int(e)) for q,e in factorint(n-1).items())
        for q,e in factors: cert(q)
        a = 2
        while not (pow(a,n-1,n)==1 and all(math.gcd(pow(a,(n-1)//q,n)-1,n)==1 for q,e in factors)):
            a += 1
        certs[str(n)] = {'factors': factors, 'base': a}
    def classify(n):
        if isprime(n): cert(n); return {'prime': True}
        if n < 2: return {'nonpositive_or_one': True}
        for d in range(2,100):
            if n%d == 0 and d<n: return {'divisor': d}
        for a in range(2,100):
            if mr_failure(n,a): return {'mr_base': a}
        factors = factorint(n)
        d = int(min(factors))
        require(1<d<n and n%d==0,'missing composite factor')
        return {'divisor':d}
    rows=[]; P=1
    for i,p in enumerate(q for q in range(2,71) if small_prime(q)):
        P*=p; row={'n':i+1,'p':p,'P':P}
        for label,sign in [('plus',1),('minus',-1)]:
            if i==0 and sign==-1: row[label]=None; continue
            witnesses=[]
            for d in range(2,P+1000):
                w=classify(P+sign*d)
                witnesses.append(w)
                if w.get('prime'): break
            row[label]={'offset':d,'witnesses':witnesses}
        rows.append(row)
        print('ROW',row['n'],p,row['plus']['offset'],row['minus'] and row['minus']['offset'],flush=True)
    special={str(n):classify(n) for n in [30013,30047,30029,30031,29741,30319,30049,30011]}
    for p in range(2,72):
        if small_prime(p): cert(p)
    obj={'schema':1,'cutoff':71,'rows':rows,'special':special,'primality_certificates':certs}
    (ROOT/'certificates.json').write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')

def reduced_forms(D):
    """Primitive positive definite reduced forms; boundary uses b>=0."""
    require(D<0 and D%4 in (0,1),'negative quadratic discriminant')
    forms=[]
    for a in range(1,math.isqrt((-D)//3)+1):
        for b in range(-a,a+1):
            if (b*b-D)%(4*a): continue
            c=(b*b-D)//(4*a)
            if a>c or math.gcd(a,math.gcd(b,c))!=1: continue
            if (abs(b)==a or a==c) and b<0: continue
            forms.append((a,b,c))
    return forms

def supersingular_j_by_count(p):
    out=[]
    for j in range(p):
        if j==0: a,b=0,1
        elif j==1728%p: a,b=1,0
        else:
            k=j*pow((1728-j)%p,-1,p)%p; a,b=3*k%p,2*k%p
        require((4*a*a*a+27*b*b)%p!=0,'nonsingular point-count model')
        trace=0
        for x in range(p):
            v=(x*x*x+a*x+b)%p
            if v:
                z=pow(v,(p-1)//2,p); trace-= 1 if z==1 else -1
        if trace==0: out.append(j)
    return out

def ogg_row(p):
    if p in (2,3): return {'p':p,'ogg':True,'g_plus':0,'low_characteristic':True}
    delta={1:0,5:1,7:1,11:2}[p%12]
    degree=p//12+delta
    h4=len(reduced_forms(-4*p)); h1=len(reduced_forms(-p)) if p%4==3 else 0
    js=supersingular_j_by_count(p)
    require(2*len(js)==h4+h1,'class-number and point-count agreement')
    numerator=2*degree-h4-h1
    require(numerator>=0 and numerator%4==0,'invalid class number genus')
    return {'p':p,'ogg':numerator==0,'g_plus':numerator//4,'ss_degree':degree,'rational_supersingular_j':js,'h_minus_4p':h4,'h_minus_p_if_applicable':h1}

def verify():
    obj=json.loads((ROOT/'certificates.json').read_text()); certs=obj['primality_certificates']
    @functools.lru_cache(None)
    def prime(n):
        require(n>=2 and str(n) in certs,'prime certificate missing '+str(n))
        c=certs[str(n)]
        if n==2: return True
        fs=[(int(q),int(e)) for q,e in c['factors']]
        require(len({q for q,e in fs})==len(fs),'duplicate factor')
        require(all(2<=q<n and e>=1 for q,e in fs),'invalid factor range')
        require(math.prod(q**e for q,e in fs)==n-1,'incomplete n-1 factorization')
        for q,e in fs: prime(q)
        a=c['base']; require(1<a<n,'Lucas base range')
        require(pow(a,n-1,n)==1,'Lucas Fermat equality')
        require(all(math.gcd(pow(a,(n-1)//q,n)-1,n)==1 for q,e in fs),'Lucas order test')
        return True
    def check(n,w):
        if w.get('prime'): return prime(n)
        if w.get('nonpositive_or_one'): require(n<2,'small composite witness'); return False
        if 'divisor' in w:
            d=w['divisor']; require(1<d<n and n%d==0,'proper divisor witness'); return False
        require('mr_base' in w and 1<w['mr_base']<n and mr_failure(n,w['mr_base']),'MR compositeness witness')
        return False
    ps=[p for p in range(2,71) if small_prime(p)]
    rows=obj['rows']; require(len(rows)==len(ps)==19,'complete finite cutoff rows')
    P=1; F=[]; L=[]
    for n,(p,row) in enumerate(zip(ps,rows),1):
        P*=p; require((row['n'],row['p'],row['P'])==(n,p,P),'primorial indexing')
        for label,sign,output in [('plus',1,F),('minus',-1,L)]:
            data=row[label]
            if n==1 and sign==-1: require(data is None,'undefined lower offset at P1=2'); output.append(None); continue
            d=data['offset']; ws=data['witnesses']
            require(len(ws)==d-1,'every tested offset present')
            for k,w in enumerate(ws,2): require(check(P+sign*k,w)==(k==d),'minimality or primality failure')
            require(d>p or (n==2 and sign==-1 and d==3),'offset lower bound')
            require(d<next(q for q in range(p+1,p*2+2) if small_prime(q))**2,'finite offset square-bound control')
            output.append(d)
    require(F==EXPECTED_F and L==EXPECTED_L,'sequence mismatch')
    for key,w in obj['special'].items(): check(int(key),w)
    for key in certs: prime(int(key))
    require(all(math.gcd(d,30030)>1 for d in range(2,17)),'small-offset wheel')
    require(all(not small_prime(30030+sign*d) for d in range(2,17) for sign in [-1,1]),'minimal common offset')
    require(small_prime(30013) and small_prime(30047),'symmetric endpoints independent trial division')
    require(small_prime(29741) and small_prime(30319) and 17**2==289,'sharp square threshold')
    require(not small_prime(30049) and 30049==151*199,'wheel-survival hostile')
    require(small_prime(23) and math.gcd(23,30030)==1 and 30030-23==37*811 and 30030+23==41*733,'two-sided false positive')
    arithmetic=[ogg_row(p) for p in range(2,72) if small_prime(p)]
    require({x['p'] for x in arithmetic if x['ogg']}==OGG,'Ogg class number predicate')
    fortunate={x for x in F if x<=71}; lesser={x for x in L if x and x<=71}
    for row in arithmetic:
        p=row['p']; row.update({'fortunate_indices':[i+1 for i,x in enumerate(F) if x==p], 'lesser_indices':[i+1 for i,x in enumerate(L) if x==p]})
    controls=[]
    for d in range(2,290):
        for sign in [-1,1]:
            endpoint=30030+sign*d
            if small_prime(endpoint):
                require(math.gcd(d,30030)==1,'endpoint forces unit offset')
                if d<289: require(small_prime(d),'square-threshold primality')
                controls.append({'offset':d,'sign':sign,'endpoint':endpoint,'offset_prime':small_prime(d)})
    # Logical quantifiers: local factor tests intersect; alternative witnesses unite.
    I=set(range(2,290)); Aplus={d for d in I if small_prime(30030+d)}; Aminus={d for d in I if small_prime(30030-d)}
    Hplus=I-Aplus; Hminus=I-Aminus
    require(I-(Aplus|Aminus)==Hplus&Hminus,'exact existential hole identity')
    require(I-(Aplus&Aminus)==Hplus|Hminus,'exact conjunction hole identity')
    wheel={d for d in I if math.gcd(d,30030)==1}
    require(Aplus|Aminus<=wheel and 19 in wheel and 19 not in Aplus,'wheel is necessary only')
    require(sum(math.gcd(d,30030)==1 for d in range(30030))==5760,'exact wheel cardinality')
    # Euclid--Mullin is a least-factor recurrence, not an offset sequence.
    em=[]; product=1
    for i in range(7):
        number=product+1
        least=next(d for d in range(2,number+1) if number%d==0)
        require(small_prime(least),'least factor primality')
        em.append(least); product*=least
    require(em==[2,3,7,43,13,53,5],'Euclid-Mullin prefix')
    summary={
        'status':'PASS','euclid_mullin_first_7':em,'wheel_units_mod_30030':5760,'cutoff':71,'complete_for_all_indices_below_cutoff':True,
        'prime_certificates':len(certs),'max_certified_prime':str(max(map(int,certs))),
        'F_1_through_19':F,'L_1_through_19':L,
        'fortunate_values_le_71':sorted(fortunate), 'lesser_values_le_71':sorted(lesser),
        'unfortunate_primes_le_71':sorted(set(p for p in range(2,72) if small_prime(p))-fortunate),
        'ogg_intersection_fortunate':sorted(OGG&fortunate),'ogg_minus_fortunate':sorted(OGG-fortunate),
        'fortunate_le_71_minus_ogg':sorted(fortunate-OGG),
        'ogg_intersection_lesser':sorted(OGG&lesser),'ogg_minus_lesser':sorted(OGG-lesser),
        'lesser_le_71_minus_ogg':sorted(lesser-OGG),
        'ogg_minus_union_of_offsets':sorted(OGG-(fortunate|lesser)),
        'arithmetic_predicates':arithmetic,'prime_endpoints_offset_2_to_289':controls,
        'certificate_sha256':hashlib.sha256((ROOT/'certificates.json').read_bytes()).hexdigest()
    }
    (ROOT/'results.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['arithmetic_predicates','prime_endpoints_offset_2_to_289']},indent=2,sort_keys=True))

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--generate',action='store_true'); args=ap.parse_args()
    if args.generate: generate()
    verify()
