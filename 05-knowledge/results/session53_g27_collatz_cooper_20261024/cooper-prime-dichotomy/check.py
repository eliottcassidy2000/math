import json,sys
from math import comb,isqrt

def require(x,msg):
    if not x: raise RuntimeError(msg)
def coefficients(N):
    c=[1,2]
    for n in range(1,N):
        a=(26*n*n+13*n+2)*c[n]+3*(3*n-1)*(3*n-2)*c[n-1]
        q,r=divmod(a,(n+1)**2)
        require(r==0,'nonintegral recurrence')
        c.append(q)
    return c
c=coefficients(8000)
require(c[:5]==[1,2,22,336,6006],'normalization')
for n in range(101):
    s=sum(comb(n,k)**2*comb(n+k,k)*comb(2*k,n) for k in range((n+1)//2,n+1))
    require(s==sum(c[k]*c[n-k] for k in range(n+1)),'independent square')
primes=[p for p in range(2,200) if all(p%d for d in range(2,isqrt(p)+1))]
rows=[]
for p in primes:
    good=p==7 or p%7 in (1,2,4)
    fails=[n for n in range(41) if (c[p*n]-c[n])%(p*p)]
    modp=[n for n in range(41) if (c[p*n]-c[n])%p]
    require((not fails)==good,'finite dichotomy')
    if not good: require(modp and modp[0]==1,'inert n=1 hostile')
    rows.append(dict(p=p,good=good,n_max=40,first_mod_p_failure=modp[0] if modp else None,first_mod_p2_failure=fails[0] if fails else None,delta_at_1_mod_p=(c[p]-2)%p))
# Independent base9 digit prediction and all finite power9 returns.
def digit_value(n):
    sign=1
    while n:
        n,d=divmod(n,9)
        if d>2:return 0
        if d==1:sign=-sign
    return sign%3
require(all(v%3==digit_value(n) for n,v in enumerate(c)), 'base9 digit law')
require(all((c[9*n]-c[n])%3==0 for n in range(8000//9+1)), 'Cartier period2')
require((c[3]+c[1])%3!=0 and (c[27]+c[1])%3!=0,'signed law hostility')
# Hostile guard: coefficientwise square root is not square root of the series.
require(c[2]**2 != 48,'coefficientwise-root confusion')
# Inert claim means failure of universal n, not every n (n=0 always succeeds).
require(all((c[0]-c[0])%p==0 for p in primes),'n=0')
require((c[3]-c[1])%3==1,'p3 explicit witness')
# Explicit failure of blind exponent strengthening.
strong=[]
for p in primes:
    if p==7 or p%7 in (1,2,4):
        fail=next((n for n in range(1,41) if (c[p*n]-c[n])%(p**3)),None)
        strong.append([p,fail])
# Corrupt one coefficient; checker must reject in both ordinary and -O modes.
corrupt=c[:];corrupt[22]+=1
caught=False
try: require((corrupt[22]-corrupt[2])%121==0,'intentional corruption')
except RuntimeError: caught=True
require(caught,'hostile integrity test disabled')
print(json.dumps({'scope':'FINITE-EXACT; does not prove all n or all primes','base9_digit_law_through':8000,'coefficient_max':8000,'independent_binomial_square_through':100,'prime_count':len(primes),'pair_count':len(primes)*41,'rows':rows,'p3_strengthening_first_failures':strong,'hostile_corruption_rejected':caught},sort_keys=True,indent=2))
