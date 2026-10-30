"""Exact random-selector phase checks. No third-party dependencies."""
from fractions import Fraction
from itertools import product
from math import comb, factorial
import json
import random

CHECKS=0
def need(ok, why):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(why)

def stirling(n,k):
    a=[1]+[0]*k
    for _ in range(n):
        a=[0]+[a[j-1]+j*a[j] for j in range(1,k+1)]
    return a[k]

def cover(p,r):
    return Fraction(factorial(p)*stirling(r,p),p**r)

def base(p,k,r):
    if r<p: return Fraction(0)
    return cover(p,r)*(1-Fraction(p**(2*k+p-2),p**(k*r)))

def expected_extra(p,k):
    ans=Fraction(factorial(p)*p**(2*k+p-2),p**((k+1)*p))
    for j in range(1,p+1): ans /= 1-Fraction(j,p**(k+1))
    return ans

# Literal ring reflections through 0 and 1: x+y*d -> x-y*d + 2*p^(k-1)*h(x mod p)*d.
# Their fixed point sets are independently compared with bent-line membership.
p,k=5,2
q=p**k
points=list(product(range(q),repeat=2))
hs=[(0,0,*tail) for tail in product(range(p),repeat=p-2)]
masks=[]
for x,y in points:
    mask=0
    for i,h in enumerate(hs):
        reflected=(x,(-y+2*(q//p)*h[x%p])%q)
        on_line=(y-(q//p)*h[x%p])%q==0
        need((reflected==(x,y))==on_line,'literal reflection versus line')
        if on_line: mask |= 1<<i
    masks.append(mask)

# Exhaust all 625 third points and 625^2 ordered fourth-point choices.
strata={}
for r in (3,4):
    counts={}
    for extra in product(range(q*q),repeat=r-2):
        j=len({0,1,*[points[e][0]%p for e in extra]})
        m=(1<<len(hs))-1
        for e in extra: m &= masks[e]
        all_n,bad=counts.get(j,(0,0))
        counts[j]=(all_n+1,bad+bool(m))
        if m: need(m.bit_count()==p**(p-j),'all surviving reflection cosets')
    for j,(total,bad) in counts.items():
        pred=Fraction(p**(2*k+j-2),p**(k*r))
        need(Fraction(bad,total)==pred,'exact conditional phase probability')
    strata[r]=counts

# Minimal transversals: all low delta digits at x=0,...,4.
low_total=low_bad=0
for ys in product(range(p),repeat=p):
    low_total+=1
    if any(all((ys[x]-s*x-b)%p==0 for x in range(p)) for s,b in product(range(p),repeat=2)):
        low_bad+=1
need(Fraction(low_bad,low_total)==Fraction(1,p**((k-1)*(p-2))),'minimum transversal phase fraction')

# Independent generic clock enumeration versus the two-pin bent-line test.
def phase_test(p,k,pins):
    q,Q=p**k,p**(k-1)
    a=pins[0]
    b=next(x for x in pins if x[0]%p!=a[0]%p)
    slope=(b[1]-a[1])*pow((b[0]-a[0])%Q,-1,Q)%Q
    intercept=(a[1]-slope*a[0])%Q
    seen={}
    for x,y in pins:
        d=(y-slope*x-intercept)%q
        if d%Q: return False
        key=x%p
        if key in seen and seen[key]!=d//Q: return False
        seen[key]=d//Q
    return True

def reflection_count(p,k,D,pins):
    q,Q=p**k,p**(k-1)
    def mul(a,b): return ((a[0]*b[0]+D*a[1]*b[1])%q,(a[0]*b[1]+a[1]*b[0])%q)
    inv=pow(1-D,-1,q)
    g=((1+D)*inv%q,2*inv%q)
    root=pins[0]
    shifted=[((x-root[0])%q,(y-root[1])%q) for x,y in pins]
    total=0
    power=(1,0)
    for _ in range(q):
        forced={0:0}
        for x,y in shifted:
            u,v=mul(power,(x,-y))
            dx,dy=(x-u)%q,(y-v)%q
            key=x%p
            if dx or dy%Q: break
            value=dy//Q
            if key in forced and forced[key]!=value: break
            forced[key]=value
        else: total+=p**(p-len(forced))
        power=mul(power,g)
    need(total%p==0,'root-fixed parameter duplication')
    return total//p

rng=random.Random(112)
generic=0
for p,k,D in ((5,2,5),(5,2,10),(5,3,5),(7,2,7)):
    q=p**k
    for r in range(2,p+3):
        for _ in range(60):
            pins=[(rng.randrange(q),rng.randrange(q)) for __ in range(r)]
            j=len({x%p for x,y in pins})
            if j<2: continue
            n=reflection_count(p,k,D,pins)
            need(n==(p**(p-j) if phase_test(p,k,pins) else 0),'clock enumeration versus phase test')
            generic+=1

# Every generated bad transversal is repaired by changing one pin by delta.
repairs=0
for p,k in ((5,2),(5,3),(7,2)):
    q,Q=p**k,p**(k-1)
    for _ in range(100):
        slope,intercept=rng.randrange(Q),rng.randrange(Q)
        h=[rng.randrange(p) for __ in range(p)]
        pins=[]
        for residue in range(p):
            x=residue+p*rng.randrange(q//p)
            pins.append((x,(slope*x+intercept+Q*h[residue])%q))
        need(phase_test(p,k,pins),'generated bad transversal')
        for index in range(p):
            changed=list(pins)
            x,y=changed[index]
            changed[index]=(x,(y+1)%q)
            need(not phase_test(p,k,changed),'one delta move repairs phase')
            repairs+=1

# General exact rational identities: occupancy inclusion-exclusion and geometric-tail sum.
records=[]
for p in (5,7,11):
  for k in (2,3,4):
    for r in range(0,3*p+1):
        c=sum(Fraction((-1)**j*comb(p,j)*(p-j)**r,p**r) for j in range(p+1))
        need(c==cover(p,r),'coupon inclusion exclusion')
        need(0<=base(p,k,r)<=cover(p,r)<=1,'base probability bounds')
    e=sum(Fraction(1,j) for j in range(1,p+1))*p
    # Rational summation of extra survival tail using inclusion-exclusion.
    delta=p**(2*k+p-2)*sum(Fraction((-1)**j*comb(p,j),1)/
                           (1-Fraction(p-j,p**(k+1))) for j in range(p+1))
    need(delta==expected_extra(p,k),'expected extra tail formula')
    records.append({'p':p,'k':k,'directed_expected':str(e),'unsigned_extra':str(delta),
                    'unsigned_extra_decimal':float(delta),
                    'minimal_transversal_bad':str(Fraction(1,p**((k-1)*(p-2))))})

# k=1 hostile: independent reflections force >=2 distinct vertices in EACH fibre.
def shallow_formula(p,r):
    ans=Fraction(0)
    for a in range(p+1):
      for b in range(p-a+1):
        c=p-a-b
        coeff=factorial(p)//(factorial(a)*factorial(b)*factorial(c))*(-p)**b*(p-1)**c
        ans += Fraction(coeff*(a*p+b)**r,p**(2*r))
    return ans
# Smaller p=3 checks are a pure wreath-action identity, not a ramified-clock theorem.
for p in (2,3):
  for r in range(0,7):
    good=0
    for tup in product(range(p*p),repeat=r):
        occupied=[set() for _ in range(p)]
        for x in tup: occupied[x//p].add(x%p)
        good+=all(len(s)>=2 for s in occupied)
    need(shallow_formula(p,r)==Fraction(good,p**(2*r)),'wreath occupancy exact enumeration')
need(shallow_formula(5,5)==0,'k1 p pins never suffice')
need(base(5,2,5)>0,'k2 p pins can suffice')
print(json.dumps({'status':'PASS','checks':CHECKS,'generic_clock_cases':generic,'repairs_checked':repairs,'pair_anchored_exhaustive_strata':strata,
                  'transversal_low_patterns':{'total':low_total,'bad':low_bad},'families':records},indent=2))
