"""Direct odd-symbol census and explicit locally character-covariant construction."""
from itertools import product
from math import gcd
import json

def require(ok,why):
    if not ok:raise RuntimeError(why)
def units(n):return [u for u in range(1,n) if gcd(u,n)==1]
def divisors(n):return [d for d in range(2,n+1) if n%d==0]
def factors(n):
    out=[];p=2
    while p*p<=n:
        if n%p==0:
            a=0
            while n%p==0:n//=p;a+=1
            out.append((p,a))
        p+=1
    if n>1:out.append((n,1))
    return out
def chi(u,p):return 1 if pow(u%p,(p-1)//2,p)==1 else -1

def predicted_count(n):
    fac=factors(n)
    if any(p%4!=3 for p,a in fac):return 0
    exponent=0
    for p,a in fac:
        z=a
        for q,b in fac:
            if p!=q:z*=b+1
        exponent+=z
    return 2**exponent

def locally_covariant(sig):
    n=len(sig)
    for m in divisors(n):
        U=units(m);rep=n//m;eps=sig[rep]
        # Recover the only possible local character directly from the symbol.
        values={u:sig[rep*u]*eps for u in U}
        if values[1]!=1 or values[m-1]!=-1:return False
        if not all(values[u*v%m]==values[u]*values[v] for u in U for v in U):return False
    return True

def exhaustive_odd(n):
    ans=[]
    for frame in product((-1,1),repeat=(n-1)//2):
        sig=[0]*n
        for d,s in enumerate(frame,1):sig[d]=s;sig[n-d]=-s
        if locally_covariant(sig):ans.append(tuple(sig))
    return ans

def constructed(n):
    fac=factors(n)
    if any(p%4!=3 for p,a in fac):return []
    layers=[]
    for m in divisors(n):
        primes=[p for p,a in factors(m)];U=units(m);choices=[]
        for on in product((False,True),repeat=len(primes)):
            if sum(on)%2!=1:continue
            values=[]
            for u in U:
                v=1
                for p,b in zip(primes,on):
                    if b:v*=chi(u,p)
                values.append(v)
            for eps in (-1,1):choices.append([(n//m*u,eps*v) for u,v in zip(U,values)])
        layers.append(choices)
    ans=[]
    for choices in product(*layers):
        sig=[0]*n
        for choice in choices:
            for d,v in choice:sig[d]=v
        require(all(sig[-d%n]==-sig[d] for d in range(n)),'constructed oddness')
        ans.append(tuple(sig))
    require(len(ans)==len(set(ans)),'unique parametrization')
    return ans

def invariants(sig):
    n=len(sig)
    conv=[sum(sig[t]*sig[(d-t)%n] for t in range(n)) for d in range(n)]
    trace4=n*sum(z*z for z in conv)
    # Literally enumerate simple directed five-cycles containing0 with0 first.
    out=[sum(1<<y for y in range(n) if sig[(y-x)%n]==1) for x in range(n)]
    into0=sum(1<<x for x in range(n) if sig[-x%n]==1)
    rooted=0
    for a in range(1,n):
        if not(out[0]>>a&1):continue
        for b in range(1,n):
            if b==a or not(out[a]>>b&1):continue
            for c in range(1,n):
                if c in (a,b) or not(out[b]>>c&1):continue
                rooted+=(out[c]&into0&~((1<<0)|(1<<a)|(1<<b)|(1<<c))).bit_count()
    require(n*rooted%5==0,'cycle orbit divisibility')
    return trace4,n*rooted//5

census=[];N21=[]
for n in range(3,22,2):
    direct=exhaustive_odd(n);made=constructed(n)
    require(set(direct)==set(made),('literal/construction agreement',n))
    require(len(direct)==predicted_count(n),('count formula',n))
    census.append({'N':n,'all_odd_symbols_examined':2**((n-1)//2),'local_character_symbols':len(direct)})
    if n==21:N21=direct
large=[]
for n in [27,49,63,77,99,147]:
    made=constructed(n)
    require(len(made)==predicted_count(n),('larger formula',n))
    require(all(locally_covariant(sig) for sig in made),('larger direct equations',n))
    large.append({'N':n,'constructed_and_checked':len(made)})
profiles={}
for sig in N21:
    inv=invariants(sig);key=str(inv)
    if key not in profiles:profiles[key]={'trace_S4':inv[0],'directed_5cycles':inv[1],'multiplicity':0,'positive_differences':[d for d in range(1,21) if sig[d]==1]}
    profiles[key]['multiplicity']+=1
require(len(profiles)>=2,'distinct invariant profiles at21')
remaining=set(N21);unit_orbits=[]
while remaining:
    sig=min(remaining)
    orbit={tuple(sig[u*d%21] for d in range(21)) for u in units(21)}
    require(orbit<=set(N21),'unit orbit closure')
    invs={invariants(s) for s in orbit}
    require(len(invs)==1,'unit orbit invariants')
    unit_orbits.append({'size':len(orbit),'trace_S4':next(iter(invs))[0],'directed_5cycles':next(iter(invs))[1]})
    remaining-=orbit
require(len(unit_orbits)==4 and all(o['size']==4 for o in unit_orbits),'four unit orbits')
require(len({o['trace_S4'] for o in unit_orbits})==4,'no isomorphism merging')

def partial_covariant(sig):
    n=len(sig)
    for m in divisors(n):
        U=units(m);rep=n//m;data=[sig[rep*u] for u in U]
        if not any(data):continue
        if 0 in data:return False
        eps=sig[rep];values={u:sig[rep*u]*eps for u in U}
        if values[1]!=1 or values[m-1]!=-1:return False
        if not all(values[u*v%m]==values[u]*values[v] for u in U for v in U):return False
    return True
partial_rows=[]
for n in range(3,16,2):
    expected={0:1}
    for m in divisors(n):
        ps=[p for p,a in factors(m)]
        if not any(p%4==3 for p in ps):continue
        coefficient=2**len(ps);weight=len(units(m));new=expected.copy()
        for degree,count in expected.items():new[degree+weight]=new.get(degree+weight,0)+coefficient*count
        expected=new
    actual={}
    for frame in product((-1,0,1),repeat=(n-1)//2):
        sig=[0]*n
        for d,s in enumerate(frame,1):sig[d]=s;sig[n-d]=-s
        if partial_covariant(sig):
            weight=sum(v!=0 for v in sig);actual[weight]=actual.get(weight,0)+1
    require(actual==expected,('partial support polynomial',n))
    partial_rows.append({'N':n,'partial_symbols':sum(actual.values()),'support_polynomial':actual})
print(json.dumps({'status':'PASS','literal_census':census,'partial_literal_census':partial_rows,'larger_constructed_controls':large,'N21_profiles':list(profiles.values()),'N21_unit_orbits':unit_orbits,'N21_exact_isomorphism_classes_within_family':4,'scope':'specified stratumwise-character family; four literal unit orbits have different profiles, giving exact within-family isomorphism classes'},indent=2,sort_keys=True))
