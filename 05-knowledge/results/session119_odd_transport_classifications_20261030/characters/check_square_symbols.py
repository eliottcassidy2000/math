"""Bounded quotient-Boolean classification and one literal modulus231 symbol."""
from itertools import product
from math import gcd
import json

def require(ok,why):
    if not ok:raise RuntimeError(why)
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
def units(n):return [u for u in range(1,n) if gcd(u,n)==1]
def legendre(x,p):return 1 if pow(x%p,(p-1)//2,p)==1 else -1
def char(x,S):return -1 if (x&S).bit_count()%2 else 1

def odd_tables(r):
    size=1<<r;mask=size-1;representatives=[x for x in range(size) if x<(x^mask)]
    ans=[]
    for signs in product((-1,1),repeat=len(representatives)):
        F=[0]*size
        for x,s in zip(representatives,signs):F[x]=s;F[x^mask]=-s
        ans.append(tuple(F))
    return ans
rows=[];r3=None
for r in range(1,5):
    tables=set(odd_tables(r));size=1<<r
    characters={tuple(e*char(x,S) for x in range(size)) for S in range(size) if S.bit_count()%2 for e in(-1,1)}
    require(len(tables)==2**(2**(r-1)),'all odd quotient tables')
    require(len(characters)==2**r and characters<=tables,'character subfamily')
    require((tables==characters)==(r<=2),'dimension threshold')
    row={'r':r,'odd_tables':len(tables),'signed_odd_characters':len(characters)}
    if r==3:
        M=tuple(1 if sum(char(x,1<<i) for i in range(3))>0 else -1 for x in range(8))
        family={tuple(e*char(x,S)*M[x] for x in range(8)) for S in range(8) if S.bit_count()%2==0 for e in(-1,1)}
        require(len(family)==8 and family==tables-characters,'complete majority family')
        for F in family:
            spectrum=[sum(F[x]*char(x,S) for x in range(8)) for S in range(8)]
            require(all(spectrum[S]==0 if S.bit_count()%2==0 else abs(spectrum[S])==4 for S in range(8)),'majority spectrum')
        r3={'character_type':8,'majority_type':8,'nonzero_unnormalized_Walsh_magnitude':4}
    rows.append(row)
# Exact count identities over a bounded arithmetic universe.
moduli=[]
for n in range(3,501,2):
    fac=factors(n)
    if any(p%4!=3 for p,a in fac):continue
    ds=[m for m in range(2,n+1) if n%m==0]
    C=sum(2**(len(factors(m))-1) for m in ds)
    L=sum(len(factors(m)) for m in ds)
    prod=1
    for p,a in fac:prod*=1+2*a
    require(2*C==prod-1,'square family exponent')
    closedL=0
    for p,a in fac:
        value=a
        for q,b in fac:
            if p!=q:value*=b+1
        closedL+=value
    require(L==closedL and ((C==L)==(len(fac)<=2)),'family dimension threshold')
    moduli.append({'N':n,'square_exponent':C,'local_character_exponent':L})
# Literal231 construction: majority on units, first-prime character on lower strata.
n=231;sig=[0]*n
for d in range(1,n):
    m=n//gcd(n,d);t=d//(n//m);ps=[p for p,a in factors(m)]
    sig[d]=(1 if sum(legendre(t,p) for p in ps)>0 else -1) if m==n else legendre(t,ps[0])
require(all(sig[-d%n]==-sig[d] for d in range(n)),'literal oddness')
U=units(n);squares={u*u%n for u in U}
require(len(squares)==15,'square subgroup order')
checks=0
for q in squares:
    for d in range(n):
        require(sig[q*d%n]==sig[d],'square invariance');checks+=1
require(sig[1]==1 and sig[2]==sig[5]==sig[10]==-1,'explicit noncharacter hostile')
require(sig[2*5%n]!=sig[2]*sig[5],'not a unit character')
require(sum(v==1 for v in sig)==115,'full tournament orientation')
print(json.dumps({'status':'PASS','quotient_rows':rows,'r3_classification':r3,'arithmetic_count_controls':moduli,'N231':{'square_invariant_symbols':8192,'local_character_symbols':4096,'square_invariance_equations':checks,'square_subgroup_order':len(squares),'positive_differences':[d for d in range(1,n) if sig[d]==1],'noncharacter_witness':{'u':2,'v':5,'uv':10,'sigma_u':sig[2],'sigma_v':sig[5],'sigma_uv':sig[10]}}},indent=2,sort_keys=True))
