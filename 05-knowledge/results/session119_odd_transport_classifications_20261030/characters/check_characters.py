"""Independent finite orbit-constraint and sign-character census. Standard library."""
from itertools import product
from math import gcd
import json

def require(ok, why):
    if not ok: raise RuntimeError(why)

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

def units(n):return [u for u in range(n) if gcd(u,n)==1]
def divisors(n):return [m for m in range(1,n+1) if n%m==0]
def legendre(u,p):return 1 if pow(u%p,(p-1)//2,p)==1 else -1

def all_characters_bruteforce(n,U):
    # Every sign character kills squares. Enumerate ALL assignments on the
    # literal square cosets, then test the quotient multiplication table.
    squares={u*u%n for u in U};remaining=set(U);cosets=[];index={}
    while remaining:
        r=min(remaining); C={r*s%n for s in squares}; j=len(cosets)
        for u in C:index[u]=j
        cosets.append((r,C));remaining-=C
    require(index[1]==0,'identity coset')
    chars=[]
    for signs in product((-1,1),repeat=len(cosets)-1):
        values=(1,)+signs
        if all(values[index[a*b%n]]==values[i]*values[j]
               for i,(a,_) in enumerate(cosets) for j,(b,_) in enumerate(cosets)):
            chars.append({u:values[index[u]] for u in U})
    return chars,len(squares),len(cosets)

def signed_components(n,U,chi):
    graph={d:[] for d in range(1,n)}
    for d in graph:
        for u in U:graph[d].append((u*d%n,chi[u]))
        graph[d].append(((-d)%n,-1))
    seen={}; free=[];forced=[]
    for start in graph:
        if start in seen:continue
        seen[start]=1;stack=[start];vertices=[];bad=False
        while stack:
            d=stack.pop();vertices.append(d)
            for e,s in graph[d]:
                value=seen[d]*s
                if e not in seen:seen[e]=value;stack.append(e)
                elif seen[e]!=value:bad=True
        (forced if bad else free).append(set(vertices))
    return free,forced

def polynomial(sizes):
    coeff={0:1}
    for s in sizes:
        nxt=coeff.copy()
        for k,v in coeff.items():nxt[k+s]=nxt.get(k+s,0)+2*v
        coeff=nxt
    return coeff

def literal_symbols(n,U,chi):
    # Oddness imposed directly by the actual pairs, not the orbit theorem.
    count=full=0;weights={}
    for vals in product((-1,0,1),repeat=(n-1)//2):
        sig=[0]*n
        for d,s in enumerate(vals,1):sig[d]=s;sig[n-d]=-s
        if all(sig[u*d%n]==chi[u]*sig[d] for u in U for d in range(1,n)):
            count+=1;w=sum(s!=0 for s in sig);weights[w]=weights.get(w,0)+1
            full+= w==n-1
    return count,full,weights

rows=[];characters=literal_count=stabilizer_checks=0
for n in range(3,152,2):
    U=units(n);fac=factors(n); chars,square_size,quotient_size=all_characters_bruteforce(n,U)
    expected=[]
    for selected in product((False,True),repeat=len(fac)):
        S=[p for (p,_),on in zip(fac,selected) if on]
        c={u:1 for u in U}
        for u in U:
            for p in S:c[u]*=legendre(u,p)
        f=1
        for p in S:f*=p
        expected.append((c,f,S))
    require({tuple(c[u] for u in U) for c in chars}=={tuple(c[u] for u in U) for c,_,_ in expected},('complete characters',n))
    require(quotient_size==2**len(fac),('square quotient',n))
    full_total=0
    for chi,f,S in expected:
        characters+=1;odd=chi[n-1]==-1
        free,forced=signed_components(n,U,chi)
        eligible=[m for m in divisors(n) if m>1 and m%f==0] if odd else []
        free_orders=sorted({n//gcd(n,d) for d in C}.pop() for C in free)
        require(free_orders==eligible,('strata',n,f))
        D=1
        for p,a in fac:D*=a if p in S else a+1
        if odd:require(len(free)==D,('divisor product',n,f))
        sizes=[len(units(m)) for m in eligible]
        require(polynomial([len(C) for C in free])==polynomial(sizes),('support polynomial',n,f))
        full=0 if forced else 2**len(free)
        full_total+=full
        if n<=15:
            c,actual_full,weights=literal_symbols(n,U,chi)
            require(c==3**len(free) and actual_full==full and weights==polynomial(sizes),('literal census',n,f))
            literal_count+=c
        # Build one nonzero symbol from each permitted stratum and verify all
        # covariance, oddness and the EXACT multiplier stabilizer directly.
        for m in eligible:
            sig=[0]*n;rep=n//m
            for u in U:
                d=u*rep%n
                require(sig[d] in (0,chi[u]),('well-defined',n,f,m))
                sig[d]=chi[u]
            require(all(sig[-d%n]==-sig[d] for d in range(n)),('odd symbol',n,f,m))
            stab={u for u in U if all(sig[u*d%n]==sig[d] for d in range(n))}
            require(stab=={u for u in U if chi[u]==1},('stabilizer',n,f,m))
            stabilizer_checks+=1
    expected_full=2**fac[0][1] if len(fac)==1 and fac[0][0]%4==3 else 0
    require(full_total==expected_full,('complete tournament family',n))
    rows.append({'N':n,'characters':len(chars),'odd_characters':sum(c[n-1]==-1 for c in chars),'full_tournament_symbols':full_total})
# Smallest useful CRT hostile to naive all-stratum extension.
require(11*3%15==3 and legendre(11,3)==-1,'N15 hostile')
print(json.dumps({'status':'PASS','odd_moduli':[3,151],'moduli_count':len(rows),'character_cases':characters,'literal_valid_symbols_counted_N_le_15':literal_count,'exact_unit_stabilizer_cases':stabilizer_checks,'hostile':{'N':15,'conductor':3,'fixed_residue':3,'multiplier':11,'character_value':-1},'rows':rows},indent=2,sort_keys=True))
