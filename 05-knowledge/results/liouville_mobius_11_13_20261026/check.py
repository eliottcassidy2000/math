"""Exact arithmetic and tournament controls; no inherited census filters."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import isqrt, prod
from pathlib import Path
import hashlib
import json
import sys

sys.stdout.reconfigure(encoding='utf-8',newline='\n')

def need(test,message):
    if not test: raise RuntimeError(message)

def factors(n):
    result={}
    p=2
    while p*p<=n:
        while n%p==0:
            result[p]=result.get(p,0)+1
            n//=p
        p+=1
    if n>1: result[n]=result.get(n,0)+1
    return result

def lam(n): return (-1)**sum(factors(n).values())
def mu(n):
    f=factors(n)
    return 0 if any(a>1 for a in f.values()) else (-1)**len(f)

def divisors(n): return [d for d in range(1,n+1) if n%d==0]

def matmul(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]

def rank(a):
    a=[[Fraction(x) for x in row] for row in a]; i=0
    for j in range(len(a[0])):
        k=next((k for k in range(i,len(a)) if a[k][j]),None)
        if k is None: continue
        a[i],a[k]=a[k],a[i]; z=a[i][j]; a[i]=[x/z for x in a[i]]
        for k in range(len(a)):
            if k!=i:
                z=a[k][j]; a[k]=[x-z*y for x,y in zip(a[k],a[i])]
        i+=1
        if i==len(a): break
    return i

def kernel_checks(r):
    m=1<<r
    a=[[0 if x&y else (-1)**(x.bit_count()+y.bit_count()) for y in range(m)] for x in range(m)]
    inv=[[(-1)**r if x|y==m-1 else 0 for y in range(m)] for x in range(m)]
    need(matmul(a,inv)==[[int(x==y) for y in range(m)] for x in range(m)],'tensor inverse')
    l=[[(-1)**(x.bit_count()+y.bit_count()) for y in range(m)] for x in range(m)]
    need(rank(l)==1,'Liouville rank one')
    return {'primes':r,'mobius_rank':m,'liouville_rank':1,'integer_inverse_verified':True}

def mul(a,b,N):
    out=[0]*(N+1)
    for i,x in enumerate(a):
        if x:
            for j in range(min(N-i,len(b)-1)+1): out[i+j]+=x*b[j]
    return out

def invseries(a,N):
    need(a[0]==1,'unit series')
    out=[1]+[0]*N
    for n in range(1,N+1): out[n]=-sum(a[k]*out[n-k] for k in range(1,n+1))
    return out

def product_coeff(fun,N):
    out=[1]+[0]*N
    for n in range(1,N+1):
        e=fun(n)
        if e==1:
            for k in range(N,n-1,-1): out[k]-=out[k-n]
        elif e==-1:
            for k in range(n,N+1): out[k]+=out[k-n]
    return out

def ghost_coeff(fun,N):
    ghost=[0]+[sum(d*fun(d) for d in divisors(n)) for n in range(1,N+1)]
    out=[1]+[0]*N
    for n in range(1,N+1):
        num=-sum(ghost[k]*out[n-k] for k in range(1,n+1))
        need(num%n==0,'integral ghost recurrence'); out[n]=num//n
    return ghost,out

def converse(rows):
    return tuple(((1<<len(rows))-1)^row^(1<<i) for i,row in enumerate(rows))

def frommask(n,value):
    out=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        a,b=(i,j) if value>>k&1 else (j,i); out[a]|=1<<b
    return tuple(out)

def relabel(rows,p):
    return tuple(sum(1<<p[j] for j in range(len(rows)) if rows[i]>>j&1)
                 for i in sorted(range(len(rows)),key=lambda i:p[i]))

def canonical(rows):
    return min(relabel(rows,p) for p in permutations(range(len(rows))))

def join(a,b):
    cross=((1<<len(b))-1)<<len(a)
    return tuple(x|cross for x in a)+tuple(x<<len(a) for x in b)

def c3sub(children):
    out=[]; sizes=list(map(len,children)); offsets=[0,sizes[0],sizes[0]+sizes[1]]
    for i,t in enumerate(children):
        j=(i+1)%3; cross=((1<<sizes[j])-1)<<offsets[j]
        out.extend((row<<offsets[i])|cross for row in t)
    return tuple(out)

def substitute(q,children):
    offsets=[]; size=0
    for child in children: offsets.append(size); size+=len(child)
    return tuple((row<<offsets[i]) | sum(((1<<len(children[j]))-1)<<offsets[j]
                 for j in range(len(q)) if q[i]>>j&1)
                 for i,child in enumerate(children) for row in child)

def cyclic(rows):
    return frozenset(t for t in combinations(range(len(rows)),3)
                     if all(((rows[t[i]]>>t[(i+1)%3])&1)==((rows[t[0]]>>t[1])&1) for i in range(3)))

def aut(rows):
    n=len(rows); deg=list(map(int.bit_count,rows)); mapping={}; total=0
    def walk(remaining,unused):
        nonlocal total
        if not remaining: total+=1; return
        options={i:[j for j in unused if deg[i]==deg[j] and all(
            ((rows[i]>>k)&1)==((rows[j]>>l)&1) for k,l in mapping.items())] for i in remaining}
        i=min(remaining,key=lambda i:len(options[i]))
        for j in options[i]:
            mapping[i]=j; walk(remaining-{i},unused-{j}); del mapping[i]
    walk(set(range(n)),set(range(n)))
    return total

def hp(rows):
    n=len(rows); dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i,count in enumerate(dp[mask]):
            choices=rows[i]&~mask
            if count:
                while choices:
                    bit=choices&-choices; choices-=bit
                    dp[mask|bit][bit.bit_length()-1]+=count
    return sum(dp[-1])

def cycles(rows,k):
    answer=[]
    for vertices in combinations(range(len(rows)),k):
        for tail in permutations(vertices[1:]):
            cycle=(vertices[0],)+tail
            if all(rows[cycle[i]]>>cycle[(i+1)%k]&1 for i in range(k)): answer.append(cycle)
    return answer

def presentation(n,bits):
    out=[0]*n; k=0
    for i,j in combinations(range(n),2):
        forward=j==i+1 or bool(bits>>k&1)
        if j>i+1: k+=1
        a,b=(i,j) if forward else (j,i); out[a]|=1<<b
    return tuple(out)

def tau(rows):
    out=list(rows)
    for i,j in combinations(range(len(rows)),2):
        if j>i+1: out[i]^=1<<j; out[j]^=1<<i
    return tuple(out)

def main():
    N=1000
    for n in range(1,N+1):
        ds=divisors(n); eps=int(n==1); square=int(isqrt(n)**2==n)
        need(sum(mu(d) for d in ds)==eps,'Mobius divisor inverse')
        need(sum(lam(d) for d in ds)==square,'Liouville square projector')
        need(sum(lam(d)*mu(n//d)**2 for d in ds)==eps,'Liouville inverse')
        need(lam(n)==sum(mu(n//(d*d)) for d in range(1,isqrt(n)+1) if n%(d*d)==0),'square lift')
        need(mu(n)==sum(mu(d)*lam(n//(d*d)) for d in range(1,isqrt(n)+1) if n%(d*d)==0),'square inversion')
        need(mu(n)==lam(n)*mu(n)**2,'pointwise mask')
    for a in range(1,101):
        for b in range(1,101):
            need(lam(a*b)==lam(a)*lam(b),'complete character')
            need(mu(a*b)==(mu(a)*mu(b) if not set(factors(a))&set(factors(b)) else 0),'overlap annihilation')
    M=[0]; L=[0]
    for n in range(1,N+1): M.append(M[-1]+mu(n)); L.append(L[-1]+lam(n))
    for n in range(1,N+1):
        need(L[n]==sum(M[n//(d*d)] for d in range(1,isqrt(n)+1)),'summatory lift')
        need(M[n]==sum(mu(d)*L[n//(d*d)] for d in range(1,isqrt(n)+1)),'summatory inverse')
    limit=80
    gm,pm=ghost_coeff(mu,limit); gl,pl=ghost_coeff(lam,limit)
    need(pm==product_coeff(mu,limit) and pl==product_coeff(lam,limit),'independent Euler coefficient paths')
    lifted=[1]+[0]*limit; lowered=[1]+[0]*limit
    for d in range(1,isqrt(limit)+1):
        a=[pm[i//(d*d)] if i%(d*d)==0 else 0 for i in range(limit+1)]
        lifted=mul(lifted,a,limit)
        b=[pl[i//(d*d)] if i%(d*d)==0 else 0 for i in range(limit+1)]
        if mu(d)==-1: b=invseries(b,limit)
        if mu(d): lowered=mul(lowered,b,limit)
    need(lifted==pl and lowered==pm,'ordinary product square lift')
    for n in range(1,limit+1):
        need(gm[n]==prod(1-p for p in factors(n)),'Mobius ghost')
        need(gl[n]==prod(sum((-p)**j for j in range(a+1)) for p,a in factors(n).items()),'Liouville ghost')
        need((1 if gl[n]>0 else -1)==lam(n),'ghost sign')

    old11=(24,1,3,6,14); old13=(18,8,3,5,14)
    pair=[tau(join((0,),presentation(4,bits))) for bits in (2,5)]
    old=[old11,old13]; examples=[]
    for t,old_t,H in zip(pair,old,(11,13)):
        need(canonical(t)==canonical(old_t),'historical isomorphism merge')
        words=sum(all(t[p[i]]>>p[i+1]&1 for i in range(4)) for p in permutations(range(5)))
        c=[len(cycles(t,k)) for k in (3,4,5)]
        need(hp(t)==words==H and c[0]==4 and c[2]==(H-9)//2,'odd cycle mechanism')
        examples.append({'rows':t,'H':H,'cycles_3_4_5':c,'aut':aut(t),'lambda':lam(H),'mu':mu(H)})
    join_examples=[]
    for i,a in enumerate(pair):
        for j,b in enumerate(pair):
            t=join(a,b); H=hp(t)
            need(H==(11,13)[i]*(11,13)[j],'literal tournament multiplication')
            join_examples.append({'factors':[(11,13)[i],(11,13)[j]],'H':H,'lambda':lam(H),'mu':mu(H)})

    four={}
    for bits in range(64):
        t=frommask(4,bits); four.setdefault(cyclic(t),{})[canonical(t)]=aut(t)
    ambiguous=[v for v in four.values() if len(v)>1]
    need(len(ambiguous)==4 and all(sorted(x.values())==[3,3] for x in ambiguous),'only ambiguous order4 fibre has Aut3')
    s4=presentation(4,2); d=join(s4,(0,)); dc=converse(d)
    need(aut(d)==aut(dc)==1 and canonical(d)!=canonical(dc) and cyclic(d)==cyclic(dc),'rigid size5 ambiguity')
    threshold=[]
    representatives=[c3sub([d]*(3-k)+[dc]*k) for k in range(4)]
    for t in representatives:
        threshold.append({'rows':t,'Aut':aut(t),'H':hp(t),'cyclic_triples':len(cyclic(t))})
    need([x['Aut'] for x in threshold]==[3,1,1,3],'order15 witness automorphisms')
    need(cyclic(tuple(threshold[0]['rows']))==cyclic(tuple(threshold[1]['rows'])),'order15 identical labelled triples')
    need([mu(x['Aut']) for x in threshold]==[-1,1,1,-1],'Mobius first failure control')
    need(len({x['H'] for x in threshold})==1,'order15 H survivor')

    # The C3 realization theorem supplies completeness; generation is independent
    # of ordinary tournament isomorphism tests and predicted orbit sizes.
    child_lifts=[frommask(5,b) for b in range(1024) if cyclic(frommask(5,b))==cyclic(d)]
    need(len(child_lifts)==8,'eight labelled rigid child lifts')
    lifts={t for children in product(child_lifts,repeat=3)
           for base in [c3sub(children)] for t in (base,converse(base))}
    need(len(lifts)==1024 and all(cyclic(t)==cyclic(representatives[0]) for t in lifts),'full fifteen fibre')
    inner=[p for p in permutations(range(5)) if cyclic(relabel(d,p))==cyclic(d)]
    group=[tuple(5*outer[i]+inside[i][j] for i in range(3) for j in range(5))
           for outer in permutations(range(3)) for inside in product(inner,repeat=3)]
    need(len(group)==384 and all(cyclic(relabel(representatives[0],p))==cyclic(representatives[0]) for p in group),'hypergraph group controls')
    orbits=[{relabel(t,p) for p in group} for t in representatives]
    need([len(o) for o in orbits]==[128,384,384,128] and set().union(*orbits)==lifts,'four full fifteen orbits')
    need(sum(Fraction(1,x['Aut']) for x in threshold)==Fraction(len(lifts),len(group))==Fraction(8,3),'fifteen groupoid mass')

    # Restricted history universe: every labelled order5 tournament, no filtering.
    census=Counter(); five={}; aut_cache={}
    for bits in range(1024):
        t=frommask(5,bits); H=hp(t); c3=len(cyclic(t)); c5=len(cycles(t,5))
        need(H==1+2*(c3+c5),'full order5 OCF control')
        census[H]+=1
        code=canonical(t)
        if code not in aut_cache:
            aut_cache[code]=sum(relabel(code,p)==code for p in permutations(range(5)))
            need(aut_cache[code]==aut(code),'independent all-permutation Aut5')
        five.setdefault(cyclic(t),set()).add(code)
    rigid_ambiguous=[(h,v) for h,v in five.items() if len(v)>1 and any(aut_cache[t]==1 for t in v)]
    need(len(rigid_ambiguous)==30 and all(len(v)==2 and all(aut_cache[t]==1 for t in v) for h,v in rigid_ambiguous),'rigid ambiguous fibre5 classification')
    need(all({canonical(d),canonical(dc)}==v for h,v in rigid_ambiguous),'unique rigid ambiguous C3 type at5')
    h9=next(frommask(5,b) for b in range(1024) if hp(frommask(5,b))==9)
    c3=(2,4,1); same_h=[]
    for t in (join(c3,c3),join(h9,(0,))):
        odd=[set(c) for k in (3,5) for c in cycles(t,k)]
        alpha=[1,len(odd),sum(not (a&b) for a,b in combinations(odd,2))]
        need(hp(t)==9==sum(a*2**k for k,a in enumerate(alpha)),'H9 weighted cycle control')
        same_h.append({'rows':t,'H':9,'alpha':alpha})
    need([x['alpha'] for x in same_h]==[[1,2,1],[1,4,0]],'H does not recover cycle collection counts')
    sharp=[]
    for p in (3,5,7,11,13):
        q=tuple(sum(1<<((i+j)%p) for j in range(1,(p+1)//2)) for i in range(p))
        need(aut(q)==p,'carousel quotient exact automorphisms')
        u=substitute(q,[d]*p); v=substitute(q,[dc]+[d]*(p-1))
        need(cyclic(u)==cyclic(v),'sharp squarefree 5p labelled C3 equality')
        au=aut(u); av=aut(v)
        need((au,av)==(p,1),'independent expanded 5p automorphisms')
        sharp.append({'p':p,'vertices':5*p,'Aut':[au,av],'cyclic_triples':len(cyclic(u))})
    necklace=[]
    for n in range(1,14):
        primitive=sum(all(word!=word[k:]+word[:k] for k in range(1,n))
                      for value in range(1<<n) for word in [tuple((value>>i)&1 for i in range(n))])
        exact=sum(mu(d)*2**(n//d) for d in divisors(n))
        need(primitive==exact and exact%n==0,'binary primitive-word census')
        wrong=sum(lam(d)*2**(n//d) for d in divisors(n))
        square_lift=sum(sum(mu(e)*2**((n//(d*d))//e) for e in divisors(n//(d*d)))
                        for d in range(1,isqrt(n)+1) if n%(d*d)==0)
        need(wrong==square_lift,'Liouville inversion adds square period layers')
        necklace.append({'length':n,'primitive_words':exact,'primitive_necklaces':exact//n})
    prime_square_hostile=[]
    for p in (11,13):
        n=p*p; exact=sum(mu(d)*2**(n//d) for d in divisors(n)); wrong=sum(lam(d)*2**(n//d) for d in divisors(n))
        need(exact%n==0 and wrong%n==2,'Liouville does not invert period repetition')
        prime_square_hostile.append({'length':n,'correct_numerator_remainder':0,'Liouville_numerator_remainder':wrong%n})
    report={'status':'FINITE-EXACT, controls of separately proved statements',
            'arithmetic_universe':N,'multiplicative_pairs':10000,'Euler_degree':limit,
            'Euler_coefficients_at_11_13':{'mu':[pm[11],pm[13]],'lambda':[pl[11],pl[13]]},
            'kernel_controls':[kernel_checks(r) for r in range(1,6)],
            'historical_11_13':examples,'actual_join_controls':join_examples,
            'H_digits_do_not_invert_cycle_collections':same_h,
            'mobius_Aut_order15':threshold,'order4_ambiguous_labelled_fibres':len(ambiguous),
            'full_order15_fibre':{'labelled_lifts':len(lifts),'hypergraph_Aut':len(group),'orbit_sizes':[len(o) for o in orbits],'mass':'8/3'},
            'sharp_squarefree_5p_controls':sharp,
            'rigid_ambiguous_labelled_fibres_order5':len(rigid_ambiguous),
            'binary_necklaces':necklace,'prime_square_necklace_hostiles':prime_square_hostile,
            'full_labelled_order5_H_census':dict(sorted(census.items())),
            'arithmetic_values':{str(n):{'factorization':factors(n),'lambda':lam(n),'mu':mu(n)}
                                 for n in [11,13,27,81,121,143,169,5715,17145,462915]}}
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=='__main__': main()
