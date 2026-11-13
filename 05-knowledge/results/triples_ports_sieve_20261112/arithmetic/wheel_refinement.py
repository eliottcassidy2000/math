"""Exact exponent-wheel refinement; no prime inference from a surviving phase.

Universe: every integer exponent k>=1, the first25 primes, and every
three-element subset of those primes. Infinite statements are proved in
WHEEL.md. The finite controls and minimal-period certificates are here.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd, isqrt, lcm, prod
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
RESULTS=HERE.parents[1]

def need(ok, label):
    if not ok: raise RuntimeError(label)

def prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def prime_factors(n):
    out=[]
    for d in range(2,isqrt(n)+1):
        if n%d==0:
            out.append(d)
            while n%d==0:n//=d
    if n>1:out.append(n)
    return out

def order2(p):
    need(p>2 and prime(p),'odd prime order')
    x=1
    for t in range(1,p):
        x=2*x%p
        if x==1:return t
    raise RuntimeError('Fermat order bound')

def clause(p):
    if p==2 or 99%p==0 or p==23:return None
    m=order2(p)
    forbidden=[r for r in range(m) if (99*pow(2,r,p)+23)%p==0]
    need(len(forbidden)<=1,'unique exponent class in cyclic subgroup')
    return dict(prime=p,modulus=m,residue=forbidden[0]) if forbidden else None

def survives(k, clauses):
    return all(k%c['modulus']!=c['residue'] for c in clauses)

def crt(a,m,b,n):
    g=gcd(m,n)
    if (b-a)%g:return None
    q=n//g
    z=0 if q==1 else ((b-a)//g)*pow(m//g,-1,q)%q
    return (a+m*z)%lcm(m,n),lcm(m,n)

def exact_density(clauses):
    value=Fraction(1)
    intersections=0
    for size in range(1,len(clauses)+1):
        for subset in combinations(clauses,size):
            state=(0,1)
            for c in subset:
                state=crt(*state,c['residue'],c['modulus'])
                if state is None:break
            if state is not None:
                intersections+=1
                value+=Fraction((-1)**size,state[1])
    return value,intersections

def minimal_word(clauses):
    period=lcm(*(c['modulus'] for c in clauses)) if clauses else 1
    word=bytes(survives(k,clauses) for k in range(period))
    for q in prime_factors(period):
        while period%q==0:
            trial=period//q
            if word!=word[trial:]+word[:trial]:break
            word=word[:trial];period=trial
    return period,word

def main():
    primes=[p for p in range(2,100) if prime(p)]
    need(len(primes)==25,'first25 primes')
    active=[c for p in primes if (c:=clause(p)) is not None]
    inactive=[p for p in primes if clause(p) is None]
    reduced=[c for c in active if c['prime'] not in (29,37)]
    need(len(active)==11 and len(reduced)==9,'declared clause counts')
    shadowed=[]
    for c in active:
        if c['prime'] in (29,37):
            need(c['modulus']%4==0 and c['residue']%4==3,'complete containment in p5 class')
            shadowed.append(dict(prime=c['prime'],covered_by=5))
    L=lcm(*(c['modulus'] for c in reduced))
    need(L==61209720,'reduced period')
    witnesses=[]
    for q in prime_factors(L):
        k=next(k for k in range(1,10001) if survives(k,reduced)!=survives(k+L//q,reduced))
        witnesses.append(dict(prime_factor=q,shift=L//q,k=k,
                              at_k=survives(k,reduced),at_shift=survives(k+L//q,reduced)))
    private=[]
    for c in reduced:
        other=[d for d in reduced if d!=c]
        k=next(k for k in range(1,10001) if not survives(k,[c]) and survives(k,other))
        private.append(dict(prime=c['prime'],k=k))
    density,intersections=exact_density(reduced)
    need(density==Fraction(5392213,10201620),'independent CRT density')
    survivors=int(density*L)
    need(Fraction(survivors,L)==density,'integral survivor count')

    base=[c for c in reduced if c['prime'] in (5,13,17,19,61)]
    base_phases=[k for k in range(360) if survives(k,base)]
    groups={'free_even':[],'restricted_even':[],'odd':[]}
    for k in base_phases:
        key='odd' if k%2 else ('restricted_even' if k%6==2 else 'free_even')
        groups[key].append(k)
    fibre_sizes={'free_even':11*13*29*41,
                 'restricted_even':10*13*29*41,
                 'odd':11*12*28*40}
    need({key:len(v) for key,v in groups.items()}==
         {'free_even':114,'restricted_even':60,'odd':25},'complete360 phase partition')
    need(sum(len(groups[key])*fibre_sizes[key] for key in groups)==survivors,'factorized count')
    # Each extension prime affects only its declared base phase and coordinate.
    phase_rules={53:(4,1,13,12),59:(2,1,29,6),
                 67:(6,2,11,1),83:(2,1,41,15)}
    for c in reduced:
        if c['prime'] not in phase_rules:continue
        g,r,p,t=phase_rules[c['prime']]
        need(lcm(g,p)==c['modulus'],'split clock modulus')
        for k in range(c['modulus']):
            need((k%c['modulus']==c['residue'])==(k%g==r and k%p==t),
                 'split forbidden cylinder')
    for k in range(1,10001):
        literal=all((99*pow(2,k,p)+23)%p for p in primes)
        need(literal==survives(k,reduced),'literal prime wheel vs compressed wheel')

    # Complete triple universe; canonicalize exact Boolean words, retaining phase.
    canonical={};triple_rows=[];class_members=defaultdict(list);cache={}
    for triple in combinations(primes,3):
        key=tuple(p for p in triple if clause(p) is not None)
        if key not in cache:
            cs=[clause(p) for p in key]
            period,word=minimal_word(cs)
            rational,_=exact_density(cs)
            need(Fraction(sum(word),period)==rational,'triple CRT vs complete word')
            signature=(period,word)
            if signature not in canonical:canonical[signature]=len(canonical)
            cache[key]=(period,sum(word),rational,canonical[signature])
        period,count,rational,index=cache[key]
        triple_rows.append(dict(primes=triple,class_index=index,minimal_period=period,
                                survivors=count,density=[rational.numerator,rational.denominator]))
        class_members[index].append(triple)
    classes=[]
    for (period,word),index in sorted(canonical.items(),key=lambda item:item[1]):
        classes.append(dict(index=index,minimal_period=period,survivors=sum(word),
                            word_sha256=sha256(word).hexdigest(),
                            triple_count=len(class_members[index]),representative=class_members[index][0]))
    need(len(triple_rows)==2300,'all triples')
    fully_free=[row for row in triple_rows if row['density']==[1,1]]
    need(len(fully_free)==364,'choose14 inactive primes,3')

    # Explicit composite progression uses the smaller proven minimal wheel period.
    q=419;order=order2(q);M=lcm(L,order)
    need(M==1162984680 and 99*2**2+23==q,'419 progression')
    need(survives(2,reduced),'419 seed phase survives')
    for t in range(1,101):
        k=2+M*t
        need(survives(k,reduced) and (99*pow(2,k,q)+23)%q==0,'infinite progression controls')
    # Every accepted phase congruent2 mod22 has one419-divisible lift among19.
    need(gcd(L,order)==22 and order//22==19,'conditional419 clock')
    captured=174*13*29*41
    conditional_419_density=Fraction(captured,M)
    extended=reduced+[dict(prime=419,modulus=418,residue=2)]
    extended_density,_=exact_density(extended)
    need(extended_density==density-conditional_419_density,'419 fibre count vs independent CRT')
    extended_witnesses=[]
    for q0 in prime_factors(M):
        k=next(k for k in range(1,10001) if survives(k,extended)!=survives(k+M//q0,extended))
        extended_witnesses.append(dict(prime_factor=q0,shift=M//q0,k=k,
                                      at_k=survives(k,extended),at_shift=survives(k+M//q0,extended)))

    # Exact totient comparisons are finite arithmetic, not an RH verification.
    triple_comparison=[]
    for triple in[(2,3,5),(2,3,7),(5,29,37),(5,13,17),(19,61,67)]:
        row=next(row for row in triple_rows if tuple(row['primes'])==triple)
        n=prod(triple);phi=prod(p-1 for p in triple)
        triple_comparison.append(dict(**row,primorial_like_product=n,phi=phi,
                                      n_over_phi=[Fraction(n,phi).numerator,Fraction(n,phi).denominator]))
    W=prod(primes)
    phiW=prod(p-1 for p in primes)
    ambient=prod(Fraction(p-2,p-1) for p in primes if p not in (2,3,11,23))
    # Independently count primes through the primary-source threshold.
    bound=1591883
    sieve=bytearray(b'\x01')*(bound+1);sieve[0:2]=b'\x00\x00'
    for p in range(2,isqrt(bound)+1):
        if sieve[p]:sieve[p*p:bound+1:p]=b'\x00'*(((bound-p*p)//p)+1)
    need(sum(sieve)==120569 and sieve[bound]==1,'Nicolas120569-th prime')
    inherited=RESULTS/'endpoint_coverage_compiler_20261111/wheel.json'
    report=dict(status='PASS',scope='All integer exponents via proved periodicity; all2300 triples amongfirst25 primes; all360 base phases; explicit finite minimality and irredundancy witnesses.',
        active=active,inactive=inactive,shadowed=shadowed,reduced=reduced,
        minimal_period=L,minimality_witnesses=witnesses,irredundancy_witnesses=private,
        survivors=survivors,density=[density.numerator,density.denominator],
        compatible_CRT_intersections=intersections,base_phase_groups=groups,
        fibre_sizes=fibre_sizes,split_rules=phase_rules,
        triple_count=len(triple_rows),distinct_triple_predicates=len(classes),
        free_triples=len(fully_free),triple_classes=classes,triple_rows=triple_rows,
        triple_comparison=triple_comparison,
        composite_progression=dict(q=q,order=order,exponents=f'2+{M}t, t>=1',
                                  accepted_classes_captured=captured,conditional_lift_denominator=19,
                                  global_density=[conditional_419_density.numerator,conditional_419_density.denominator]),
        with_419=dict(minimal_period=M,minimality_witnesses=extended_witnesses,
                      density=[extended_density.numerator,extended_density.denominator],
                      prime_exception_k=2),
        primorial_25=W,phi_primorial_25=phiW,ambient_affine_unit_density=[ambient.numerator,ambient.denominator],
        nicolas_threshold=dict(index=120569,prime=bound,prime_count=sum(sieve)),
        direct_exponent_controls=10000,
        inherited_wheel_sha256=sha256(inherited.read_bytes()).hexdigest())
    (HERE/'wheel_refinement.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='PASS',minimal_period=L,survivors=survivors,
                          triple_predicates=len(classes),triple_count=len(triple_rows),
                          nicolas_prime=bound,composite_progression=M)))

if __name__=='__main__':main()
