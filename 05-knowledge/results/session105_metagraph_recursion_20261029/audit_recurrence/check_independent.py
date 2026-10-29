#!/usr/bin/env python3
"""Independent pair-orbit partition oracle; imports targets only for comparison."""
from collections import Counter
from fractions import Fraction as Q
from math import factorial, gcd, comb
from pathlib import Path
import hashlib,json,sys
sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'enumeration'))
import integer_response as target
import baseline87, baseline98
sys.path.insert(0,str(ROOT/'updown'))
import compressed_smith

def check(ok,msg):
    if not ok: raise RuntimeError(msg)

def parts(n,cap=None):
    if n==0:
        yield (); return
    cap=min(n,n if cap is None else cap)
    for p in range(cap,0,-1):
        if p%2:
            for rest in parts(n-p,p): yield (p,)+rest

def orbit_energy(mu):
    """Literal permutation action on unordered pairs, no gcd energy formula."""
    perm=[]; start=0
    for p in mu:
        perm.extend(range(start+1,start+p));perm.append(start);start+=p
    unseen={(i,j) for j in range(start) for i in range(j)}
    count=0
    while unseen:
        pair=next(iter(unseen)); count+=1
        while pair in unseen:
            unseen.remove(pair);pair=tuple(sorted((perm[pair[0]],perm[pair[1]])))
    return count

def central(mu):
    return __import__('functools').reduce(lambda x,pk:x*pk[0]**pk[1]*factorial(pk[1]),Counter(mu).items(),1)
def oddfact(n):
    f=factorial(n)
    return f//(f&-f)
def pow2(e):return Q(2**e) if e>=0 else Q(1,2**-e)
def weight(mu,r,odd):
    d=sum(mu); w=pow2(orbit_energy(mu)+r*len(mu))/central(mu)
    return w*(oddfact(d) if odd else factorial(d+r)//factorial(r))

def literal_fixed_sourcefree(mu):
    n=sum(mu); perm=[];start=0
    for p in mu:
        perm.extend(range(start+1,start+p));perm.append(start);start+=p
    edges=[(i,j) for j in range(n) for i in range(j)]
    index={pair:k for k,pair in enumerate(edges)}
    fixed=0; sourcefree=0
    for mask in range(1<<len(edges)):
        def beats(i,j):
            bit=(mask>>index[tuple(sorted((i,j)))])&1
            return bool(bit) if i<j else not bit
        if any(beats(i,j)!=beats(perm[i],perm[j]) for i,j in edges):continue
        fixed+=1
        if not any(all(i==j or beats(i,j) for j in range(n)) for i in range(n)):
            sourcefree+=1
    return fixed,sourcefree

def main():
    values=[]; partition_count=0; gauges=0; negatives=0; carries=0
    for n in range(31):
        total=Q(0)
        for mu in parts(n):
            partition_count+=1
            total+=Q(2**orbit_energy(mu),central(mu))
            if n>22:continue
            for p in range(n if n%2 else n-1,2,-2):
                past=tuple(a for a in mu if a>p); k=mu.count(p)
                r=n-sum(past); rr=r-k*p
                after=past+(p,)*k
                s=sum(gcd(a,p)-p for a in past)
                E=k*(p-1)//2+p*k*(k-1)//2+k*s+rr*k
                for odd in (False,True):
                    old=weight(past,r,odd); new=weight(after,rr,odd)
                    if odd:
                        v=factorial(k); twos=(v&-v).bit_length()-1
                        factor=Q(oddfact(n-rr),oddfact(n-r)*p**k*oddfact(k)); exp=E-twos
                    else:factor=Q(factorial(r),factorial(rr)*p**k*factorial(k));exp=E
                    check(old.denominator==new.denominator==factor.denominator==1,'integer gauge')
                    check(old*factor*pow2(exp)==new,'history transfer')
                    check(target.exact_shift(int(old*factor),exp)==new,'checked shift')
                    gauges+=1;negatives+=exp<0
                width=max(1,n.bit_length()); js=list(range(3,min(p-2,rr)+1,2))
                ds=[sum(j-gcd(a,j) for a in past) for j in js]
                inc=[k*(j-gcd(p,j)) for j in js]
                packed=sum(d<<(width*i) for i,d in enumerate(ds))+sum(d<<(width*i) for i,d in enumerate(inc))
                expected=[sum(j-gcd(a,j) for a in after) for j in js]
                check(all(0<=d<=n for d in expected),'carry bound')
                check([(packed>>(width*i))&((1<<width)-1) for i in range(len(js))]==expected,'packed addition')
                carries+=1
        check(total.denominator==1,'Burnside integer');values.append(int(total))
        for tail in (1,3):
            for cache in (False,True):check(target.count(n,tail,cache)[0]==total,'tuple target')
        for odd in (False,True):check(target.count_packed(n,odd)[0]==total,'packed target')
    check(negatives>0,'negative shifts exercised')
    try:target.exact_shift(25,-3)
    except ArithmeticError:pass
    else:raise RuntimeError('bad division accepted')
    check(orbit_energy((9,))==4 and central((9,))==9,'composite-nine control')
    prefix,_=baseline87.degree_bucket_dp(200)
    check(prefix[:31]==values,'actual Session87 prefix')
    large=[]
    for n in (40,60,100,150,200):
        a,stats=target.count_packed(n)
        check(a==prefix[n]==target.count_packed(n,False)[0],'large packed gauges')
        if n<=100:
            old,rows=baseline98.fixed_degree(n)
            check(a==old and stats['stages']==rows[:-1],'actual Session98 state rows')
        large.append({'n':n,'sha256_decimal':hashlib.sha256(str(a).encode()).hexdigest()})
    # Independent finite core coefficient extraction for B and transformed harmonic series.
    B=[Q(factorial(n)*values[n],2**comb(n,2)) for n in range(31)]
    eta=[Q(factorial(n)*(values[n]-(values[n-1] if n else 0)),2**comb(n,2)) for n in range(31)]
    coreB=[Q(0) for _ in B]; coreEta=[Q(0) for _ in B]
    cores=0; banks={}
    for m in range(31):
        for mu in parts(m):
            if 1 in mu:continue
            cores+=1;c=len(mu);q=pow2(c-m)
            banks.setdefault((m-c)//2,[]).append(mu)
            C=Q(factorial(m),central(mu))*pow2(orbit_energy(mu)-comb(m,2))
            for n in range(m,31):
                term=C*comb(n,m)*q**(n-m);coreB[n]+=term;coreEta[n]+=term
                if n<30:coreEta[n+1]-=(n+1)*pow2(-n)*term
                f=n-m
                check(f==0 or f*pow2(1-c-f)<=1,'paired core Taylor positivity')
    check(coreB==B and coreEta==eta,'harmonic core-pair coefficients')
    for k in range(11):
        top=max(map(sum,banks[k]))
        check(top==3*k and [mu for mu in banks[k] if sum(mu)==top]==[(3,)*k], 'unique top pole core')
        A=Q(factorial(3*k),3**k*factorial(k))*pow2(3*k*k+k)
        check(A>0 and -(3*k+1)*(2*4**k)*A<0,'both pole amplitudes nonzero')
    for n in range(1,31):check(eta[n]==B[n]-n*pow2(1-n)*B[n-1],'harmonic transform')
    sourcefree_counts=[]
    for n in range(7):
        total=Q(0)
        for mu in parts(n):
            fixed,sf=literal_fixed_sourcefree(mu)
            f=mu.count(1); c=len(mu)-f
            ratio=Q(0) if f==0 else f*pow2(1-c-f)
            check(fixed==2**orbit_energy(mu),'literal fixed orientations')
            check(sf==fixed*(1-ratio),'literal source-free Burnside filter')
            total+=Q(sf,central(mu))
        check(total==values[n]-(values[n-1] if n else 0),'source-free class count')
        sourcefree_counts.append(int(total))
    # Inherited integer diagonal, independently converted using integer gcd/lcm
    # for small ranks and factorized valuation arrays through n=10.
    def factorize(d):
        result={};p=2
        while p*p<=d:
            while d%p==0:result[p]=result.get(p,0)+1;d//=p
            p+=1
        if d>1:result[d]=result.get(d,0)+1
        return result
    factors={d:factorize(d) for d in range(1,201)}
    primes=sorted({p for fs in factors.values() for p in fs})
    smith_summary=[]
    for n in range(1,201):
        h=[prefix[r]-(prefix[r-1] if r else 0) for r in range(n)]
        histogram,bound=compressed_smith.compressed_smith(n,prefix)
        check(sum(mult for d,mult in histogram)==prefix[n-1],'compressed Smith rank')
        check(all(mult>0 for d,mult in histogram),'positive interval lengths')
        check(all(b%a==0 for (a,_),(b,_) in zip(histogram,histogram[1:])),'Smith divisibility')
        hist_factors=[(factorize(d),mult) for d,mult in histogram]
        independent_bound=1+sum(max(factors[d].get(p,0) for d in range(1,n+1)) for p in primes if p<=n)
        check(len(histogram)<=bound==independent_bound,'Smith interval bound')
        for p in primes:
            original_order=sum(h[r]*factors[n-r].get(p,0) for r in range(n))
            smith_order=sum(mult*fs.get(p,0) for fs,mult in hist_factors)
            check(original_order==smith_order,'Smith determinant-order exponent')
            original_rank=sum(h[r] for r in range(n) if (n-r)%p)
            smith_rank=sum(mult for d,mult in histogram if d%p)
            check(original_rank==smith_rank,'Smith finite-field rank')
        if n<=10:
            diagonal=[n-r for r in range(n) for _ in range(h[r])]
            expanded=[1]*len(diagonal)
            for p in primes:
                if p>n:break
                for i,v in enumerate(sorted(factors[d].get(p,0) for d in diagonal)):expanded[i]*=p**v
            check(sorted(Counter(expanded).items())==histogram,'expanded Smith arrays')
            if n<=8:
                for i in range(len(diagonal)):
                    for j in range(i+1,len(diagonal)):
                        a,b=diagonal[i],diagonal[j];g=gcd(a,b)
                        diagonal[i],diagonal[j]=g,a//g*b
                check(diagonal==expanded,'independent gcd-lcm Smith reduction')
        if n in (5,10,25,50,100,200):smith_summary.append({'n':n,'intervals':len(histogram),'bound':bound,'rank_digits':len(str(prefix[n-1]))})
    atlas_path=ROOT/'updown/smith_atlas_200.json'
    check(atlas_path.read_bytes()==(ROOT/'updown/smith_atlas_200.optimized.json').read_bytes(),'atlas normal optimized bytes')
    atlas=json.loads(atlas_path.read_text())
    check(atlas['count_prefix_sha256']==hashlib.sha256('\n'.join(map(str,prefix)).encode()).hexdigest(),'atlas count prefix')
    for row in atlas['atlas']:
        n=row['n']; expected,bound=compressed_smith.compressed_smith(n,prefix)
        actual=[(int(entry['factor']),int(entry['multiplicity'])) for entry in row['invariant_factors']]
        check(actual==expected and int(row['rank'])==prefix[n-1] and row['interval_bound']==bound,'stored atlas factors')
        exponents={str(p):str(sum((prefix[r]-(prefix[r-1] if r else 0))*factors[n-r].get(p,0) for r in range(n))) for p in primes if p<=n}
        check(row['torsion_order_prime_factorization']=={p:e for p,e in exponents.items() if e!='0'},'stored atlas order')
    pins=json.loads((ROOT/'enumeration/SOURCE_PINS.json').read_text())
    for copied,meta in pins.items():
        copy=ROOT/'enumeration'/copied
        check(hashlib.sha256(copy.read_bytes()).hexdigest()==meta['sha256'],'baseline hash')
    files=['enumeration/integer_response.py','enumeration/INTEGER_PROJECTIVE.md','enumeration/SOURCE_PINS.json','harmonic/HARMONIC_GF.md','updown/UPDOWN_HARMONIC_AUDIT.md','updown/compressed_smith.py','updown/smith_atlas_200.json']
    hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files}
    print(json.dumps(dict(status='PASS',oracle='literal unordered-pair orbit counts per odd partition',partition_through=30,partitions=partition_count,gauge_history_through=22,gauge_checks=gauges,negative_history_transitions=negatives,packed_field_checks=carries,large=large,harmonic_coefficient_through=30,cores=cores,smith_all_degrees_through=200,smith_summary=smith_summary,sourcefree_class_counts=sourcefree_counts,values=values,source_hashes=hashes),indent=2))
if __name__=='__main__':main()
