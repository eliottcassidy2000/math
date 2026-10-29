#!/usr/bin/env python3
"""Exact guards; no numerical/asymptotic theorem is inferred from samples."""
import itertools, json, math
from collections import Counter
from fractions import Fraction
from pathlib import Path

def primes(n):
    sieve=bytearray(b'\1')*(n+1); sieve[:2]=b'\0\0'
    for p in range(2,math.isqrt(n)+1):
        if sieve[p]: sieve[p*p:n+1:p]=b'\0'*((n-p*p)//p+1)
    return [i for i in range(2,n+1) if sieve[i]]

def chi(a,p):
    v=pow(a%p,(p-1)//2,p)
    return -1 if v==p-1 else v

def factor(n):
    n=abs(n); out={}; p=2
    while p*p<=n:
        while n%p==0: out[p]=out.get(p,0)+1; n//=p
        p+=1
    if n>1: out[n]=out.get(n,0)+1
    return out

def squareclass_code(bases):
    fac=[factor(a) for a in bases]
    cols=sorted(set().union(*(f.keys() for f in fac)))
    rows=[sum((f.get(p,0)%2)<<i for i,p in enumerate(cols)) for f in fac]
    allowed={tuple((-1)**((a<0)+(r&z).bit_count()) for a,r in zip(bases,rows)) for z in range(1<<len(cols))}
    return cols,allowed

def jacobi(a,n):
    assert n>0 and n%2
    ans=1; a%=n
    while a:
        while a%2==0:
            a//=2
            if n%8 in (3,5): ans=-ans
        a,n=n,a
        if a%4==n%4==3: ans=-ans
        a%=n
    return ans if n==1 else 0

def order(a,p):
    d=p-1
    for q in factor(d):
        while d%q==0 and pow(a,d//q,p)==1: d//=q
    return d

def v2(n):
    return (n&-n).bit_length()-1

def main():
    report={}; ps=primes(100000); checked=0
    for p in primes(300)[1:]:
        for a in range(1,p):
            d=order(a,p)
            assert (chi(a,p)==1)==((p-1)//2%d==0)
            assert (chi(a,p)==-1)==(v2(d)==v2(p-1))
            checked+=1
    report['order_gate_checks']=checked
    report['wieferich_examples']=[]
    for p,d in [(1093,364),(3511,1755)]:
        assert p in ps and order(2,p)==d and pow(2,p-1,p*p)==1
        report['wieferich_examples'].append(dict(p=p,order_2=d,character_2=chi(2,p),fermat_residue_mod_p_squared=1))
    cases=[([2,3,6],4),([4,2,5,3,7,13],32),([2,3,7,5,2],16),([-1,2,-2],2),([1093,3511,13,139,19],32)]
    report['codes']=[]
    for bases,size in cases:
        cols,allowed=squareclass_code(bases); assert len(allowed)==size
        counts=Counter(tuple(chi(a,p) for a in bases) for p in ps if p%4==3 and all(a%p for a in bases))
        assert set(counts)<=allowed
        assert set(counts)==allowed # finite coverage only
        Q=8*math.prod(p for p in cols if p!=2)
        crt=None
        if Q<100000:
            residue_counts=Counter(tuple(jacobi(a,q) for a in bases) for q in range(3,Q,4) if math.gcd(q,Q)==1)
            # For composite odd q, Jacobi(a/q) is the correct Dirichlet character evaluation.
            assert set(residue_counts)==allowed
            assert len(set(residue_counts.values()))==1
            crt={'modulus':Q,'classes_per_color':next(iter(residue_counts.values()))}
        report['codes'].append({'bases':bases,'allowed_colors':size,'crt_exact':crt,'prime_bound':100000,'prime_count_min':min(counts.values()),'prime_count_max':max(counts.values())})
    assert chi(1,7)==chi(4,7)==1 and chi(4,7)==1 and chi(13,7)==-1
    report['additive_quotient_hostile']={'p':7,'n':[1,4],'chi_n':[1,1],'chi_3n_plus_1':[1,-1]}
    # Three-vertex labelled tournament universe: six transitive, two cycles.
    universe=Counter()
    for bits in itertools.product((0,1),repeat=3):
        edges={(i,j) if b else (j,i) for (i,j),b in zip([(0,1),(0,2),(1,2)],bits)}
        degrees=[sum(x==i for x,y in edges) for i in range(3)]
        universe['cycle' if degrees==[1,1,1] else 'transitive']+=1
    assert universe=={'cycle':2,'transitive':6}
    assert Fraction(universe['cycle'],8)==Fraction(math.factorial(3),3*8)
    assert Fraction(universe['transitive'],8)==Fraction(math.factorial(3),1*8)
    report['labelled_three_vertex_universe']=dict(universe)
    roots=[3,7,43]
    assert 1807==13*139
    for p in roots:
        assert chi(13,p)==chi(139,p)
        assert chi(p,13)==-chi(p,139)
    report['sylvester_split_port_checks']=[{'root':p,'incoming_13':chi(13,p),'incoming_139':chi(139,p),'outgoing_13':chi(p,13),'outgoing_139':chi(p,139)} for p in roots]
    report['status']='all exact assertions passed; prime sample counts are diagnostics only'
    Path(__file__).with_name('verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__': main()
