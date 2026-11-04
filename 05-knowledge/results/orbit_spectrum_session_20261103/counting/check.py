#!/usr/bin/env python3
"""Independent finite controls for the moved-mass certified-head calculator."""
import hashlib,importlib.util,json,sys
from collections import Counter
from fractions import Fraction
from itertools import permutations, combinations
from math import factorial,gcd
from pathlib import Path
from heads import head,normalized_prefix,normalized_tail,reconstruct,need,cores
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
BASE=ROOT/'05-knowledge/results/session87_four_vertices_enumeration_20261028/fast_counts.py'
need(hashlib.sha256(BASE.read_bytes()).hexdigest()=='3844a4a590daaf01410cf99d3584441bfa2f8e5e95371ccaafe534fd152b0ece', 'inherited baseline hash')
spec=importlib.util.spec_from_file_location('s87',BASE)
baseline=importlib.util.module_from_spec(spec);spec.loader.exec_module(baseline)


def partitions(n,minimum=1):
    if n==0:
        yield ()
    for k in range(minimum,n+1,2):
        for rest in partitions(n-k,k):
            yield (k,)+rest


def davis(n):
    value=Fraction(0)
    for part in partitions(n):
        mult=Counter(part)
        exponent=sum((k-1)//2 for k in part)+sum(gcd(a,b) for a,b in combinations(part,2))
        denominator=1
        for k,m in mult.items():
            denominator*=k**m*factorial(m)
        value+=Fraction(1<<exponent,denominator)
    need(value.denominator==1,'independent Davis integrality')
    return value.numerator


def unpack(record):
    return (Fraction(int(record['lower_numerator_hex'],16),int(record['lower_denominator_hex'],16)),
            Fraction(int(record['upper_numerator_hex'],16),int(record['upper_denominator_hex'],16)))


def summary(record):
    return {k:v for k,v in record.items() if k not in (
        'lower_numerator_hex','lower_denominator_hex','upper_numerator_hex','upper_denominator_hex')}


def pair_orbit_controls(nmax=8):
    checks=0
    # Enumerate each odd cycle type, form its literal permutation and count
    # actual unordered pair orbits; no gcd-energy formula is used here.
    for n in range(nmax+1):
        for part in partitions(n):
            perm=list(range(n));offset=0
            for cycle in part:
                for i in range(cycle):
                    perm[offset+i]=offset+(i+1)%cycle
                offset+=cycle
            unseen=set(combinations(range(n),2)); sizes=[]
            while unseen:
                a=next(iter(unseen));orb=set();b=a
                while b not in orb:
                    orb.add(b);b=tuple(sorted((perm[b[0]],perm[b[1]])))
                unseen-=orb;sizes.append(len(orb))
            expected=sum((k-1)//2 for k in part)+sum(gcd(a,b) for a,b in combinations(part,2))
            need(len(sizes)==expected,'literal pair orbit count')
            if 1 not in part:
                need(all(s>=3 for s in sizes),'non-singleton odd core pair orbits >=3')
                need(len(sizes)<=n*(n-1)//6,'core energy bound')
            checks+=1
    return checks


def main():
    prefix,stats=baseline.degree_bucket_dp(200)
    need(all(prefix[n]==davis(n) for n in range(41)),'independent Fraction counts through40')
    head_checks=0;table=[]
    for n in range(201):
        # Full queries at small orders and several real prefix resolutions.
        for bits in (1,16,64,128):
            rec=head(n,bits,max_core_mass=40)
            expected_shift=max(0,prefix[n].bit_length()-bits)
            need(rec['bit_length']==prefix[n].bit_length(),'exact bit length')
            need(rec['discarded_bits']==expected_shift,'exact truncation shift')
            need(int(rec['head_hex'],16)==prefix[n]>>expected_shift,'exact leading bits')
            lower,upper=unpack(rec)
            scale=1 << (max(0,n*(n-1)//2-factorial(n).bit_length()+1-bits))
            need(lower<=Fraction(prefix[n],scale)<=upper,'actual count in interval')
            head_checks+=1
            if bits==128 and n in (20,32,40,64,100,128,160,200):
                table.append(summary(rec))
    tail_checks=0
    for n in range(13,65):
        actual=Fraction(prefix[n]*factorial(n),1 << (n*(n-1)//2))
        for cutoff in range(min(n,18)+1):
            known,count=normalized_prefix(n,cutoff)
            bound=normalized_tail(n,cutoff)
            need(0<=actual-known<=bound,'uniform explicit remainder bound')
            tail_checks+=1
    pairchecks=pair_orbit_controls()
    # Hostile: a tiny normalized error does not authorize full-integer rounding.
    identity=Fraction(1 << (100*99//2),factorial(100))
    need(prefix[100]-identity>1,'identity-only full rounding is false')
    insufficient=False
    try:
        head(40,128,max_core_mass=0)
    except ArithmeticError:
        insufficient=True
    need(insufficient,'insufficient precision cap must reject')
    # Exact high/low reconstruction is valid only after checking interval width.
    reconstructions=0
    for n in range(13,81):
        rec=head(n,64,max_core_mass=40)
        lo,hi=unpack(rec)
        shift=max(0,n*(n-1)//2-factorial(n).bit_length()+1-64)
        lo*=1<<shift;hi*=1<<shift
        width=hi-lo
        modulus=1
        while modulus<=width:
            modulus*=2
        need(reconstruct(lo,hi,prefix[n]%modulus,modulus)==prefix[n], 'independent residue reconstruction')
        reconstructions+=1
    bad_width=False
    try:
        reconstruct(0,4,1,4)
    except ValueError:
        bad_width=True
    need(bad_width,'reject nonunique high/low interval')
    large=[]
    for n in (1000,10000):
        rec=head(n,256,max_core_mass=0)
        lo,hi=unpack(rec)
        need(lo.numerator//lo.denominator==hi.numerator//hi.denominator,'large rational certificate')
        need(rec['core_types']==1,'identity suffices at large n')
        # A second algebraic path: division of 2^(factorial_bits+255) by n!,
        # with the independent simple error bound; no recursion/cores involved.
        f=factorial(n);q=n*(n-1)//2
        shift=q-f.bit_length()+1-256
        num=1 << (q-shift)
        direct,rem=divmod(num,f)
        rr=Fraction(n,2**((n-1)//3))
        tail=rr**3/(1-rr)
        need(Fraction(rem,f)+Fraction(num,f)*tail<1,'distance to next integer dominates tail')
        need(direct==int(rec['head_hex'],16),'independent identity quotient')
        large.append(summary(rec))
    output={'status':'PASS','universe':{
        'exact_prefix':200,'independent_Davis':40,'head_queries':head_checks,
        'tail_checks':tail_checks,'literal_pair_cycle_types':pairchecks,
        'hybrid_reconstruction_controls':reconstructions},
        'scope':'Certified leading bits; full-integer spectrum coverage remains open.',
        'head_table':table,'large_heads':large,
        'controls':{'identity_full_rounding_rejected':True,'insufficient_cap_rejected':True,'nonunique_modular_interval_rejected':True},
        'baseline_work':{k:v for k,v in stats.items() if k!='rows'},
        'source_hashes':{str(BASE.relative_to(ROOT)):hashlib.sha256(BASE.read_bytes()).hexdigest(),
                         'heads.py':hashlib.sha256((HERE/'heads.py').read_bytes()).hexdigest()}}
    print(json.dumps(output,indent=2,sort_keys=True))
if __name__=='__main__':
    main()
