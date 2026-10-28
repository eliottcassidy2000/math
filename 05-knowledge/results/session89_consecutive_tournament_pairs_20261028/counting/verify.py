#!/usr/bin/env python3
"""Independent exact partitions, literal orbit universes, inherited baselines."""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import permutations
import json
from math import factorial,gcd
from pathlib import Path
import sys

from baseline_sources import load, SOURCES, COMMIT
from channels import scalar_prefix,fused_prefix


def check(condition,message):
    if not condition:
        raise AssertionError(message)


def partitions(n, minimum=1):
    if n == 0:
        yield ()
    for k in range(minimum,n+1):
        for tail in partitions(n-k,k):
            yield (k,)+tail


def partition_data(parts):
    length=len(parts)
    e=sum((k-1)//2 for k in parts)+sum(gcd(a,b) for i,a in enumerate(parts) for b in parts[i+1:])
    z=1
    for k,m in Counter(parts).items():
        z*=k**m*factorial(m)
    return e,length,z


def partition_count(n,r,s):
    total=Fraction()
    for parts in partitions(n):
        if all(k%2 for k in parts):
            e,t,z=partition_data(parts)
            total+=Fraction(2**(r*e+s*t),z)
    check(total.denominator == 1,'partition count nonintegral')
    return int(total)


def matrix_count(n):
    total=Fraction()
    for parts in partitions(n):
        z=partition_data(parts)[2]
        fixedbits=sum(gcd(a,b) for a in parts for b in parts)
        total+=Fraction(2**fixedbits,z)
    check(total.denominator == 1,'matrix count nonintegral')
    return int(total)


def literal_maps(n):
    edges=[(i,j) for i in range(n) for j in range(i+1,n)]
    pos={pair:i for i,pair in enumerate(edges)}
    maps=[]
    for perm in permutations(range(n)):
        tm=[]
        for code in range(1<<len(edges)):
            image=0
            for index,(i,j) in enumerate(edges):
                a,b=perm[i],perm[j]
                bit=((code>>index)&1) ^ (a>b)
                if bit:
                    image|=1<<pos[min(a,b),max(a,b)]
            tm.append(image)
        vm=[]
        for code in range(1<<n):
            vm.append(sum(((code>>i)&1)<<perm[i] for i in range(n)))
        maps.append((tm,vm))
    return edges,maps


def literal_orbits(n,r,s):
    edges,maps=literal_maps(n)
    edge_bits=len(edges)
    shift=edge_bits*r
    universe=1<<(shift+n*s)
    seen=bytearray(universe)
    orbits=0
    sizes=Counter()
    for code in range(universe):
        if seen[code]:
            continue
        orbits+=1
        orbit=set()
        for tm,vm in maps:
            image=0
            for layer in range(r):
                image|=tm[(code>>(edge_bits*layer))&((1<<edge_bits)-1)]<<(edge_bits*layer)
            for bit in range(s):
                image|=vm[(code>>(shift+n*bit))&((1<<n)-1)]<<(shift+n*bit)
            seen[image]=1
            orbit.add(image)
        sizes[len(orbit)]+=1
    check(sum(k*v for k,v in sizes.items()) == universe,'literal universe incomplete')
    return {'n':n,'r':r,'s':s,'labelled_universe':universe,'orbits':orbits,
            'orbit_size_histogram':dict(sorted(sizes.items()))}


def literal_matrix_orbits(n):
    perms=list(permutations(range(n)))
    universe=1<<(n*n)
    seen=bytearray(universe)
    count=0
    for code in range(universe):
        if seen[code]:
            continue
        count+=1
        for perm in perms:
            image=0
            for i in range(n):
                for j in range(n):
                    image|=((code>>(n*i+j))&1)<<(n*perm[i]+perm[j])
            seen[image]=1
    return count


def hostile_merge():
    # At the length-one boundary, both cores have degree 10 and response t=2.
    cores=[(3,7),(5,5)]
    entries=[partition_data(parts) for parts in cores]
    check(entries == [(5,2,21),(9,2,50)],'unexpected hostile data')
    m0=sum((Fraction(1,z) for e,t,z in entries),Fraction())
    m1=sum((Fraction(2**e,z) for e,t,z in entries),Fraction())
    m2=sum((Fraction(2**(2*e),z) for e,t,z in entries),Fraction())
    defect=m0*m2-m1*m1
    check(defect == Fraction(1536,7),'wrong weighted Frobenius defect')
    check(m2 != m1*m1/m0,'squaring merged first moment unexpectedly works')
    check(4*m2 != 4*m1*m1/m0,'root attachment factor repairs nothing')
    return {'cores':[list(c) for c in cores], 'degree':10,'terminal_response':2,
            'energy':[5,9],'z':[21,50],
            'zeroth_moment':str(m0),'first_moment':str(m1),'second_moment':str(m2),
            'weighted_frobenius_defect':str(defect)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-root',type=Path)
    parser.add_argument('--receipt',type=Path)
    args=parser.parse_args()
    old=load('session22',args.baseline_root)
    s87=load('session87',args.baseline_root)
    channels=((1,0),(1,1),(2,0),(2,1),(3,2))
    values,stats=fused_prefix(30,channels)
    for i,(r,s) in enumerate(channels):
        independent=[partition_count(n,r,s) for n in range(31)]
        scalar,scalar_stats=scalar_prefix(30,r,s)
        historical,old_stats=old.channel_dp(30,r,s)
        check(values[i] == independent == scalar == historical,f'channel {(r,s)} disagreement')
        check(stats['rows'] == scalar_stats['rows'] == old_stats['rows'],'state rows differ')
    # The strongest scalar baseline, not just the old Session 22 channel code.
    for n in (0,1,2,3,7,10,11,30,80,120,160):
        pair,pstats=fused_prefix(n)
        a,astats=s87.degree_bucket_dp(n)
        b,bstats=scalar_prefix(n,2,1)
        check(pair == [a,b],f'fused pair mismatch {n}')
        check(pstats['rows']==astats['rows']==bstats['rows'],'pair state graph changed')
    for s in (0,1,2,3):
        value,st=scalar_prefix(80,1,s)
        baseline,bst=s87.degree_bucket_dp(80,colors=2**s)
        check(value==baseline and st['rows']==bst['rows'],'generalization changed S87')
    # Terminal-fold ablation and ordinary/no-fold state rows agree.
    for r,s in channels:
        folded,st=scalar_prefix(30,r,s)
        unf,ust=scalar_prefix(30,r,s,terminal_fold=False)
        check(folded==unf and st['rows']==ust['rows'],'terminal fold mismatch')
    # Empty prefix, one channel, arbitrary pair, duplicates, and generic batch.
    for batch in (((2,1),),((3,2),(1,1)),((2,1),(2,1)),channels):
        for n in (0,1,12):
            fused,_=fused_prefix(n,batch)
            check(fused==[scalar_prefix(n,*c)[0] for c in batch],'batch edge case')
    rejected=0
    for bound,batch in ((True,[(1,0)]),(-1,[(1,0)]),(2,[]),(2,[(0,1)]),
                        (2,[(1,-1)]),(2,[(True,0)]),(2,[(1,False)]),
                        (2,[(1,)]),(2,[('1',0)])):
        try:
            fused_prefix(bound,batch)
        except ValueError:
            rejected+=1
        else:
            raise AssertionError(f'invalid input accepted {bound} {batch}')
    literal=[]
    for n in range(5):
        for r,s in ((1,0),(1,1),(2,1)):
            row=literal_orbits(n,r,s)
            check(row['orbits']==partition_count(n,r,s),'literal tournament orbit mismatch')
            literal.append(row)
    for n in range(4):
        row=literal_orbits(n,3,2)
        check(row['orbits']==partition_count(n,3,2),'literal multilayer mismatch')
        literal.append(row)
    matrices=[]
    for n in range(5):
        value=literal_matrix_orbits(n)
        check(value==matrix_count(n),'literal matrix mismatch')
        matrices.append(value)
    check(matrices[2]==10 and partition_count(2,2,1)==8,'matrix hostile boundary')
    # For odd cycle types only: ordered matrix bits = 2E + length.
    oddtypes=0
    for n in range(16):
        for parts in partitions(n):
            if all(k%2 for k in parts):
                e,t,z=partition_data(parts)
                check(sum(gcd(a,b) for a in parts for b in parts)==2*e+t,'odd matrix identity')
                oddtypes+=1
    pair=fused_prefix(12)[0][1]
    rooted=scalar_prefix(12,1,1)[0]
    a=scalar_prefix(12)[0]
    check(pair[2]==8 != a[2]*rooted[2]==4,'independent class-product hostile')
    receipt={'status':'PASS','python':sys.version,'optimization_level':sys.flags.optimize,
             'verification_driver_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
             'baseline_commit':COMMIT,
             'baseline_sources':SOURCES,
             'code_sha256':sha256((Path(__file__).parent/'channels.py').read_bytes()).hexdigest(),
             'checks':{'partition_channels':list(channels),'partition_prefix_bound':30,
                       'fused_vs_strongest_prefix_bound':160,'s87_colored_agreement_bound':80,
                       'invalid_inputs_rejected':rejected,'odd_matrix_identity_types':oddtypes},
             'literal_tournament_universes':literal,'literal_matrix_counts_0_to_4':matrices,
             'paired_counts_0_to_12':pair,'rooted_next_counts_0_to_12':rooted,
             'merge_hostile':hostile_merge()}
    output=json.dumps(receipt,indent=2)+'\n'
    if args.receipt:
        args.receipt.write_text(output)
    print(output)


if __name__ == '__main__':
    main()
