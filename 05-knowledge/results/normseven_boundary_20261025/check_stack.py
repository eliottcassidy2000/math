#!/usr/bin/env python3
"""Literal finite action-groupoid decoder and its coarse-fibre loss."""
from collections import Counter
from itertools import combinations, product
from pathlib import Path
import hashlib
import json

GATES=0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(why)


def span(vectors):
    result={0}
    for v in vectors:
        result |= {x ^ v for x in tuple(result)}
    return frozenset(result)


def subspaces(n):
    out={frozenset({0})}
    for size in range(1,n+1):
        out |= {span(v) for v in combinations(range(1,1<<n),size)}
    return out


def coset(x,I):
    return min(x ^ i for i in I)


def check(n,H,K,I):
    G=range(1<<n)
    decompose={h ^ k:(h,k) for h,k in product(H,K)}
    need(len(decompose) == 1<<n and len(H)*len(K) == 1<<n, 'complementary kernels')
    X={coset(g,I) for g in G}
    for x,g in product(X,G):
        decoded=coset(x ^ decompose[g][0],I)
        for h,k in product(H,K):
            changed=coset(x ^ h ^ decompose[g ^ h ^ k][0],I)
            need(decoded == changed, 'decoder is invariant under both representative changes')
        h,k=decompose[g]
        need(g ^ h ^ k == 0 and decoded == coset(x ^ h,I), 'canonical representative exists')
    need({coset(x ^ decompose[0][0],I) for x in X} == X, 'all sources recover')
    if H & I == {0} and K & I == {0}:
        HI={h ^ i for h,i in product(H,I)}
        KI={k ^ i for k,i in product(K,I)}
        counts=Counter((coset(x,HI),coset(x,KI)) for x in X)
        need(set(counts.values()) == {len(I)}, 'coarse product loses exactly inertia-sized fibres')
        need(len(counts)*len(I) == len(X), 'stack completion restores actual sources')
    return len(X)


def main():
    cases=0
    for n in (1,2,3):
        spaces=subspaces(n)
        for H,K in product(spaces,repeat=2):
            if len(H)*len(K) != 1<<n or H & K != {0}:
                continue
            for I in spaces:
                check(n,H,K,I)
                cases += 1
    def kernel(labels):
        chars=[sum(((a & x).bit_count()%2)<<i for i,x in enumerate(labels[1:])) for a in range(8)]
        return {g for g in range(64) if all((g & c).bit_count()%2 == 0 for c in chars)}
    H=kernel((1,2,4,5,6,3,7))
    K=kernel((1,2,4,6,3,7,5))
    triples=[{0,1,2},{0,3,4},{0,5,6}]
    strata=[set()]+[{i} for i in range(7)]+[set(t) for t in combinations(range(7),2)
             if not any(set(t)<=z for z in triples)]+triples
    counts=Counter()
    for zero in strata:
        I=span(63 if i==0 else 1<<(i-1) for i in zero)
        need(H & I == {0} and K & I == {0}, 'actual norm-seven kernels are free')
        counts[len(zero),check(6,H,K,I)] += 1
    # If H=K={0,1} in F2, the map to both quotients has a kernel.
    noncomplementary={0,1}
    need(len({(coset(g,noncomplementary),coset(g,noncomplementary)) for g in range(2)}) != 2,
         'noncomplementary group map hostile')
    print(json.dumps({'status':'FINITE-EXACT','gates':GATES,
                     'all_complement_and_inertia_cases_dimensions_1_to_3':cases,
                     'normseven_strata_including_open':len(strata),
                     'recovered_source_counts':{str(k):v for k,v in sorted(counts.items())},
                     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
