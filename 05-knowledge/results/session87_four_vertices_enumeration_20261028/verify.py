#!/usr/bin/env python3
"""Independent exact controls, no timings or known terms used by production."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations, product
from math import factorial, gcd, isqrt
from pathlib import Path
import argparse
import json
import sys

from baseline_sources import configure, source_path, SOURCES

from fast_counts import (count_direct, count_prefix, degree_bucket_dp,
                         packed_integer_dp, rooted_prefix, odd_partition_numbers,
                         palette_prefixes)

HERE = Path(__file__).resolve().parent
CHECKS = Counter()
sys.set_int_max_str_digits(0)


def need(ok, label):
    CHECKS[label] += 1
    if not ok:
        raise RuntimeError(label)


def load(name, filename):
    spec = spec_from_file_location(name, source_path(filename))
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parts(n, low=1, odd=True):
    if n == 0:
        yield ()
    for k in range(low,n+1,2 if odd else 1):
        for tail in parts(n-k,k,odd):
            yield (k,)+tail


def exponent(p):
    return sum((k-1)//2 for k in p)+sum(gcd(a,b) for a,b in combinations(p,2))


def centralizer(p):
    z = 1
    for k,m in Counter(p).items():
        z *= k**m*factorial(m)
    return z


def fraction_count(n, colors=1, graph=False):
    total = Fraction(0)
    for p in parts(n,odd=not graph):
        e = (sum(k//2 for k in p)+sum(gcd(a,b) for a,b in combinations(p,2))
             if graph else exponent(p))
        total += Fraction((1 << e)*colors**len(p),centralizer(p))
    need(total.denominator == 1,'independent Fraction integrality')
    return int(total)


def literal_permutation_burnside(n, colors=1):
    """Enumerate all labelled tournament bits and all vertex permutations."""
    pairs = list(combinations(range(n),2))
    index = {e:i for i,e in enumerate(pairs)}
    total = 0
    for p in permutations(range(n)):
        action = [(index[tuple(sorted((p[u],p[v])))],p[u]>p[v]) for u,v in pairs]
        fixed = 0
        for mask in range(1 << len(pairs)):
            image = sum((((mask >> i)&1)^flip) << j for i,(j,flip) in enumerate(action))
            fixed += image == mask
        # Literal color assignments, not the cycle-count formula.
        color_fixed = sum(all(c[i] == c[p[i]] for i in range(n))
                          for c in product(range(colors),repeat=n))
        total += fixed*color_fixed
    value,rem = divmod(total,factorial(n))
    need(not rem,'literal permutation Burnside integrality')
    return value


def literal_graph_classes(n, royle=False, degree_even=False):
    """All simple graph masks, literal relabeling; tiny control only."""
    edges = list(combinations(range(n),2))
    index = {e:i for i,e in enumerate(edges)}
    actions = []
    for p in permutations(range(n)):
        mapping = [index[tuple(sorted((p[u],p[v])))] for u,v in edges]
        inversions = sum((p[u]>p[v]) << i for i,(u,v) in enumerate(edges))
        actions.append((mapping,inversions))
    representatives = set()
    for mask in range(1 << len(edges)):
        if degree_even and any(sum((mask>>i)&1 for i,e in enumerate(edges) if v in e)%2
                               for v in range(n)):
            continue
        images = []
        even = True
        for mapping,inversions in actions:
            image = sum(((mask>>i)&1) << j for i,j in enumerate(mapping))
            images.append(image)
            if image == mask and (mask & inversions).bit_count()%2:
                even = False
        if not royle or even:
            representatives.add(min(images))
    return len(representatives)


def isprime(p):
    if p < 2:
        return False
    return all(p%d for d in range(2,isqrt(p)+1))


def certified_crt(n, colors=1):
    """Small independent modular check, not a claimed faster method.

    Exact bound C_n(c)<=c^n 2^C(n,2). Product of primes strictly exceeding
    that bound gives unique reconstruction, not heuristic bit estimates.
    """
    bound = colors**n << (n*(n-1)//2)
    x,modulus = 0,1
    candidate = max(n+1,1000003)
    prime_count = 0
    data = [(exponent(p),centralizer(p),len(p)) for p in parts(n)]
    while modulus <= bound:
        while not isprime(candidate):
            candidate += 1
        prime = candidate
        residue = sum(pow(2,e,prime)*pow(colors,t,prime)*pow(z,-1,prime)
                      for e,z,t in data)%prime
        step = ((residue-x)%prime)*pow(modulus,-1,prime)%prime
        x += modulus*step
        modulus *= prime
        candidate += 1
        prime_count += 1
    need(0 <= x <= bound,'CRT reconstruction within exact labelled upper bound')
    return x,{'prime_count':prime_count,'modulus_bits':modulus.bit_length(),
              'bound_bits':bound.bit_length(),'unique':modulus>bound}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--large',type=int,default=200)
    ap.add_argument('--baseline-root',type=Path)
    ap.add_argument('--receipt',type=Path)
    args = ap.parse_args()
    configure(args.baseline_root)
    old = load('s22','local_scale.py')
    channels = load('channels','channels.py')
    s86 = load('s86','recreate_a000568.py')
    for N in list(range(31))+[50,80,100,args.large]:
        baseline,stats = old.local_dp(N)
        for fold in (False,True):
            values,newstats = degree_bucket_dp(N,terminal_fold=fold)
            need(values == baseline,'degree-bucket exact full prefix')
            need(newstats['rows'] == stats['rows'],'unchanged state and transition rows')
            need(newstats['precomputed_transfers'] == stats['precomputed_transfers'],
                 'unchanged degree-transfer count')
        if N <= 100:
            need(packed_integer_dp(N)[0] == baseline,'full-key packing exploratory control')
            need(count_prefix(N)[0] == baseline,'integer 1-and-3-tail kernel')
            need(count_direct(N)[0] == baseline[N],'incremental single-order bucket sum')
    values,_ = degree_bucket_dp(max(50,args.large))
    for n in range(51):
        need(values[n] == fraction_count(n),'independent odd-partition Fraction formula')
        need(values[n] == s86.OEIS_A000568_0_50[n],'frozen OEIS reference audit only')
    literal_types = 0
    for n in range(15):
        for p in parts(n,odd=False):
            expected = 0 if any(k%2==0 for k in p) else 1 << exponent(p)
            need(s86.direct_fixed_orientations(p) == expected,'literal ordered-pair orbits including even-cycle hostiles')
            literal_types += 1
    for n in range(6):
        need(literal_permutation_burnside(n) == values[n],
             'literal labelled tournaments and all permutations')
    palettes = {}
    for c in (1,2,3,4,8):
        sequence,stats = degree_bucket_dp(35,colors=c)
        for n in range(21):
            need(sequence[n] == fraction_count(n,colors=c),'colored partition formula')
        if c & (c-1) == 0:
            inherited,st = channels.channel_dp(35,a=1,b=c.bit_length()-1)
            need(sequence == inherited,'inherited power-of-two channel')
            need(stats['rows'] == st['rows'],'color-independent DP state rows')
        if c <= 3:
            for n in range(5):
                need(sequence[n] == literal_permutation_burnside(n,c),
                     'literal vertex-color assignment Burnside')
        palettes[str(c)] = sequence[:11]
    pointed = rooted_prefix(20)
    colored,_ = degree_bucket_dp(19,colors=2)
    need(pointed[1:] == colored,'rooted equals two-color attachment channel')
    shared,_ = palette_prefixes(100,(1,2,3,4,8,16))
    for c,sequence in shared.items():
        need(sequence == degree_bucket_dp(100,c)[0],'reused integer kernel color evaluation')
    for k in range(4):
        colored,_ = degree_bucket_dp(18,colors=2**k)
        for n in range(k,19):
            direct = Fraction(0)
            for p in parts(n):
                ones = p.count(1)
                falling = (factorial(ones)//factorial(ones-k)) if ones >= k else 0
                direct += Fraction(falling*(1 << exponent(p)),centralizer(p))
            need(direct == (colored[n-k] << (k*(k-1)//2)),
                 'ordered distinct roots equals attachment-color channel')
    graph_counts,royle_counts,degree_even_counts = [],[],[]
    for n in range(6):
        graph_counts.append(literal_graph_classes(n))
        royle_counts.append(literal_graph_classes(n,royle=True))
        degree_even_counts.append(literal_graph_classes(n,degree_even=True))
        need(graph_counts[-1] == fraction_count(n,graph=True),'literal ordinary graph classes')
        need(royle_counts[-1] == values[n],'literal Royle-even classes equal tournaments')
    need(degree_even_counts[4] == 3 and royle_counts[4] == 4,
         'degree-even is NOT Royle-even hostile n4')
    crt = {}
    for n,c in ((0,1),(20,1),(30,1),(20,3)):
        reconstructed,certificate = certified_crt(n,c)
        expected = degree_bucket_dp(n,c)[0][n]
        need(reconstructed == expected,'independent certified CRT reconstruction')
        crt[f'{n},{c}'] = certificate
    for bad in (-1,True,1.5,'5'):
        try:
            degree_bucket_dp(bad)
        except (ValueError,TypeError):
            need(True,'reject invalid bound')
        else:
            need(False,'reject invalid bound')
    for bad in (0,-1,True,1.5):
        try:
            degree_bucket_dp(4,bad)
        except (ValueError,TypeError):
            need(True,'reject invalid palette')
        else:
            need(False,'reject invalid palette')
    receipt = {'status':'PASS','checks':dict(CHECKS),'large_prefix_through':args.large,
        'literal_cycle_types_through14':literal_types,'palettes_through10':palettes,
        'pointed_through10':pointed[:11],'ordinary_graphs_through5':graph_counts,
        'Royle_even_graphs_through5':royle_counts,'degree_even_graphs_through5':degree_even_counts,
        'CRT_certificates':crt,'count_at_large_sha256':sha256(str(values[args.large]).encode()).hexdigest(),
        'source_sha256':dict({p.name:sha256(p.read_bytes()).hexdigest()
                         for p in [HERE/'fast_counts.py',HERE/'verify.py',HERE/'baseline_sources.py']},
                         **{SOURCES[name][0]:sha256(source_path(name).read_bytes()).hexdigest()
                            for name in ('local_scale.py','channels.py','recreate_a000568.py')})}
    output = json.dumps(receipt,indent=2)+'\n'
    if args.receipt:
        args.receipt.write_text(output)
    print(output)


if __name__ == '__main__':
    main()
