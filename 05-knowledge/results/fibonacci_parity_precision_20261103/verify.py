"""Exact low-bit and golden-SCC complement checks; run from any directory."""
from collections import Counter
from fractions import Fraction
from importlib.util import spec_from_file_location, module_from_spec
from itertools import combinations, permutations
from math import factorial, gcd, prod
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def load(name, path):
    spec = spec_from_file_location(name, ROOT / path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


PRECISION = '05-knowledge/results/port_residual_integral_20261025/precision/check.py'
CHANNELS = '05-knowledge/results/session89_consecutive_tournament_pairs_20261028/counting/channels.py'
SOURCE_TEXTS = [
    '05-knowledge/results/scc_basins_20261020d/arithmetic.md',
    '05-knowledge/results/session81_counterexample_families_20261028/CONFLICT_IMAGES.md',
    '05-knowledge/results/port_residual_integral_20261025/precision/README.md',
    '05-knowledge/results/metagraph_scheduled_recursion_20261029/PRECISION.md',
    '05-knowledge/results/overnight_rank_return_20261029/DEFECT_TRANSFER.md',
    '05-knowledge/results/merged_metagraph_purity_and_strong_core_20261020/STRONG_CORE_COUNT_RECURSIONS.md',
    '05-knowledge/results/session149_tournament_spectrum_20261103/toothpick/CORRECTION.md',
    '01-canon/theorems/THM-1420-no-linear-tournament-invariants-and-the-chirality-that-survives.md',
    '01-canon/theorems/THM-646-line-metagraph-score-complement-skeleton.md',
    '01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md',
]


def prime_factors(n):
    answer = []
    p = 2
    while p*p <= n:
        if n % p == 0:
            answer.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        answer.append(n)
    return answer


def phi(n):
    primes = prime_factors(n)
    return n // prod(primes) * prod(p-1 for p in primes)


def v2(n):
    need(n > 0, 'positive valuation argument')
    return (n & -n).bit_length()-1


def coprime_triples(n):
    return sum(gcd(a,b) == gcd(a,c) == gcd(b,c) == 1
               for a in range(1, n//3+1, 2)
               for b in range(a+2, (n-a)//2+1, 2)
               if b < (c := n-a-b))


def residue8(n):
    need(type(n) is int and n >= 0, 'nonnegative integer order')
    if n <= 4:
        return 1
    if n % 2:
        return (pow(n,-1,8)+2*pow(n-2,-1,8)+4*coprime_triples(n)) % 8
    if n % 4 == 0:
        main = (n-1)*phi(n)//2
    else:
        m = n//2
        main = m*prod(1-p for p in prime_factors(m))//2
    return (main + (2*phi(n//3) if n % 3 == 0 else 0)) % 8


def odd_partitions(n, minimum=1):
    if n == 0:
        yield ()
    for a in range(minimum, n+1, 2):
        for tail in odd_partitions(n-a, a):
            yield (a,) + tail


def davis(n):
    answer = Fraction(0)
    for parts in odd_partitions(n):
        mult = Counter(parts)
        exponent = (n-len(parts))//2 + sum(gcd(a,b) for a,b in combinations(parts,2))
        z = prod(d**m*factorial(m) for d,m in mult.items())
        answer += Fraction(1 << exponent, z)
    need(answer.denominator == 1, ('integral Davis count',n))
    return int(answer)


def chart(n, bits):
    rows = [0]*n
    j = 0
    for a in range(n):
        for b in range(a+1,n):
            forward = b == a+1 or bool(bits >> j & 1)
            if b != a+1:
                j += 1
            if forward:
                rows[a] |= 1 << b
            else:
                rows[b] |= 1 << a
    return tuple(rows)


def strong(rows):
    n = len(rows)
    for transpose in (False,True):
        seen, todo = {0}, [0]
        while todo:
            u = todo.pop()
            for v in range(n):
                edge = rows[v] >> u & 1 if transpose else rows[u] >> v & 1
                if edge and v not in seen:
                    seen.add(v)
                    todo.append(v)
        if len(seen) != n:
            return False
    return True


def complement(rows):
    n = len(rows)
    all_bits = (1 << n)-1
    target = [(all_bits ^ row) & ~(1 << a) for a,row in enumerate(rows)]
    for a in range(n-1):
        need(rows[a] >> (a+1) & 1, 'physical path required')
        target[a] |= 1 << (a+1)
        target[a+1] &= ~(1 << a)
    return tuple(target)


def reflected_converse(rows):
    n = len(rows)
    return tuple(sum(1 << b for b in range(n) if rows[n-1-b] >> (n-1-a) & 1)
                 for a in range(n))


def join(left, right):
    k,n = len(left),len(right)
    mask = ((1 << n)-1) << k
    return tuple(x | mask for x in left) + tuple(x << k for x in right)


def reverse_order_width(rows):
    n = len(rows)
    # In the reversed path order, a backward edge was forward in the old order.
    return max((sum(bool(rows[a] >> b & 1) for a in range(t) for b in range(t,n))
                for t in range(n+1)), default=0)


def path_count(rows):
    n = len(rows)
    table = [[0]*n for _ in range(1 << n)]
    for a in range(n):
        table[1 << a][a] = 1
    for mask in range(1,1 << n):
        for a,number in enumerate(table[mask]):
            if not number:
                continue
            available = rows[a] & ~mask
            while available:
                bit = available & -available
                available -= bit
                table[mask | bit][bit.bit_length()-1] += number
    return sum(table[-1])


def canonical(rows):
    n = len(rows)
    pairs = list(combinations(range(n),2))
    return min(sum((rows[p[a]] >> p[b] & 1) << j for j,(a,b) in enumerate(pairs))
               for p in permutations(range(n)))


def check_arithmetic():
    fast = load('golden_channels',CHANNELS)
    precision = load('golden_precision',PRECISION)
    counts, stats = fast.scalar_prefix(100)
    published = [1,1,1,2,4,12,56,456,6880,191536,9733056,903753248,
                 154108311168,48542114686912,28401423719122304,
                 31021002160355166848,63530415842308265100288,
                 244912778438520759443245824,1783398846284777975419600287232,
                 24605641171260376770598003978281472,
                 645022068557873570931850526424042500096]
    need(counts[:21] == published, 'published A000568 prefix')
    hostile = []
    valuation_checks = 0
    for n,value in enumerate(counts):
        B = n//2 if n >= 3 else 0
        need(value % (1 << B) == 0, ('baseline',n))
        need((value >> B) % 8 == residue8(n), ('third-bit formula',n))
        if n <= 30:
            need(davis(n) == value, ('independent rational count',n))
        if n >= 4 and n % 2 == 0:
            guess = n//2+v2(phi(n))-1
            if v2(phi(n)) <= 3:
                need(v2(value) == guess, ('proved valuation lane',n))
                valuation_checks += 1
            elif v2(value) != guess:
                hostile.append([n,v2(value),guess])
    checked = list(range(101)) + [127,128,139,149,192,199,200,257,300]
    for n in checked:
        found,_ = precision.normalized_residue(n,3)
        need(found == residue8(n), ('bounded-cycle independent route',n))
    need(hostile[0] == [40,26,23], 'first valuation extrapolation failure')
    # Resolve the first hostile by cycle-length, retaining precision 7.
    contributions = Counter()
    for length in precision.cycle_counts(40,7):
        for parts in precision.fixed_length_partitions(40,length):
            exponent,z,twos = precision.term(parts)
            relative = exponent-twos-20
            if relative < 7:
                contributions[length] += (1 << relative)*pow(z >> twos,-1,128)
    contributions = {str(k):v % 128 for k,v in sorted(contributions.items())}
    need(sum(contributions.values()) % 128 == (counts[40] >> 20) % 128, 'hostile layers')
    return {'full_integer_orders':[0,100], 'independent_Fraction_orders':[0,30],
            'bounded_cycle_orders':checked, 'valuation_lane_cases':valuation_checks,
            'counter_stats':{k:v for k,v in stats.items() if k != 'rows'},
            'valuation_extrapolation_failures':hostile,
            'n40_normalized_mod128_by_cycle_count':contributions,
            'even_normalized_mod8_orders4_through40':[residue8(n) for n in range(4,41,2)]}


def check_diagrams():
    atom3 = [chart(3,b) for b in range(2) if strong(chart(3,b))]
    atom4 = [chart(4,b) for b in range(8) if strong(chart(4,b))]
    need(len(atom3) == 1 and len(atom4) == 5, 'seven-symbol alphabet')
    atoms = [(0,)] + atom3 + atom4
    atom_widths = [reverse_order_width(complement(a)) for a in atoms]
    banks = [[((),0,0)]] + [[] for _ in range(10)]
    rows = []
    maximum_width = 0
    witness = None
    for n in range(1,11):
        for atom in atoms:
            if len(atom) > n:
                continue
            for prefix,a,b in banks[n-len(atom)]:
                banks[n].append((join(prefix,atom),a+(len(atom)==3),b+(len(atom)==4)))
        unique = {T for T,_,_ in banks[n]}
        need(len(unique) == len(banks[n]), ('SCC word injectivity',n))
        blue = targets_blue = 0
        spectrum = Counter()
        for T,a,b in banks[n]:
            K = complement(T)
            need(complement(K) == T, ('literal involution',n))
            need((reflected_converse(T) == T) == (reflected_converse(K) == K),
                 ('colour transport',n))
            blue += reflected_converse(T) == T
            targets_blue += reflected_converse(K) == K
            width = reverse_order_width(K)
            maximum_width = max(maximum_width,width)
            if width >= 4 and witness is None:
                witness = {'n':n,'source_rows':T,'target_rows':K,'width':width}
            need(width <= 4, ('bounded query width',n,T))
            if n >= 5:
                need(strong(K), ('strong actual target',n))
                need(K not in unique, ('forest locus exit',n))
            if n <= 8:
                need(path_count(T) == 3**a*5**b, ('source spectrum',n))
                spectrum[path_count(K)] += 1
        rows.append({'n':n,'marked_sources':len(banks[n]),'blue_sources':blue,
                     'target_H_histogram' : dict(sorted(spectrum.items()))})
        need(blue == targets_blue, 'blue transport mass')
    marked = [1] + [r['marked_sources'] for r in rows]
    for n in range(1,11):
        need(marked[n] == sum(c*marked[n-d] for d,c in [(1,1),(3,1),(4,5)] if d <= n),
             ('marked recurrence',n))
    # B(x)=(1+x+x^3+3x^4)/(1-x^2-x^6-5x^8).
    blues = [1] + [r['blue_sources'] for r in rows]
    for n in range(1,11):
        need(blues[n] == {1:1,3:1,4:3}.get(n,0)+
             sum(c*blues[n-d] for d,c in [(2,1),(6,1),(8,5)] if d <= n), ('blue recurrence',n))
    four = [chart(4,b) for b in range(8)]
    keys = [canonical(T) for T in four]
    names = {}
    for key,T in zip(keys,four):
        degree = sorted(x.bit_count() for x in T)
        names[key] = { (0,1,2,3):'TT4', (0,2,2,2):'C3_to_1',
                       (1,1,1,3):'1_to_C3',(1,1,2,2):'Q4'}[tuple(degree)]
    fibres = Counter(names[k] for k in keys)
    edges = Counter(tuple(sorted([names[keys[b]],names[keys[b^7]]])) for b in range(4))
    need(fibres == Counter({'TT4':1,'C3_to_1':1,'1_to_C3':1,'Q4':5}), 'n4 fibres')
    need(edges == Counter({('Q4','TT4'):1,('C3_to_1','Q4'):1,('1_to_C3','Q4'):1,('Q4','Q4'):1}),
         ('star plus loop',edges))
    # Incoming odd-fibre cancellation: local re-marking pairs all but one
    # chart in each fibre. Alternate those edges with actual complements.
    rho, representatives = {}, set()
    for key in set(keys):
        fibre = [b for b,k in enumerate(keys) if k == key]
        need(len(fibre) % 2 == 1, 'odd class fibre')
        representatives.add(fibre[0])
        for a,b in zip(fibre[1::2],fibre[2::2]):
            rho[a],rho[b] = b,a
    trails, used = [],set()
    for start in sorted(representatives):
        if start in used:
            continue
        trail = [start]
        current = start
        while True:
            current ^= 7
            need(current not in trail, 'trail cannot close at a degree-one endpoint')
            trail.append(current)
            if current in representatives:
                break
            current = rho[current]
            need(current not in trail, 'no internal repeated vertex')
            trail.append(current)
        need(keys[start] != keys[current], 'distinct endpoint classes')
        used.update([start,current])
        trails.append(trail)
    need(used == representatives and trails == [[0,7],[3,4,5,2,1,6]], 'explicit class pairing trails')
    return {'alphabet_sizes':[1,1,5], 'atom_widths':atom_widths, 'word_universe':rows,
            'maximum_retained_order_width':maximum_width,'width4_witness':witness,
            'n4_fibres':dict(fibres),'n4_chart_classes':[names[k] for k in keys],
            'n4_blue_charts':[b for b,T in enumerate(four) if reflected_converse(T) == T],
            'n4_class_pairing_trails':trails,
            'n4_lines':[{'ends':ends,'multiplicity':m} for ends,m in sorted(edges.items())]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args = parser.parse_args()
    data = {'arithmetic':check_arithmetic(),'diagrams':check_diagrams()}
    data['source_hashes_lf'] = {p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').encode('utf-8')).hexdigest()
                                for p in [PRECISION,CHANNELS]+SOURCE_TEXTS}
    encoded = json.dumps(data,indent=2,sort_keys=True)+'\n'
    if args.write:
        (HERE/'results.json').write_text(encoded,encoding='utf-8',newline='\n')
    else:
        need(json.loads((HERE/'results.json').read_text(encoding='utf-8')) == json.loads(encoded), 'saved exact receipt')
    print(json.dumps({'status':'PASS','checked_full_counts':101,
                      'marked_family_sources':sum(r['marked_sources'] for r in data['diagrams']['word_universe']),
                      'maximum_width':data['diagrams']['maximum_retained_order_width'],
                      'n40_layers':data['arithmetic']['n40_normalized_mod128_by_cycle_count']}))


if __name__ == '__main__':
    main()
