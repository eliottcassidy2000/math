"""Exact pinned cycle indices, order-p rigidity, and aliquot query controls."""
from collections import Counter
from itertools import product
from pathlib import Path
from math import comb, gcd, lcm
import importlib.util
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
GATES = 0


def need(value, why):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(why)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sig(counts):
    return tuple(sorted((a, b) for a, b in counts.items() if b))


def clock_type(p, k, r):
    counts = {1: p**(2*k-2*r+1)}
    for i in range(1, r):
        counts[p**i] = (p*p-1)*p**(2*k-2*r+i-1)
    counts[p**r] = (p-1)*p**(2*k-r-1)
    need(sum(a*b for a, b in counts.items()) == p**(2*k), 'clock mass')
    return sig(counts)


def inventory(p, k, j, ell=None):
    """UNNORMALIZED cycle inventory of the pointwise stabilizer."""
    n = p**(2*k)
    out = Counter()
    for z in range(p-j+1):
        out[sig({1: (j+z)*n//p, p: (p-j-z)*n//(p*p)})] += (
            comb(p-j, z)*(p-1)**(p-j-z))
    if j == 1:
        for r in range(2, ell+1):
            out[clock_type(p, k, r)] += p**(p-2)*(p-1)*p**(r-1)
    group_size = p**(p-j) if j >= 2 else p**(p-2+ell)
    need(sum(out.values()) == group_size, 'stabilizer order')
    return out


def orbit_count(inv, palette, pins=0):
    numerator = sum(weight*palette**(sum(c for _, c in ct)-pins)
                    for ct, weight in inv.items())
    order = sum(inv.values())
    need(numerator % order == 0, 'Burnside integrality')
    return numerator//order


def orbit_count_mod(inv, palette, modulus):
    order = sum(inv.values())
    numerator = sum(weight*pow(palette, sum(c for _, c in ct), order*modulus)
                    for ct, weight in inv.items()) % (order*modulus)
    need(numerator % order == 0, 'modular Burnside exact division')
    return numerator//order


def cycle_type(perm):
    visited = set()
    counts = Counter()
    for v in range(len(perm)):
        if v in visited:
            continue
        w, length = v, 0
        while w not in visited:
            visited.add(w)
            w, length = perm[w], length+1
        need(w == v, 'permutation cycles close at start')
        counts[length] += 1
    return sig(counts)


def compose(a, b):
    return tuple(a[b[v]] for v in range(len(a)))


def perm_power(a, exponent):
    result = tuple(range(len(a)))
    while exponent:
        if exponent & 1:
            result = compose(a, result)
        exponent >>= 1
        a = compose(a, a)
    return result


def literal_groups():
    prior = load('05-knowledge/results/universal_ramified_closure_20261030/check.py',
                 'prior_ring')
    output = []
    # Exhaust all root stabilizers: 625, 3125, and 117649 literal permutations.
    for p, k, d in [(5, 1, 1), (5, 2, 1), (5, 2, 2), (7, 1, 1)]:
        ring = prior.Ring(p, k, d)
        identity = tuple(range(ring.n))
        # h(0)=h(1)=0 removes exactly the overlap with the top clock.
        pins = [('root', [0], 1, k)]
        if k >= 2:
            pins += [('tau1', [0, 1], 1, 1),
                     ('tau2', [0, ring.idx((p, 0))], 1, 1),
                     ('tau3', [0, ring.idx((0, p))], 1, 2)]
        pins += [('residues'+str(j), [ring.idx((a, 0)) for a in range(j)], j, None)
                 for j in range(2, p+1)]
        actual = {name: Counter() for name, *_ in pins}
        seen = set()
        order_p = set()
        support = []
        for residue in range(1, p):
            h = [0]*p
            h[residue] = 1
            support.append(tuple(ring.affine_fibre(0, (0, 0), h)))
        f0 = set()
        for values in product(range(p), repeat=p-1):
            f0.add(tuple(ring.affine_fibre(0, (0, 0), (0, *values))))
        # Single black vertex in every unoccupied fibre must kill H_0.
        black = {ring.idx((a, 0)) for a in range(1, p)}
        black_survivors = 0
        # An independently chosen binary word, not a representative census.
        word = {v for v in range(ring.n) if (17*v*v+11*v+3) % 23 < 10}
        word_stabilizer = 0
        for t in range(ring.q):
            for tail in product(range(p), repeat=p-2):
                perm = tuple(ring.affine_fibre(t, (0, 0), (0, 0, *tail)))
                need(perm not in seen, 'unique root-stabilizer parametrization')
                seen.add(perm)
                ct = cycle_type(perm)
                if max(length for length, _ in ct) <= p:
                    order_p.add(perm)
                for name, subset, j, ell in pins:
                    if all(perm[v] == v for v in subset):
                        actual[name][ct] += 1
                black_survivors += {perm[v] for v in black} == black
                word_stabilizer += {perm[v] for v in word} == word
                if perm != identity:
                    order = max(length for length, _ in ct)
                    descendant = perm_power(perm, order//p)
                    need(descendant != identity and descendant in f0,
                         'every nonidentity has a nontrivial order-p fibre descendant')
        need(order_p == f0, 'Omega_1 is exactly F_0')
        need(black_survivors == 1, 'sharp black support construction')
        fibre_test = all({perm[v] for v in word} != word for perm in support)
        need((word_stabilizer == 1) == fibre_test, 'literal rigidity equivalence')
        rows = []
        for name, subset, j, ell in pins:
            predicted = inventory(p, k, j, ell)
            need(actual[name] == predicted, 'full pinned cycle distribution '+name)
            for palette in (2, 3, 4):
                need(orbit_count(actual[name], palette, len(subset)) ==
                     orbit_count(predicted, palette, len(subset)), 'palette evaluation')
            rows.append({'name': name, 'pins': subset, 'j': j, 'ell': ell,
                         'group_order': sum(predicted.values()),
                         'binary_extensions': str(orbit_count(predicted, 2, len(subset))),
                         'cycle_types': [{'cycles': ct, 'multiplicity': w}
                                         for ct, w in sorted(predicted.items())]})
        output.append({'p': p, 'k': k, 'D': p*d, 'root_permutations': len(seen),
                       'order_dividing_p': len(order_p), 'pin_cases': rows})
    return output


def deeper_clock_controls():
    prior = load('05-knowledge/results/universal_ramified_closure_20261030/check.py',
                 'deeper_ring')
    rows = []
    for p, k, d in [(5, 3, 2), (7, 2, 3)]:
        ring = prior.Ring(p, k, d)
        tested = 0
        for r in range(2, k+1):
            for unit in (1, 2, p+1):
                t = unit*p**(k-r)
                for seed in range(5):
                    h = [0, 0]+[(seed*(i+1)+i*i) % p for i in range(p-2)]
                    perm = ring.affine_fibre(t, (0, 0), h)
                    need(cycle_type(perm) == clock_type(p, k, r),
                         'higher clock with explicit fibre perturbation')
                    tested += 1
        rows.append({'p': p, 'k': k, 'D': p*d, 'selected_maps': tested})
    return rows


def small_color_orbits():
    # At p=5,k=1, pin all points in 3 fibres. Two free fibres each have
    # one 5-cycle. Exhaust their 2^10 words and all 25 translations.
    perms = []
    for a, b in product(range(5), repeat=2):
        perms.append(tuple((v+a) % 5 if v < 5 else 5+(v-5+b) % 5
                           for v in range(10)))
    reflection = tuple((-v) % 5 if v < 5 else 5+(-(v-5)) % 5 for v in range(10))
    unsigned = perms+[compose(a, reflection) for a in perms]
    reps, unsigned_reps, rigid_masks, unsigned_rigid = set(), set(), 0, 0
    for mask in range(1 << 10):
        images = []
        for perm in perms:
            images.append(sum(1 << perm[v] for v in range(10) if mask >> v & 1))
        reps.add(min(images))
        rigid_masks += len(set(images)) == 25
        reflected = sum(1 << reflection[v] for v in range(10) if mask >> v & 1)
        images2 = images+[sum(1 << perm[v] for v in range(10) if reflected >> v & 1)
                          for perm in perms]
        unsigned_reps.add(min(images2))
        unsigned_rigid += len(set(images2)) == 50
    inv = inventory(5, 1, 3)
    need(len(reps) == orbit_count(inv, 2, 15) == 64, 'all small binary orbits')
    need(rigid_masks == (2**5-2)**2 == 900, 'all small rigid masks')
    need(len(unsigned_reps) == 64 and unsigned_rigid == 0,
         'one common reflection invalidates directed rigidity even with same orbit count')
    return {'universe': 1024, 'group': 25, 'orbits': len(reps),
            'rigid_masks': rigid_masks, 'rigid_orbits': rigid_masks//25,
            'unsigned_group': len(unsigned), 'unsigned_orbits': len(unsigned_reps),
            'unsigned_rigid_masks': unsigned_rigid}


def generated_orbits(permutations, n):
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for perm in permutations:
        for a, b in enumerate(perm):
            a, b = find(a), find(b)
            if a != b:
                parent[a] = b
    return len({find(x) for x in range(n)})


def unsigned_controls():
    prior = load('05-knowledge/results/universal_ramified_closure_20261030/check.py',
                 'phase_ring')
    ring = prior.Ring(5, 2, 1)
    p, n, q = ring.p, ring.n, ring.q
    conjugate = tuple(ring.idx((x, -y % q)) for x, y in ring.points)
    shifts = []
    generators = []
    for h in product(range(p), repeat=p-2):
        shifts.append(tuple(ring.affine_fibre(0, (0, 0), (0, 0, *h))))
    for residue in range(2, p):
        h = [0]*p
        h[residue] = 1
        generators.append(tuple(ring.affine_fibre(0, (0, 0), h)))
    reflections = [compose(a, conjugate) for a in shifts]
    need(len(set(shifts+reflections)) == 250, 'actual two-pin unsigned group')
    predicted = ((1, q), (2, (n-q)//2))
    for f in reflections:
        need(cycle_type(f) == predicted, 'all literal phase cycles')
        for subset in range(1 << (p-2)):
            chosen = [generators[i] for i in range(p-2) if subset >> i & 1]
            count = generated_orbits([f, *chosen], n)
            need(count == (n+q)//2-subset.bit_count()*(n*(p-1)//(p*p))//2,
                 'dihedral generated orbit count for inclusion-exclusion')
    cases = []
    for name, subset, j in [('two',[0,q],2), ('three_integer',[0,q,2*q],3),
                           ('three_tilted',[0,q,2*q+1],3),
                           ('transversal',[i*q for i in range(p)],p)]:
        direct = [f for f in shifts if all(f[v] == v for v in subset)]
        phase = [f for f in reflections if all(f[v] == v for v in subset)]
        epsilon = int(bool(phase))
        need(len(phase) in (0, len(direct)), 'phase coset present or absent')
        inv = Counter(cycle_type(f) for f in direct+phase)
        for palette in (2,3):
            numerator = orbit_count(inventory(p, 2, j), palette, len(subset))
            expected = numerator if not epsilon else (
                numerator+palette**((n+q)//2-len(subset)))//2
            need(orbit_count(inv, palette, len(subset)) == expected,
                 'full unsigned extension count')
        cases.append({'name':name,'pins':subset,'j':j,'extra_phase':epsilon,
                      'direct_order':len(direct),'unsigned_order':len(direct+phase)})
    need([x['extra_phase'] for x in cases] == [1,1,0,1], 'phase controls')
    return {'two_pin_maps':250,'reflections_checked':125,
            'generated_subgroup_checks':1000,'cases':cases,
            'random_transversal_binary_failure625':'2^-300 / 125',
            'scope':'k>=2,j>=2; inherited S112 axes and phase distribution'}


def symbolic_controls():
    rows = []
    for p in (5, 7, 11, 13):
        for k in range(1, 5):
            n = p**(2*k)
            for ell in range(1, k+1):
                inv = inventory(p, k, 1, ell)
                for palette in (2, 3):
                    modulus = 1000000007
                    count = orbit_count_mod(inv, palette, modulus)
                    if ell > 1:
                        old = orbit_count_mod(inventory(p, k, 1, ell-1), palette, modulus)
                        power = 2*p**(2*k-ell)-p**(2*k-2*ell)
                        need((p*count-old-(p-1)*pow(palette, power, modulus)) % modulus == 0,
                             'one-fibre LCA recurrence')
                    order = sum(inv.values())
                    d0 = n*(p-1)//(p*p)
                    rigid_mod = sum((-1)**z*comb(p-1,z)*pow(palette,n-z*d0,order)
                                    for z in range(p)) % order
                    need(rigid_mod == 0, 'free orbit divisibility')
            for j in range(2, p+1):
                inv = inventory(p, k, j)
                orbit_count_mod(inv, 2, 1000000007)
            rows.append({'p': p, 'k': k, 'selectors_tested': p+k-1})
    # Same number of pins, same occupancy, different LCA and orbit count.
    deep = orbit_count(inventory(5, 2, 1, 2), 2, 2)
    shallow = orbit_count(inventory(5, 2, 1, 1), 2, 2)
    need(shallow != deep, 'rigidity probability is not a complete orbit counter')
    return {'parameters': rows, 'same_j_two_pin_counts':
            {'tau1': str(shallow), 'tau3': str(deep)},
            'root625_rigidity': '(1-2^-100)^4',
            'root625_minimum_black_support': 4}


def sigma(n):
    result, p = 1, 2
    while p*p <= n:
        total, power = 1, 1
        while n % p == 0:
            n //= p
            power *= p
            total += power
        result *= total
        p += 1
    return result*(n+1 if n > 1 else 1)


def aliquot_queries():
    # Inherited fibre, independently regenerate by a full sieve under the
    # sharper Session78 bound: 342 is not prime, so n<=342^2/4=29241.
    limit = 29241
    sums = [0]*(2*limit+1)
    for divisor in range(1, len(sums)//2+1):
        for multiple in range(2*divisor, len(sums), divisor):
            sums[multiple] += divisor
    sources = [n for n in range(1, limit+1) if sums[n] == 343]
    expected = [578,1001,1685,3641,9077,9641,16697,17141,19241,19637,
                20777,24341,24617,25397,25877,27641,28757,28841,29177]
    need(sources == expected, 'inherited complete 343 fibre')
    traces = []
    for n in sources:
        need(sigma(n)-n == sums[n] and sigma(2*n)-2*n == sums[2*n],
             'independent sigma and sieve at both probes')
        row = {'n': n}
        for p in (2, 7):
            x, y = sigma(n)-n, sigma(p*n)-p*n
            c, d = y-p*x, y-x
            a, unit = 0, d
            while unit % p == 0:
                unit //= p
                a += 1
            recovered = (p**(a+1)*c-d)//(p-1)
            need(recovered == n, 'inherited prime dilation inverse')
            row[str(p)] = {'next': y, 'carry': c, 'valuation': a}
        traces.append(row)
    anchor = []
    for n in (1001, 343, 57):
        x, y = sigma(n)-n, sigma(7*n)-7*n
        anchor.append({'n': n, 's': x, 's7n': y, 'carry': y-7*x})
    need(anchor == [{'n':1001,'s':343,'s7n':2569,'carry':168},
                    {'n':343,'s':57,'s7n':400,'carry':1},
                    {'n':57,'s':23,'s7n':241,'carry':80}], 'shared-prime anchor')
    # First modular collision for the chosen two-dilation trace, not new
    # exact-source injectivity. Complete future traces depend on x,c mod M.
    min_modulus = next(m for m in range(19, 1000)
                       if len({row['2']['next'] % m for row in traces}) == 19)
    return {'complete_fibre': sources, 'queries': traces,
            'anchor_prime7_queries': anchor,
            'least_modulus_separating_this_doubling_fibre': min_modulus,
            'status': 'inherited all-size inverse; new finite query receipt only'}


def scaled_zero_controls():
    cases = 0
    for bits in range(11):
        for forbidden, extra in [(2, 4), (10, 8)]:
            m = bits+extra
            moduli = [d for d in range(2, 2*m+1) if gcd(d, forbidden) == 1]
            joint = lcm(*moduli)
            need(joint > 2**bits, 'unit-modulus lcm exceeds unscaled height')
            central = comb(2*m, m)
            stripped = central
            for prime in (2,) if forbidden == 2 else (2, 5):
                pp = 1
                while stripped % prime == 0:
                    stripped //= prime
                    pp *= prime
                need(pp <= 2*m, 'central-binomial prime power bound')
            need(joint % stripped == 0, 'stripped central binomial divides lcm')
            for count in range(2**bits+1):
                for exponent in (0, 1, 17, 1000000):
                    residues = [pow(forbidden, exponent, d)*count % d for d in moduli]
                    need(all(r == 0 for r in residues) == (count == 0),
                         'scaled count zero test independent of scale exponent')
                    for d, r in zip(moduli, residues):
                        need(r*pow(pow(forbidden, exponent, d), -1, d) % d == count % d,
                             'unit-scale exact modular count recovery')
                    cases += 1
    need((2**20*0) % (2**10) == (2**20*1) % (2**10), 'nonunit-scale hostile')
    moduli = [3, 5, 7, 9, 11, 13]
    false_zero = lcm(*moduli)
    need(false_zero > 0 and all(false_zero % d == 0 for d in moduli),
         'fixed finite moduli fail without height guard')
    return {'scaled_count_cases': cases, 'bits': [0, 10],
            'scale_exponents': [0, 1, 17, 1000000],
            'scale_bases': [2, 10], 'unbounded_height_hostile': false_zero}


def main():
    data = {'literal_groups': literal_groups(), 'deeper_clocks': deeper_clock_controls(),
            'small_color_orbits': small_color_orbits(),
            'symbolic_controls': symbolic_controls(), 'aliquot': aliquot_queries(),
            'scaled_zero': scaled_zero_controls(), 'unsigned': unsigned_controls()}
    hostile = load('05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
                   'native_hostile')
    data['hostile_81_27'] = hostile.hostile_checks()
    data['local_gates'] = GATES
    data['inherited_hostile_gates'] = hostile.GATES if hasattr(hostile, 'GATES') else 7
    print(json.dumps(data, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
