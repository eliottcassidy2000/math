"""Exact controls for arithmetic codes, affine holonomy and bounded lifts.

Run from the repository: python -B <this path> [--output receipt.json].
No third-party dependencies and no removable assertions.
"""
import argparse
import itertools as it
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def prime(n):
    return n >= 2 and all(n % d for d in range(2, math.isqrt(n) + 1))


def xor(values):
    z = 0
    for v in values:
        z ^= v
    return z


def rank2(rows):
    pivots = {}
    for x in rows:
        while x:
            p = x.bit_length()
            if p in pivots:
                x ^= pivots[p]
            else:
                pivots[p] = x
                break
    return len(pivots)


def rankq(matrix):
    a = [[Fraction(v) for v in row] for row in matrix]
    if not a:
        return 0
    j = 0
    for c in range(len(a[0])):
        p = next((p for p in range(j, len(a)) if a[p][c]), None)
        if p is None:
            continue
        a[j], a[p] = a[p], a[j]
        v = a[j][c]
        a[j] = [x / v for x in a[j]]
        for p in range(len(a)):
            if p != j:
                v = a[p][c]
                a[p] = [x - v*y for x, y in zip(a[p], a[j])]
        j += 1
        if j == len(a):
            break
    return j


def codes():
    out = {}
    for r, basis in [(3, [2, 3, 5]), (4, [2, 3, 5, 7])]:
        n = 2**r - 1
        bank = []
        for v in range(1, n + 1):
            d = math.prod(basis[j] for j in range(r) if v >> j & 1)
            t = 1
            while not prime(d*t*t - 1):
                t += 1
            bank.append((v, d, d*t*t - 1, t))
        require(len({b[2] for b in bank}) == n, "distinct primes")
        syndromes = [0]*(2**n)
        sigmas = [1]*(2**n)
        weights = Counter()
        formula_checks = 0
        for x in range(2**n):
            if x:
                low = x & -x
                j = low.bit_length() - 1
                syndromes[x] = syndromes[x ^ low] ^ (j + 1)
                sigmas[x] = sigmas[x ^ low] * (bank[j][2] + 1)
            require((syndromes[x] == 0) == (math.isqrt(sigmas[x])**2 == sigmas[x]), "direct sigma/code agreement")
            if syndromes[x] == 0:
                weights[x.bit_count()] += 1
            parity = x.bit_count() % 2
            for c in range(1, 2**r):
                hyperplane_sum = sum((x >> (v-1)) & 1 for v in range(1, n+1) if (v & c).bit_count() % 2 == 0) % 2
                require(((c & syndromes[x]).bit_count() % 2) == (parity ^ hyperplane_sum), "parity plus hyperplane")
                formula_checks += 1
        require(sum(weights.values()) == 2**(n-r), "code dimension")
        rows = [sum(1 << (v-1) for v in range(1, n+1) if (v & c).bit_count() % 2) for c in range(1, 2**r)]
        allbits = (1 << n) - 1
        distances = sorted({min(h.bit_count(), (h ^ allbits).bit_count()) for h in rows})
        require(distances == [2**(r-1)-1], "sharp one-check query cost")
        # Exact joint laws for the original hostile, not separate sign marginals.
        joint_checks = 0
        witness = None
        if r == 3:
            c0 = [x for x in range(128) if syndromes[x] == 0]
            c1 = [x for x in range(128) if syndromes[x] == 1]
            def law(coset, ids):
                return Counter((x.bit_count() % 2, tuple((x >> j) & 1 for j in ids)) for x in coset)
            for k in range(3):
                for ids in it.combinations(range(n), k):
                    require(law(c0, ids) == law(c1, ids), "joint parity and <=2 bits")
                    joint_checks += 1
            ids = (1, 3, 5)  # column vectors 2,4,6; primes 2,19,59
            require(law(c0, ids) != law(c1, ids), "three-bit joint witness")
            witness = [bank[j][2] for j in ids]
        out[str(r)] = dict(bank=bank, code_size=sum(weights.values()), weight_enumerator=dict(sorted(weights.items())),
                           subsets=2**n, joint_formula_checks=formula_checks,
                           one_check_queries_without_parity=2**(r-1), one_check_queries_with_parity=distances[0],
                           full_membership_queries_without_parity=n, full_membership_queries_with_parity=n-1,
                           parity_joint_law_comparisons=joint_checks, joint_witness_primes=witness)
    return out


def signed_graphs():
    totals = Counter()
    # All simple signed graphs through four vertices; 0=absent, +/-1=gain.
    for n in range(1, 5):
        pairs = list(it.combinations(range(n), 2))
        for signs in it.product((0, 1, -1), repeat=len(pairs)):
            edges = [(u, v, s) for (u, v), s in zip(pairs, signs) if s]
            adj = [[] for _ in range(n)]
            matrix = []
            for u, v, s in edges:
                adj[u].append((v, s)); adj[v].append((u, s))
                row = [0]*n; row[u] = -s; row[v] = 1
                matrix.append(row)
            rho = {}; balanced = components = 0
            for root in range(n):
                if root in rho:
                    continue
                components += 1; rho[root] = 1; stack = [root]; good = True
                while stack:
                    u = stack.pop()
                    for v, s in adj[u]:
                        if v not in rho:
                            rho[v] = s*rho[u]; stack.append(v)
                        elif rho[v] != s*rho[u]:
                            good = False
                balanced += good
            require(rankq(matrix) == n-balanced, "signed rational rank")
            require(rank2([sum((abs(x) % 2) << j for j, x in enumerate(row)) for row in matrix]) == n-components, "binary loss of signs")
            totals[n] += 1
    a, b, c = 5, 20, 44
    require([a+b, b+c, c+a] == [25, 64, 49], "square triangle")
    require(Fraction(25-64+49, 2) == a, "root fixed by odd holonomy")
    require(Fraction(4-4+1, 2) == Fraction(1, 2), "rational but nonintegral hostile")
    # A positive, distinct, even cycle has a true one-dimensional translation family.
    labels = [1, 8, 136, 440]
    require([labels[i]+labels[(i+1)%4] for i in range(4)] == [9, 144, 576, 441], "even-cycle square control")
    shifted = [labels[i] + (-1)**i for i in range(4)]
    require([labels[i]+labels[(i+1)%4] for i in range(4)] == [shifted[i]+shifted[(i+1)%4] for i in range(4)], "even affine gauge")
    return dict(graph_counts=dict(totals), total=sum(totals.values()), square_triangle=[a,b,c],
                nonintegral_square_edge_cycle=[4,4,1], even_affine_labels=labels,
                even_control_all_edges_square=True)


def frobenius():
    result = {}
    for r, modulus in [(3, 0b1011), (4, 0b10011)]:
        def mul(a, b):
            z = 0
            while b:
                if b & 1:
                    z ^= a
                a <<= 1; b >>= 1
                if a >> r:
                    a ^= modulus
            return z
        def orbit(a):
            values = [a]; b = mul(a,a)
            while b != a:
                values.append(b); b = mul(b,b)
                require(len(values) <= r, "Frobenius closes")
            return values
        orbits = []; seen = set()
        for a in range(1,2**r):
            if a not in seen:
                o = orbit(a); orbits.append(o); seen.update(o)
        normal = next(orbit(a) for a in range(1,2**r) if len(orbit(a)) == r and rank2(orbit(a)) == r)
        for word in range(2**r):
            a = xor(normal[i] for i in range(r) if word >> i & 1)
            rotated = ((word << 1) & ((1<<r)-1)) | (word >> (r-1))
            b = xor(normal[i] for i in range(r) if rotated >> i & 1)
            require(mul(a,a) == b, "normal coordinate rotation")
        result[str(r)] = dict(polynomial_bits=modulus, orbits=orbits, normal_basis=normal,
                              orbit_length_counts=dict(Counter(len(o) for o in orbits)))
    return result


def geometry15(bank):
    primes = {row[0]:row[2] for row in bank}
    points = set(range(1,16))
    lines = sorted({tuple(sorted((a,b,a^b))) for a,b in it.combinations(points,2)})
    require(len(lines) == 35, "PG(3,2) lines")
    def spreads(allowed):
        found = []
        def step(remaining, chosen):
            if not remaining:
                found.append(tuple(chosen)); return
            first = min(remaining)
            for line in allowed:
                if first in line and set(line) <= remaining:
                    step(remaining-set(line), chosen+[line])
        step(points, [])
        return found
    pairs = list(it.combinations(range(4),2))
    def form(a,b,mask):
        return sum(((a>>i&1)*(b>>j&1) ^ (a>>j&1)*(b>>i&1)) for k,(i,j) in enumerate(pairs) if mask>>k&1) % 2
    forms = [mask for mask in range(64) if rank2([sum(form(1<<i,1<<j,mask)<<j for j in range(4)) for i in range(4)])==4]
    require(len(forms)==28, "nondegenerate alternating forms")
    # Frobenius in F16 = F2[a]/(a^4+a+1), computed as reduction of squares.
    def square(a):
        value = sum(((a>>i)&1) << (2*i) for i in range(4))
        while value.bit_length()>4:
            value ^= 0b10011 << (value.bit_length()-5)
        return value
    invariant = [m for m in forms if all(form(square(a),square(b),m)==form(a,b,m) for a in points for b in points)]
    require(len(invariant)==2, "two Frobenius-compatible symplectic forms")
    # Standard pairing of coordinate0 with2, coordinate1 with3.
    standard = (1 << pairs.index((0,2))) | (1 << pairs.index((1,3)))
    isotropic = [line for line in lines if form(line[0],line[1],standard)==0]
    require(len(isotropic)==15, "doily lines")
    require(all(sum(v in line for line in lines)==7 for v in points), "seven lines per point")
    require(all(sum(v in line for line in isotropic)==3 for v in points), "three isotropic lines per point")
    allspreads, isospreads = spreads(lines), spreads(isotropic)
    require(len(allspreads)==56 and len(isospreads)==6, "spread counts")
    marked=[]
    for m in invariant:
        iso=[line for line in lines if form(line[0],line[1],m)==0]
        spr=spreads(iso)
        fixed=sum(set(s)=={tuple(sorted(square(v) for v in line)) for line in s} for s in spr)
        marked.append(dict(form_mask=m, isotropic_line_count=len(iso), spreads=len(spr), Frobenius_fixed_spreads=fixed))
    def mul(a,b):
        z=0
        while b:
            if b&1: z ^= a
            a <<= 1; b >>= 1
            if a>>4: a ^= 19
        return z
    def trace(a):
        return a ^ square(a) ^ square(square(a)) ^ square(square(square(a)))
    subfield = [a for a in points if square(square(a))==a]
    fieldspread = sorted({tuple(sorted(mul(a,x) for a in subfield)) for x in points})
    compatible = [m for m in forms if all(form(line[0],line[1],m)==0 for line in fieldspread)]
    def traceform(c,x,y):
        return trace(mul(c,mul(x,square(square(y)))))
    traceforms = [sum(traceform(c,1<<i,1<<j)<<k for k,(i,j) in enumerate(pairs)) for c in subfield]
    require(set(traceforms)==set(compatible) and len(compatible)==3, "subfield-spread trace classification")
    require(set(compatible)&set(invariant)=={60}, "unique transport and spread compatible form")
    require([m for m in invariant if form(2,12,m)==0]==[60], "one marked line suffices")
    sample = [[primes[v] for v in line] for line in isospreads[0]]
    require(all(math.isqrt(math.prod(p+1 for p in line))**2==math.prod(p+1 for p in line) for line in sample), "five arithmetic square triples")
    return dict(lines=35, lines_through_point=7, isotropic_lines=15, all_spreads=56,
                symplectic_spreads=6, lines_per_spread=5, alternating_forms=28,
                sample_prime_partition=sample, standard_form_mask=standard,
                Frobenius_compatible_forms=marked, subfield_nonzero_vectors=subfield,
                subfield_spread=fieldspread, subfield_spread_compatible_forms=compatible,
                unique_spread_and_Frobenius_form=60, single_line_separator=[2,12,14])


def collatz():
    checks = 0
    for sign in (-1, 1):
        for k in range(1, 11):
            seen = set()
            for bits in it.product((0, 1), repeat=k):
                b = r = 0
                for i, e in enumerate(bits):
                    b = (3**e)*b + e*(2**i); r += e
                residue = (-sign*b*pow(3**r, -1, 2**k)) % (2**k)
                require(residue not in seen, "parity bijection"); seen.add(residue)
                next_bits = []
                for lift in (0, 1):
                    start = residue + (lift+1)*2**k
                    x = start
                    got = []
                    for _ in bits:
                        e = x % 2; got.append(e)
                        x = (3*x+sign)//2 if e else x//2
                    require(tuple(got) == bits, "word realized")
                    require(x*2**k == 3**r*start+sign*b, "affine endpoint")
                    next_bits.append(x % 2)
                require(next_bits[0] != next_bits[1], "next bit not free")
                checks += 1
    return dict(signs=[-1,1], max_depth=10, words_checked=checks,
                claim="parity/residue bijection and two opposite next-bit lifts; no convergence claim")


def degree_and_geometry():
    checks = 0; rows = []
    for moduli in [(5,), (7,), (15,), (5,7), (5,15), (7,15), (8,9)]:
        m = math.lcm(*moduli)
        for d in (m-2, m-1):
            aliases = [(a,b) for a in range(d+1) for b in range(d+1-a)
                       if all((a-b-1) % q == 0 for q in moduli) and a-b != 1]
            require((len(aliases) == 0) == (d == m-2), "sharp charge bandwidth")
            if d == m-1:
                require((0,m-1) in aliases, "boundary shear")
            checks += 1
        rows.append(dict(moduli=moduli, effective_modulus=m, linear_through_degree=m-2, first_nonlinear_equivariant_shear_degree=m-1))
    gram = [[Fraction(1) if i==j else Fraction(1,2) for j in range(3)] for i in range(3)]
    require(rankq(gram)==3, "K4 unit Gram cannot be planar")
    # Two anchors do not resolve reflection; a third noncollinear anchor does.
    anchors = [(0,0),(2,0),(0,2)]
    dist = lambda p,q: sum((x-y)**2 for x,y in zip(p,q))
    p, q = (1,1), (1,-1)
    require([dist(p,a) for a in anchors[:2]] == [dist(q,a) for a in anchors[:2]], "mirror hostile")
    require(dist(p,anchors[2]) != dist(q,anchors[2]), "third anchor resolves bit")
    return dict(charge_tests=checks, degree_table=rows, K4_gram_rank=3,
                mirror_points=[p,q], anchors=anchors)


def main():
    arithmetic = codes()
    result = dict(arithmetic_codes=arithmetic, geometry15=geometry15(arithmetic['4']['bank']), signed_incidence=signed_graphs(),
                  frobenius=frobenius(), parity_lifts=collatz(), bounded_lifts=degree_and_geometry())
    text = json.dumps(result, indent=2, sort_keys=True)+"\n"
    parser = argparse.ArgumentParser(); parser.add_argument('--output'); args = parser.parse_args()
    if args.output:
        Path(args.output).write_text(text, encoding='utf-8', newline='\n')
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
