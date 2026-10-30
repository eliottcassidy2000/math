"""Exact checks; run with normal Python and -O. Standard library only."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import gcd, prod
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import random
import sys

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(100000)
HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent/'regular_parent_symmetry_20261029/verify.py'
spec = importlib.util.spec_from_file_location('native_parent', SOURCE)
native = importlib.util.module_from_spec(spec)
spec.loader.exec_module(native)
GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(why)


def sha(x):
    return hashlib.sha256(json.dumps(x, sort_keys=True).encode()).hexdigest()


def mul(a, b):
    c = defaultdict(int)
    for i, x in a.items():
        for j, y in b.items():
            c[i+j] += x*y
    return dict(c)


def polypow(a, n):
    c = {0: 1}
    for _ in range(n):
        c = mul(c, a)
    return c


def inventory_step(poly, p):
    c = defaultdict(int, polypow(poly, p))
    for i, v in poly.items():
        c[p*i] -= v
    need(all(v % p == 0 and v >= 0 for v in c.values()), 'nonconstant orbit divisibility')
    c = {i: v//p for i, v in c.items() if v}
    for i, v in poly.items():
        c[1+p*i] = c.get(1+p*i, 0)+v
    return c


def inventories():
    rows = []
    for p in (3, 5, 7, 11, 13):
        poly = {0: 2}
        free = 2
        for k in range(4):
            n, g = p**k, (p**k-1)//(p-1)
            need(p**g * sum(F(v, p**i) for i, v in poly.items()) == 2**n, 'labelled mass')
            need(poly[0] == free, 'free recurrence')
            child = poly.copy()
            child[g-k] -= 1
            need(min(child.values()) >= 0, 'clone subtraction')
            need(p**g * sum(F(v, p**i) for i, v in child.items()) == 2**n-n, 'unmarked mass')
            if k:
                flat = {j: (2**(p**(k-j))-2**(p**(k-j-1)))//p**(k-j) for j in range(k)}
                flat[k] = 2
                flat[0] -= 1
                need(p**k * sum(F(v, p**i) for i, v in flat.items()) == 2**n-n, 'flat mass')
            else:
                flat = {0: 1}
            rows.append({'p': p, 'k': k, 'child_order': n+1,
                         'parent_group_exponent': g,
                         'hierarchical_children': str(sum(child.values())),
                         'hierarchical_rigid': str(child[0]),
                         'flat_children': str(sum(flat.values())),
                         'flat_rigid': str(flat[0]),
                         'inventory_sha256': sha(sorted(child.items()))})
            if k < 3:
                poly = inventory_step(poly, p)
                free = (free**p-free)//p
    return rows


def color_type(bits, p):
    if len(bits) == 1:
        return bits[0], 1
    m = len(bits)//p
    children = [color_type(bits[i*m:(i+1)*m], p) for i in range(p)]
    types = tuple(x[0] for x in children)
    key = min(types[i:]+types[:i] for i in range(p))
    return key, prod(x[1] for x in children)*(p if len(set(types)) == 1 else 1)


def native_checks():
    reports = []
    for p, k in ((3, 1), (5, 1), (7, 1), (3, 2)):
        base = native.carousel(p)
        parent = [0]
        for _ in range(k):
            parent = native.lex_product(base, parent)
        n = len(parent)
        count = Counter()
        keys = {}
        for mask in range(1 << n):
            key, predicted = color_type(tuple(mask >> i & 1 for i in range(n)), p)
            child = native.add_vertex(parent, mask)
            maps = native.iso_maps(child, child)
            need(len(maps) == predicted, 'native full child group')
            need(all(g[-1] == n for g in maps), 'intrinsic new root')
            count[predicted] += 1
            keys[key] = predicted
        poly = {0: 2}
        for _ in range(k):
            poly = inventory_step(poly, p)
        expected = Counter({p**e: v for e, v in poly.items()})
        need(Counter(keys.values()) == expected, 'color orbit inventory')
        reports.append({'p': p, 'k': k, 'all_binary_cuts': 1 << n,
                        'rooted_orbits': len(keys), 'labelled_group_histogram': dict(count)})
    # Independent flat-nine group checks, all cuts.
    parent = native.carousel(9)
    hist = Counter()
    for mask in range(512):
        predicted = sum(native.permute_mask(mask, tuple((i+t) % 9 for i in range(9))) == mask
                        for t in range(9))
        need(len(native.iso_maps(native.add_vertex(parent, mask), native.add_vertex(parent, mask))) == predicted,
             'flat9 actual graph group')
        hist[predicted] += 1
    reports.append({'flat9_all_cuts': 512, 'labelled_group_histogram': dict(hist)})
    # Larger actual targets use deterministic cuts, not a full census.
    rng = random.Random(310507)
    for p in (5, 7):
        parent = native.lex_product(native.carousel(p), native.carousel(p))
        n = p*p
        groups = []
        for _ in range(8):
            mask = rng.randrange(1 << n)
            _, predicted = color_type(tuple(mask >> i & 1 for i in range(n)), p)
            child = native.add_vertex(parent, mask)
            need(len(native.iso_maps(child, child)) == predicted, 'larger native target')
            groups.append((str(mask), predicted))
        reports.append({'p': p, 'k': 2, 'sample_cuts': groups})
    return reports


AFFINE = {'paley7': (7, 3), 'paley11': (11, 5), 'local13': (13, 3)}


@lru_cache(None)
def moment(kind, k, q, s):
    if k == 0:
        return F(1)
    b, l = AFFINE.get(kind, (int(kind) if kind.isdigit() else 0, 1))
    h = b*l
    m = b**(k-1)
    a = F(h**((m-1)//(b-1)), q**m)
    v = moment(kind, k-1, q, s)
    if kind in AFFINE:
        t = (b-1)//l
        return (v**b + b*(l**s-1)*a**(t*(l-1))*v*moment(kind, k-1, q, l*(s-1)+1)**t
                +(h**s-b*l**s+b-1)*a**(b-1)*moment(kind, k-1, q, b*(s-1)+1))
    p = b
    return v**p + (p**s-1)*a**(p-1)*moment(kind, k-1, q, p*(s-1)+1)


def moment_checks():
    rec = []
    for kind in ('3', '5', '7', *AFFINE):
        b, l = AFFINE.get(kind, (int(kind) if kind.isdigit() else 0, 1))
        h, d = b*l, (b-1)*(l-1)//l if l > 1 else b-1
        for q in (2, 4, 8):
            poly = {0: q}
            orbit_count = q
            for k in range(1, 4):
                if kind not in AFFINE:
                    poly = inventory_step(poly, b)
                else:
                    orbit_count = (orbit_count**b+b*(l-1)*orbit_count**(1+(b-1)//l)
                                   +(b-1)*orbit_count)//h
                for s in (1, 2, 3):
                    val = moment(kind, k, q, s)
                    if kind not in AFFINE:
                        other = F(b**((b**k-1)//(b-1)), q**(b**k))*sum(v*b**((s-1)*i) for i, v in poly.items())
                        need(val == other, 'moment recursion versus orbit polynomial')
                    elif s == 1:
                        need(val == F(h**((b**k-1)//(b-1))*orbit_count, q**(b**k)),
                             'affine first moment versus independent Burnside orbit recurrence')
                    low = moment(kind, 1, q, s)**(b**(k-1))
                    upper = F(1)
                    for j in range(2, k+1):
                        m = b**(j-1)
                        gm = h**((m-1)//(b-1))
                        bump = (h**s-1)*(h-1)*F(gm**s, q**m)**d
                        upper *= (1+bump)**(b**(k-j))
                    need(low <= val <= low*upper, 'finite weighted ratio bound')
                    rec.append((kind, k, q, s, str(val.numerator), str(val.denominator)))
    # Literal nonuniform binary distributions on all 3^2 leaves.
    nonuniform = []
    for weights in ((F(1, 2), F(1, 2)), (F(1, 3), F(2, 3)), (F(1, 10), F(9, 10))):
        vals = {s: F(0) for s in (1, 2, 3)}
        for mask in range(512):
            bits = tuple(mask >> i & 1 for i in range(9))
            st = color_type(bits, 3)[1]
            prob = prod(weights[x] for x in bits)
            for s in vals:
                vals[s] += prob*st**s
        cp = sum(w**3 for w in weights)
        for s, val in vals.items():
            low = (1+(3**s-1)*cp)**3
            bump = (3**s-1)*2*(3**s*max(weights)**3)**2
            need(low <= val <= low*(1+bump), 'nonuniform weighted bound')
        nonuniform.append({'weights': list(map(str, weights)), 'c3': str(cp),
                           'moments': {s: str(v) for s, v in vals.items()}})
    # Independent native full groups and literal nonuniform binary color laws.
    affine = []
    for kind, (p, l) in AFFINE.items():
        residues = ({x*x % p for x in range(1, p)} if p in (7, 11)
                    else {1, 2, 3, 5, 6, 9})
        parent = [sum(1 << y for y in range(p) if (y-x) % p in residues) for x in range(p)]
        g = native.iso_maps(parent, parent)
        need(len(g) == p*l, 'full affine base group')
        cycle_types = Counter()
        for perm in g:
            unseen, sizes = set(range(p)), []
            while unseen:
                v = min(unseen)
                cycle = set()
                while v not in cycle:
                    cycle.add(v)
                    v = perm[v]
                unseen -= cycle
                sizes.append(len(cycle))
            cycle_types[tuple(sorted(sizes))] += 1
        need(cycle_types == Counter({(1,)*p: 1, (p,): p-1, (1,)+(l,)*((p-1)//l): p*(l-1)}),
             'affine full cycle law')
        hist = Counter()
        for bits in product(range(2), repeat=p):
            size = sum(all(bits[i] == bits[perm[i]] for i in range(p)) for perm in g)
            hist[(sum(bits), size)] += 1
        for weights in ((F(1, 2), F(1, 2)), (F(1, 3), F(2, 3))):
            actual = Counter()
            for (ones, size), count in hist.items():
                actual[size] += count*weights[0]**(p-ones)*weights[1]**ones
            cl, cp = sum(w**l for w in weights), sum(w**p for w in weights)
            need(actual[p*l] == cp and actual[l] == p*(cl**((p-1)//l)-cp), 'affine collision law')
            need(set(actual) <= {1, l, p*l} and sum(actual.values()) == 1, 'complete base color law')
        affine.append({'base': kind, 'group_order': len(g), 'all_binary_cuts': 2**p,
                       'histogram_sha256': sha(sorted(hist.items()))})
    return {'exact_moment_queries': len(rec), 'sha256': sha(rec), 'nonuniform': nonuniform,
            'affine_base_checks': affine}


def factor(n):
    out = {}
    d = 2
    while d*d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0)+1
            n //= d
        d = 3 if d == 2 else d+2
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out


def order(a, n):
    need(gcd(a, n) == 1, 'unit order')
    phi = n
    for p in factor(n):
        phi = phi//p*(p-1)
    d = phi
    for p in factor(phi):
        while d % p == 0 and pow(a, d//p, n) == 1:
            d //= p
    return d


def valuation(n, p):
    e = 0
    while n % p == 0:
        e += 1
        n //= p
    return e


def arithmetic_checks():
    rows = []
    for p in (3, 5, 7, 11, 13):
        for k in range(1, 11):
            n = p**k+1
            fs = factor(n)
            need(prod(ell**e for ell, e in fs.items()) == n, 'complete factorization')
            need(fs[2] == (valuation(p+1, 2) if k % 2 else 1), 'two-adic split')
            data = []
            for ell, e in fs.items():
                if ell == 2:
                    continue
                d = order(p, ell)
                need(valuation(d, 2) == valuation(k, 2)+1 and (2*k) % d == 0, 'negative order')
                need(k % 2 or ell % 4 == 1, 'even exponent Paley barrier')
                need(pow(p, k, ell**e) == ell**e-1, 'full prime-power source')
                data.append((ell, e, d))
            need(k == 1 or any(d == 2*k for _, _, d in data), 'primitive factor finite control')
            rows.append({'p': p, 'k': k, 'N': n, 'factors': fs, 'odd_factor_orders': data})
    return rows


def bank_checks():
    records = []
    for p, ell in ((3, 7), (5, 7), (3, 5), (7, 5), (5, 13), (7, 43)):
        h = order(p, ell)
        a = valuation(h, 2)-1
        u = h//2**(a+1)
        delta = valuation(p**h-1, ell)
        depths = []
        for e in range(1, 5):
            n = ell**e
            ue = u*ell**max(0, e-delta)
            need(order(p, n) == 2**(a+1)*ue, 'all-depth order formula')
            bits = (ell**min(e, delta)-1)//(2*u) + max(0, e-delta)*ell**(delta-1)*(ell-1)//(2*u)
            independently = sum((ell**j-ell**(j-1))//(2*(order(p, ell**j)//2**(a+1))) for j in range(1, e+1))
            need(bits == independently, 'valuation-layer bit sum')
            depths.append({'e': e, 'u': ue, 'orientation_bits': bits})
            if n > 169 or bits > 12:
                continue
            r = pow(p, 2**a, n)
            generator = r*r % n
            unseen = set(range(1, n))
            pairs = []
            while unseen:
                x = min(unseen)
                orbit, y = set(), x
                while y not in orbit:
                    orbit.add(y)
                    y = y*generator % n
                negative = {(-v) % n for v in orbit}
                need(not orbit & negative, 'skew orbit guard')
                unseen -= orbit | negative
                pairs.append((orbit, negative))
            need(len(pairs) == bits, 'actual multiplier orbits')
            for choices in range(1 << bits):
                s = set().union(*(pair[choices >> i & 1] for i, pair in enumerate(pairs)))
                negative = {(-v) % n for v in s}
                need(len(s) == (n-1)//2 and not s & negative, 'actual tournament connection set')
                need({r*v % n for v in s} == negative, 'literal converse multiplier')
                need({generator*v % n for v in s} == s, 'literal symmetry multiplier')
                if choices < 2:
                    for x in range(n):
                        for y in range(n):
                            if x != y:
                                need(((y-x) % n in s) != ((r*y-r*x) % n in s), 'full arc reversal')
            depths[-1]['exhausted_states'] = 1 << bits
        records.append({'p': p, 'ell': ell, 'a': a, 'delta': delta, 'depths': depths})
    return records


def frobenius_roots():
    # F3[X]/Phi7(X), irreducible because ord_7(3)=6. Vectors have length6.
    def multiply(a, b):
        c = [0]*11
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                c[i+j] += x*y
        for i in range(10, 5, -1):
            for j in range(6):
                c[i-6+j] -= c[i]
        return tuple(v % 3 for v in c[:6])
    one, zeta = (1, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0)
    roots = [one]
    for _ in range(6):
        roots.append(multiply(roots[-1], zeta))
    need(multiply(roots[-1], zeta) == one and len(set(roots)) == 7, 'primitive root literal')
    f = [roots.index(multiply(multiply(x, x), x)) for x in roots]
    need(f == [3*i % 7 for i in range(7)], 'Frobenius exponent transport')
    for x, y in product(range(7), repeat=2):
        if x != y:
            arc = (y-x) % 7 in (1, 2, 4)
            need(arc != ((f[y]-f[x]) % 7 in (1, 2, 4)), 'Frobenius is converse on roots')
            need(arc == ((f[f[y]]-f[f[x]]) % 7 in (1, 2, 4)), 'squared Frobenius preserves')
    return {'field': 'F3[X]/(1+X+...+X^6)', 'root_vectors': roots, 'Frobenius_on_exponents': f}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=HERE/'certificate.json')
    args = parser.parse_args()
    out = {'status': 'FINITE-EXACT, with all-size proofs in PROOFS.md',
           'inventories': inventories(), 'native_graphs': native_checks(),
           'moments': moment_checks(), 'arithmetic': arithmetic_checks(),
           'signed_banks': bank_checks(), 'Frobenius': frobenius_roots(),
           'inherited_native_engine_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest()}
    out['gates'] = GATES
    text = json.dumps(out, indent=2, sort_keys=True)+'\n'
    args.output.write_text(text, encoding='utf-8')
    print(json.dumps({'gates': GATES, 'receipt_lf_sha256': hashlib.sha256(text.encode()).hexdigest(),
                      'output': str(args.output)}))


if __name__ == '__main__':
    main()
