"""Exact source disintegration and a finite kernel for infinite sibling flows.

Run with Python (standard library only). The finite kernel is constructed by
bounded explicit searches, then verified separately. No universal search bound
or positivity outside the certified sibling families is asserted.
"""
from fractions import Fraction as F
from collections import defaultdict
import json
from pathlib import Path

CHECKS = 0


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise AssertionError(why)


def v2(n):
    if n <= 0:
        raise ValueError("positive integer required")
    return (n & -n).bit_length() - 1


def odd(n):
    if type(n) is not int or n < 1 or not n & 1:
        raise ValueError("positive odd integer required")


def U(n):
    odd(n)
    return (3*n+1) >> v2(3*n+1)


def S(n, k=1):
    odd(n)
    if type(k) is not int or k < 0:
        raise ValueError("nonnegative depth required")
    return 4**k*n + (4**k-1)//3


def gamma(m):
    if type(m) is not int or m < 1:
        raise ValueError("positive source index required")
    return F(2, 4**m.bit_length())


def mu(n):
    odd(n)
    return gamma((n+1)//2)


def nu(n):
    odd(n)
    return F(8, 3*4**n.bit_length())


def q(h):
    return F(3, 4**h)


def base(n):
    """Independent inverse-sibling reader; no Collatz orbit search."""
    odd(n)
    k = 0
    while n % 8 == 5:
        n = (n-1)//4
        k += 1
    return n, k


def first_predecessor(y):
    odd(y)
    if y % 3 == 0:
        return None
    return ((4 if y % 3 == 1 else 2)*y-1)//3


def path_to_base_root(b, limit=10000):
    """Discovery only: supplied finite universe, explicit cap, fails closed."""
    path = [b]
    seen = {b}
    while path[-1] != 1:
        if len(path) > limit:
            raise RuntimeError("uncertified base path")
        c, _ = base(U(path[-1]))
        if c in seen:
            raise RuntimeError("nonroot quotient cycle")
        seen.add(c)
        path.append(c)
    return path


def verify_kernel(entries):
    """Finite semantic verifier, including exact edges and root reachability."""
    for row in entries:
        if set(row) != {'b', 'parent', 'depth', 'valuation'}:
            raise ValueError("unexpected kernel fields")
        odd(row['b'])
        if (type(row['depth']) is not int or row['depth'] < 0 or
                type(row['valuation']) is not int or row['valuation'] < 1):
            raise ValueError("invalid integer fields")
        if row['parent'] is not None:
            odd(row['parent'])
    nodes = {row['b']: row for row in entries}
    if len(nodes) != len(entries) or 1 not in nodes:
        raise ValueError("duplicate bases or missing ROOT")
    for b, row in nodes.items():
        odd(b)
        if v2(3*b+1) not in (1, 2):
            raise ValueError("nonprimitive base")
        if b == 1:
            if row != {'b': 1, 'parent': None, 'depth': 0, 'valuation': 2}:
                raise ValueError("bad ROOT record")
        else:
            c, k, a = row['parent'], row['depth'], row['valuation']
            if c not in nodes or type(k) is not int or k < 0:
                raise ValueError("ungrounded parent")
            if a != v2(3*b+1) or (3*b+1) != 2**a*S(c, k):
                raise ValueError("false literal edge")
    for b in nodes:
        seen = set()
        while b != 1:
            if b in seen:
                raise ValueError("cyclic dependency")
            seen.add(b)
            b = nodes[b]['parent']
    return nodes


def compile_root_word(entries, n):
    """Decode a strict first-hit word from a checked kernel, without search."""
    nodes = verify_kernel(entries)
    b, depth = base(n)
    if b not in nodes:
        raise ValueError("source outside certified sibling families")
    words = {1: ()}
    path = []
    c = b
    while c not in words:
        path.append(c)
        c = nodes[c]['parent']
    for c in reversed(path):
        row = nodes[c]
        parent, k = row['parent'], row['depth']
        suffix = words[parent]
        if k:
            suffix = ((2+2*k,) if parent == 1 else
                      (suffix[0]+2*k,)+suffix[1:])
        words[c] = (row['valuation'],)+suffix
    result = words[b]
    if depth:
        result = ((2+2*depth,) if b == 1 else
                  (result[0]+2*depth,)+result[1:])
    return result


def main():
    # Normalization and entropy are proved symbolically in the companion note.
    for cutoff in range(2, 21):
        height_head = F(2, 3) + sum(F(2*l, 3*2**l) for l in range(2, cutoff))
        height_tail = F(2*(cutoff+1), 3*2**(cutoff-1))
        check(height_head+height_tail == F(5, 3), 'nu ordinary height mean')
        climb_head = sum(h*q(h) for h in range(1, cutoff))
        climb_tail = 4*(cutoff+F(1, 3))/4**cutoff
        check(climb_head+climb_tail == F(4, 3), 'q climb mean')
    for n in range(1, 16384, 2):
        h = v2(n+1)
        t = (n+1) >> h
        check(mu(n) == q(h)*nu(t), ('product', n))
        m = (n+1)//2
        dyadic = m & (m-1) == 0
        check(nu(n)/mu(n) == F(4 if dyadic else 1, 3), ('comparison', n))
        b, k = base(n)
        check(S(b, k) == n and v2(3*b+1) in (1, 2), ('base', n))
        check(k == (v2(3*n+1)-1)//2, ('depth', n))
        for j in range(1, 7):
            z = S(n, j)
            check(U(z) == U(n), ('sibling target', n, j))
            check(nu(z) == nu(n)/16**j, ('nu scaling', n, j))
            check(mu(z)/mu(n) == F(4 if dyadic else 1, 16**j), ('mu scaling', n, j))
    for m in range(1, 513):
        for d in range(4):
            check(gamma(4*m+d) == gamma(m)/16, ('radix regrouping', m, d))
    check(sum(gamma(m) for m in (1, 2, 3)) == F(3, 4), 'three heads')

    incoming_classes = defaultdict(int)
    for y in range(1, 1024, 2):
        b = first_predecessor(y)
        if b is None:
            check(y % 3 == 0, 'empty inverse fibre')
            incoming_classes['0'] += 1
            continue
        check(U(b) == y and base(b) == (b, 0), ('minimal inverse', y))
        # An independent direct valuation enumeration versus sibling powers.
        direct = []
        for a in range(1, 26):
            num = 2**a*y-1
            if num % 3 == 0 and (num//3) & 1:
                direct.append(num//3)
        check(direct == [S(b, k) for k in range(len(direct))], ('inverse enumeration', y))
        j = 12
        lhs = sum(nu(S(b, k)) for k in range(j))
        check(lhs + nu(b)*F(16, 15*16**j) == F(16, 15)*nu(b), ('exact tail', y))
        ratio = F(16, 15)*nu(b)/nu(y)
        check(ratio in (F(4, 15), F(16, 15), F(64, 15)), ('incoming ratio', y))
        incoming_classes[str(ratio)] += 1

    # Compile a finite grounded base forest, then close each base under all S^k.
    seeds = [b for b in range(1, 128, 2) if v2(3*b+1) in (1, 2)]
    paths = [path_to_base_root(b) for b in seeds]
    nodes = sorted({b for path in paths for b in path})
    entries = []
    for b in nodes:
        c, k = base(U(b))
        entries.append({'b': b, 'parent': None if b == 1 else c,
                        'depth': 0 if b == 1 else k, 'valuation': v2(3*b+1)})
    verify_kernel(entries)
    check(True, 'finite semantic verifier')
    longest = 0
    for b in nodes:
        for k in range(7):
            n = S(b, k)
            word = compile_root_word(entries, n)
            longest = max(longest, len(word))
            for a in word:
                check(n > 1 and a == v2(3*n+1), ('strict compiled edge', b, k))
                n = (3*n+1) >> a
            check(n == 1, ('compiled ROOT', b, k))
    missing = next(b for b in range(1, 10000, 2) if b not in nodes and base(b) == (b, 0))
    try:
        compile_root_word(entries, missing)
    except ValueError:
        check(True, 'unlisted base is not silently certified')
    else:
        raise AssertionError('accepted unlisted base')
    for mutate in ('edge', 'root', 'duplicate', 'boolean', 'float'):
        bad = [dict(row) for row in entries]
        if mutate == 'edge':
            next(row for row in bad if row['b'] != 1)['depth'] += 1
        elif mutate == 'root':
            bad[0]['parent'] = 1
        elif mutate == 'duplicate':
            bad.append(dict(bad[0]))
        elif mutate == 'boolean':
            bad[0]['depth'] = False
        else:
            bad[0]['valuation'] = 2.0
        try:
            verify_kernel(bad)
        except ValueError:
            check(True, ('rejected', mutate))
        else:
            raise AssertionError(('accepted malformed kernel', mutate))

    r, rho = F(1, 16), F(1, 2)
    g = defaultdict(F)
    for j, path in enumerate(paths, 1):
        mass = F(1, 2**j*len(path))
        g[1] += mass
        for b in reversed(path[:-1]):
            c, k = base(U(b))
            mass *= rho*(1-r)*r**k
            g[b] += mass
    for b in nodes:
        check(g[b] > 0, ('positive supported base', b))
        if b != 1:
            c, k = base(U(b))
            check(g[b] <= rho*(1-r)*r**k*g[c], ('finite flow inequality', b))
    check(sum(g.values()) < 1, 'summable base bound')
    mass_nu = F(16, 15)*sum(nu(b) for b in nodes)
    mass_mu = sum(mu(b)*F(19 if ((b+1)//2)&(((b+1)//2)-1)==0 else 16, 15) for b in nodes)
    check(0 < mass_mu < 1 and 0 < mass_nu < 1, 'proper infinite certified subsets')

    payload = {'status': 'FINITE-EXACT kernel; its sibling closure is infinite and certified',
               'seed_bound': 127, 'seeds': seeds, 'r': str(r), 'rho': str(rho),
               'bases': entries, 'coverage_mu': str(mass_mu), 'coverage_nu': str(mass_nu)}
    target = Path(__file__).resolve().parents[2]/'05-knowledge/results/collatz_three_bit_sibling_flow_20261005.json'
    target.write_text(json.dumps(payload, indent=2)+'\n', encoding='utf-8')
    print('Exact checks:', CHECKS)
    print('Source disintegration: mu(2^h*t-1) = (3/4^h)*nu(t)')
    print('Entropy: (8/3-log2(3)) + (log2(3)+1/3) = 3 bits')
    print('Incoming ratios through odd target1023 (full, not killed at1):', dict(sorted(incoming_classes.items())))
    print('Kernel seed bases:', len(seeds), 'closed bases:', len(nodes), 'largest base:', max(nodes))
    print('Maximum quotient path length:', max(len(p)-1 for p in paths))
    print('Longest compiled first-hit ROOT word (bases and sibling depths0..6):', longest)
    print('Infinite sibling-closure coverage mu:', mass_mu)
    print('Infinite sibling-closure coverage nu:', mass_nu)
    print('All supported base flow inequalities passed; total base weight < 1')
    print('Kernel JSON:', target.name)
    print('No positive weight or ROOT certificate asserted outside these sibling families.')


if __name__ == '__main__':
    main()
