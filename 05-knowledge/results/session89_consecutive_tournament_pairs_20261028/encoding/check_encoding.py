#!/usr/bin/env python3
"""Independent stdlib checks for the fixed-frame matrix/pair correspondence.
Run: python3 check_encoding.py
No external files, packages, randomization, or repository writes.
"""
from collections import Counter
from fractions import Fraction
from itertools import permutations
from math import factorial, gcd
import json

CHECKS = 0

def need(condition, message="verification failed"):
    """An active check under normal and optimized Python alike."""
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


def pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def encode(m, n):
    """Return actual adjacency rows for T_n and U_(n+1), root n."""
    t = [0] * n
    u = [0] * (n + 1)
    for i, j in pairs(n):
        if (m >> (i * n + j)) & 1:
            t[i] |= 1 << j
        else:
            t[j] |= 1 << i
        if (m >> (j * n + i)) & 1:
            u[i] |= 1 << j
        else:
            u[j] |= 1 << i
    for i in range(n):
        if (m >> (i * n + i)) & 1:
            u[n] |= 1 << i
        else:
            u[i] |= 1 << n
    return tuple(t), tuple(u)


def decode(t, u, n):
    m = 0
    for i, j in pairs(n):
        m |= ((t[i] >> j) & 1) << (i * n + j)
        m |= ((u[i] >> j) & 1) << (j * n + i)
    for i in range(n):
        m |= ((u[n] >> i) & 1) << (i * n + i)
    return m


def relabel_adj(rows, p):
    out = [0] * len(rows)
    for i, row in enumerate(rows):
        for j in range(len(rows)):
            out[p[i]] |= ((row >> j) & 1) << p[j]
    return tuple(out)


def direct_pair_action(m, p):
    n = len(p)
    t, u = encode(m, n)
    return decode(relabel_adj(t, p), relabel_adj(u, p + (n,)), n)


def alpha(m, p):
    """Natural aligned-pair action transported to fixed-frame matrices."""
    n = len(p)
    out = 0
    for i, j in pairs(n):
        a, b = sorted((p[i], p[j]))
        flip = int(p[i] > p[j])
        out |= (((m >> (i * n + j)) & 1) ^ flip) << (a * n + b)
        out |= (((m >> (j * n + i)) & 1) ^ flip) << (b * n + a)
    for i in range(n):
        out |= ((m >> (i * n + i)) & 1) << (p[i] * n + p[i])
    return out


def ordinary(m, p):
    n = len(p)
    out = 0
    for i in range(n):
        for j in range(n):
            out |= ((m >> (i * n + j)) & 1) << (p[i] * n + p[j])
    return out


def rowcol(m, p, q):
    n = len(p)
    return sum(((m >> (i * n + j)) & 1) << (p[i] * n + q[j])
               for i in range(n) for j in range(n))


def transpose(m, n):
    return sum(((m >> (i * n + j)) & 1) << (j * n + i)
               for i in range(n) for j in range(n))


def cycle_type(p):
    unused = set(range(len(p)))
    lengths = []
    while unused:
        i = next(iter(unused))
        length = 0
        while i in unused:
            unused.remove(i)
            length += 1
            i = p[i]
        lengths.append(length)
    return tuple(sorted(lengths))


def partitions(n, low=1):
    if not n:
        yield ()
    for a in range(low, n + 1):
        for tail in partitions(n - a, a):
            yield (a,) + tail


def representative(lam):
    p = []
    start = 0
    for size in lam:
        p.extend(range(start + 1, start + size))
        p.append(start)
        start += size
    return tuple(p)


def z(lam):
    out = 1
    for size, multiplicity in Counter(lam).items():
        out *= size ** multiplicity * factorial(multiplicity)
    return out


def edge_orbits(lam):
    """Used only when every part is odd."""
    need(all(a % 2 for a in lam), "check at source line 146")
    return sum((a - 1) // 2 for a in lam) + sum(
        gcd(a, b) for i, a in enumerate(lam) for b in lam[i + 1:])


def ordered_orbits(lam):
    return sum(gcd(a, b) for a in lam for b in lam)


def fixed_formula(lam, natural):
    if natural and any(a % 2 == 0 for a in lam):
        return 0
    return 1 << ordered_orbits(lam)


def swapped_fixed_formula(lam):
    """Fix(alpha_pi composed with transpose/core exchange)."""
    if any(a % 4 == 0 for a in lam):
        return 0
    square = tuple(a // gcd(a, 2) for a in lam for _ in range(gcd(a, 2)))
    return 1 << (edge_orbits(square) + len(lam))


def burnside(n, natural):
    ans = sum((Fraction(fixed_formula(lam, natural), z(lam))
               for lam in partitions(n)), Fraction())
    need(ans.denominator == 1, "check at source line 172")
    return ans.numerator


def affine_fixed_count(fn, dimension):
    """Solve (linear part + I)x=translation; no orbit-count formulas."""
    translation = fn(0)
    columns = [fn(1 << i) ^ translation ^ (1 << i) for i in range(dimension)]
    basis = {}
    for j in range(dimension):
        row = sum(((col >> j) & 1) << i for i, col in enumerate(columns))
        rhs = (translation >> j) & 1
        while row:
            pivot = row.bit_length() - 1
            if pivot in basis:
                oldrow, oldrhs = basis[pivot]
                row ^= oldrow
                rhs ^= oldrhs
            else:
                basis[pivot] = (row, rhs)
                break
        if not row and rhs:
            return 0
    return 1 << (dimension - len(basis))


def exhaustive_orbits(n, actions):
    seen = set()
    count = 0
    for m in range(1 << (n * n)):
        if m not in seen:
            orbit = {fn(m) for fn in actions}
            need(m in orbit, "check at source line 204")
            need(not (seen & orbit), "check at source line 205")
            seen.update(orbit)
            count += 1
    need(len(seen) == 1 << (n * n), "check at source line 208")
    return count


def run():
    report = {"schema": 1, "universe": "all binary n by n matrices; no filters"}
    # Independent adjacency round trips cover every labelled object through n=4.
    for n in range(5):
        for m in range(1 << (n * n)):
            need(decode(*encode(m, n), n) == m, "check at source line 217")
    report["exhaustive_bijection_n_max"] = 4
    # Independent direct adjacency relabelling covers all n<=3 objects/permutations.
    for n in range(4):
        for p in permutations(range(n)):
            for m in range(1 << (n * n)):
                need(alpha(m, p) == direct_pair_action(m, p), "check at source line 223")
    report["exhaustive_action_equivariance_n_max"] = 3
    # Affine maps are equal iff equal on zero and the coordinate basis.
    for n in range(5):
        ps = list(permutations(range(n)))
        probes = [0] + [1 << i for i in range(n * n)]
        for p in ps:
            for q in ps:
                pq = tuple(p[q[i]] for i in range(n))
                for m in probes:
                    need(alpha(alpha(m, q), p) == alpha(m, pq), "check at source line 233")
            for m in probes:
                need(alpha(m, p) == direct_pair_action(m, p), "check at source line 235")
                need(transpose(alpha(m, p), n) == alpha(transpose(m, n), p), "check at source line 236")
    report["all_affine_action_laws_n_max"] = 4
    # Completely independent Gaussian fixed-point solver on every cycle type to 12.
    type_count = 0
    for n in range(13):
        for lam in partitions(n):
            p = representative(lam)
            need(cycle_type(p) == lam, "check at source line 243")
            for natural, fn in ((True, alpha), (False, ordinary)):
                need(affine_fixed_count(lambda m: fn(m, p), n*n) == fixed_formula(lam, natural), "check at source line 245")
            need(affine_fixed_count(lambda m: alpha(transpose(m, n), p), n*n) == swapped_fixed_formula(lam), "check at source line 246")
            if all(a % 2 for a in lam):
                need(ordered_orbits(lam) == 2 * edge_orbits(lam) + len(lam), "check at source line 248")
            type_count += 1
    report["gaussian_cycle_type_checks"] = {"max_n": 12, "partitions": type_count, "actions_per_partition": 3}
    # Orbit enumeration uses direct action, not any partition or Gaussian solver.
    counts = []
    for n in range(5):
        ps = list(permutations(range(n)))
        natural = exhaustive_orbits(n, [lambda m, p=p: direct_pair_action(m, p) for p in ps])
        matrix = exhaustive_orbits(n, [lambda m, p=p: ordinary(m, p) for p in ps])
        swapped = exhaustive_orbits(n, [lambda m, p=p, swap=swap: direct_pair_action(transpose(m, n) if swap else m, p)
                                       for p in ps for swap in (False, True)])
        swapped_formula = (Fraction(natural) + sum((Fraction(swapped_fixed_formula(lam), z(lam))
                                                    for lam in partitions(n)), Fraction())) / 2
        need(natural == burnside(n, True), "check at source line 261")
        need(matrix == burnside(n, False), "check at source line 262")
        need(swapped == swapped_formula, "check at source line 263")
        counts.append({"n": n, "natural_aligned_pair": natural, "ordinary_matrix": matrix,
                       "natural_plus_core_exchange": swapped})
    need(counts[2]["natural_aligned_pair"] == 8, "check at source line 266")
    need(counts[2]["ordinary_matrix"] == 10, "check at source line 267")
    report["exhaustive_orbit_counts"] = counts
    ps = list(permutations(range(2)))
    need(exhaustive_orbits(2, [lambda m,p=p,q=q: rowcol(m,p,q) for p in ps for q in ps]) == 7, "check at source line 270")
    need(exhaustive_orbits(2, [lambda m,p=p,s=s: ordinary(transpose(m,2) if s else m,p)
                                for p in ps for s in (False,True)]) == 9, "check at source line 271")
    report["n2_additional_counts"] = {"independent_row_column": 7, "ordinary_plus_transpose": 9}
    # Exact equivariance locus: complementary opposite entries, arbitrary diagonal.
    for n in range(5):
        ps = list(permutations(range(n)))
        es = pairs(n)
        for packed in range(1 << (len(es) + n)):
            m = sum(((packed >> k) & 1) << (i*n+j) | (1 ^ ((packed >> k) & 1)) << (j*n+i)
                    for k, (i,j) in enumerate(es))
            m |= sum(((packed >> (len(es)+i)) & 1) << (i*n+i) for i in range(n))
            for p in ps:
                need(alpha(m, p) == ordinary(m, p), "check at source line 283")
    report["complementary_core_equivariance_exhaustive_n_max"] = 4
    report["partition_counts"] = [{"n": n, "natural_aligned_pair": burnside(n,True),
                                    "ordinary_matrix": burnside(n,False)} for n in range(11)]
    report["active_checks"] = CHECKS
    report["status"] = "PASS"
    return report


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
