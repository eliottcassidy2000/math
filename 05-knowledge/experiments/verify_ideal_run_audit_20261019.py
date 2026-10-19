#!/usr/bin/env python3
"""Independent exact regression for the ideal-run rank/moment reduction.

No imports from earlier packets; no numerical optimization or floating point.
A theorem is not inferred from the finite census: this tests its bijections,
identities, and implementation against actual small posets.
"""
from functools import lru_cache
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path


def closure(pred):
    pred = list(pred)
    for i in range(len(pred)):
        for j in range(i):
            if pred[i] >> j & 1:
                pred[i] |= pred[j]
    return tuple(pred)


def natural_posets(n):
    pairs = list(combinations(range(n), 2))
    seen = set()
    for edges in range(1 << len(pairs)):
        pred = [0] * n
        for k, (i, j) in enumerate(pairs):
            if edges >> k & 1:
                pred[j] |= 1 << i
        seen.add(closure(pred))
    return sorted(seen)


def linear_extensions(pred):
    n = len(pred)
    full = (1 << n) - 1
    out = []
    def rec(mask, word):
        if mask == full:
            out.append(tuple(word))
            return
        for i in range(n):
            if not mask >> i & 1 and pred[i] & mask == pred[i]:
                rec(mask | 1 << i, word + [i])
    rec(0, [])
    return out


def count_extensions(pred):
    n = len(pred)
    full = (1 << n) - 1
    @lru_cache(None)
    def rec(mask):
        if mask == full:
            return 1
        return sum(rec(mask | 1 << i) for i in range(n)
                   if not mask >> i & 1 and pred[i] & mask == pred[i])
    return rec(0)


def count_before(pred, i, j):
    pred = list(pred)
    pred[j] |= 1 << i
    return count_extensions(tuple(pred))


def initial_run(word, ideal):
    for k, i in enumerate(word):
        if not ideal >> i & 1:
            return k
    return len(word)


def small_collar_descriptors(max_n=3):
    descriptors = []
    for n in range(max_n + 1):
        for pred in natural_posets(n):
            words = linear_extensions(pred)
            e = len(words)
            ideals = [M for M in range(1 << n)
                      if all(pred[i] & M == pred[i] for i in range(n) if M >> i & 1)]
            for K in ideals:
                for H in ideals:
                    if H & K != H:
                        continue
                    r = [initial_run(w, K) for w in words]
                    d = Fraction(e + sum(r), e)
                    v = Fraction(sum((ss + 1) * (2 * rr - ss + 2) // 2
                                     for w, rr in zip(words, r)
                                     for ss in [initial_run(w, H)]), e)
                    descriptors.append((pred, H, K, e, d, v))
    return descriptors


def audit_both_sides():
    descriptors = small_collar_descriptors()
    total = above_seven = 0
    for pa, Ha, Ka, ea, c, u in descriptors:
        na = len(pa)
        allA = ((1 << na) - 1) << 5
        for pb, Hb, Kb, eb, d, v in descriptors:
            nb = len(pb)
            offB = 5 + na
            core = [allA & ~(Ka << 5), allA,
                    1 | (allA & ~(Ha << 5)), 1 | allA,
                    1 | 2 | 4 | allA]
            # A is the order dual of the descriptor poset.
            for i in range(na):
                core.append(sum(1 << (5 + j) for j in range(na) if pa[j] >> i & 1))
            for i, p in enumerate(pb):
                q = (p << offB) | allA | 1 | 2 | 8
                if not Hb >> i & 1:
                    q |= 4
                if not Kb >> i & 1:
                    q |= 16
                core.append(q)
            core = tuple(core)
            E = count_extensions(core)
            Nxy = count_before(core, 1, 3)
            Nxz = count_before(core, 1, 2)
            Nzy = count_before(core, 2, 3)
            factor = ea * eb
            E0 = 2*c*d+c+2*c*v+2*u*d+u+d+1+v
            X0 = c*d+c+c*v+u*d+u+d+1+v
            Z0 = c*d+c+2*c*v+d+1+v
            Y0 = c*d+c+2*u*d+u+d+1
            assert (E, Nxy, Nxz, Nzy) == tuple(factor * t for t in (E0, X0, Z0, Y0))
            assert any(E <= 3 * t <= 2 * E for t in (Nxy, Nxz, Nzy))
            C, D, U, V = c-1, d-1, u-1, v-1
            assert c <= u <= c*c and d <= v <= d*d
            assert 3*X0-2*E0 == 2-C-D-C*D-C*V-D*U
            assert E0-3*Z0 == -C*D-4*C*V-7*C+2*D*U-D+3*U-6*V-10
            assert 14*E0-11*(Z0+Y0) == 6*C*D+6*C*V+6*D*U+4*C+4*D+9*U+9*V
            if 3*Nxy > 2*E:
                assert C < 2 and D < 2
                upper = (3*C+5)*(C-2)+D*(2*C*C-C-7)
                assert E0-3*Z0 <= upper < 0
                assert E0-3*Y0 < 0
            total += 1
            # z has exactly two core incomparabilities plus |Ha|+|Hb|.
            if 2 + Ha.bit_count() + Hb.bit_count() > 7:
                above_seven += 1
    # The old positive mixture lies outside the proved realizability region.
    fake_d, fake_v = Fraction(19, 7), Fraction(59, 7)
    assert fake_v > fake_d*fake_d
    return {'collar_descriptors_per_side': len(descriptors),
            'actual_two_sided_posets': total,
            'actual_tests_with_z_degree_above_seven': above_seven,
            'old_mixture_moment_violation': str(fake_v-fake_d*fake_d)}


def audit():
    totals = {'posets': 0, 'ideal_rank_tests': 0, 'nested_ideal_core_tests': 0,
              'degree_seven_core_tests': 0}
    by_n = {}
    for n in range(6):
        ps = natural_posets(n)
        by_n[n] = len(ps)
        for pred in ps:
            totals['posets'] += 1
            words = linear_extensions(pred)
            e = len(words)
            full = (1 << n) - 1
            ideals = [M for M in range(1 << n)
                      if all(pred[i] & M == pred[i] for i in range(n) if M >> i & 1)]
            for K in ideals:
                runs = [initial_run(w, K) for w in words]
                tail = [sum(r >= j for r in runs) for j in range(n + 1)]
                # New marker t=n is minimal and below precisely B\K.
                augmented = tuple(p | ((1 << n) if not K >> i & 1 else 0)
                                  for i, p in enumerate(pred)) + (0,)
                rank = [0] * (n + 1)
                for w in linear_extensions(augmented):
                    rank[w.index(n)] += 1
                assert tail == rank
                assert all(tail[j] ** 2 >= tail[j - 1] * tail[j + 1]
                           for j in range(1, n))
                for i in range(n + 1):
                    for j in range(n + 1):
                        qsum = tail[i + j] if i + j <= n else 0
                        assert e * qsum <= tail[i] * tail[j]
                msum = sum(runs)
                second_sum = sum(r * (r - 1) // 2 for r in runs)
                assert e * second_sum <= msum ** 2
                d = Fraction(e + msum, e)
                totals['ideal_rank_tests'] += 1
                for H in ideals:
                    if H & K != H:
                        continue
                    vsum = sum((s + 1) * (2 * r - s + 2) // 2
                               for w, r in zip(words, runs)
                               for s in [initial_run(w, H)])
                    v = Fraction(vsum, e)
                    assert d <= v <= d ** 2
                    # Construct the actual rooted-five poset, independently.
                    # Core b=0,x=1,z=2,y=3,a=4; B starts at label 5.
                    core = [0, 0, 1, 1, 1 | 2 | 4]
                    for i, p in enumerate(pred):
                        q = (p << 5) | 1 | 2 | 8
                        if not H >> i & 1:
                            q |= 4
                        if not K >> i & 1:
                            q |= 16
                        core.append(q)
                    core = tuple(core)
                    E = count_extensions(core)
                    Nxy = count_before(core, 1, 3)
                    Nzy = count_before(core, 2, 3)
                    assert E == e * (5 * d + 3 * v + 3)
                    assert Nxy == e * (3 * d + 2 * v + 3)
                    assert Nzy == e * (4 * d + 3)
                    if d >= 3:
                        assert E <= 3 * Nxy <= 2 * E
                    else:
                        assert E <= 3 * Nzy <= 2 * E
                    totals['nested_ideal_core_tests'] += 1
                    # Degrees are counted with full transitive closure.
                    cp = list(core)
                    for _ in range(len(cp)):
                        for i in range(len(cp)):
                            for j in range(len(cp)):
                                if cp[i] >> j & 1:
                                    cp[i] |= cp[j]
                    degrees = [sum(i != j and not(cp[i] >> j & 1)
                                   and not(cp[j] >> i & 1) for j in range(len(cp)))
                               for i in range(len(cp))]
                    if max(degrees, default=0) <= 7:
                        totals['degree_seven_core_tests'] += 1
    return {'status': 'PASS', 'method': 'exact integers and fractions',
            'scope': 'all naturally labeled B of size at most five; all nested ideals H subset K',
            'posets_by_size': by_n, **totals, 'both_sided_tests': audit_both_sides()}

if __name__ == '__main__':
    result = audit()
    path = Path(__file__).with_name('ideal_run_independent_audit_20261019.json')
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
