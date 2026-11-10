#!/usr/bin/env python3
"""Self-contained exact physical graft / smooth-center certificate. Stdlib only."""
from itertools import permutations
from math import gcd, isqrt, lcm
import json
import sys


def need(test, message):
    if not test:
        raise RuntimeError(message)


def validate(a):
    n = len(a)
    need(all(0 <= x < 1 << n and not (x >> i & 1) for i, x in enumerate(a)), 'invalid row')
    need(all(((a[i] >> j) ^ (a[j] >> i)) & 1 for i in range(n) for j in range(i + 1, n)), 'not tournament')


def tables(a):
    validate(a)
    n = len(a)
    incoming = [sum(1 << u for u in range(n) if a[u] >> v & 1) for v in range(n)]
    end = [[0] * n for _ in range(1 << n)]
    h = [0] * (1 << n)
    h[0] = 1
    for s in range(1, 1 << n):
        vs = s
        while vs:
            bit = vs & -vs
            v = bit.bit_length() - 1
            t = s ^ bit
            if t == 0:
                end[s][v] = 1
            else:
                choices = t & incoming[v]
                total = 0
                while choices:
                    bitu = choices & -choices
                    total += end[t][bitu.bit_length() - 1]
                    choices ^= bitu
                end[s][v] = total
            vs ^= bit
        h[s] = sum(end[s])
    return h, end


def reverse(a):
    full = (1 << len(a)) - 1
    return [full ^ x ^ (1 << i) for i, x in enumerate(a)]


def pair(a):
    h, _ = tables(a)
    full = (1 << len(a)) - 1
    return h[full], sum(h[s] * h[full ^ s] for s in range(full + 1))


def state(a, root):
    h, end = tables(a)
    _, start = tables(reverse(a))
    full = (1 << len(a)) - 1
    return (h[full], end[full][root], sum(start[s][root] * h[full ^ s] for s in range(1, full + 1)), h[full ^ (1 << root)])


def ab(z):
    D, E, U, d = z
    return D + 2 * d - E - U, U - d


def tt(n):
    return [sum(1 << j for j in range(i + 1, n)) for i in range(n)]


def endpoint(k):
    a = tt(k + 2)
    a[0] ^= 1 << (k + 1)
    a[k + 1] |= 1
    return a


def join(a, b):
    n, m = len(a), len(b)
    return [x | (((1 << m) - 1) << n) for x in a] + [x << n for x in b]


def graft(a, root, b):
    n, m = len(a), len(b)
    new = [x if i == root else x | (((1 << m) - 1) << n) for i, x in enumerate(a)]
    return new + [(x << n) | (1 << root) for x in b]


def strong(a):
    full = (1 << len(a)) - 1
    for r in range(len(a)):
        seen = 1 << r
        while True:
            new = seen
            for u in range(len(a)):
                if seen >> u & 1:
                    new |= a[u]
            if new == seen:
                break
            seen = new
        if seen != full:
            return False
    return True


def brute_tables(a):
    n = len(a)
    h = [0] * (1 << n)
    starts = [[0] * n for _ in h]
    ends = [[0] * n for _ in h]
    h[0] = 1
    for s in range(1, 1 << n):
        vertices = [i for i in range(n) if s >> i & 1]
        for p in permutations(vertices):
            if all(a[u] >> v & 1 for u, v in zip(p, p[1:])):
                h[s] += 1
                starts[s][p[0]] += 1
                ends[s][p[-1]] += 1
    return h, starts, ends


def brute_state(a, root):
    h, starts, ends = brute_tables(a)
    full = (1 << len(a)) - 1
    return h[full], ends[full][root], sum(starts[s][root] * h[full ^ s] for s in range(1, full + 1)), h[full ^ (1 << root)]


def labelled(n, code):
    a = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if code & 1:
                a[i] |= 1 << j
            else:
                a[j] |= 1 << i
            code >>= 1
    return a


def prime(x):
    if x < 2:
        return False
    if x % 2 == 0:
        return x == 2
    return all(x % d for d in range(3, isqrt(x) + 1, 2))


def factor(x):
    factors = []
    d = 2
    while d * d <= x:
        while x % d == 0:
            factors.append(d)
            x //= d
        d = 3 if d == 2 else d + 2
    if x > 1:
        factors.append(x)
    return factors


def order2(p):
    d = 1
    while pow(2, d, p) != 1:
        d += 1
    return d


def zero_classes(A, B, p):
    d = order2(p)
    if B % p == 0:
        return d, list(range(d)) if A % p == 0 else []
    if A % p == 0:
        return d, []
    target = -A * pow(B, -1, p) % p
    return d, [r for r in range(d) if pow(2, r, p) == target]


def word(w):
    a, root = [2, 0], 0
    for letter in w:
        n = len(a)
        if letter == 'O':
            a = join(a, [0])
        else:
            a = graft(a, root, [0])
            if letter == 'J':
                root = n
    return a, root


def main():
    if '--fail' in sys.argv:
        need(False, 'intentional failure gate')
    A, B = [8, 1, 3, 6], [20, 1, 2, 7, 14]
    sources = [(A, 3, -1, (5, 1, 7, 1)), (A, 2, 1, (5, 1, 9, 3)),
               (B, 1, 3, (9, 3, 9, 3)), (B, 2, 5, (9, 3, 11, 5))]
    source_records = []
    for rows, root, alpha, expected in sources:
        need(state(rows, root) == brute_state(rows, root) == expected, 'source state')
        need(ab(expected) == (alpha, 6) and strong(rows), 'source interface/strongness')
        source_records.append(dict(rows=rows, root=root, state=expected, a=alpha, b=6))

    # Complete labelled census through order five; merge only the full physical state.
    source_bank, gadget_bank = {}, {}
    tournaments = rooted = 0
    for n in range(1, 6):
        for code in range(1 << (n * (n - 1) // 2)):
            rows = labelled(n, code)
            h, end = tables(rows)
            _, start = tables(reverse(rows))
            full = (1 << n) - 1
            tournaments += 1
            if n <= 3:
                hq = h[full], sum(h[s] * h[full ^ s] for s in range(full + 1))
                gadget_bank.setdefault(hq, rows)
            for root in range(n):
                if rows[root] == 0:
                    continue
                z = (h[full], end[full][root], sum(start[s][root] * h[full ^ s] for s in range(1, full + 1)), h[full ^ (1 << root)])
                aa, bb = ab(z)
                need(aa % 2 == 1 and bb > 0, 'physical parity')
                source_bank.setdefault(z, (rows, root))
                rooted += 1
    graft_checks = 0
    for z, (rows, root) in sorted(source_bank.items()):
        D, E, U, d = z
        aa, bb = ab(z)
        for (H, Q), gadget in sorted(gadget_bank.items()):
            out = graft(rows, root, gadget)
            expected = H * aa + Q * bb, H * d, Q * bb + H * d, H * d
            actual = state(out, root)
            need(actual == expected and ab(actual) == (H * aa, Q * bb), 'full graft law')
            graft_checks += 1

    windows = []
    graph_checks = 0
    for m in range(1, 9):
        counts = []
        for rows, root, alpha, z in sources:
            out = graft(rows, root, tt(m))
            value = tables(out)[0][-1]
            need(value == 6 * 2**m + alpha and strong(out), 'dyadic window')
            counts.append(value)
            graph_checks += 1
        windows.append(dict(m=m, values=counts))
    arbitrary = []
    for gadget in ([0], [2, 4, 1], endpoint(4), join([2, 4, 1], [0])):
        H, Q = pair(gadget)
        vals = []
        for rows, root, alpha, z in sources:
            out = graft(rows, root, gadget)
            value = tables(out)[0][-1]
            need(value == 6 * Q + alpha * H and strong(out), 'arbitrary gadget window')
            vals.append(value)
            graph_checks += 1
        arbitrary.append(dict(gadget_rows=gadget, H=H, Q=Q, values=vals, offset_center_gcd=gcd(H, 6*Q)))
    # The iterated physical construction is literally the single ordinal-gadget graft.
    iterated_checks = 0
    for rows, root, alpha, z in sources:
        for left in ([0], [2, 4, 1], tt(2)):
            for right in ([0], [2, 4, 1], tt(2)):
                out = graft(graft(rows, root, left), root, right)
                flat = graft(rows, root, join(left, right))
                need(out == flat, 'graft word flattening')
                h1, q1 = pair(left)
                h2, q2 = pair(right)
                need(pair(join(left, right)) == (h1*h2, q1*q2), 'join pair')
                need(ab(state(out, root)) == (alpha*h1*h2, 6*q1*q2), 'word interface')
                iterated_checks += 1

    # Every fixed-source pair in this census, under every gadget pair, for t=0..5.
    # Three equal differences suffice algebraically; the proof supplies all orders.
    profiles = sorted(set(ab(z) for z in source_bank))
    gap_tests = 0
    for i, (a1, b1) in enumerate(profiles):
        for a2, b2 in profiles[i+1:]:
            for H, Q in gadget_bank:
                differences = [H*(a2-a1) + Q*2**t*(b2-b1) for t in range(6)]
                constant = len(set(differences)) == 1
                need(constant == (b1 == b2), 'constant gap criterion')
                if constant and abs(differences[0]) == 2:
                    need(H == 1 and abs(a2-a1) == 2, 'twin criterion')
                gap_tests += 1

    primes = [p for p in range(3, 32, 2) if prime(p)]
    sieve_checks = 0
    for aa in range(-11, 12, 2):
        for bb in range(2, 18, 2):
            for p in primes:
                d, zeros = zero_classes(aa, bb, p)
                for t in range(60):
                    need(((aa + bb*2**t) % p == 0) == (t % d in zeros), 'general residue sieve')
                    sieve_checks += 1
    residue_records = []
    full_period = 1
    for p in primes:
        d = order2(p)
        low = zero_classes(-1, 6, p)[1]
        high = zero_classes(1, 6, p)[1]
        residue_records.append(dict(p=p, order=d, lower_zero_residues=low, upper_zero_residues=high))
        full_period = lcm(full_period, d)
    survivors = []
    for m in range(full_period):
        modular = all((6*pow(2,m,p)-1) % p and (6*pow(2,m,p)+1) % p for p in primes)
        simple = m % 6 == 5 and m % 5 != 1 and m % 11 != 2
        need(bool(modular) == simple, 'reduced wheel')
        if modular:
            survivors.append(m)
    reduced = [m for m in range(330) if m % 6 == 5 and m % 5 != 1 and m % 11 != 2]
    need(len(reduced) == 40 and len(survivors) == 3360, 'sieve counts')
    prime_pairs = {}
    for gap, offsets in [(2, (-1, 1)), (4, (1, 5)), (6, (-1, 5))]:
        records = []
        for m in range(1, 31):
            low, high = [6*2**m + x for x in offsets]
            if prime(low) and prime(high):
                consecutive = not any(prime(x) for x in range(low+2, high, 2))
                records.append(dict(m=m, lower=low, upper=high, consecutive=consecutive))
        prime_pairs[str(gap)] = records
    hostile = [6*2**23-1, 6*2**23+1]
    need(23 in reduced and all(not prime(x) for x in hostile), 'sieve not primality')

    # Brown roots are independently supplied, never extracted from square counts.
    root5 = endpoint(2)
    root11 = graft(A, 3, [0])
    root71, root = word('OOZJOZ')
    brown = []
    for h, rows in [(5, root5), (11, root11), (71, root71)]:
        need(tables(rows)[0][-1] == h, 'Brown root count')
        brute = sum(all(rows[u] >> v & 1 for u, v in zip(p,p[1:])) for p in permutations(range(len(rows))))
        need(brute == h, 'Brown permutation root count')
        square = join(rows, rows)
        need(tables(square)[0][-1] == h*h, 'Brown square literal count')
        brown.append(dict(H=h, rows=rows, square_H=h*h, square_order=len(square), square_rows=square))

    # Nonunit intercepts may produce Euclid values; this guards the exact scope.
    z31=state(B,3)
    graph31=graft(B,3,[0])
    pv31=tables(graph31)[0][-1]
    need(ab(z31)==(-5,18) and pv31==31 and strong(graph31), 'nonunit Euclid hostile')

    # A bounded consequence control for an all-order primorial-disjointness proof.
    euclid = []
    P = 1
    for p in [2,3,5,7,11,13,17,19,23,29,31,37]:
        P *= p
        euclid.append(P+1)
    need(not set(euclid) & {6*2**m+a for m in range(1,101) for a in (-1,1,3,5)}, 'primorial hostile')
    print(json.dumps(dict(
        status='PASS; all-order statements proved in PROOF.md; computations are finite-exact',
        source_bank=source_records,
        exhaustive_controls=dict(labelled_tournaments_orders_1_to_5=tournaments, nonsink_rooted_instances=rooted,
            full_source_states=len(source_bank), ab_profiles=len(profiles), gadget_pairs_orders_1_to_3=len(gadget_bank),
            full_state_grafts=graft_checks, displayed_window_graphs=graph_checks, iterated_grafts=iterated_checks,
            constant_gap_cases=gap_tests, general_sieve_cases=sieve_checks),
        transitive_windows=windows, arbitrary_gadget_windows=arbitrary,
        prime_sieve=dict(primes=primes, residues=residue_records, order_lcm=full_period,
            survivors_in_lcm=len(survivors), reduced_period=330, reduced_residues=reduced, density='4/33'),
        prime_pairs_m_1_to_30=prime_pairs,
        false_positive_m23=[dict(value=x, prime_factors=factor(x)) for x in hostile],
        Brown_root_and_square_graphs=brown,
        primorial_control=dict(euclid_values=euclid, tested_m_range=[1,100], intersections=[],
            nonunit_intercept_hostile=dict(source_rows=B,root=3,state=z31,a=-5,b=18,gadget_rows=[0],output_rows=graph31,H=pv31)),
        limitations=['No prime infinitude or Polignac conjecture is proved.',
          'The four-window does not intersect positive-index prime-product-plus-one values.',
          'No universal +2 operation or arbitrary-target placement is supplied.',
          'Fixed-seed unit-intercept obstruction is construction-relative.']
    ), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
