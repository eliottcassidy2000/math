#!/usr/bin/env python3
"""Exact scratch probe for cyclic finite-difference fibre obstructions.

The full W=0 Hom lattice uses the THM-4241 O-basis [u,f,g,h].  Here O is
Z[omega], omega^2+omega+1=0.  Precomposition by the attachment generator tau
acts by

    T(u)=-omega*u, T(f)=g, T(g)=-omega*f,
    T(h)=omega^2*h-omega*f.

If m takes one value on the twelve-point tau orbit, then (T^k-1)m vanishes
on all twelve points.  Every nonzero such difference must therefore have
curve degree at least twelve.  This script enumerates the exact degree-34 and
degree-42 shells and records their minimum positive cyclic-difference degree.
"""

from collections import Counter, defaultdict


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def sub(x, y):
    return x[0] - y[0], x[1] - y[1]


def neg(x):
    return -x[0], -x[1]


def mul(x, y):
    a, b = x
    c, d = y
    return a * c - b * d, a * d + b * c - b * d


def conj(x):
    a, b = x
    return a - b, -b


def norm(x):
    return mul(x, conj(x))[0]


ZERO = (0, 0)
ONE = (1, 0)
OMEGA = (0, 1)
OMEGA2 = (-1, -1)
MINUS_OMEGA = (0, -1)

HIDDEN_H = (
    ((6, 0), (-4, -2)),
    ((-2, 2), (6, 0)),
)


def hidden_degree(x, y):
    coeffs = (x, y)
    total = ZERO
    for i in range(2):
        for j in range(2):
            total = add(total, mul(mul(coeffs[i], HIDDEN_H[i][j]), conj(coeffs[j])))
    assert total[1] == 0
    return total[0]


def projections(m):
    """Return visible coefficient data and hidden projection coefficients."""
    a, b, c, d = m
    ell_a = add(mul((2, 0), b), mul(OMEGA2, d))
    ell_b = add(mul((2, 0), c), d)
    return a, d, ell_a, ell_b


def degree(m):
    a, d, ell_a, ell_b = projections(m)
    qell = hidden_degree(ell_a, ell_b)
    assert qell % 4 == 0
    return 4 * norm(a) + norm(d) + qell // 4


def t_action(m):
    a, b, c, d = m
    return (
        mul(MINUS_OMEGA, a),
        mul(MINUS_OMEGA, add(c, d)),
        b,
        mul(OMEGA2, d),
    )


def t_power(m, k):
    ans = m
    for _ in range(k):
        ans = t_action(ans)
    return ans


def vector_sub(x, y):
    return tuple(sub(a, b) for a, b in zip(x, y))


def enumerate_shells():
    eis = []
    for x0 in range(-8, 9):
        for x1 in range(-8, 9):
            x = (x0, x1)
            if norm(x) <= 42:
                eis.append(x)
    hidden_by_parity = defaultdict(list)
    for la0 in range(-12, 13):
        for la1 in range(-12, 13):
            ell_a = (la0, la1)
            for lb0 in range(-12, 13):
                for lb1 in range(-12, 13):
                    ell_b = (lb0, lb1)
                    qell = hidden_degree(ell_a, ell_b)
                    if qell <= 168:
                        key = (la0 % 2, la1 % 2, lb0 % 2, lb1 % 2)
                        hidden_by_parity[key].append((ell_a, ell_b, qell))

    shells = {34: [], 42: []}
    # Reuse the projection enumeration: d and hidden projection ell determine
    # b,c uniquely when the parity class matches.
    for d in eis:
        omega2d = mul(OMEGA2, d)
        parity = (omega2d[0] % 2, omega2d[1] % 2, d[0] % 2, d[1] % 2)
        for ell_a, ell_b, qell in hidden_by_parity[parity]:
            delta_a = sub(ell_a, omega2d)
            delta_b = sub(ell_b, d)
            b = (delta_a[0] // 2, delta_a[1] // 2)
            c = (delta_b[0] // 2, delta_b[1] // 2)
            base = norm(d) + qell // 4
            if qell % 4 or base > 42:
                continue
            for a in eis:
                q = base + 4 * norm(a)
                if q in shells:
                    m = (a, b, c, d)
                    assert degree(m) == q
                    assert degree(t_action(m)) == q
                    shells[q].append(m)
    assert len(shells[34]) == 36288
    assert len(shells[42]) == 16992
    return shells


def main():
    shells = enumerate_shells()
    print("W=0 CYCLIC DIFFERENCE PROBE")
    for q, vectors in shells.items():
        min_hist = Counter()
        first_hist = Counter()
        sub12 = Counter()
        sub12_projection = Counter()
        sub12_vectors = []
        zero_counts = Counter()
        hidden_degree12_difference = Counter()
        hidden_degree12_vectors = set()
        survivors = []
        for m in vectors:
            diffs = []
            for k in range(1, 12):
                dq = degree(vector_sub(t_power(m, k), m))
                if dq == 0:
                    zero_counts[k] += 1
                else:
                    diffs.append((dq, k))
                    if dq < 12:
                        sub12[(k, dq)] += 1
                        a0, d0, la0, lb0 = projections(m)
                        qell0 = hidden_degree(la0, lb0)
                        qvis0 = 4 * q - qell0
                        sub12_projection[(k, dq, qvis0, qell0)] += 1
                        if k == 3:
                            sub12_vectors.append(m)
                    # T^3 and T^9 fix the visible v-eigenline.  At degree 12
                    # their visible u contribution must vanish, so the
                    # difference is in the hidden lattice and THM-4247 applies.
                    if k in (3, 9) and dq == 12:
                        assert m[0] == ZERO
                        _, _, la0, lb0 = projections(m)
                        qell0 = hidden_degree(la0, lb0)
                        hidden_degree12_difference[(k, qell0)] += 1
                        hidden_degree12_vectors.add(m)
            first_hist[degree(vector_sub(t_action(m), m))] += 1
            mind = min(dq for dq, _ in diffs)
            min_hist[mind] += 1
            if mind >= 12:
                survivors.append(m)
        print(f"full_degree={q} count={len(vectors)}")
        print(f"first_difference_hist={dict(sorted(first_hist.items()))}")
        print(f"minimum_positive_difference_hist={dict(sorted(min_hist.items()))}")
        print(f"sub12_by_power_degree={dict(sorted(sub12.items()))}")
        print(f"sub12_projection_hist={dict(sorted(sub12_projection.items()))}")
        print(f"zero_difference_by_power={dict(sorted(zero_counts.items()))}")
        print(f"hidden_degree12_difference_hist={dict(sorted(hidden_degree12_difference.items()))}")
        print(f"hidden_degree12_difference_union={len(hidden_degree12_vectors)}")
        print(f"cyclic_fibre_bound_survivors={len(survivors)}")
        if survivors:
            print(f"first_survivor={survivors[0]}")
        if sub12_vectors:
            print(f"first_sub12_vector={sub12_vectors[0]}")


if __name__ == "__main__":
    main()
