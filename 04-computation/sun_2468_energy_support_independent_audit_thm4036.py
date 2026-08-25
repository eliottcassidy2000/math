#!/usr/bin/env python3
"""Independent hostile finite audit of the proposed Sun energy peeling lemma."""

from collections import Counter
from itertools import product
from math import comb


LOWER = {2: 2, 4: 3, 6: 5, 8: 7}
DEGREES = (2, 4, 6, 8)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def f(d, t):
    return comb(t, d)


def tau(n):
    out = 1
    p = 2
    while p * p <= n:
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        if exponent:
            out *= exponent + 1
        p += 1
    if n > 1:
        out *= 2
    return out


def atom_values(d, X):
    out = []
    t = LOWER[d]
    while f(d, t) <= X:
        out.append((t, f(d, t)))
        t += 1
    return tuple(out)


def tail_sum_counts(degrees, X):
    values = [tuple(value for _index, value in atom_values(d, X)) for d in degrees]
    counts = Counter()
    for row in product(*values):
        total = sum(row)
        if total <= X:
            counts[total] += 1
    return counts


def energy(counts):
    return sum(value * value for value in counts.values())


def audit_strictness_and_divisors(H=5000):
    print("STRICTNESS_AND_DIVISOR_BOUND")
    for d in DEGREES:
        L = LOWER[d]
        for q in range(1, 31):
            old = None
            for v in range(L, L + 200):
                value = f(d, v + q) - f(d, v)
                if old is not None:
                    require(value > old, f"g_q not strict at {(d, q, v)}")
                old = value
        boundary = tuple(f(d, L + q) - f(d, L) for q in range(1, 6))

        # Once the q=1 difference exceeds H, every later v and every q>=1
        # exceeds H. This makes the enumeration exhaustive for |h|<=H.
        vmax = L
        while f(d, vmax + 1) - f(d, vmax) <= H:
            vmax += 1
        reps = Counter()
        for v in range(L, vmax + 1):
            q = 1
            while True:
                h = f(d, v + q) - f(d, v)
                if h > H:
                    break
                reps[h] += 1
                q += 1
        worst = max((count / tau(comb(d, d) * 1) for count in (1,)), default=0)
        require(all(count <= tau(factorial(d) * h) for h, count in reps.items()),
                f"difference-divisor bound failed for degree {d}")
        max_slack_row = max(
            ((h, count, tau(factorial(d) * h)) for h, count in reps.items()),
            key=lambda row: row[1] / row[2],
        )
        print(
            f"d={d};L={L};zero_boundary_gq={boundary};strict=PASS;"
            f"R_bound_h1_{H}=PASS;tightest={max_slack_row}"
        )


def factorial(n):
    out = 1
    for k in range(2, n + 1):
        out *= k
    return out


def audit_peel(d, tail, X):
    outer = atom_values(d, X)
    tail_counts = tail_sum_counts(tail, X)
    full_counts = Counter()
    for g, multiplicity in tail_counts.items():
        for _u, value in outer:
            if value + g <= X:
                full_counts[value + g] += multiplicity
    direct_energy = energy(full_counts)

    exact_diagonal = 0
    exact_offdiagonal = 0
    ordered_tail_offdiag = 0
    for gs, cs in tail_counts.items():
        for gt, ct in tail_counts.items():
            multiplicity = cs * ct
            h = gt - gs
            if h:
                ordered_tail_offdiag += multiplicity
            outer_pairs = 0
            for _u, fu in outer:
                if fu + gs > X:
                    continue
                for _v, fv in outer:
                    if fu - fv == h:
                        outer_pairs += 1
            contribution = multiplicity * outer_pairs
            if h == 0:
                exact_diagonal += contribution
            else:
                exact_offdiagonal += contribution

    require(direct_energy == exact_diagonal + exact_offdiagonal,
            f"energy decomposition failed at {(d, tail, X)}")
    tail_B = sum(tail_counts.values())
    tail_E = energy(tail_counts)
    require(ordered_tail_offdiag == tail_B * tail_B - tail_E,
            f"ordered tail accounting failed at {(d, tail, X)}")

    max_outer = len(outer)
    max_R = max(
        (tau(factorial(d) * abs(gt - gs))
         for gs in tail_counts for gt in tail_counts if gs != gt),
        default=0,
    )
    require(exact_diagonal <= max_outer * tail_E,
            f"diagonal estimate failed at {(d, tail, X)}")
    require(exact_offdiagonal <= max_R * tail_B * tail_B,
            f"off-diagonal estimate failed at {(d, tail, X)}")
    print(
        f"peel={(d,) + tail};X={X};E={direct_energy};diag={exact_diagonal};"
        f"offdiag={exact_offdiagonal};tail_B={tail_B};tail_E={tail_E};"
        f"ordered_tail_offdiag={ordered_tail_offdiag};decomposition=PASS"
    )


def audit_fixed_AP(X=5000):
    counts = tail_sum_counts(DEGREES, X)
    global_E = energy(counts)
    print("FIXED_AP_CAUCHY_SCHWARZ")
    for q in range(1, 9):
        for r in range(q):
            fibres = tuple((n, a) for n, a in counts.items() if n % q == r and 1 <= n <= X)
            S = sum(a for _n, a in fibres)
            R = sum(a > 0 for _n, a in fibres)
            Eclass = sum(a * a for _n, a in fibres)
            require(S * S <= R * Eclass <= R * global_E,
                    f"fixed-AP Cauchy-Schwarz failed at {(q, r)}")
        print(f"q={q};all_classes_exact_CS=PASS")


def main():
    audit_strictness_and_divisors()
    for X in (50, 200, 1000):
        audit_peel(6, (8,), X)
        audit_peel(4, (6, 8), X)
        audit_peel(2, (4, 6, 8), X)
    audit_fixed_AP()
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
