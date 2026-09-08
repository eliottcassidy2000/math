"""Exact controls for the pure-shear nuclear-moment theorem.

Universe: rational polar certificates built from 2x2 neutral blocks,
weighted directed cycles in dimensions 3..6, and rational orthogonal
conjugations. These controls do not prove the all-real theorem: its proof
uses polar decomposition and the elementary hollow-basis induction.
Run with Python 3, normally and with -O. No assertions can be disabled.
"""
from fractions import Fraction as F
from itertools import product


CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def zero(n):
    return [[F(0) for _ in range(n)] for _ in range(n)]


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def transpose(a):
    return list(map(list, zip(*a)))


def add(a, b):
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(c, a):
    return [[c * x for x in row] for row in a]


def mul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def inner(a, b):
    return sum(x * y for ar, br in zip(a, b) for x, y in zip(ar, br))


def outer(x, y):
    return [[a * b for b in y] for a in x]


def rank(a):
    a = [[F(x) for x in row] for row in a]
    i = 0
    for j in range(len(a[0])):
        k = next((k for k in range(i, len(a)) if a[k][j]), None)
        if k is None:
            continue
        a[i], a[k] = a[k], a[i]
        pivot = a[i][j]
        a[i] = [x / pivot for x in a[i]]
        for k in range(i + 1, len(a)):
            c = a[k][j]
            a[k] = [x - c * y for x, y in zip(a[k], a[i])]
        i += 1
        if i == len(a):
            break
    return i


def rotation(n, axis, c=F(3, 5), s=F(4, 5)):
    a = eye(n)
    i, j = axis
    a[i][i] = a[j][j] = c
    a[i][j], a[j][i] = -s, s
    return a


def conjugate(u, a):
    return mul(mul(u, a), transpose(u))


def audit(m, q, h, summands, lengths):
    """h is independently known PSD from the explicit construction."""
    n = len(m)
    check(mul(transpose(q), q) == eye(n), "orthogonal polar factor")
    check(h == transpose(h), "symmetric positive polar factor")
    check(mul(q, h) == m, "polar reconstruction")
    check(mul(transpose(m), m) == mul(h, h), "squared singular-value identity")
    check(trace(m) == 0, "trace-free mean")
    check(sum(lengths) == trace(h), "nuclear norm by positive polar factor")
    check(rank(m) == len(summands), "minimal nonzero atom count")
    total = zero(n)
    for s, length in zip(summands, lengths):
        check(rank(s) == 1 and trace(s) == 0, "pure shear")
        check(mul(s, s) == zero(n), "nilpotence independent check")
        check(inner(s, s) == length * length and length > 0, "Frobenius norm")
        total = add(total, s)
    check(total == m, "sum reconstruction")
    c = trace(h)
    mean = zero(n)
    second = F(0)
    variance = F(0)
    for s, length in zip(summands, lengths):
        weight = length / c
        atom = scale(c / length, s)
        mean = add(mean, scale(weight, atom))
        second += weight * inner(atom, atom)
        delta = add(atom, scale(-1, m))
        variance += weight * inner(delta, delta)
        check(inner(atom, atom) == c * c, "equal-radius equality condition")
    check(mean == m and second == c * c, "optimal moment certificate")
    check(variance == c * c - inner(m, m), "variance identity")


def curl(a):
    return [a[2][1] - a[1][2], a[0][2] - a[2][0], a[1][0] - a[0][1]]


def stretching(a):
    w = curl(a)
    norm = sum(x * x for x in w)
    if not norm:
        raise ValueError("alpha undefined at zero vorticity")
    return sum(w[i] * a[i][j] * w[j] for i in range(3) for j in range(3)) / norm


def main():
    cycle_cases = 0
    for n in range(3, 7):
        # 81 positive weighted cycles for each n; repeated weights do not
        # affect invertibility. Their diagonal H is manifestly positive.
        for first in product((F(1, 2), F(1), F(2)), repeat=4):
            weights = [first[i % 4] for i in range(n)]
            q, h, summands = zero(n), zero(n), []
            for j, w in enumerate(weights):
                q[(j + 1) % n][j] = 1
                h[j][j] = w
                s = zero(n)
                s[(j + 1) % n][j] = w
                summands.append(s)
            m = mul(q, h)
            u = mul(rotation(n, (0, 1)), rotation(n, (1, 2), F(5, 13), F(12, 13)))
            audit(conjugate(u, m), conjugate(u, q), conjugate(u, h),
                  [conjugate(u, s) for s in summands], weights)
            cycle_cases += 1

    neutral_cases = 0
    for a, b in product((F(1, 3), F(1, 2), F(1), F(2), F(3)), repeat=2):
        # Noncommuting Q,H controls. H=sum p p^T is positive definite,
        # with eigenvalues 2a^2,2b^2. Q reflects the second coordinate.
        for n in (2, 3, 4):
            q = eye(n)
            q[1][1] = -1
            p1, p2 = [a, a] + [F(0)] * (n - 2), [b, -b] + [F(0)] * (n - 2)
            h = add(outer(p1, p1), outer(p2, p2))
            summands = []
            for p in (p1, p2):
                qp = [sum(q[i][j] * p[j] for j in range(n)) for i in range(n)]
                summands.append(outer(qp, p))
            u = rotation(n, (0, n - 1))
            audit(conjugate(u, mul(q, h)), conjugate(u, q), conjugate(u, h),
                  [conjugate(u, s) for s in summands], [2 * a * a, 2 * b * b])
            neutral_cases += 1

    # Hostile: average zero-stretch shears gives positive mean stretching.
    cycle = [[F(0), F(1), F(0)], [F(0), F(0), F(1)], [F(1), F(0), F(0)]]
    atoms = []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        atom = zero(3)
        atom[i][j] = 3
        atoms.append(atom)
        check(stretching(atom) == 0, "atom self-stretching zero")
    check(stretching(cycle) == 1, "strict positive mean stretching")
    check(mul(transpose(cycle), cycle) == eye(3), "cycle nuclear norm three")
    # The distance formula gives d^2=3-1=2: |Mp|=1 for all unit p,
    # and p=e1 makes p^T M p=0. This checks the two exact certificates.
    check(cycle[0][0] == 0 and inner(cycle, cycle) == 3, "distance squared two")
    check(sum(inner(x, x) for x in atoms) / 3 == 9, "cycle second moment nine")

    # A minimum-cardinality compositional hostile from the catalysis lane.
    s = zero(3)
    s[0][1] = 1
    t = outer([F(1)] * 3, [F(1), F(0), F(-1)])
    check(trace(t) == 0 and rank(t) == 1, "two-atom hostile is legal")
    check(stretching(s) == stretching(t) == 0, "two atoms no self stretch")
    check(stretching(add(s, t)) == F(-3, 5), "two-atom composition detects stretch")

    # Boundaries: trace-nonzero means are impossible; PSD tracezero means
    # are zero; zero atom handles M=0 without dividing by the nuclear norm.
    check(trace(eye(3)) != 0, "trace constraint cannot be discarded")
    check(rank(zero(3)) == 0 and inner(zero(3), zero(3)) == 0, "zero boundary")
    print("FINITE-EXACT rational controls; no PDE realization claim")
    print(f"positive weighted-cycle polar certificates: {cycle_cases}")
    print(f"neutral-block polar certificates including rank-deficient cases: {neutral_cases}")
    print("cycle mean: alpha=1, shear distance squared=2, nuclear norm=3")
    print("cycle optimal second moment=9, optimal variance=6")
    print("two-shear sum: alpha=-3/5 (each atom has alpha=0)")
    print(f"active exact gates: {CHECKS}")


if __name__ == "__main__":
    main()
