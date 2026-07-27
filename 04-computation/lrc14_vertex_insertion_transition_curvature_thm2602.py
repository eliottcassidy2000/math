#!/usr/bin/env python3
"""Exact finite controls for THM-2602.

This companion distinguishes three objects:

* diagonal vertex insertions and terminal Fourier readouts;
* joint edge correspondences with identical one-vertex marginals;
* ordered, composable transition kernels and their twisted return.

Only integer arithmetic and arithmetic modulo the split prime 547 are used.
"""

from itertools import product


P = 13
N = 7


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def zero(n):
    return [[0 for _ in range(n)] for _ in range(n)]


def matmul(a, b, modulus=None):
    n = len(a)
    c = zero(n)
    for i in range(n):
        for k in range(n):
            if a[i][k] == 0:
                continue
            for j in range(n):
                c[i][j] += a[i][k] * b[k][j]
        if modulus is not None:
            c[i] = [x % modulus for x in c[i]]
    return c


def matsub(a, b, modulus=None):
    c = [[x - y for x, y in zip(ar, br)] for ar, br in zip(a, b)]
    if modulus is not None:
        c = [[x % modulus for x in row] for row in c]
    return c


def matpow(a, exponent, modulus=None):
    r = eye(len(a))
    b = a
    e = exponent
    while e:
        if e & 1:
            r = matmul(r, b, modulus)
        b = matmul(b, b, modulus)
        e //= 2
    return r


def shift(step):
    """Row=source, column=target: q maps to q+step."""
    return [[int(qp == (q + step) % P) for qp in range(P)]
            for q in range(P)]


def diagonal(values):
    return [[values[i] if i == j else 0 for j in range(P)]
            for i in range(P)]


def row_sums(a):
    return [sum(row) for row in a]


def col_sums(a):
    return [sum(a[i][j] for i in range(len(a)))
            for j in range(len(a))]


def exact_order(x, order, modulus):
    if pow(x, order, modulus) != 1:
        return False
    primes = []
    n = order
    d = 2
    while d * d <= n:
        if n % d == 0:
            primes.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        primes.append(n)
    return all(pow(x, order // p, modulus) != 1 for p in primes)


def selector_bank():
    z0 = {(7, 4), (7, 5), (7, 6)}
    z6 = {(6, 4), (6, 5), (6, 6)}
    selector_count = 0
    marker_checks = 0
    last_choice_checks = 0
    for s in range(1, P):
        choices = []
        for ell in range(N):
            choices.append(tuple(v for v, zeros in ((0, z0), (6, z6))
                                 if (s, ell) not in zeros))
        last_choice = tuple(values[-1] for values in choices)
        require(sum(last_choice[ell] - last_choice[(ell - 1) % N]
                    for ell in range(N)) % P == 0,
                "last-choice rail selector failed to telescope")
        last_choice_checks += 1
        for h in product(*choices):
            selector_count += 1
            coboundary = sum(h[ell] - h[(ell - 1) % N]
                              for ell in range(N)) % P
            require(coboundary == 0, "selector coboundary did not telescope")
            for a in range(1, P):
                require((N * a + coboundary) % P == (N * a) % P != 0,
                        "vertex selector changed nonzero holonomy")
                marker_checks += 1
    require(selector_count == 1312, "wrong selector count")
    require(marker_checks == 15744, "wrong marker/selector count")
    return selector_count, marker_checks, last_choice_checks


def diagonal_algebra_checks():
    basis = [diagonal([int(q == i) for q in range(P)]) for i in range(P)]
    commutators = 0
    for x in basis:
        for y in basis:
            require(matmul(x, y) == matmul(y, x),
                    "diagonal insertions did not commute")
            commutators += 1

    # A simultaneous exact Fourier conjugation still gives a commutative
    # algebra.  It is a representation change, not an inserted physical edge.
    modulus = 547
    omega = next(x for x in range(2, modulus)
                 if exact_order(x, P, modulus))
    fourier = [[pow(omega, i * j, modulus) for j in range(P)]
               for i in range(P)]
    inv_p = pow(P, -1, modulus)
    inverse = [[inv_p * pow(omega, -i * j, modulus) % modulus
                for j in range(P)] for i in range(P)]
    require(matmul(fourier, inverse, modulus) == eye(P),
            "Fourier inverse failed")
    conjugates = [matmul(matmul(fourier, d, modulus), inverse, modulus)
                  for d in basis]
    conjugate_commutators = 0
    for x in conjugates:
        for y in conjugates:
            require(matmul(x, y, modulus) == matmul(y, x, modulus),
                    "simultaneous Fourier conjugates did not commute")
            conjugate_commutators += 1
    return commutators, conjugate_commutators, modulus, omega


def four_atom_hostile():
    # Atom order: 00,01,10,11.  The mixed Boolean difference is nonzero,
    # while the two coordinate multiplication operators commute exactly.
    p = [1, 1, 0, 0]
    q = [1, 0, 1, 0]
    f = [x * y for x, y in zip(p, q)]
    mixed = f[0] - f[1] - f[2] + f[3]
    mp = [[p[i] if i == j else 0 for j in range(4)] for i in range(4)]
    mq = [[q[i] if i == j else 0 for j in range(4)] for i in range(4)]
    require(mixed == 1, "four-atom mixed difference should be one")
    require(matmul(mp, mq) == matmul(mq, mp),
            "four-atom insertion commutator should vanish")
    return mixed


def same_marginal_holonomy_hostile():
    identity = eye(P)
    marginal_checks = 0
    monodromy_checks = 0
    for a in range(1, P):
        translated = shift(a)
        require(row_sums(identity) == row_sums(translated) == [1] * P,
                "source marginals differ")
        require(col_sums(identity) == col_sums(translated) == [1] * P,
                "target marginals differ")
        marginal_checks += 2 * P
        require(matpow(identity, N) == identity,
                "identity monodromy failed")
        require(matpow(translated, N) == shift(N * a),
                "translation monodromy failed")
        require(shift(N * a) != identity,
                "nonzero seven-edge holonomy disappeared")
        monodromy_checks += 1
    return marginal_checks, monodromy_checks


def formal_label_translation_control():
    # Work at a split prime exactly as in the THM-2594 controls.  A delta
    # table has every primitive Fourier coefficient equal to one.  Under
    # T_a A(ell,theta)=A(ell+1,theta+a), its mode multiplier is
    # xi^beta*zeta^(alpha*a); after seven edges it is zeta^(7*alpha*a).
    modulus = 547
    root91 = next(x for x in range(2, modulus)
                  if exact_order(x, 91, modulus))
    xi = pow(root91, 13, modulus)
    zeta = pow(root91, 7, modulus)
    require(exact_order(xi, 7, modulus), "xi does not have order seven")
    require(exact_order(zeta, 13, modulus), "zeta does not have order thirteen")
    edge_checks = 0
    seam_checks = 0
    for a in range(1, P):
        for alpha in range(1, P):
            for beta in range(1, N):
                edge = (pow(xi, beta, modulus) *
                        pow(zeta, alpha * a, modulus)) % modulus
                seam = pow(zeta, N * alpha * a, modulus)
                require(edge != 1,
                        "primitive edge multiplier unexpectedly trivial")
                require((edge - 1) % modulus != 0,
                        "edge difference unexpectedly zero")
                require(seam != 1,
                        "seven-edge label seam unexpectedly trivial")
                require((seam - 1) % modulus != 0,
                        "seven-edge label difference unexpectedly zero")
                edge_checks += 1
                seam_checks += 1
    require(edge_checks == 864 and seam_checks == 864,
            "wrong primitive label-control count")
    return modulus, root91, edge_checks, seam_checks


def twisted_return_checks():
    partial_fill_checks = 0
    diagonal_failures = 0
    full_returns = 0
    for a in range(1, P):
        base_return = shift(N * a)

        # Every product of physical-root diagonal vertex insertions remains
        # diagonal, so it has no q -> q-7a entry.
        d = diagonal([q + 1 for q in range(P)])
        diagonal_product = matpow(d, N)
        twisted = matmul(diagonal_product, base_return)
        require(all(twisted[q][q] == 0 for q in range(P)),
                "diagonal vertex product passed twisted-return test")
        diagonal_failures += P

        # A correction -a on t of the seven ordered edges leaves (7-t)a.
        for t in range(N + 1):
            correction = matpow(shift(-a), t)
            combined = matmul(correction, base_return)
            require(combined == shift((N - t) * a),
                    "partial-fill displacement mismatch")
            require((combined == eye(P)) == (t == N),
                    "partial-fill closure boundary mismatch")
            partial_fill_checks += 1

        full_correction = matpow(shift(-a), N)
        require(full_correction == shift(-N * a),
                "full correction has wrong monodromy")
        require(matmul(full_correction, base_return) == eye(P),
                "full twisted return did not close")
        require(all(full_correction[q][(q - N * a) % P] == 1
                    for q in range(P)),
                "full correction lacks the required paths")
        full_returns += P
    return partial_fill_checks, diagonal_failures, full_returns


def main():
    selectors, marker_checks, last_choice = selector_bank()
    comm, conjugate_comm, p_fourier, omega = diagonal_algebra_checks()
    mixed = four_atom_hostile()
    marginal_checks, monodromy_checks = same_marginal_holonomy_hostile()
    p91, root91, edge_checks, seam_checks = formal_label_translation_control()
    partial, diagonal_failures, full_returns = twisted_return_checks()

    print("THM-2602 exact vertex/transition curvature controls")
    print(f"selector_bank={selectors} marker_selector_checks={marker_checks} "
          f"last_choice_selectors={last_choice}")
    print(f"diagonal_basis_commutators={comm} "
          f"fourier_conjugate_commutators={conjugate_comm} "
          f"split_prime={p_fourier} root13={omega}")
    print(f"four_atom_mixed_difference={mixed} insertion_commutator=0")
    print(f"same_vertex_marginal_entries_checked={marginal_checks} "
          f"distinct_monodromies={monodromy_checks}")
    print(f"label_translation_prime={p91} root91={root91} "
          f"primitive_edge_checks={edge_checks} seam_checks={seam_checks}")
    print(f"partial_fill_checks={partial} diagonal_twisted_zeros={diagonal_failures} "
          f"full_twisted_returns={full_returns}")
    print("verdict=PASS: vertex calculus telescopes; ordered common-carrier "
          "transition is the minimal missing operation")


if __name__ == "__main__":
    main()
