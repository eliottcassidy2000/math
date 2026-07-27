#!/usr/bin/env python3
"""Exact referee for THM-2548.

The computation is dependency-free.  It verifies the integral transfer
sequence, all C7 x C13 characters in F_547, the full-norm factorization,
the Hall aligned/swap boundary, and the neutral-gain dynamic-program controls.
"""

from __future__ import annotations


P_ROOT = 13
Q_CLOCK = 7
N = P_ROOT * Q_CLOCK
FIELD = 547


def prime_factors(n: int) -> list[int]:
    ans: list[int] = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            ans.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        ans.append(n)
    return ans


def primitive_root(p: int) -> int:
    factors = prime_factors(p - 1)
    for g in range(2, p):
        if all(pow(g, (p - 1) // ell, p) != 1 for ell in factors):
            return g
    raise AssertionError("no primitive root")


def transfer(x: list[int]) -> list[int]:
    """D=sum_{j=0}^6 S^j, with (Sx)_i=x_{i-1}."""
    return [sum(x[(i - j) % N] for j in range(Q_CLOCK)) for i in range(N)]


def shift(x: list[int], steps: int = 1) -> list[int]:
    return [x[(i - steps) % N] for i in range(N)]


def root_fibre_sums(y: list[int]) -> list[int]:
    """The orbit indices i=k+7t have fixed clock k and run over all roots."""
    return [sum(y[k + Q_CLOCK * t] for t in range(P_ROOT))
            for k in range(Q_CLOCK)]


def integral_preimage(y: list[int]) -> list[int]:
    """Invert D integrally on tables with clock-independent root sum."""
    sums = root_fibre_sums(y)
    assert len(set(sums)) == 1
    x = [0] * N
    for k in range(Q_CLOCK):
        for t in range(1, P_ROOT):
            i = k + Q_CLOCK * t
            prev = k + Q_CLOCK * (t - 1)
            x[i] = x[prev] + y[i] - y[(i - 1) % N]
        # The omitted closing recurrence is exactly Y_k-Y_(k-1)=0.
        assert (x[k] - x[k + Q_CLOCK * (P_ROOT - 1)]
                == y[k] - y[(k - 1) % N])
    dx = transfer(x)
    errors = [y[i] - dx[i] for i in range(N)]
    assert len(set(errors)) == 1
    error = errors[0]
    # Adding one constant root fibre adds that constant to every 7-window.
    for t in range(P_ROOT):
        x[Q_CLOCK * t] += error
    assert transfer(x) == y
    return x


def rank_mod(matrix: list[list[int]], p: int) -> int:
    a = [[v % p for v in row] for row in matrix]
    rows = len(a)
    cols = len(a[0]) if rows else 0
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][c], p - 2, p)
        a[r] = [(v * inv) % p for v in a[r]]
        for i in range(rows):
            if i != r and a[i][c]:
                q = a[i][c]
                a[i] = [(a[i][j] - q * a[r][j]) % p
                        for j in range(cols)]
        r += 1
        if r == rows:
            break
    return r


def transfer_matrix() -> list[list[int]]:
    matrix = [[0] * N for _ in range(N)]
    for col in range(N):
        e = [0] * N
        e[col] = 1
        image = transfer(e)
        for row, value in enumerate(image):
            matrix[row][col] = value
    return matrix


def image_basis() -> list[list[int]]:
    basis: list[list[int]] = []
    # Sum-zero basis inside each of the seven 13-point root fibres.
    for k in range(Q_CLOCK):
        anchor = k + Q_CLOCK * (P_ROOT - 1)
        for t in range(P_ROOT - 1):
            y = [0] * N
            y[k + Q_CLOCK * t] = 1
            y[anchor] = -1
            basis.append(y)
    # One vector of common root-fibre sum one.
    y = [0] * N
    for k in range(Q_CLOCK):
        y[k] = 1
    basis.append(y)
    assert len(basis) == N - (Q_CLOCK - 1)
    return basis


def kernel_basis() -> list[list[int]]:
    basis: list[list[int]] = []
    for k in range(Q_CLOCK - 1):
        x = [0] * N
        for t in range(P_ROOT):
            x[k + Q_CLOCK * t] = 1
            x[Q_CLOCK - 1 + Q_CLOCK * t] = -1
        basis.append(x)
    return basis


def sumset(sets: list[set[int]]) -> set[int]:
    out = {0}
    for values in sets:
        out = {(x + y) % P_ROOT for x in out for y in values}
    return out


def main() -> None:
    matrix = transfer_matrix()
    rank = rank_mod(matrix, 1_000_003)
    assert rank == 85

    operator_identity_checks = 0
    for i in range(N):
        e = [0] * N
        e[i] = 1
        de = transfer(e)
        lhs = [u - v for u, v in zip(shift(de), de)]
        rhs = [u - v for u, v in zip(shift(e, Q_CLOCK), e)]
        assert lhs == rhs
        operator_identity_checks += 1

    kb = kernel_basis()
    assert len(kb) == 6
    assert all(transfer(v) == [0] * N for v in kb)
    assert rank_mod(kb, 1_000_003) == 6

    ib = image_basis()
    assert len(ib) == 85
    assert rank_mod(ib, 1_000_003) == 85
    preimages = [integral_preimage(v) for v in ib]
    assert all(transfer(x) == y for x, y in zip(preimages, ib))
    # The six root-sum differences are a primitive quotient map: e_k maps
    # to its kth standard basis vector.  Hence im(D) is primitive and the
    # Smith form has 85 ones and six zeros.
    constraint_controls = []
    for k in range(Q_CLOCK - 1):
        e = [0] * N
        e[k] = 1
        sums = root_fibre_sums(e)
        constraint_controls.append(tuple(sums[j] - sums[-1]
                                         for j in range(Q_CLOCK - 1)))
    assert constraint_controls == [tuple(1 if i == j else 0 for j in range(6))
                                   for i in range(6)]

    # Absolute clock sheets collide positively: the transfer remembers the
    # root-charged quotient, not a preferred clock origin.
    sheet_images = []
    for k in (0, 1):
        sheet = [0] * N
        for t in range(P_ROOT):
            sheet[k + Q_CLOCK * t] = 1
        sheet_images.append(transfer(sheet))
    assert sheet_images[0] == sheet_images[1] == [1] * N

    generator = primitive_root(FIELD)
    xi = pow(generator, (FIELD - 1) // Q_CLOCK, FIELD)
    zeta = pow(generator, (FIELD - 1) // P_ROOT, FIELD)
    assert pow(xi, Q_CLOCK, FIELD) == 1 and xi != 1
    assert pow(zeta, P_ROOT, FIELD) == 1 and zeta != 1

    character_checks = 0
    factor_checks = 0
    root_charged_survive = 0
    clock_only_killed = 0
    full_nonconstant_killed = 0
    for a in range(1, P_ROOT):
        for alpha in range(P_ROOT):
            for beta in range(Q_CLOCK):
                # Fourier multiplier for the pullback convention; changing
                # every sign merely conjugates the calculation.
                q = (pow(xi, (-beta) % Q_CLOCK, FIELD)
                     * pow(zeta, (-alpha * a) % P_ROOT, FIELD)) % FIELD
                d = sum(pow(q, j, FIELD) for j in range(Q_CLOCK)) % FIELD
                expected_zero = alpha == 0 and beta != 0
                assert (d == 0) == expected_zero
                if alpha != 0:
                    root_charged_survive += 1
                if expected_zero:
                    clock_only_killed += 1

                full = sum(pow(q, j, FIELD) for j in range(N)) % FIELD
                assert (full == 0) == (alpha != 0 or beta != 0)
                if alpha != 0 or beta != 0:
                    full_nonconstant_killed += 1
                root_norm = sum(pow(q, Q_CLOCK * t, FIELD)
                                for t in range(P_ROOT)) % FIELD
                assert full == d * root_norm % FIELD
                character_checks += 1
                factor_checks += 1

    assert root_charged_survive == 12 * 12 * 7
    assert clock_only_killed == 12 * 6
    assert full_nonconstant_killed == 12 * 90

    # If the root translation is zero, the seven-step transfer becomes a
    # separate clock norm on every root and kills all 13*6 beta!=0 modes.
    zero_holonomy_killed = 0
    for alpha in range(P_ROOT):
        for beta in range(1, Q_CLOCK):
            q = pow(xi, (-beta) % Q_CLOCK, FIELD)
            d = sum(pow(q, j, FIELD) for j in range(Q_CLOCK)) % FIELD
            assert d == 0
            zero_holonomy_killed += 1
    assert zero_holonomy_killed == 78

    # THM-2539 -> THM-2535 coefficient-level cut amplification.  The old
    # nonzero coefficient is factored out; this checks the new geometric
    # multiplier for 3 slopes, 6 clock colours, 12 target colours, 6 scales.
    lambdas = (1, 2, 3)
    cut_amplification = 0
    for lam in lambdas:
        tau = lam
        for kappa in range(1, Q_CLOCK):
            for b in range(1, P_ROOT):
                alpha = lam * b % P_ROOT
                for scale in range(1, Q_CLOCK):
                    scale_inv = pow(scale, -1, Q_CLOCK)
                    base = (pow(zeta, (-alpha * tau) % P_ROOT, FIELD)
                            * pow(xi, (kappa * scale_inv) % Q_CLOCK,
                                  FIELD)) % FIELD
                    geometric = sum(pow(base, j, FIELD)
                                    for j in range(Q_CLOCK)) % FIELD
                    assert pow(base, Q_CLOCK, FIELD) != 1
                    assert geometric != 0
                    cut_amplification += 1
    assert cut_amplification == 1296

    # THM-2545 aligned/swap survives tensoring with the positive transfer of
    # one mapping-torus atom.  One-point margins agree; the diagonal does not.
    point = [0] * N
    point[0] = 1
    packet = transfer(point)
    assert sum(packet) == 7 and all(v in (0, 1) for v in packet)
    aligned = ((1, 0), (0, 1))
    swapped = ((0, 1), (1, 0))
    margins = []
    diagonals = []
    for coupling in (aligned, swapped):
        left = [sum(packet) * sum(coupling[i][j] for j in range(2))
                for i in range(2)]
        right = [sum(packet) * sum(coupling[i][j] for i in range(2))
                 for j in range(2)]
        margins.append((left, right))
        diagonals.append(sum(packet) * sum(coupling[i][i] for i in range(2)))
    assert margins[0] == margins[1] == ([7, 7], [7, 7])
    assert diagonals == [14, 0]

    # Directed neutral-gain controls for the conditional holotopy/Hall bridge.
    invoice = (-7) % P_ROOT
    neutral_flat = sumset([{0}] * Q_CLOCK)
    neutral_hostile = sumset([{invoice}] + [{0}] * (Q_CLOCK - 1))
    neutral_abundant = sumset([{0, 1, 2}] * Q_CLOCK)
    assert invoice not in neutral_flat
    assert invoice in neutral_hostile
    assert neutral_abundant == set(range(P_ROOT))
    # Even after a later active selector exists everywhere, a cyclic root
    # shift b=h+1 has no Hall diagonal and uniform one-point margins.
    hall_active_diagonal = sum(1 for h in range(P_ROOT)
                               if (h + 1) % P_ROOT == h)
    assert hall_active_diagonal == 0

    print(f"field=F_{FIELD} generator={generator} xi={xi} zeta={zeta}")
    print(f"transfer_rank={rank} kernel_rank={len(kb)} image_rank={len(ib)}")
    print("smith_form=1^85,0^6 primitive_image=1")
    print(f"operator_identity_checks={operator_identity_checks} kernel_basis_checks={len(kb)} integral_image_preimages={len(preimages)}")
    print("positive_clock_sheet_collision=1")
    print(f"character_checks={character_checks} factor_checks={factor_checks}")
    print(f"root_charged_survive={root_charged_survive}/{12*12*7} clock_only_killed={clock_only_killed}/{12*6}")
    print(f"full_norm_nonconstant_killed={full_nonconstant_killed}/{12*90}")
    print(f"zero_holonomy_clock_modes_killed={zero_holonomy_killed}/78")
    print(f"signed_cut_amplification={cut_amplification}/1296")
    print(f"hall_transfer_margins={margins[0][0]} diagonal_aligned={diagonals[0]} diagonal_swap={diagonals[1]}")
    print(f"neutral_gain_invoice={invoice} flat={sorted(neutral_flat)} hostile={sorted(neutral_hostile)} abundant={len(neutral_abundant)}")
    print("active_shift_hall_diagonal=0")


if __name__ == "__main__":
    main()
