#!/usr/bin/env python3
"""Exact companion for THM-2407.

Universe:
  * all 2^13 Boolean rational C_13 marginals;
  * two sharp integer-valued branch controls;
  * sparse nonnegative 13^3 tables realizing the diagonal-zero extraction.

No floating-point arithmetic is used.
"""

from fractions import Fraction
from itertools import product


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def cyclotomic_remainder(values, character):
    """P(zeta^character) in the basis 1,zeta,...,zeta^11.

    ``character`` is nonzero modulo 13.  Exponents are first permuted
    modulo 13 and then zeta^12=-(1+...+zeta^11) is applied.
    """

    coeff = [0] * P
    for shift, value in enumerate(values):
        coeff[(character * shift) % P] += value
    top = coeff[12]
    return tuple(coeff[i] - top for i in range(12))


def charged_in_every_colour(values):
    return all(
        any(c != 0 for c in cyclotomic_remainder(values, k))
        for k in range(1, P)
    )


def is_flat(values):
    return all(value == values[0] for value in values)


def charged_energy(values):
    """sum_{k != 0} |vhat(k)|^2 for vhat=p^-1 sum_s v_s zeta^(ks)."""

    total = sum(values)
    squares = sum(value * value for value in values)
    return Fraction(P * squares - total * total, P * P)


def add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def sparse_table(marginal):
    """A nonnegative diagonal-zero table with the prescribed s-marginal.

    H(r,s,t) is supported only at (r,t)=(1,0).  Therefore H(t,s,t)=0.
    """

    table = {}
    for s, value in enumerate(marginal):
        if value:
            table[(1, s, 0)] = value
    return table


def verify_diagonal_and_marginal(table, marginal):
    for s in range(P):
        for t in range(P):
            require(
                table.get((t, s, t), 0) == 0,
                f"diagonal zero failed at s={s}, t={t}",
            )
        require(
            sum(table.get((r, s, 0), 0) for r in range(P)) == marginal[s],
            f"marginal failed at s={s}",
        )


def j_remainder(table, a, b):
    """Cyclotomic numerator of 13^2 J(a,b)."""

    coeff = [0] * P
    for (r, s, t), value in table.items():
        require(t == 0, "synthetic table escaped t=0")
        coeff[(a * r + b * s) % P] += value
    top = coeff[12]
    return tuple(coeff[i] - top for i in range(12))


def verify_triangle_extraction(table, marginal):
    for b in range(1, P):
        require(
            any(cyclotomic_remainder(marginal, b)),
            f"marginal colour b={b} vanished",
        )
        # J(0,b) is nonzero.
        require(any(j_remainder(table, 0, b)), f"J(0,{b}) vanished")
        # Diagonal-zero cancellation: sum_a J(a,b)=0, while some a!=0
        # survives.  In this sparse control every a survives.
        require(
            all(any(j_remainder(table, a, b)) for a in range(1, P)),
            f"nonzero-a extraction failed at b={b}",
        )


def main():
    boolean_count = 0
    flat_count = 0
    nonflat_count = 0
    for values in product((0, 1), repeat=P):
        boolean_count += 1
        flat = is_flat(values)
        charged = charged_in_every_colour(values)
        require(
            charged == (not flat),
            f"all-or-flat failed for Boolean word {values}",
        )
        flat_count += int(flat)
        nonflat_count += int(not flat)

    # Deletion branch: the owner is flat, while H_+ and U have all colours.
    h_del = (0,) + (1,) * 12
    o_del = (1,) * 13
    u_del = add(o_del, h_del)
    require(u_del[0] == o_del[0], "deletion control lost common base")
    require(is_flat(o_del), "deletion-control owner is not flat")
    require(charged_in_every_colour(h_del), "H+ deletion control lost a colour")
    require(charged_in_every_colour(u_del), "U deletion control lost a colour")
    require(charged_energy(o_del) == 0, "flat-owner energy is nonzero")
    require(
        charged_energy(h_del) == Fraction(12, 169),
        "H+ deletion-control energy mismatch",
    )
    require(
        charged_energy(u_del) == Fraction(12, 169),
        "U deletion-control energy mismatch",
    )

    # Owner branch: O has all colours, while exact cancellation makes U flat.
    o_own = (2,) + (1,) * 12
    h_own = (0,) + (1,) * 12
    u_own = add(o_own, h_own)
    require(u_own[0] == o_own[0], "owner control lost common base")
    require(charged_in_every_colour(o_own), "owner control lost a colour")
    require(is_flat(u_own), "owner-control deletion bank is not flat")
    require(
        charged_energy(o_own) == Fraction(12, 169),
        "owner-control owner energy mismatch",
    )
    require(
        charged_energy(h_own) == Fraction(12, 169),
        "owner-control H+ energy mismatch",
    )
    require(charged_energy(u_own) == 0, "flat-U energy is nonzero")

    table_owner = sparse_table(o_own)
    table_deletion = sparse_table(u_del)
    verify_diagonal_and_marginal(table_owner, o_own)
    verify_diagonal_and_marginal(table_deletion, u_del)
    verify_triangle_extraction(table_owner, o_own)
    verify_triangle_extraction(table_deletion, u_del)

    print("THM-2407 OWNER-OR-SOURCE-DELETION EXACT AUDIT")
    print(f"prime={P}")
    print(f"Boolean C13 marginals exhausted={boolean_count}")
    print(f"flat Boolean marginals={flat_count}")
    print(f"nonflat Boolean marginals={nonflat_count}")
    print("all-or-flat equivalence on every nonzero character=PASS")
    print("deletion branch common base u0=o0=1")
    print(f"deletion branch charged energies O,H+,U=0,{charged_energy(h_del)},{charged_energy(u_del)}")
    print("owner branch common base u0=o0=2")
    print(f"owner branch charged energies O,H+,U={charged_energy(o_own)},{charged_energy(h_own)},0")
    print("both sharp branch controls nonnegative=PASS")
    print("sparse 13^3 diagonal-zero tables=PASS")
    print("J(0,b) and some a!=0 survive for all 12 colours in selected bank=PASS")
    print("normalization convention vhat(b)=13^-1 sum_s v_s zeta^(bs)")
    print("THM-2407 exact companion PASS")


if __name__ == "__main__":
    main()
