#!/usr/bin/env python3
"""Exact Cohn-core repair ladder and odd/even rectangle audit.

This companion does not assert the external Weihrauch reduction described in
the motivating draft.  It verifies the finite combinatorial lemma used there
and an exact algebraic analogue in the first transport equation for the
determinant-one matrix

    [[1+xy, x^2], [-y^2, 1-xy]].

The local algebra is over QQ.  The abstract unweighted cycle is deliberately
kept separate from the factorial-weighted Cohn seam: parity supplies a support
mode, while the edge-gain product decides whether that mode actually closes.
All checks use ``require`` so ordinary and optimized Python runs have
identical validity gates.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path

import sympy as sp


Q = Fraction


def require(label: str, condition: bool, checks: list[str]) -> None:
    if not condition:
        raise RuntimeError(f"FAIL: {label}")
    checks.append(label)


def subsets(items: tuple[int, ...]):
    for size in range(1, len(items) + 1):
        yield from combinations(items, size)


def rectangle_profile(n: int) -> tuple[int, int]:
    """Count full-row and balanced 2-by-(n/2) rectangles of size n."""

    rows = (0, 1)
    cols = tuple(range(n))
    full_rows = 0
    balanced = 0
    for a in subsets(rows):
        for b in subsets(cols):
            if len(a) * len(b) != n:
                continue
            if len(a) == 1 and len(b) == n:
                full_rows += 1
            elif len(a) == 2 and 2 * len(b) == n:
                balanced += 1
            else:
                raise RuntimeError("unexpected size-n rectangle geometry")
    return full_rows, balanced


def rectangle_coboundary_count(n: int) -> int:
    """Enumerate 2-by-n binary tables with zero 2-by-2 rectangle parity."""

    count = 0
    for bits in product((0, 1), repeat=2 * n):
        top = bits[:n]
        bottom = bits[n:]
        if all(top[0] ^ top[j] ^ bottom[0] ^ bottom[j] == 0 for j in range(1, n)):
            # Pin x_0=0.  Then y_j=e_0j and x_1=e_10+y_0 recover
            # the table.  The unpinned potentials have one global-bit gauge.
            y_bits = top
            x_one = bottom[0] ^ y_bits[0]
            reconstructed = top + tuple(x_one ^ bit for bit in y_bits)
            if reconstructed != bits:
                raise RuntimeError("rectangle-potential reconstruction failed")
            count += 1
    return count


def orthogonal_grid_count(n: int) -> int:
    """Construct and validate every nonstandard orthogonal grid geometry.

    The standard row/column grid contributes one.  For n=2m, a nonstandard
    grid is determined by an unordered half-set {S,S^c} and independent
    perfect matchings from S to S^c in the two target rows.
    """

    omega = {(row, col) for row in (0, 1) for col in range(n)}
    count = 1
    if n % 2:
        return count
    m = n // 2
    all_cols = set(range(n))
    for s_tuple in combinations(range(n), m):
        s = set(s_tuple)
        complement = all_cols - s
        # Tooth labels are forgotten, so retain one of S,S^c.
        if tuple(sorted(s)) > tuple(sorted(complement)):
            continue
        for matching_zero in permutations(sorted(complement)):
            for matching_one in permutations(sorted(complement)):
                tooth_zero = {(row, col) for row in (0, 1) for col in s}
                tooth_one = omega - tooth_zero
                stripes = []
                for row, matching in ((0, matching_zero), (1, matching_one)):
                    stripes.extend(
                        ({(row, source), (row, target)})
                        for source, target in zip(sorted(s), matching)
                    )
                if set().union(*stripes) != omega:
                    raise RuntimeError("stripe partition does not cover omega")
                if sum(len(stripe) for stripe in stripes) != len(omega):
                    raise RuntimeError("stripe partition overlaps")
                for tooth in (tooth_zero, tooth_one):
                    if any(len(tooth & stripe) != 1 for stripe in stripes):
                        raise RuntimeError("tooth/stripe intersection is not singleton")
                count += 1
    return count


def cyclic_matrix(n: int) -> sp.Matrix:
    return sp.Matrix(
        n,
        n,
        lambda i, j: int(j == i) + int(j == ((i - 1) % n)),
    )


def weighted_cyclic_matrix(weights: tuple[sp.Rational, ...]) -> sp.Matrix:
    n = len(weights)
    return sp.Matrix(
        n,
        n,
        lambda i, j: weights[i] * int(j == i) + int(j == ((i - 1) % n)),
    )


def main() -> None:
    checks: list[str] = []
    x, y = sp.symbols("x y")

    core = sp.Matrix([[1 + x * y, x**2], [-y**2, 1 - x * y]])
    require("core determinant one", sp.expand(core.det()) == 1, checks)
    top_curl = sp.diff(core[0, 0], y) - sp.diff(core[0, 1], x)
    bottom_curl = sp.diff(core[1, 0], y) - sp.diff(core[1, 1], x)
    require("core top curl", top_curl == -x, checks)
    require("core bottom curl", bottom_curl == -y, checks)

    def transport(z_expr: sp.Expr) -> sp.Expr:
        return sp.expand(
            (1 + x * y) * sp.diff(z_expr, y)
            - x**2 * sp.diff(z_expr, x)
            - x * z_expr
        )

    # The degree-D truncation is forced term by term.  It solves the repair
    # equation up to one terminal monomial and can never solve it exactly.
    for degree in range(13):
        z_degree = sum(
            sp.Rational((-1) ** i, sp.factorial(i + 2)) * x**i * y ** (i + 2)
            for i in range(degree + 1)
        )
        residual = sp.Rational((-1) ** degree, sp.factorial(degree + 2))
        residual *= x ** (degree + 1) * y ** (degree + 2)
        require(
            f"forced factorial truncation degree {degree}",
            sp.expand(transport(z_degree) - y - residual) == 0,
            checks,
        )

    formal_closed = (sp.exp(-x * y) - 1 + x * y) / x**2
    require(
        "formal exponential solves the transport equation",
        sp.simplify(transport(formal_closed) - y) == 0,
        checks,
    )
    upper_formal = -(sp.exp(x * y) - 1 - x * y) / y**2
    require(
        "symmetric upper formal transport",
        sp.simplify(
            (1 - x * y) * sp.diff(upper_formal, x)
            + y**2 * sp.diff(upper_formal, y)
            + y * upper_formal
            + x
        )
        == 0,
        checks,
    )

    # Generic finite polynomial identity behind the all-W one-left/one-right
    # obstruction.  The proof uses its highest-degree term; this check guards
    # every sign and derivative in the identity itself.
    w_coeffs = sp.symbols("w0:6")
    z_coeffs = sp.symbols("z0:6")
    monomials_two = (1, x, y, x**2, x * y, y**2)
    w_generic = sum(c * monomial for c, monomial in zip(w_coeffs, monomials_two))
    z_generic = sum(c * monomial for c, monomial in zip(z_coeffs, monomials_two))
    p_lower = -y**2 + (1 + x * y) * z_generic
    q_lower = 1 - x * y + x**2 * z_generic + w_generic * p_lower
    require(
        "generic right-shear lower-curl identity",
        sp.expand(
            sp.diff(p_lower, y)
            - sp.diff(q_lower, x)
            - (-y + transport(z_generic) - sp.diff(w_generic * p_lower, x))
        )
        == 0,
        checks,
    )

    # After b_i=(i+2)! a_i, the homogeneous coefficient recurrence is
    # b_i+b_(i-1)=0.  Closing it on an n-cycle gives I+shift.
    cyclic_rows: list[str] = []
    for n in range(3, 13):
        matrix = cyclic_matrix(n)
        expected_rank = n if n % 2 else n - 1
        expected_det = 2 if n % 2 else 0
        require(f"cycle determinant n={n}", matrix.det() == expected_det, checks)
        require(f"cycle rank n={n}", matrix.rank() == expected_rank, checks)
        if n % 2:
            solution = matrix.inv() * sp.eye(n)[:, 0]
            require(
                f"odd forced cycle half-integral n={n}",
                all(2 * value in (-1, 1) for value in solution),
                checks,
            )
            ambiguity = 0
        else:
            alternating = sp.Matrix([(-1) ** i for i in range(n)])
            require(
                f"even alternating kernel n={n}",
                matrix * alternating == sp.zeros(n, 1),
                checks,
            )
            augmented = matrix.row_join(sp.eye(n)[:, 0])
            require(
                f"even unit forcing obstruction n={n}",
                augmented.rank() == n,
                checks,
            )
            ambiguity = 1
        full_rows, balanced = rectangle_profile(n)
        require(f"two full rows n={n}", full_rows == 2, checks)
        if n % 2:
            require(f"odd rectangle rigidity n={n}", balanced == 0, checks)
        else:
            require(
                f"even balanced rectangles n={n}",
                balanced == sp.binomial(n, n // 2),
                checks,
            )
        cyclic_rows.append(
            f"n={n}:det={expected_det},rank={expected_rank},"
            f"alternating_kernel={ambiguity},balanced_rectangles={balanced}"
        )

        # The actual Cohn ladder has weights alpha_i=i+2.  Flattening its
        # interior recurrence by factorials does not flatten the wrap seam;
        # its weighted holonomy remains (n+1)!, so no cyclic kernel survives.
        cohn_weights = tuple(sp.Rational(i + 2) for i in range(n))
        cohn_cycle = weighted_cyclic_matrix(cohn_weights)
        cohn_expected_det = sp.factorial(n + 1) - (-1) ** n
        require(
            f"weighted Cohn cycle determinant n={n}",
            cohn_cycle.det() == cohn_expected_det,
            checks,
        )
        require(
            f"weighted Cohn cycle full rank n={n}",
            cohn_cycle.rank() == n,
            checks,
        )

    reciprocal_cycle = weighted_cyclic_matrix((sp.Rational(2), sp.Rational(1, 2)))
    require("reciprocal two-cycle has even kernel", reciprocal_cycle.det() == 0, checks)
    require(
        "reciprocal two-cycle kernel witness",
        reciprocal_cycle * sp.Matrix([1, -2]) == sp.zeros(2, 1),
        checks,
    )

    for n in range(2, 7):
        count = rectangle_coboundary_count(n)
        require(
            f"rectangle coboundaries n={n}",
            count == 2 ** (n + 1),
            checks,
        )
        expected_grids = 1 if n % 2 else 1 + sp.factorial(n) // 2
        require(
            f"orthogonal grid classification n={n}",
            orthogonal_grid_count(n) == expected_grids,
            checks,
        )

    semantic_lines = [
        "status=PROVED-SYMBOLIC+FINITE-EXACT; no JC(2) or Weihrauch conclusion",
        "core=[[1+xy,x^2],[-y^2,1-xy]];det=1;curls=(-x,-y)",
        ("lower_repair=(1+xy)Z_y-x^2 Z_x-xZ=y;"
         "a_i=(-1)^i/(i+2)!;no finite polynomial truncation"),
        ("formal_solutions=lower:(exp(-xy)-1+xy)/x^2;"
         "upper:-(exp(xy)-1-xy)/y^2"),
        ("all_W_identity=curl_lower(E_-(Z) C E_+(W))="
         "-y+L(Z)-d_x(W[-y^2+(1+xy)Z]);proof in reflection"),
        "scaled_cycle_operator=I+backward_shift",
        *cyclic_rows,
        ("rectangle_size_n=full row when gcd(2,n)=1;"
         "even n also has 2-by-(n/2)"),
        ("orthogonal_grids=1 for odd n;1+n!/2 for even n;"
         "extra datum is relative matching holonomy"),
        ("weighted_cycle_det=product(alpha_i)-(-1)^n;"
         "Cohn product=(n+1)! so every actual cyclic closure is nonsingular"),
        ("first_even_weighted_kernel=two-cycle weights (2,1/2);"
         "support parity without multiplier holonomy is insufficient"),
        "zero_rectangle_parity_tables=2^(n+1);potentials unique modulo one global bit",
        ("search_router=couple terminal residues on support cycles and retain "
         "weighted holonomy;single acyclic factorial ladders are exact no-go controls"),
    ]
    semantic_digest = sha256(("\n".join(semantic_lines) + "\n").encode()).hexdigest()
    source_digest = sha256(Path(__file__).read_bytes()).hexdigest()
    for line in semantic_lines:
        print(line)
    print(f"checks={len(checks)}")
    print(f"semantic_sha256={semantic_digest}")
    print(f"source_sha256={source_digest}")


if __name__ == "__main__":
    main()
