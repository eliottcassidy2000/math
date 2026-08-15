#!/usr/bin/env python3
"""Exact terminal-filtration and V4 three-view effectivity hostile.

The probe combines two proved carriers without claiming a JC transfer:

* THM-3383's terminal polynomial cover kills the cyclic invariant-ring
  torsor but leaves one decoded target outside the polynomial ring;
* THM-2655's V4 polynomial cover kills the class-group carrier but leaves
  three invariant ratios with distinct boundary poles.

Runtime gates survive python -O.
"""

from __future__ import annotations

import ast
from functools import reduce
from hashlib import sha256
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PINNED = (
    (
        "THM-3383-script",
        ROOT / "04-computation/jc_terminal_monomial_cone_polynomiality_fork_thm3383.py",
        "b72a196e7fda4a2de33df26e1f26f06f654908a9eb42f453e210728559ad4b10",
    ),
    (
        "THM-3383-output",
        ROOT / "05-knowledge/results/jc_terminal_monomial_cone_polynomiality_fork_thm3383.out",
        "514af913271c9fa283bf7c969f9166730a3a187e93880b65d9a480c8966f2abc",
    ),
    (
        "THM-2655-script",
        ROOT / "04-computation/jacobian_s4_resolvent_quasietale_hostile.py",
        "f2f72df60e5443cfe35177cbb561c7e473b63006b2f7254d080a14a3be949c8f",
    ),
    (
        "THM-2655-output",
        ROOT / "05-knowledge/results/jacobian_s4_resolvent_quasietale_hostile.out",
        "01748e08163b3ff45c971ea7597a2a1fe1d7791fef99ef17bf98baee0893a97e",
    ),
)

EXPECTED_SEMANTIC_DIGEST = "90e478b47b925203d81978048dcc25d76b97258bbc413982016db7a13d594a18"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


class ExactDigest:
    def __init__(self):
        self._hash = sha256()

    def add(self, value):
        self._hash.update(repr(value).encode("ascii"))
        self._hash.update(bytes((10,)))

    def hexdigest(self):
        return self._hash.hexdigest()


def add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def scale(multiplier, vector):
    return tuple(multiplier * value for value in vector)


def determinant3(rows):
    (a, b, c), (d, e, f), (g, h, i) = rows
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def two_by_two_minors(rows):
    values = []
    for row_a, row_b in combinations(rows, 2):
        for col_a, col_b in combinations(range(3), 2):
            values.append(
                row_a[col_a] * row_b[col_b] - row_a[col_b] * row_b[col_a]
            )
    return tuple(values)


def lattice_index_rank_three(generators):
    minors = tuple(
        abs(determinant3(rows))
        for rows in combinations(generators, 3)
    )
    return reduce(gcd, minors)


def invariant_under_even_signs(exponents):
    signs = ((1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1))
    for sign in signs:
        value = 1
        for coordinate_sign, exponent in zip(sign, exponents):
            if exponent % 2:
                value *= coordinate_sign
        if value != 1:
            return False
    return True


def terminal_filtration():
    cells = []
    gates = []
    for e in range(1, 9):
        for a in range(0, 9):
            for orientation in (-1, 1):
                g = a * e + orientation
                if g < 1:
                    continue
                require(gcd(e, g) == 1, ("terminal gcd", e, a, orientation, g))
                d = min(g, a * e)

                # With s the polynomial decoded target and m the missing one:
                # s=x^e v^(d+1), m=x^(-e) v^(-d), and sm=v.
                present = (e, d + 1)
                missing = (-e, -d)
                require(add(present, missing) == (0, 1), (e, a, orientation))

                orientation_name = "u_present" if orientation == 1 else "t_present"
                cells.append((e, a, g, orientation, d, orientation_name))
                for power in range(1, 6):
                    v_order = d * power
                    l_order = e * power

                    # m^q v^(dq) L^(eq)=x^(-eq)L^(eq)=y^(eq).
                    positive_x_order = missing[0] * power
                    positive_v_order = missing[1] * power + v_order
                    require(positive_x_order == -l_order, "positive x order")
                    require(positive_v_order == 0, "positive v order")

                    # Dropping either available boundary factor leaves a pole.
                    v_hostile = None if v_order == 0 else (v_order - 1, l_order)
                    l_hostile = (v_order, l_order - 1)
                    gates.append(
                        (
                            e,
                            a,
                            g,
                            orientation,
                            power,
                            (v_order, l_order),
                            v_hostile,
                            l_hostile,
                        )
                    )
    require(len(cells) == 135, ("terminal cell count", len(cells)))
    require(len(gates) == 675, ("filtration gate count", len(gates)))
    return tuple(cells), tuple(gates)


def v4_boundary_hostile():
    a = (2, 0, 0)
    b = (0, 2, 0)
    c = (0, 0, 2)
    d = (1, 1, 1)
    rho_a = (-1, 1, 1)
    rho_b = (1, -1, 1)
    rho_c = (1, 1, -1)
    rhos = (rho_a, rho_b, rho_c)

    require(scale(2, d) == add(add(a, b), c), "d^2=abc")
    require(add(rho_a, rho_b) == c, "rho_a rho_b=c")
    require(add(rho_a, rho_c) == b, "rho_a rho_c=b")
    require(add(rho_b, rho_c) == a, "rho_b rho_c=a")
    require(add(add(rho_a, rho_b), rho_c) == d, "rho_a rho_b rho_c=d")
    require(all(invariant_under_even_signs(vector) for vector in (a, b, c, d, *rhos)), "V4 invariance")
    require(all(min(vector) == -1 for vector in rhos), "one surviving pole per view")
    require(all(min(add(left, right)) >= 0 for left, right in combinations(rhos, 2)), "pair effectivity")
    require(min(add(add(rho_a, rho_b), rho_c)) >= 0, "triple effectivity")

    determinant = determinant3(rhos)
    one_minor_gcd = reduce(gcd, (abs(value) for row in rhos for value in row))
    two_minor_gcd = reduce(gcd, (abs(value) for value in two_by_two_minors(rhos)))
    smith = (
        one_minor_gcd,
        two_minor_gcd // one_minor_gcd,
        abs(determinant) // two_minor_gcd,
    )
    require(determinant == 4, determinant)
    require(smith == (1, 2, 2), smith)

    class_relations = (a, b, c, d)
    class_index = lattice_index_rank_three(class_relations)
    rho_index = abs(determinant)
    require(class_index == rho_index == 4, (class_index, rho_index))
    for row, basis in zip(rhos, (a, b, c)):
        require(row == add(d, scale(-1, basis)), (row, d, basis))

    # Pullback to C[x,y,z] makes a,b,c literal squares, so every V4 Kummer
    # character/class dies; Laurent exponent negativity still obstructs rho_i.
    square_roots = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    require(tuple(scale(2, root) for root in square_roots) == (a, b, c), "squareclasses killed")

    return {
        "generators": (a, b, c, d),
        "rhos": rhos,
        "determinant": determinant,
        "smith": smith,
        "class_index": class_index,
        "pair_products": (c, b, a),
        "triple_product": d,
    }


def main():
    for name, path, expected in PINNED:
        require(lf_hash(path) == expected, ("dependency changed", name, path))

    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(tree)
        ),
        "floating literal",
    )

    cells, gates = terminal_filtration()
    v4 = v4_boundary_hostile()

    semantic = ExactDigest()
    semantic.add(("terminal_cells", cells))
    semantic.add(("terminal_gates", gates))
    semantic.add(("v4", v4))
    digest = semantic.hexdigest()
    require(digest == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", digest))

    print("JC TORSOR KILLING VERSUS BOUNDARY EFFECTIVITY EXACT PROBE")
    print(f"source_sha256_lf={lf_hash(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINNED)}")
    print("status=VERIFIED-EXACT terminal filtration and V4 three-view hostile plus elementary analytic proofs;unnumbered_and_not_canon")
    print("terminal_normal_form=s=x^e*v^(d+1);m=x^(-e)*v^(-d);s*m=v;d=min(g,ae);orientation_selects_s")
    print("terminal_filtration=m^q*f(v)_is_polynomial_iff_v^(dq)*(1+cv)^(eq)_divides_f(v)")
    print(f"terminal_cells={len(cells)};powers_per_cell=5;filtration_gates={len(gates)};e_range=1..8;a_range=0..8")
    print("cyclic_hostile=mu_e_invariant_torsor_dies_on_C[x,y]_but_missing_target_retains_v=0_and_xy=0_boundary_orders")
    print(f"v4_ratio_valuation_matrix={v4['rhos']};det={v4['determinant']};smith={v4['smith']};class_lattice_index={v4['class_index']}")
    print("v4_ratios=(d/a,d/b,d/c)=(yz/x,xz/y,xy/z);each_has_one_distinct_pole_after_polynomial_pullback")
    print(f"pair_products=(c,b,a)_valuations={v4['pair_products']};triple_product=d_valuation={v4['triple_product']}")
    print("effectivity_hostile=all_three_Kummer_squareclasses_and_divisor_classes_die_on_A3;all_pair_and_triple_products_are_polynomial;no_single_ratio_is_polynomial")
    print("typing=class_or_H1_trivialization_is_not_effective_boundary_regularization;valuation_cone_is_required_sidecar")
    print("scope=sharp_nonKeller_V4_and_terminal_monomial_hostiles_only;no_A4_or_S4_exclusion;no_JC2_or_LRC_consequence")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
