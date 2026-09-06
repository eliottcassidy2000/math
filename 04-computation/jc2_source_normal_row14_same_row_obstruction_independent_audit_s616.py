#!/usr/bin/env python3
"""Independent row-14 cokernel audit relative to proved THM-4410.

This certificate does not import the row-14 scout.  It instantiates the
proved THM-4410 terminal state, builds the row-14 coefficient matrix directly,
and uses its left nullspace rather than the primary elimination routine to
derive the source conditions and every exact-valuation-14 response.
"""

from __future__ import annotations

import contextlib
import hashlib
import importlib.util
import io
from pathlib import Path
import subprocess
import sys
import types

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


Q13_PATH = Path(
    "04-computation/jc2_source_normal_weight18_row13_continuation_scout_s616.py"
)
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def load_row13(path: Path):
    """Load the proved dependency under ordinary or sparse checkout."""

    if path.is_file():
        spec = importlib.util.spec_from_file_location("row14_audit_thm4410", path)
        module = importlib.util.module_from_spec(spec)
        runner = lambda: spec.loader.exec_module(module)
    else:
        source = subprocess.check_output(
            ["git", "show", f"HEAD:{path.as_posix()}"], text=True
        )
        module = types.ModuleType("row14_audit_thm4410")
        module.__file__ = str(path.resolve())
        runner = lambda: exec(
            compile(source, path.as_posix(), "exec"), module.__dict__
        )
    with contextlib.redirect_stdout(io.StringIO()):
        runner()
    return module


Q13 = load_row13(Q13_PATH)
Q11 = Q13.Q11
R9 = Q13.R9
R8 = Q13.R8
x = Q13.x


def exact(value: sp.Expr) -> sp.Expr:
    return R9.exact(value)


def primitive(value: sp.Expr) -> sp.Expr:
    return R9.primitive(value)


def proportional(left: sp.Expr, right: sp.Expr) -> bool:
    return R9.proportional(left, right)


def resolve(value: sp.Expr, mapping: dict[sp.Symbol, sp.Expr]) -> sp.Expr:
    answer = exact(value)
    for _ in range(16):
        updated = exact(answer.subs(mapping, simultaneous=False))
        if updated == answer:
            return answer
        answer = updated
    raise AssertionError("source graph did not resolve")


def unique_primitives(values: tuple[sp.Expr, ...] | list[sp.Expr]) -> list[sp.Expr]:
    answer: list[sp.Expr] = []
    for value in values:
        value = exact(value)
        if value == 0:
            continue
        candidate = primitive(value)
        if all(not proportional(candidate, old) for old in answer):
            answer.append(candidate)
    return answer


def expression_hash(value: sp.Expr) -> str:
    return hashlib.sha256(sp.srepr(sp.expand(value)).encode("ascii")).hexdigest()


def instantiate_thm4410_terminal():
    """Instantiate, but do not reaudit, the proved row-13 terminal state."""

    (
        a12,
        c12,
        free12,
        g13,
        source_graph,
        _,
        lambda18,
        nu18,
    ) = Q13.build_row_twelve_state()

    kappa20 = sp.symbols("kappa20")
    source13 = exact(g13 + kappa20 * x**6)
    select13, conditions13, matrix13, pivots13 = R9.solve_bracket_fibre(
        a12, c12, 13, source13, free12
    )
    check(
        (matrix13.shape, matrix13.rank(), pivots13, len(conditions13))
        == ((14, 9), 9, tuple(range(9)), 1),
        "THM4410 bracket instantiation",
    )
    equation13 = primitive(conditions13[0])
    check(
        exact(sp.diff(equation13, kappa20))
        == 146361421124462229507072000000000,
        "THM4410 constant pivot",
    )
    kappa_graph = exact(sp.solve(equation13, kappa20)[0])

    a12_selected = R9.impose(a12, select13, {kappa20: kappa_graph})
    c12_selected = R9.impose(c12, select13, {kappa20: kappa_graph})
    a13_trial, c13_trial, theta13 = R9.append_tangent(
        a12_selected, c12_selected, 13, "row14_audit_theta13"
    )
    depth13, _, _ = R9.projected_depth_residuals(a13_trial, c13_trial, 13)
    terminal13, depth_conditions13, depth_matrix13, depth_pivots13 = (
        R9.eliminate_linear_fibre(depth13, theta13)
    )
    check(not depth_conditions13, "THM4410 depth instantiation")
    check(
        (depth_matrix13.shape, depth_matrix13.rank(), depth_pivots13)
        == ((66, 14), 5, (9, 10, 11, 12, 13)),
        "THM4410 terminal selector",
    )

    full_graph = dict(source_graph)
    full_graph[kappa20] = kappa_graph
    full_graph = {
        symbol: resolve(value, full_graph) for symbol, value in full_graph.items()
    }
    a13 = [resolve(value, full_graph) for value in R9.impose(a13_trial, terminal13)]
    c13 = [resolve(value, full_graph) for value in R9.impose(c13_trial, terminal13)]
    free13 = tuple(
        symbol for index, symbol in enumerate(theta13)
        if index not in depth_pivots13
    )
    check(len(free13) == 9, "THM4410 terminal fibre dimension")

    _, _, _, _, grows, bracket_subs, gate, graph10 = Q11.row_ten_state()
    g14_base = exact(grows[14].subs(bracket_subs).subs(gate).subs(graph10))
    g14 = exact(
        g14_base
        + lambda18 * R9.tcoeff(R8.p**3 * R8.y**4, 14)
        + nu18 * R9.tcoeff(R8.y**6, 14)
        + kappa20 * R9.tcoeff(R8.p * R8.y**6, 14)
    )
    return a13, c13, free13, resolve(g14, full_graph)


def main() -> None:
    a13, c13, free13, g14 = instantiate_thm4410_terminal()

    difference = exact(g14 - R8.predicted_G(14, a13, c13))
    equations = tuple(R9.xcoeff(difference, degree) for degree in range(15))
    check(
        exact(difference - sum(equations[degree] * x**degree for degree in range(15)))
        == 0,
        "row14 coefficient universe",
    )

    zero_fibre = {symbol: 0 for symbol in free13}
    base = sp.Matrix([exact(value.subs(zero_fibre)) for value in equations])
    matrix = sp.Matrix(
        [[exact(sp.diff(value, symbol)) for symbol in free13] for value in equations]
    )
    reconstructed = matrix * sp.Matrix(free13) + base
    check(
        all(exact(equations[index] - reconstructed[index]) == 0 for index in range(15)),
        "row14 equations affine-linear in terminal fibre",
    )
    check(matrix.shape == (15, 9), "row14 matrix shape")
    check(matrix.rank() == 9, "row14 matrix rank")
    check(matrix.rref()[1] == tuple(range(9)), "row14 column pivots")

    left_kernel = tuple(matrix.T.nullspace())
    check(len(left_kernel) == 6, "row14 raw cokernel dimension")
    raw_conditions = tuple(exact((vector.T * base)[0]) for vector in left_kernel)
    source_conditions = unique_primitives(raw_conditions)
    check(len(source_conditions) == 2, "two active row14 source conditions")

    # Exhaustiveness of the Hamiltonian valuation filtration: at each
    # valuation the leading monomials have distinct powers of x.
    for valuation in range(15):
        diagonal = tuple(
            R8.p ** (valuation - 2*b) * R8.y**b
            for b in range(valuation // 2 + 1)
        )
        check(
            tuple(exact(R9.tcoeff(value, valuation)) for value in diagonal)
            == tuple(x**b for b in range(valuation // 2 + 1)),
            f"valuation-{valuation} leading diagonal",
        )

    rho = sp.symbols("rho14_0:8")
    channel_vector = sp.Matrix(
        [rho[degree] if degree < 8 else 0 for degree in range(15)]
    )
    channel_conditions = tuple(
        exact((vector.T * (base + channel_vector))[0]) for vector in left_kernel
    )
    response = sp.Matrix(
        [
            [exact(sp.diff(value, rho[index])) for index in range(8)]
            for value in channel_conditions
        ]
    )
    check(response.shape == (6, 8), "direct cokernel response shape")
    check(response.rank() == 1, "direct cokernel response rank one")
    check(
        tuple(response[:, index].is_zero_matrix for index in range(8))
        == (False, True, False, True, False, True, False, True),
        "direct parity response",
    )

    pivot_row = next(index for index in range(6) if response[index, 6] != 0)
    ratios = tuple(
        exact(response[pivot_row, index] / response[pivot_row, 6])
        for index in range(8)
    )
    expected_ratios = (
        sp.Rational(115, 24),
        0,
        sp.Rational(115, 108),
        0,
        sp.Rational(23, 30),
        0,
        1,
        0,
    )
    check(ratios == expected_ratios, "direct response ratios")
    for index in range(8):
        check(
            response[:, index] == ratios[index] * response[:, 6],
            f"direct proportional response {index}",
        )

    primitive_channels = unique_primitives(channel_conditions)
    check(len(primitive_channels) == 2, "two conditions after all channels")
    rho_set = set(rho)
    unpaid = tuple(
        value for value in primitive_channels
        if value.free_symbols.isdisjoint(rho_set)
    )
    paid = tuple(
        value for value in primitive_channels
        if not value.free_symbols.isdisjoint(rho_set)
    )
    check(len(unpaid) == len(paid) == 1, "one untouched cokernel direction")

    base_variables = (R8.Phi, R8.eta, R8.alpha11, R9.c51)
    origin = {symbol: 0 for symbol in base_variables}
    unpaid_origin = exact(unpaid[0].subs(origin))
    check(unpaid_origin != 0, "unpaid condition proper at origin")
    paid_pivot = exact(sp.diff(paid[0], rho[6]))
    check(paid_pivot != 0 and not paid_pivot.free_symbols, "constant p2y6 pivot")

    dense = dict(zip(base_variables, map(sp.Integer, (1, 2, 3, 5))))
    dense_unpaid = exact(unpaid[0].subs(dense))
    check(dense_unpaid != 0, "dense off-hypersurface hostile")

    print("source_normal_row14_same_row_obstruction_independent_audit")
    print("imports=proved_THM4410_state;row14_scout_import=no")
    print("method=direct_15x9_coefficient_matrix_and_six_dimensional_left_nullspace")
    print("row14_matrix=15x9;rank=9;column_pivots=0,1,2,3,4,5,6,7,8;raw_cokernel=6")
    print("active_source_conditions=2;degrees=4,4")
    print("valuation14_weights=28,27,26,25,24,23,22,21")
    print("response_rank=1;odd_channels_zero=True;ratios_to_p2y6=115/24,0,115/108,0,23/30,0,1,0")
    print(f"paid_condition_hash={expression_hash(paid[0])};p2y6_pivot={paid_pivot}")
    print(f"unpaid_condition_hash={expression_hash(unpaid[0])};origin_value={unpaid_origin};dense_value={dense_unpaid}")
    print("filtration_law=distinct_x^b_leading_diagonals_for_valuations_0_through_14")
    print("same_row_nogo=all_H_of_t_valuation_at_least_14_miss_one_active_condition")
    print("least_weight_visible=p^2*y^6_weight22;y^7_weight21_is_odd_and_invisible")
    print("scope=restricted_THM4410_tail_row14_bracket_only;no_depth_or_complete_weight22_or_JC2")
    print(f"checks={CHECKS};result=PASS")


if __name__ == "__main__":
    main()
