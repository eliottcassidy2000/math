#!/usr/bin/env python3
"""Exact row-fourteen same-row obstruction scout beyond THM-4410.

Assume the restricted row-thirteen continuation certified by the companion
session scout: start with the complete fixed residual-weight-at-most-fourteen
source and the three late channels

    lambda18*p**3*y**4 + nu18*y**6 + kappa20*p*y**6.

Rebuild its global row-eleven through row-thirteen source graph, compute the
complete row-fourteen bracket quotient, and audit every exact-valuation-
fourteen monomial p**(14-2*b)*y**b, b=0,...,7.  The exact same-row response
has rank one against two active conditions, so no projected-depth claim is
made.  This is relative to proved THM-4410 and is not a complete
weight-at-most-22 face.
"""

from __future__ import annotations

import contextlib
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


def load_row13_module(path: Path):
    """Load the proved THM-4410 companion even under sparse checkout."""

    if path.is_file():
        spec = importlib.util.spec_from_file_location("row14_thm4410", path)
        module = importlib.util.module_from_spec(spec)
        runner = lambda: spec.loader.exec_module(module)
    else:
        source = subprocess.check_output(
            ["git", "show", f"HEAD:{path.as_posix()}"], text=True
        )
        module = types.ModuleType("row14_thm4410")
        module.__file__ = str(path.resolve())
        runner = lambda: exec(compile(source, path.as_posix(), "exec"), module.__dict__)
    with contextlib.redirect_stdout(io.StringIO()):
        runner()
    return module


Q13 = load_row13_module(Q13_PATH)


Q11 = Q13.Q11
R9 = Q13.R9
R8 = Q13.R8
x = Q13.x
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


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


def unique_primitives(values: list[sp.Expr]) -> list[sp.Expr]:
    answer: list[sp.Expr] = []
    for value in values:
        value = exact(value)
        if value == 0:
            continue
        candidate = primitive(value)
        if all(not proportional(candidate, old) for old in answer):
            answer.append(candidate)
    return answer


def row_thirteen_state():
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
    kappa_source13 = exact(g13 + kappa20 * x**6)
    select13, conditions13, matrix13, pivots13 = R9.solve_bracket_fibre(
        a12, c12, 13, kappa_source13, free12
    )
    check(
        (matrix13.shape, matrix13.rank(), pivots13, len(conditions13))
        == ((14, 9), 9, tuple(range(9)), 1),
        "conditional row-thirteen bracket geometry",
    )
    equation13 = primitive(conditions13[0])
    response13 = exact(sp.diff(equation13, kappa20))
    check(
        response13 == 146361421124462229507072000000000,
        "conditional row-thirteen kappa response",
    )
    kappa_graph = exact(sp.solve(equation13, kappa20)[0])
    check(not sp.fraction(sp.together(kappa_graph))[1].free_symbols,
          "conditional row-thirteen constant denominator")

    a12_selected = R9.impose(a12, select13, {kappa20: kappa_graph})
    c12_selected = R9.impose(c12, select13, {kappa20: kappa_graph})
    a13_trial, c13_trial, theta13 = R9.append_tangent(
        a12_selected, c12_selected, 13, "row14_scout_theta13"
    )
    depth13, amatrix13, cmatrix13 = R9.projected_depth_residuals(
        a13_trial, c13_trial, 13
    )
    terminal13, depth_conditions13, depth_matrix13, depth_pivots13 = (
        R9.eliminate_linear_fibre(depth13, theta13)
    )
    check(
        (amatrix13.shape, amatrix13.rank(), cmatrix13.shape, cmatrix13.rank())
        == ((133, 308), 97, (147, 491), 117),
        "conditional row-thirteen depth universes",
    )
    check(
        (depth_matrix13.shape, depth_matrix13.rank(), depth_pivots13)
        == ((66, 14), 5, (9, 10, 11, 12, 13)),
        "conditional row-thirteen depth selector",
    )
    check(not depth_conditions13,
          "conditional row-thirteen depth automatic")

    full_graph = dict(source_graph)
    full_graph[kappa20] = kappa_graph
    full_graph = {
        symbol: resolve(value, full_graph) for symbol, value in full_graph.items()
    }
    base = {R8.Phi, R8.eta, R8.alpha11, R9.c51}
    check(
        set().union(*(value.free_symbols for value in full_graph.values())) <= base,
        "conditional row-thirteen graph over affine four",
    )
    check(
        all(not sp.fraction(sp.together(value))[1].free_symbols
            for value in full_graph.values()),
        "conditional row-thirteen graph denominators constant",
    )
    a13 = R9.impose(a13_trial, terminal13)
    c13 = R9.impose(c13_trial, terminal13)
    a13 = [resolve(value, full_graph) for value in a13]
    c13 = [resolve(value, full_graph) for value in c13]
    free13 = tuple(
        symbol for index, symbol in enumerate(theta13)
        if index not in depth_pivots13
    )
    check(len(free13) == 9, "conditional row-thirteen terminal affine nine")

    # Recover the complete H14 source row through the already proved row-ten
    # graph, then append all three late tails and resolve the source graph.
    _, _, _, _, grows, bracket_subs, gate, graph10 = Q11.row_ten_state()
    g14_base = exact(grows[14].subs(bracket_subs).subs(gate).subs(graph10))
    g14 = exact(
        g14_base
        + lambda18 * R9.tcoeff(R8.p**3 * R8.y**4, 14)
        + nu18 * R9.tcoeff(R8.y**6, 14)
        + kappa20 * R9.tcoeff(R8.p * R8.y**6, 14)
    )
    g14 = resolve(g14, full_graph)
    return a13, c13, free13, g14, full_graph, kappa20, response13


def main() -> None:
    a13, c13, free13, g14, inherited_graph, kappa20, response13 = (
        row_thirteen_state()
    )

    difference14 = exact(g14 - R8.predicted_G(14, a13, c13))
    equations14 = [R9.xcoeff(difference14, degree) for degree in range(15)]
    check(
        exact(difference14 - sum(equations14[j] * x**j for j in range(15))) == 0,
        "row-fourteen coefficient universe exhaustive",
    )
    select14, conditions14, matrix14, pivots14 = R9.solve_bracket_fibre(
        a13, c13, 14, g14, free13
    )
    reduced14 = [exact(value.subs(select14)) for value in equations14]
    primitives14 = unique_primitives(reduced14)
    raw_cokernel14 = matrix14.T.nullspace()

    # At every lower valuation, distinct powers x**b form the leading
    # diagonal.  Hence cancellations cannot hide a lower-valuation monomial:
    # modulo valuation >14, every Hamiltonian of valuation at least 14 is a
    # linear combination of the eight channels audited below.
    for valuation in range(14):
        diagonal = tuple(
            R8.p ** (valuation - 2*b) * R8.y**b
            for b in range(valuation // 2 + 1)
        )
        check(
            tuple(exact(R9.tcoeff(value, valuation)) for value in diagonal)
            == tuple(x**b for b in range(valuation // 2 + 1)),
            f"valuation-{valuation} leading independence",
        )

    # Full exact-valuation-fourteen diagonal, including every channel through
    # the weight-21 boundary.  Its leading row is x**b and all earlier rows
    # vanish identically.
    rho = sp.symbols("rho14_0:8")
    monomials14 = tuple(R8.p ** (14 - 2*b) * R8.y**b for b in range(8))
    weights14 = tuple(28 - b for b in range(8))
    leading14 = tuple(exact(R9.tcoeff(value, 14)) for value in monomials14)
    check(leading14 == tuple(x**b for b in range(8)),
          "valuation-fourteen leading diagonal")
    check(
        all(exact(R9.tcoeff(value, row)) == 0
            for value in monomials14 for row in range(14)),
        "valuation-fourteen firewall",
    )
    perturbed_g14 = exact(g14 + sum(rho[b] * leading14[b] for b in range(8)))
    select_all, conditions_all, matrix_all, pivots_all = R9.solve_bracket_fibre(
        a13, c13, 14, perturbed_g14, free13
    )
    check(matrix_all == matrix14 and pivots_all == pivots14,
          "valuation-fourteen common tangent quotient")
    difference_all = exact(perturbed_g14 - R8.predicted_G(14, a13, c13))
    equations_all = [R9.xcoeff(difference_all, degree) for degree in range(15)]
    reduced_all = [exact(value.subs(select_all)) for value in equations_all]
    response_full = sp.Matrix([
        [exact(sp.diff(value, rho[b])) for b in range(8)]
        for value in reduced_all
    ])

    student14 = R9.R9.primitive_student_row(14)
    student_base = exact(sum(student14[j] * equations14[j] for j in range(15)))
    student_responses = tuple(exact(sum(
        student14[j] * sp.diff(equations_all[j], rho[b]) for j in range(15)
    )) for b in range(8))
    check(student_responses == (1671525, 0, 371450, 0, 267444, 0, 348840, 0),
          "valuation-fourteen Student response table")
    check(response_full.rank() == 1,
          "valuation-fourteen full response rank one")
    check(tuple(response_full[:, b].is_zero_matrix for b in range(8))
          == (False, True, False, True, False, True, False, True),
          "valuation-fourteen parity boundary")
    response_ratios = tuple(exact(value / student_responses[6])
                            for value in student_responses)
    for b in range(8):
        check(response_full[:, b] == response_ratios[b] * response_full[:, 6],
              f"valuation-fourteen full response b={b}")

    channel_primitives = unique_primitives(reduced_all)
    check(len(channel_primitives) == 2,
          "two primitive conditions remain with all same-row channels")
    rho_set = set(rho)
    unpaid = tuple(
        value for value in channel_primitives
        if value.free_symbols.isdisjoint(rho_set)
    )
    paid = tuple(
        value for value in channel_primitives
        if not value.free_symbols.isdisjoint(rho_set)
    )
    check(len(unpaid) == len(paid) == 1,
          "one paid and one untouched primitive condition")
    base_origin = {R8.Phi: 0, R8.eta: 0, R8.alpha11: 0, R9.c51: 0}
    check(exact(unpaid[0].subs(base_origin)) != 0,
          "untouched condition is a proper nonzero hypersurface")
    check(exact(sp.diff(paid[0], rho[6])) != 0,
          "least-weight visible channel has constant pivot")

    print("row14_weight22_conditional_continuation_scout")
    print("status=FINITE_EXACT_RELATIVE_TO_PROVED_THM4410")
    print("inherited_source_projection=global_A4(Phi,eta,alpha11,c51)")
    print("inherited_solved=(c23,c70,lambda18,nu18,kappa20)")
    print(f"inherited_kappa20_response={response13}")
    print(
        f"row14_bracket_matrix={matrix14.shape};rank={matrix14.rank()};"
        f"pivots={pivots14};raw_cokernel={len(raw_cokernel14)}"
    )
    print(f"row14_selected_nonzero={sum(value != 0 for value in reduced14)}")
    print(f"row14_primitive_condition_count={len(primitives14)}")
    for index, value in enumerate(primitives14):
        variables = tuple(sorted(value.free_symbols, key=str))
        poly = sp.Poly(value, *variables, domain=sp.QQ) if variables else None
        print(
            f"row14_condition_{index}_terms="
            f"{len(poly.terms()) if poly is not None else 1};"
            f"degree={poly.total_degree() if poly is not None else 0};"
            f"symbols={tuple(map(str, variables))}"
        )
        print(f"row14_condition_{index}={sp.sstr(value)}")
    print(f"valuation14_weights={weights14}")
    print(f"valuation14_full_response_rank={response_full.rank()}")
    print(f"valuation14_student_responses={sp.sstr(student_responses)}")
    print(f"valuation14_ratios_to_p2y6={sp.sstr(response_ratios)}")
    for b in range(8):
        print(
            f"valuation14_channel=b{b},monomial=p^{14-2*b}*y^{b},"
            f"weight={weights14[b]},student={student_responses[b]},"
            f"full_zero={response_full[:, b].is_zero_matrix}"
        )
    print(f"row14_student_base={sp.sstr(primitive(student_base))}")
    print(f"row14_conditions_with_channels={len(conditions_all)}")
    for index, value in enumerate(channel_primitives):
        print(f"row14_channel_condition_{index}={sp.sstr(value)}")
    print(f"row14_unpaid_same_row_condition={sp.sstr(unpaid[0])}")
    print("row14_same_row_nogo=every_Hamiltonian_of_t_valuation_at_least_14_has_response_rank_at_most_1_against_2_active_conditions")
    print("row14_least_weight_visible=p^2*y^6;residual_weight=22;odd_weight21_y^7_is_invisible")
    print("scope=restricted_tail_row14_bracket_only;not_complete_weight22;no_projected_depth;not_JC2")
    print(f"checks={CHECKS};result=PASS")


if __name__ == "__main__":
    main()
