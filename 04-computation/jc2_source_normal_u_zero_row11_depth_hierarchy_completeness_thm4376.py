#!/usr/bin/env python3
"""Primary exact certificate for proposed THM-4376.

This deliberately imports the canonized THM-4366 row-ten carrier.  It then
independently appends row eleven
and tests projected P_2/P_3 depth compatibility.  It makes no JC(2) claim.
"""

from __future__ import annotations

from pathlib import Path
import contextlib
import io
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

import jc2_source_normal_bracket_hasse_rows8_thm4308 as R8  # noqa: E402
import jc2_source_normal_student_stein_row9_thm4315 as R9  # noqa: E402
import jc2_source_normal_u_zero_row11_hierarchy_selected_extinction_thm4366 as H11  # noqa: E402


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise RuntimeError(f"check failed: {label}")
    CHECKS += 1


def capture_inherited_locals() -> dict[str, object]:
    captured: dict[str, object] = {}

    def tracer(frame, event, arg):
        if frame.f_code is H11.nonzero_branch.__code__ and event == "return":
            captured.update(frame.f_locals)
        return tracer

    arows, crows, bracket_subs = R8.build_bracket_rows()
    _, _, theta8, terminal8 = R8.hasse_audit(arows, crows, bracket_subs)
    old = sys.gettrace()
    sys.settrace(tracer)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            H11.nonzero_branch(arows, crows, bracket_subs, theta8, terminal8)
    finally:
        sys.settrace(old)
    require(bool(captured), "capture inherited THM-4366 exact construction")
    return captured


def selected_solution(matrix: sp.Matrix, rhs: sp.Matrix, columns: tuple[int, ...], variables):
    row_indices = tuple(matrix.T.rref()[1])
    square = matrix.extract(row_indices, columns)
    require(square.rows == square.cols and square.det() != 0, "selected square invertible")
    other = tuple(j for j in range(matrix.cols) if j not in columns)
    effective_rhs = rhs.extract(row_indices, (0,))
    if other:
        effective_rhs -= matrix.extract(row_indices, other) * sp.Matrix([variables[j] for j in other])
    values = square.inv() * effective_rhs
    return tuple(sp.cancel(value) for value in values)


def exact_zero(expression: sp.Expr) -> bool:
    return sp.cancel(sp.together(sp.expand(expression))) == 0


def hierarchy_value(rows, m: int, ell: int, q: int):
    s = (ell + 1) // 2
    return sp.factor(
        sp.cancel(
            sum(
                (-1) ** (n - s)
                * sp.binomial(m + q - n, q)
                * R8.xcoeff(rows[n], 2 * n - ell)
                for n in range(s, m + 1)
            )
        )
    )


def hierarchy_vector(coordinates, m: int, ell: int, q: int) -> sp.Matrix:
    s = (ell + 1) // 2
    index = {coordinate: j for j, coordinate in enumerate(coordinates)}
    result = sp.zeros(len(coordinates), 1)
    for n in range(s, m + 1):
        coordinate = (n, 2 * n - ell)
        if coordinate in index:
            result[index[coordinate]] = (-1) ** (n - s) * sp.binomial(m + q - n, q)
    return result


def hierarchy_bank(rows, m: int, depth: int):
    labels = []
    values = []
    for ell in range(2, 2 * m + 1):
        rho = (ell + 2) // 3
        for q in range(rho):
            if m + q >= ell and depth <= m + q - ell:
                labels.append((ell, q))
                values.append(hierarchy_value(rows, m, ell, q))
    return labels, values


def admissible_labels(m: int, depth: int):
    return {
        (ell, q)
        for ell in range(2, 2 * m + 1)
        for q in range((ell + 2) // 3)
        if m + q >= ell and depth <= m + q - ell
    }


def main() -> None:
    data = capture_inherited_locals()
    rows_a10 = data["row10_a_joint"]
    rows_c10 = data["row10_c_joint"]
    q10s = tuple(data["theta10_symbols"])
    free10 = q10s[:8]
    P = H11.P
    qp = H11.QP

    require(len(rows_a10) == len(rows_c10) == 11, "inherited rows through ten")
    require(free10 == q10s[:8], "inherited affine-eight fibre")
    require(qp.degree() == 6 and qp.TC() != 0 and sp.gcd(qp, qp.diff()) == sp.Poly(1, P, domain=sp.QQ), "six nonzero simple source points")
    require(sp.gcd(H11.Q, H11.R11) == sp.Poly(1, H11.Y, domain=sp.QQ), "inherited bracket coprimality")

    b11 = R8.B_row(11, rows_a10, rows_c10)
    a11base, c11base = R8.particular_row(11, b11)
    theta11s = tuple(sp.symbols("theta11_0:12"))
    theta11 = sum(theta11s[j] * R8.x**j for j in range(12))
    a11 = sp.expand(a11base + theta11 * sp.diff(R8.A0, R8.x))
    c11 = sp.expand(c11base + theta11 * sp.diff(R8.C0, R8.x))
    rows_a11 = list(rows_a10) + [a11]
    rows_c11 = list(rows_c10) + [c11]

    acoords, amat = R9.depth_matrix(2, 11)
    ccoords, cmat = R9.depth_matrix(3, 11)
    avec = sp.Matrix([R8.xcoeff(rows_a11[n], r) for n, r in acoords])
    cvec = sp.Matrix([R8.xcoeff(rows_c11[n], r) for n, r in ccoords])
    ares = [sp.expand((left.T * avec)[0]) for left in amat.T.nullspace()]
    cres = [sp.expand((left.T * cvec)[0]) for left in cmat.T.nullspace()]
    variables = free10 + theta11s
    am, ar = sp.linear_eq_to_matrix(ares, variables)
    cm, cr = sp.linear_eq_to_matrix(cres, variables)
    jm = am.col_join(cm)
    jr = ar.col_join(cr)

    require((len(acoords), amat.cols, amat.rank(), len(amat.T.nullspace())) == (102, 228, 77, 25), "A row-eleven universe")
    require((len(ccoords), cmat.cols, cmat.rank(), len(cmat.T.nullspace())) == (114, 361, 94, 20), "C row-eleven universe")
    require(len(variables) == 20, "twenty extension variables")
    require((am.rank(), am.row_join(ar).rank(), am.rref()[1]) == (4, 4, (6, 7, 18, 19)), "A extension ranks")
    require((cm.rank(), cm.row_join(cr).rank(), cm.rref()[1]) == (3, 4, (7, 18, 19)), "C generic extension ranks")
    require((jm.rank(), jm.row_join(jr).rank(), jm.rref()[1]) == (4, 5, (6, 7, 18, 19)), "joint generic extension ranks")

    print("ROW11_DEPTH_UNIVERSES")
    print("A", len(acoords), amat.cols, amat.rank(), len(amat.T.nullspace()))
    print("C", len(ccoords), cmat.cols, cmat.rank(), len(cmat.T.nullspace()))
    print("VARIABLES", len(variables), variables)
    print("A_RANKS", am.rank(), am.row_join(ar).rank(), am.rref()[1])
    print("C_RANKS", cm.rank(), cm.row_join(cr).rank(), cm.rref()[1])
    print("JOINT_RANKS", jm.rank(), jm.row_join(jr).rank(), jm.rref()[1])

    # Test whether all THM-4364 hierarchy rows already span the complete
    # terminal depth equation spaces on this specific recursive source jet.
    alabels, avals = hierarchy_bank(rows_a11, 11, 2)
    clabels, cvals = hierarchy_bank(rows_c11, 11, 3)
    for label, value in zip(alabels, avals):
        vector = hierarchy_vector(acoords, 11, *label)
        require(vector.T * amat == sp.zeros(1, amat.cols), f"A hierarchy module annihilator {label}")
        require(exact_zero((vector.T * avec)[0] - value), f"A hierarchy evaluation {label}")
    for label, value in zip(clabels, cvals):
        vector = hierarchy_vector(ccoords, 11, *label)
        require(vector.T * cmat == sp.zeros(1, cmat.cols), f"C hierarchy module annihilator {label}")
        require(exact_zero((vector.T * cvec)[0] - value), f"C hierarchy evaluation {label}")
    ahm, ahr = sp.linear_eq_to_matrix(avals, variables)
    chm, chr_ = sp.linear_eq_to_matrix(cvals, variables)
    aaug = ahm.row_join(-ahr)
    caug = chm.row_join(-chr_)
    full_aaug = am.row_join(-ar)
    full_caug = cm.row_join(-cr)
    print("A_HIERARCHY", len(alabels), ahm.rank(), aaug.rank(), "STACK", full_aaug.col_join(aaug).rank())
    print("C_HIERARCHY", len(clabels), chm.rank(), caug.rank(), "STACK", full_caug.col_join(caug).rank())
    abasis = tuple(aaug.T.rref()[1])
    cbasis = tuple(caug.T.rref()[1])
    require((len(alabels), ahm.rank(), aaug.rank(), full_aaug.col_join(aaug).rank()) == (24, 4, 4, 4), "A hierarchy terminal completeness")
    require((len(clabels), chm.rank(), caug.rank(), full_caug.col_join(caug).rank()) == (19, 3, 4, 4), "C hierarchy terminal completeness")
    require([alabels[i] for i in abasis] == [(10, 1), (11, 2), (12, 3), (13, 4)], "A hierarchy basis labels")
    require([clabels[i] for i in cbasis] == [(9, 1), (10, 2), (10, 3), (11, 3)], "C hierarchy basis labels")
    print("A_HIERARCHY_AUG_BASIS", [alabels[i] for i in abasis])
    print("C_HIERARCHY_AUG_BASIS", [clabels[i] for i in cbasis])
    anew = sorted(admissible_labels(11, 2) - admissible_labels(10, 2))
    cnew = sorted(admissible_labels(11, 3) - admissible_labels(10, 3))
    amap = dict(zip(alabels, avals))
    cmap = dict(zip(clabels, cvals))
    require(anew == [(9, 0), (10, 1), (11, 2), (12, 3), (13, 4)], "A new clock positions")
    require(cnew == [(8, 0), (9, 1), (10, 2), (11, 3)], "C new clock positions")
    require(amap[(9, 0)] == 0 and all(amap[label] != 0 for label in anew[1:]), "A unique leading-edge silence")
    require(cmap[(8, 0)] == 0 and all(cmap[label] != 0 for label in cnew[1:]), "C unique leading-edge silence")
    require(R8.xcoeff(a11, 13) == 0 and hierarchy_value(rows_a10, 10, 9, 0) == 0, "A silence mechanism")
    require(R8.xcoeff(c11, 14) == 0 and hierarchy_value(rows_c10, 10, 8, 0) == 0, "C silence mechanism")
    print("A_NEW_CLOCK_POSITIONS", anew)
    print("C_NEW_CLOCK_POSITIONS", cnew)
    print("A_AUTOMATIC_SILENCE", [(label, amap[label]) for label in anew if amap[label] == 0])
    print("C_AUTOMATIC_SILENCE", [(label, cmap[label]) for label in cnew if cmap[label] == 0])
    anew_active = [label for label in anew if amap[label] != 0]
    cnew_active = [label for label in cnew if cmap[label] != 0]
    anm, anr = sp.linear_eq_to_matrix([amap[label] for label in anew_active], variables)
    cnm, cnr = sp.linear_eq_to_matrix([cmap[label] for label in cnew_active], variables)
    anp = tuple(anm.rref()[1])
    cnp = tuple(cnm.rref()[1])
    require((anm.rank(), [variables[j] for j in anp], sp.factor(anm[:, anp].det())) == (4, [q10s[6], q10s[7], theta11s[10], theta11s[11]], sp.Rational(5, 4)), "A active-new primitive rank")
    require((cnm.rank(), [variables[j] for j in cnp], sp.factor(cnm[:, cnp].det())) == (3, [q10s[7], theta11s[10], theta11s[11]], sp.Rational(27, 128)), "C active-new primitive rank")
    require(full_aaug.col_join(anm.row_join(-anr)).rank() == full_aaug.rank(), "A new positions span full equations")
    require(cm.col_join(cnm).rank() == cm.rank(), "C new positions span coefficient equations")
    print("A_ACTIVE_NEW", anew_active, "RANK", anm.rank(), "PIVOTS", [variables[j] for j in anp], "DET", sp.factor(anm[:, anp].det()))
    print("C_ACTIVE_NEW", cnew_active, "RANK", cnm.rank(), "PIVOTS", [variables[j] for j in cnp], "DET", sp.factor(cnm[:, cnp].det()))
    print("A_NEW_SPANS_FULL", full_aaug.col_join(anm.row_join(-anr)).rank() == full_aaug.rank())
    print("C_NEW_COEFF_SPANS", cm.col_join(cnm).rank() == cm.rank())

    # Select A's pivot variables, then inspect all remaining C conditions.
    apiv = tuple(am.rref()[1])
    asol = selected_solution(am, ar, apiv, variables)
    asubs = {variables[column]: asol[j] for j, column in enumerate(apiv)}
    require(all(exact_zero(value.subs(asubs)) for value in ares), "A selected solution")
    c_after_a = [sp.factor(sp.cancel(value.subs(asubs))) for value in cres]
    remaining = tuple(v for v in variables if v not in asubs)
    crm, crr = sp.linear_eq_to_matrix(c_after_a, remaining)
    require((len(remaining), crm.rank(), crm.row_join(crr).rank(), crm.rref()[1]) == (16, 0, 1, ()), "C after A is purely scalar")
    print("C_AFTER_A_REMAINING", len(remaining), remaining)
    print("C_AFTER_A_RANKS", crm.rank(), crm.row_join(crr).rank(), crm.rref()[1])

    # If possible, solve remaining coefficient directions and list scalar
    # residual numerators; these are the exact mixed-module gluing ideal.
    rpiv = tuple(crm.rref()[1])
    rsol = selected_solution(crm, crr, rpiv, remaining)
    rsubs = {remaining[column]: rsol[j] for j, column in enumerate(rpiv)}
    residuals = [sp.factor(sp.cancel(value.subs(rsubs))) for value in c_after_a]
    nonzero = [(i, value) for i, value in enumerate(residuals) if value != 0]
    require([index for index, _ in nonzero] == [12, 17], "two C scalar residual rows")
    print("C_AFTER_A_SCALAR_RESIDUAL_COUNT", len(nonzero))
    common = qp
    for index, value in nonzero:
        numerator, _ = sp.fraction(value)
        ppoly = sp.Poly(numerator, P, domain=sp.QQ)
        common = sp.gcd(common, ppoly)
        rem = ppoly.rem(qp)
        print("RESIDUAL", index, "DEG", ppoly.degree(), "GCD_Q_DEG", sp.gcd(ppoly, qp).degree(), "MOD_Q", sp.factor(rem.as_expr()))
        print("RESIDUAL_OVER_Q", index, sp.factor(sp.cancel(value / qp.as_expr())))
    expected_ratios = {
        12: -1 / (sp.Integer(15381515239492820800781250) * P**5),
        17: -1 / (sp.Integer(5127171746497606933593750) * P**5),
    }
    require(all(exact_zero(value / qp.as_expr() - expected_ratios[index]) for index, value in nonzero), "exact Q residual ratios")
    require(common.monic() == qp.monic(), "mixed scalar ideal exactly Q")
    print("COMMON_GCD_WITH_Q", common.as_expr())

    c_hier_after_a = [sp.factor(sp.cancel(value.subs(asubs))) for value in cvals]
    c_hier_nonzero = [(label, value) for label, value in zip(clabels, c_hier_after_a) if value != 0]
    require(len(c_hier_nonzero) == 1 and c_hier_nonzero[0][0] == (10, 3), "one inherited C hierarchy obstruction")
    require(exact_zero(c_hier_nonzero[0][1] / qp.as_expr() - 1 / (sp.Integer(15381515239492820800781250) * P**5)), "C hierarchy generator ratio")
    print("C_HIERARCHY_AFTER_A")
    for label, value in zip(clabels, c_hier_after_a):
        if value != 0:
            print(label, "OVER_Q", sp.factor(sp.cancel(value / qp.as_expr())))

    # The exact first-entry hierarchy orders at row eleven.
    print("FIRST_ENTRY_A")
    for ell in range(2, 25):
        rho = (ell + 2) // 3
        m0 = ell + 2 - rho + 1
        if m0 == 11:
            q = rho - 1
            print((ell, q), hierarchy_value(rows_a11, 11, ell, q))
    print("FIRST_ENTRY_C")
    for ell in range(2, 25):
        rho = (ell + 2) // 3
        m0 = ell + 3 - rho + 1
        if m0 == 11:
            q = rho - 1
            print((ell, q), hierarchy_value(rows_c11, 11, ell, q))

    print("JOINT_COMPATIBILITY iff Q(Phi^2)=0; fibre_dimension=16")
    print("POST_BRACKET_CAUSAL_SCOPE THM-4366 already kills Q at row-eleven bracket; this is an alternative projected-depth survival/redundancy statement only")
    print(f"CHECKS={CHECKS}")
    print("PASS")


if __name__ == "__main__":
    main()
