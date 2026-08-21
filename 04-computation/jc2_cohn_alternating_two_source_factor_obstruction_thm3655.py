#!/usr/bin/env python3
"""Exact symbolic companion for THM-3655.

The universal obstruction is proved in the theorem by degree and boundary
arguments.  This companion verifies the matrix/curl reductions, the four
load-bearing jets, and a finite hostile bank in both alternating orders.
"""

from __future__ import annotations

import ast
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


EXPECTED_SEMANTIC_SHA256 = "d1f05c8841156820f7ab6573e16315995dfbd197dc2e02e403f27179dc7298e7"
CHECKS = 0


def require(condition: bool, payload: object) -> None:
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


x, y = sp.symbols("x y")
u = sp.Function("u")(x, y)
w = sp.Function("w")(x, y)
A = 1 + x * y
D = 1 - x * y
C = sp.Matrix(((A, x**2), (-y**2, D)))


def eplus(value):
    return sp.Matrix(((1, value), (0, 1)))


def eminus(value):
    return sp.Matrix(((1, 0), (value, 1)))


def curl(row):
    return sp.expand(sp.diff(row[0], y) - sp.diff(row[1], x))


def zero(expression) -> bool:
    return sp.expand(expression) == 0


def main() -> None:
    require(sp.expand(C.det()) == 1, "Cohn determinant")

    plus_minus = sp.simplify(C * eplus(w) * eminus(u))
    r_top = x**2 + A * w
    r_bottom = D - y**2 * w
    expected_plus_minus = sp.Matrix((
        (A + u * r_top, r_top),
        (-y**2 + u * r_bottom, r_bottom),
    ))
    require(all(zero(value) for value in plus_minus - expected_plus_minus),
            "plus-minus row form")
    plus_minus_curls = tuple(curl(plus_minus.row(index)) for index in range(2))
    plus_minus_reduced = (
        x + sp.diff(u * r_top, y) - sp.diff(r_top, x),
        -2 * y + sp.diff(u * r_bottom, y) - sp.diff(r_bottom, x),
    )
    require(all(zero(left - right) for left, right in
                zip(plus_minus_curls, plus_minus_reduced, strict=True)),
            "plus-minus curls")

    minus_plus = sp.simplify(C * eminus(u) * eplus(w))
    s_top = A + x**2 * u
    s_bottom = -y**2 + D * u
    expected_minus_plus = sp.Matrix((
        (s_top, x**2 + w * s_top),
        (s_bottom, D + w * s_bottom),
    ))
    require(all(zero(value) for value in minus_plus - expected_minus_plus),
            "minus-plus row form")
    minus_plus_curls = tuple(curl(minus_plus.row(index)) for index in range(2))
    minus_plus_reduced = (
        sp.diff(s_top, y) - sp.diff(w * s_top, x) - 2 * x,
        sp.diff(s_bottom, y) - sp.diff(w * s_bottom, x) + y,
    )
    require(all(zero(left - right) for left, right in
                zip(minus_plus_curls, minus_plus_reduced, strict=True)),
            "minus-plus curls")
    require(sp.expand(plus_minus.det()) == 1 == sp.expand(minus_plus.det()),
            "decorated determinants")

    # The four boundary data used by the proof.  The two curve substitutions
    # are performed only after multiplying out the vanishing factor, avoiding
    # any undefined evaluation of the arbitrary polynomial placeholder.
    require(zero(r_top - x**2 - A * w), "forward top A-congruence")
    require(zero(r_bottom.subs(y, 0) - 1), "forward bottom constant jet")
    require(zero(sp.diff(r_bottom, y).subs(y, 0) + x),
            "forward bottom linear jet")
    require(zero(s_top.subs(x, 0) - 1), "reverse top constant jet")
    require(zero(sp.diff(s_top, x).subs(x, 0) - y),
            "reverse top linear jet")
    require(zero(s_bottom + y**2 - D * u), "reverse bottom D-congruence")

    # Cheap hostile controls.  They are not the universal proof; they catch
    # sign/order mistakes in the symbolic reductions.
    bank = (1, -1, x, -x, y, -y, x**2, -x**2, x*y, -x*y, y**2, -y**2)
    hostile_counts = {"plus_minus": [0, 0], "minus_plus": [0, 0]}
    for u_value in bank:
        for w_value in bank:
            substitutions = {u: u_value, w: w_value,
                             sp.diff(u, x): sp.diff(u_value, x),
                             sp.diff(u, y): sp.diff(u_value, y),
                             sp.diff(w, x): sp.diff(w_value, x),
                             sp.diff(w, y): sp.diff(w_value, y)}
            for name, curls in (("plus_minus", plus_minus_curls),
                                ("minus_plus", minus_plus_curls)):
                vanished = tuple(zero(value.subs(substitutions)) for value in curls)
                hostile_counts[name][0] += int(vanished[0])
                hostile_counts[name][1] += int(vanished[1])
                require(not vanished[0] and not vanished[1],
                        (name, u_value, w_value, vanished))

    proof_gates = (
        "pm_top:y-leading log-derivative contradiction;w=0 gives x*u_y=1",
        "pm_bottom:x-leading product derivative contradiction;u=0 gives y=2y mod y^2",
        "mp_top:y-leading product derivative contradiction;w=0 gives x*u_y=1",
        "mp_bottom:x-leading log-derivative contradiction;u=0 gives y^2*w_x=y",
    )
    formula_digests = {
        "plus_minus": digest(tuple(sp.srepr(value) for value in plus_minus_curls)),
        "minus_plus": digest(tuple(sp.srepr(value) for value in minus_plus_curls)),
    }
    semantic_record = {
        "ring": "characteristic-zero polynomial ring k[x,y]",
        "orders": ("C E_+(w) E_-(u)", "C E_-(u) E_+(w)"),
        "formula_digests": formula_digests,
        "hostile_bank": len(bank),
        "hostile_pairs_per_order": len(bank) ** 2,
        "hostile_zero_row_counts": hostile_counts,
        "proof_gates": proof_gates,
        "scope": "two alternating right factors and at most one left elementary factor",
    }
    semantic = digest(semantic_record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic))

    source = Path(__file__).resolve()
    source_bytes = source.read_bytes()
    require(b"\r\n" not in source_bytes, "source raw LF")
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source_bytes.decode("utf-8")))),
            "Python assert node present")

    print("== THM-3655 Cohn alternating two-source-factor obstruction ==")
    print("determinants=(C:1,plus_minus:1,minus_plus:1)")
    print("orders=" + repr(semantic_record["orders"]))
    print("formula_digests=" + repr(formula_digests))
    print("boundary_gates=(A_congruence,bottom_y_jets,top_x_jets,D_congruence)=PASS")
    print(f"hostile_bank={len(bank)};pairs_per_order={len(bank)**2};zero_row_counts={hostile_counts}")
    print("proof_gates=" + repr(proof_gates))
    print("semantic_sha256=" + semantic)
    print(f"CHECKS={CHECKS}")
    print("status=SYMBOLIC-EXACT COMPANION;UNIVERSAL DEGREE PROOF IN THEOREM")
    print("scope=no longer words/all non-elementary classes/Keller pair/JC2")


if __name__ == "__main__":
    main()
