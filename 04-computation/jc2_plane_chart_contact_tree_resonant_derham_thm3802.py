#!/usr/bin/env python3
"""Exact companion for THM-3802's intrinsic plane-chart residue law."""

from __future__ import annotations

import hashlib
import itertools
import json

import sympy as sp


CHECKS = 0


def check(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise AssertionError(label)
    CHECKS += 1


def same(lhs: sp.Expr, rhs: sp.Expr, label: str) -> None:
    check(sp.cancel(sp.expand(lhs - rhs)) == 0, label)


f, w, a = sp.symbols("f w a")


def coeff(poly: sp.Expr, degree: int) -> sp.Expr:
    return sp.Poly(sp.expand(poly), f).coeff_monomial(f**degree)


# Exact chart transitions, primitive signs, and termwise Laurent primitives
# across six pole orders.  The profiles include constant leaves, lower moving
# jets, the resonant jet, and nonuniform f-adic contact orders.
for n in range(1, 7):
    low_1 = f if n >= 3 else sp.Integer(0)
    low_2 = 2 * f**2 if n >= 4 else sp.Integer(0)
    betas = (
        f ** (n - 1),
        3 + low_1 + 2 * f ** (n - 1),
        8 - 2 * low_1 + low_2 - f ** (n - 1),
        15 + 3 * low_1 - low_2 + 4 * f ** (n - 1),
    )
    for beta in betas:
        check(sp.Poly(beta, f).degree() < n, f"n={n} truncated leaf jet")

    rho = tuple(coeff(beta, n - 1) for beta in betas)
    for i, j in itertools.combinations(range(len(betas)), 2):
        ai = (w - betas[i]) / f**n
        aj = (w - betas[j]) / f**n
        same(aj - ai, (betas[i] - betas[j]) / f**n,
             f"n={n} chart transition {i},{j}")

        # eta=-a df, so eta_j-eta_i=(beta_j-beta_i)f^-n df.
        primitive_difference = sp.expand((betas[j] - betas[i]) / f**n)
        check(primitive_difference.coeff(f, -1) == rho[j] - rho[i],
              f"n={n} logarithmic residue {i},{j}")

        difference = sp.Poly(sp.expand(betas[j] - betas[i]), f)
        for (degree,), coefficient_value in difference.terms():
            if degree == n - 1:
                continue
            primitive = (
                coefficient_value
                * f ** (degree - n + 1)
                / (degree - n + 1)
            )
            same(
                sp.diff(primitive, f),
                coefficient_value * f ** (degree - n),
                f"n={n} nonresonant primitive {i},{j},degree={degree}",
            )


# The Cech q=0 row is the simplex cochain complex.  The q=1 row loses its
# vertex term because H1(A2)=0; its edge kernel is therefore the reduced leaf
# module of dimension h-1 and it has no higher cohomology.
for h in range(1, 9):
    vertices = list(range(h))
    edges = list(itertools.combinations(vertices, 2))
    triangles = list(itertools.combinations(vertices, 3))
    tetrahedra = list(itertools.combinations(vertices, 4))
    edge_index = {edge: i for i, edge in enumerate(edges)}
    triangle_index = {face: i for i, face in enumerate(triangles)}

    d0 = sp.zeros(len(edges), h)
    for row, (i, j) in enumerate(edges):
        d0[row, i] = -1
        d0[row, j] = 1

    d1 = sp.zeros(len(triangles), len(edges))
    for row, (i, j, k) in enumerate(triangles):
        d1[row, edge_index[(j, k)]] = 1
        d1[row, edge_index[(i, k)]] = -1
        d1[row, edge_index[(i, j)]] = 1

    d2 = sp.zeros(len(tetrahedra), len(triangles))
    for row, face in enumerate(tetrahedra):
        for omitted in range(4):
            triangle = tuple(face[j] for j in range(4) if j != omitted)
            d2[row, triangle_index[triangle]] = (-1) ** omitted

    check(d0.rank() == max(0, h - 1), f"h={h} reduced leaf rank")
    check(d1 * d0 == sp.zeros(len(triangles), h), f"h={h} d1d0=0")
    check(d2 * d1 == sp.zeros(len(tetrahedra), len(edges)), f"h={h} d2d1=0")
    check(len(edges) - d1.rank() == max(0, h - 1),
          f"h={h} q1 edge kernel")
    check(len(triangles) - d2.rank() - d1.rank() == 0,
          f"h={h} q1 higher exactness")
    check(h - d0.rank() == 1, f"h={h} q0 constants only")


# Equality and coordinate-invariance controls.  With at least two distinct
# leaves, exponent one can never have constant resonance; the positive exact
# controls therefore start at n=2.
for n in range(2, 7):
    r, s = sp.symbols(f"r{n} s{n}")
    base = (sp.Integer(0), sp.Integer(2), sp.Integer(7))
    common = tuple(beta + r * f ** (n - 1) for beta in base)
    shifted = tuple(beta + s * f ** (n - 1) for beta in common)
    check(len(set(coeff(beta, n - 1) for beta in common)) == 1,
          f"n={n} exact constant resonance")
    before = [coeff(beta, n - 1) for beta in common]
    after = [coeff(beta, n - 1) for beta in shifted]
    check(all(sp.expand(after[i] - before[i] - s) == 0 for i in range(3)),
          f"n={n} common target translation")
    check(all(sp.expand((after[i] - after[0]) - (before[i] - before[0])) == 0
              for i in range(3)), f"n={n} reduced class invariant")

check(len(set(fixed_roots := (sp.Integer(0), sp.Integer(2), sp.Integer(5)))) == 3,
      "n=1 distinct leaves force nonconstant resonance")


# Canonical controls from the two parent mechanisms.
check(tuple(coeff(beta, 0) for beta in fixed_roots) == fixed_roots,
      "THM3791 fixed exponent-one class")
check(tuple(coeff(beta, 2) for beta in fixed_roots) == (0, 0, 0),
      "THM3791 fixed exponent-three exactness")

lambda_1, lambda_2, ell = sp.symbols("lambda_1 lambda_2 ell", nonzero=True)
moving_axis = (f**2, lambda_1, lambda_2)
check(tuple(coeff(beta, 2) for beta in moving_axis) == (1, 0, 0),
      "THM3789 moving-axis control")

confluent = (f**2, sp.Integer(0), ell)
check(tuple(coeff(beta, 2) for beta in confluent) == (1, 0, 0),
      "THM3797 confluent contact-tree control")
check(sp.Poly(confluent[0] - confluent[1], f).as_dict() == {(2,): 1},
      "THM3797 contact order two")
check(coeff(confluent[2] - confluent[1], 0) == ell,
      "THM3797 transverse third leaf")


semantic = {
    "class": "rho_i=[f^(n-1)]beta_i;[omega]=rho_mod_constant_leaf_vector",
    "controls": "THM3791_split;THM3789_moving_axis;THM3797_confluent_tree",
    "darboux": "nonconstant_rho_implies_no_global_polynomial_Darboux_pair",
    "hypotheses": "actual_smooth_affine_A2_chart_cover;all_multi_overlaps=D(f);w=beta_i+f^n*a_i",
    "scope_guard": "leaf_jets_do_not_imply_existence_affineness_or_hypersurface_origin",
    "topology": "H1=0;H2=reduced_leaf_module=k^h/k1;all_other_positive_degrees_zero",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3802-plane-chart-contact-tree-resonant-de-rham-law")
print("hypotheses=actual_smooth_affine_A2_cover;multi_overlaps=D(f);w=beta_i+f^n*a_i")
print("cohomology=H1=0;H2=k^h/k*1;other_positive_degrees_zero")
print("class=[omega]=([f^(n-1)]beta_i)_i_mod_constants")
print("exactness=omega_exact_iff_resonant_leaf_vector_constant")
print("darboux=nonconstant_resonance_implies_no_pair")
print("controls=THM3791_split;THM3789_axis;THM3797_confluent_contact_tree")
print("guard=atlas_is_hypothesis_not_conclusion_from_leaf_labels")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
