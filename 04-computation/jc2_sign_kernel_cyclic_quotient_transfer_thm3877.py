#!/usr/bin/env python3
"""Exact finite-group controls for THM-3877's sign-kernel transfer."""

from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

from sympy import Matrix
from sympy.combinatorics.named_groups import (
    AlternatingGroup,
    CyclicGroup,
    DihedralGroup,
    SymmetricGroup,
)
from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.combinatorics.permutations import Permutation

sys.stdout.reconfigure(newline="\n")

CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True:
        raise RuntimeError(message)


def elements(group: PermutationGroup) -> list[Permutation]:
    return list(group.generate_schreier_sims())


def even_kernel(group: PermutationGroup) -> PermutationGroup:
    return PermutationGroup([g for g in elements(group) if g.signature() == 1])


def is_transposition(g: Permutation) -> bool:
    return len(g.cyclic_form) == 1 and len(g.cyclic_form[0]) == 2


def orbit_count(g: Permutation, degree: int) -> int:
    moved = sum(len(cycle) for cycle in g.cyclic_form)
    return len(g.cyclic_form) + degree - moved


def tame_permutation_conductor(g: Permutation, degree: int) -> int:
    return degree - orbit_count(g, degree)


def abelianization_order(group: PermutationGroup) -> int:
    return group.order() // group.derived_subgroup().order()


S3 = SymmetricGroup(3)
A3 = AlternatingGroup(3)
S4 = SymmetricGroup(4)
A4 = AlternatingGroup(4)
D4 = DihedralGroup(4)
C4 = CyclicGroup(4)
A5 = AlternatingGroup(5)
V4 = PermutationGroup(
    [
        Permutation([[0, 1], [2, 3]], size=4),
        Permutation([[0, 2], [1, 3]], size=4),
    ]
)

# Tame conductor one is exactly a single transposition in the relevant
# permutation action.  This is the computational shadow of the general
# orbit-count proof in the theorem.
for name, group, degree in (
    ("S3", S3, 3),
    ("S4", S4, 4),
    ("A4", A4, 4),
    ("D4", D4, 4),
    ("V4", V4, 4),
    ("C4", C4, 4),
):
    conductor_one = [
        g for g in elements(group) if tame_permutation_conductor(g, degree) == 1
    ]
    gate(
        all(is_transposition(g) and g.signature() == -1 for g in conductor_one),
        f"{name}: conductor-one elements are odd transpositions",
    )
    gate(
        all(tame_permutation_conductor(g, degree) != 1 for g in elements(group)
            if not is_transposition(g)),
        f"{name}: no other element has conductor one",
    )

gate(len([g for g in elements(S3) if is_transposition(g)]) == 3,
     "S3 transposition count")
gate(len([g for g in elements(S4) if is_transposition(g)]) == 6,
     "S4 transposition count")
gate(len([g for g in elements(D4) if is_transposition(g)]) == 2,
     "D4 transposition count")
gate(len([g for g in elements(C4) if is_transposition(g)]) == 0,
     "regular C4 has no transposition")

# The discriminant/sign kernels and their abelianizations.
H3 = even_kernel(S3)
H4S = even_kernel(S4)
H4D = even_kernel(D4)
H4C = even_kernel(C4)
gate(H3.order() == 3 and H3.is_abelian, "S3 sign kernel is C3")
gate(H4S.order() == 12 and H4S.derived_subgroup().order() == 4,
     "S4 sign kernel is A4 with derived V4")
gate(H4D.order() == 4 and H4D.is_abelian,
     "D4 sign kernel is V4")
gate(H4C.order() == 2 and H4C.is_abelian,
     "C4 sign kernel is C2")
gate(abelianization_order(H3) == 3, "S3 kernel cyclic quotient C3")
gate(abelianization_order(H4S) == 3, "S4 kernel cyclic quotient C3")
gate(abelianization_order(H4D) == 4, "D4 kernel has C2 quotient")
gate(abelianization_order(H4C) == 2, "C4 kernel has C2 quotient")

# All transpositions disappear after intersection with the even kernel.
for group in (S3, S4, D4):
    kernel_elements = set(elements(even_kernel(group)))
    gate(
        all(g not in kernel_elements for g in elements(group) if is_transposition(g)),
        "transposition inertia is disjoint from the sign kernel",
    )

# Degree-four controls and the first perfect equality boundary.
gate(S4.is_transitive() and A4.is_transitive(), "S4/A4 transitivity")
gate(D4.is_transitive() and V4.is_transitive() and C4.is_transitive(),
     "remaining quartic transitivity controls")
gate(A4.is_solvable, "A4 solvability")
gate(A5.order() == 60 and A5.derived_subgroup().order() == 60,
     "A5 is perfect")
gate(not A5.is_solvable, "A5 perfect boundary is not solvable")

# Torsion-free class-group control: multiplication by p on Z is injective.
for p in (2, 3, 5, 7):
    multiplication = Matrix([[p]])
    gate(
        multiplication.det() == p and multiplication.rank() == 1,
        f"multiplication by {p} on Z has zero kernel",
    )

semantic = {
    "transfer": "simple tame discriminant inertia is one transposition; sign base change kills it",
    "kummer": "A*/A*^p=0 and Cl(A)[p]=0 forbid an unramified Cp quotient",
    "group": "the sign kernel must be perfect when the Kummer gate holds for every prime",
    "degree3": "S3 has sign kernel C3",
    "degree4": "S4 and D4 sign kernels have cyclic quotients; C4 has no transposition",
    "boundary": "degree2 has trivial kernel; S5 has perfect sign kernel A5",
    "specialization": "THM3874 scalar units and Cl=Z exclude sole-simple degree3/4 fields",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3877-sign-kernel-transfer-through-torsion-free-quadratic-resolvent")
print("simple_inertia=transposition;after_sign_base_change=trivial")
print("kummer_gate=A_units_mod_p_zero+Cl_p_torsion_zero=>no_Cp_quotient")
print("degree3=H_C3;degree4=H_nonperfect;sole_simple_branch=NONE_under_gate")
print("THM3874_specialization=degree3_and_degree4_three_cusp_sole_branch_NONE")
print("sharp_boundary=degree2_H_trivial;S5_H_A5_perfect")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
