#!/usr/bin/env python3
"""Exact companion for THM-3786's irregular 2x3 support closure."""

from __future__ import annotations

import ast
import hashlib
from collections import defaultdict
from itertools import product
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def buckets(d: int, e: int, f: int) -> dict[int, list[tuple[int, int]]]:
    """Contribution buckets after subtracting the common bottom sum."""
    answer: dict[int, list[tuple[int, int]]] = defaultdict(list)
    left = (0, d)
    right = (0, e, e + f)
    for i, u in enumerate(left):
        for j, v in enumerate(right):
            answer[u + v].append((i, j))
    return dict(answer)


def collision_offsets(d: int, e: int, f: int) -> tuple[int, ...]:
    return tuple(sorted(k for k, values in buckets(d, e, f).items()
                        if len(values) > 1))


# Pair-sum census.  Positivity leaves exactly d=e, d=f, or d=e+f; two
# collisions occur exactly on the common-step AP d=e=f.
RANGE = 16
collision_counts = defaultdict(int)
packet_rows: list[str] = []
for d, e, f in product(range(1, RANGE + 1), repeat=3):
    actual = collision_offsets(d, e, f)
    predicted = []
    if d == e:
        predicted.append(d)
    if d == e + f:
        predicted.append(d)
    if d == f:
        predicted.append(d + e)
    predicted_tuple = tuple(sorted(set(predicted)))
    gate(actual == predicted_tuple, f"collision census d={d},e={e},f={f}")
    gate(len(actual) <= 2, f"bucket multiplicity d={d},e={e},f={f}")
    gate((len(actual) == 2) == (d == e == f),
         f"double collision iff AP d={d},e={e},f={f}")
    if not actual:
        collision_counts["none"] += 1
    elif d == e == f:
        collision_counts["common_step"] += 1
    elif d == e:
        collision_counts["lower"] += 1
    elif d == f:
        collision_counts["upper"] += 1
    else:
        collision_counts["central"] += 1
    packet_rows.append(f"{d},{e},{f}:{actual}")

gate(collision_offsets(4, 1, 2) == (), "no-collision hostile control")
gate(collision_offsets(2, 2, 1) == (2,), "lower collision control")
gate(collision_offsets(2, 1, 2) == (3,), "upper collision control")
gate(collision_offsets(3, 1, 2) == (3,), "central collision control")
gate(collision_offsets(2, 2, 2) == (2, 4), "common-step open control")


# Coarse but exact commutation atlas.  C/N records constant/nonconstant
# profiles.  This is a necessary superset: it deliberately ignores parity
# beyond the exact rule that a constant profile can occur only in an even
# nonnegative weight.  The proof closes every surviving superset family.
def profile_types(weight: int) -> tuple[str, ...]:
    kinds = ["N"]
    if weight == 0 or (weight > 0 and weight % 2 == 0):
        kinds.append("C")
    return tuple(kinds)


def commuting_possible(w1: int, k1: str, w2: int, k2: str) -> bool:
    if k1 == "C" and k2 == "C":
        return True
    if k1 == "C" and k2 == "N":
        return w1 == 0
    if k1 == "N" and k2 == "C":
        return w2 == 0
    return (w1 == 0 and w2 == 0) or w1 * w2 > 0


def feasible_statuses(weights: tuple[int, ...],
                      edges: tuple[tuple[int, int], ...]) -> list[tuple[str, ...]]:
    survivors = []
    for kinds in product(*(profile_types(weight) for weight in weights)):
        if all(commuting_possible(weights[i], kinds[i], weights[j], kinds[j])
               for i, j in edges):
            survivors.append(kinds)
    return survivors


ATLAS_WEIGHT = 64
ATLAS_GAP = 18
atlas_counts = defaultdict(int)
atlas_packet: list[str] = []

# Lower collision d=e, final gap f!=d.  The hub is G_2.
for d in range(1, ATLAS_GAP + 1):
    for f in range(1, ATLAS_GAP + 1):
        if f == d:
            continue
        for a in range(-ATLAS_WEIGHT, ATLAS_WEIGHT + 1):
            weights = (a, a + d, -3 - a - d, -3 - a, -3 - a + f)
            statuses = feasible_statuses(weights, ((0, 2), (1, 3), (0, 4), (1, 4)))
            if not statuses:
                continue
            gate((a, d, f) in ((-2, 2, 1), (-1, 1, 2)),
                 f"unexpected lower atlas survivor a={a},d={d},f={f}")
            gate(all(row == ("N", "C", "N", "N", "C") for row in statuses),
                 f"lower constant seam a={a},d={d},f={f}")
            atlas_counts["lower_boundary"] += len(statuses)
            atlas_packet.append(f"L:{a},{d},{f}:{statuses}")

# Upper collision d=f, first gap e!=d.  The hub is G_0.
for d in range(1, ATLAS_GAP + 1):
    for e in range(1, ATLAS_GAP + 1):
        if e == d:
            continue
        for a in range(-ATLAS_WEIGHT, ATLAS_WEIGHT + 1):
            weights = (a, a + d, -3 - a - d - e, -3 - a - d, -3 - a)
            statuses = feasible_statuses(weights, ((0, 2), (0, 3), (1, 2), (1, 4)))
            if not statuses:
                continue
            if a + d == 0:
                label = "upper_zero_Ftop"
            elif a == -3 and d in (1, 2):
                label = "upper_zero_Gtop"
            elif a == -2 and d == 1:
                label = "upper_synchronized"
            else:
                raise RuntimeError(
                    f"unexpected upper atlas survivor a={a},d={d},e={e}"
                )
            atlas_counts[label] += len(statuses)
            atlas_packet.append(f"U:{label}:{a},{d},{e}:{statuses}")

# Central collision d=e+f.  The hub is G_1.
for e in range(1, ATLAS_GAP + 1):
    for f in range(1, ATLAS_GAP + 1):
        d = e + f
        for a in range(-ATLAS_WEIGHT, ATLAS_WEIGHT + 1):
            weights = (a, a + d, -3 - a - d, -3 - a - f, -3 - a)
            statuses = feasible_statuses(weights, ((0, 2), (0, 3), (1, 3), (1, 4)))
            if not statuses:
                continue
            if a + d == 0:
                label = "central_zero_Ftop"
            elif a == -3 and d == 2:
                label = "central_zero_Gtop"
            elif weights[3] == 0:
                label = "central_constant_hub_to_2x2"
            else:
                raise RuntimeError(
                    f"unexpected central atlas survivor a={a},e={e},f={f}"
                )
            atlas_counts[label] += len(statuses)
            atlas_packet.append(f"C:{label}:{a},{e},{f}:{statuses}")


# Exact scalar exits for the only nontrivial boundary packets.
s = sp.symbols("s")
p = sp.Function("p")(s)
q = sp.Function("q")(s)
lam, mu = sp.symbols("lam mu", nonzero=True)


def W(a: int, first: sp.Expr, b: int, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        a * first * sp.diff(second, s)
        - b * sp.diff(first, s) * second
    )


# Both lower-boundary scalar profiles have negative weight and hence contain
# s^2-1.  Their bracket therefore vanishes at either arm s=+/-1.
arm = s**2 - 1
lower_scalar = W(-2, arm * p, -1, arm * q)
gate(sp.simplify(lower_scalar.subs(s, 1)) == 0,
     "lower boundary scalar vanishes at positive arm")
gate(sp.simplify(lower_scalar.subs(s, -1)) == 0,
     "lower boundary scalar vanishes at negative arm")

# In the upper synchronized family the hub forces F_-2 proportional to
# F_-1^2; the two other commuting edges make the scalar summands commute too.
upper_scalar = W(-2, p**2, -1, mu * p) + W(-1, p, -2, lam * p**2)
gate(sp.simplify(upper_scalar) == 0, "upper synchronized scalar is zero")
for e in range(1, 9):
    if e == 1:
        continue
    hub = p ** (e + 2)
    gate(W(-2, p**2, -2 - e, hub) == 0,
         f"upper hub commutes with low spoke e={e}")
    gate(W(-1, p, -2 - e, hub) == 0,
         f"upper hub commutes with high spoke e={e}")

# On the central constant-hub sheet, deleting that global constant leaves
# exactly an equal-gap 2x2 support: both middle sums are -3.
for e in range(3, 13):
    for f in range(1, 9):
        f_weights = (-f - 3, e - 3)
        g_weights = (-e, f)
        gate(f_weights[1] - f_weights[0] == e + f,
             f"central F gap e={e},f={f}")
        gate(g_weights[1] - g_weights[0] == e + f,
             f"central G gap e={e},f={f}")
        gate(f_weights[0] + g_weights[1] == -3,
             f"central first scalar address e={e},f={f}")
        gate(f_weights[1] + g_weights[0] == -3,
             f"central second scalar address e={e},f={f}")

semantic_rows = (
    "surface=THM3783;r^3g=z^2-r/4;Euler_bracket_degree_plus3",
    "supports=F:{a,a+d};G:{b,b+e,b+e+f};d,e,f_positive",
    "collisions=only_d=e_or_d=f_or_d=e+f;double_iff_d=e=f",
    "hub_lemma=nonconstant_hub_sign_synchronization;constant_and_weight0_typed",
    "lower_irregular=two_zero_weight_boundaries;negative_arm_vanishing",
    "upper_irregular=zero_boundaries_or_-2_-1_synchronization;scalar_zero",
    "central_irregular=single_homogeneous_exit_or_constant_hub_to_2x2",
    "conclusion=all_irregular_2x3_empty;only_common_step_AP_remains_open",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()
packet = hashlib.sha256("\n".join(packet_rows + atlas_packet).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3786-quadratic-etale-surface-irregular-two-by-three-support-no-go")
print("field=algebraically_closed_characteristic_zero")
print("surface=r^3*g=z^2-r/4;weights_r_z_g=2,1,-4;bracket_degree=3")
print("supports=F:{a,a+d};G:{b,b+e,b+e+f};d,e,f>0")
print("collision_relations=d=e_or_d=f_or_d=e+f;double_collision_iff_d=e=f")
print("collision_counts=" + ",".join(f"{k}:{collision_counts[k]}" for k in sorted(collision_counts)))
print("atlas_range=weights[-64,64];gaps[1,18]")
print("atlas_counts=" + ",".join(f"{k}:{atlas_counts[k]}" for k in sorted(atlas_counts)))
print("hub_lemma=nonconstant_sign_sync_with_zero_constant_exception;constant_hub_typed")
print("lower_exit=two_weight_zero_boundaries_then_negative_arm_vanishing")
print("upper_exit=homogeneous_boundaries_or_-2_-1_full_synchronization")
print("central_exit=homogeneous_boundaries_or_constant_hub_deletes_to_2x2")
print("open_boundary=common_step_AP_d=e=f_only")
print(f"packet_sha256={packet}")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
