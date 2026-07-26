#!/usr/bin/env python3
"""Exact referee for THM-2459.

Dependency-free; every truth-bearing check uses ``require`` so optimized
Python performs the same verification.
"""

from fractions import Fraction as F
from itertools import combinations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def vadd(*vectors):
    if not vectors:
        return ()
    return tuple(sum((v[i] for v in vectors), F(0)) for i in range(len(vectors[0])))


def vscale(a, vector):
    return tuple(a * x for x in vector)


def vsub(left, right):
    return tuple(x - y for x, y in zip(left, right))


def norm2(vector):
    return sum((x * x for x in vector), F(0))


def qpair(table):
    """Mean-zero projection on two coordinates."""
    x, y = table
    return ((x - y) / 2, (y - x) / 2)


def subsets(items):
    items = tuple(items)
    for size in range(len(items) + 1):
        for part in combinations(items, size):
            yield frozenset(part)


# The exact Hilbert identity and its sharp collinear realization.
identity_cases = 0
sharp_nonnegative_cases = 0
for r in range(1, 128):
    # A noncollinear deterministic control for the algebraic identity.
    generic_z0 = (F(r, 3), F(-2, r + 1), F(5, 7))
    generic_outside = [
        (F(j, r + 2), F(r - j, 2 * r + 1), F((-1) ** j, j + 1))
        for j in range(1, r + 1)
    ]
    generic_candidates = [generic_z0] + [
        vadd(generic_z0, qj) for qj in generic_outside
    ]
    generic_total = vadd(generic_z0, *generic_outside)
    generic_rhs = vsub(
        vadd(*generic_candidates[1:]), vscale(F(r - 1), generic_z0)
    )
    require(generic_rhs == generic_total, f"generic Hilbert identity failed at r={r}")

    e = (F(1), F(-1))
    z0 = vscale(F(-1), e)
    outside = [vscale(F(2), e) for _ in range(r)]
    candidates = [z0] + [vadd(z0, qj) for qj in outside]
    qtotal = vadd(z0, *outside)
    rhs = vsub(vadd(*candidates[1:]), vscale(F(r - 1), z0))
    require(rhs == qtotal, f"Hilbert identity failed at r={r}")
    require(qtotal == vscale(F(2 * r - 1), e), f"sharp total failed at r={r}")
    require(max(norm2(z) for z in candidates) == norm2(e), f"sharp max failed at r={r}")
    require(norm2(qtotal) == F((2 * r - 1) ** 2) * norm2(e), f"sharp ratio failed at r={r}")
    identity_cases += 1

    # A nonnegative two-entry realization.  Constant parts are killed by Q.
    eps = F(1, 8 * r)
    base = (F(0), 2 * eps)  # Q(base)=-eps*(1,-1)
    extensions = [(4 * eps, F(0)) for _ in range(r)]  # Q=2eps*(1,-1)
    require(qpair(base) == vscale(-eps, e), f"base realization failed at r={r}")
    require(all(qpair(t) == vscale(2 * eps, e) for t in extensions),
            f"extension realization failed at r={r}")
    coordinate_total = vadd(base, *extensions)
    require(all(F(0) <= x <= F(1) for x in coordinate_total),
            f"nonnegative capacity failed at r={r}")
    sharp_nonnegative_cases += 1


# Canonical N=128 constants.
N = 128
root_entries = N * N
max_outside = 126
drift_linear_denominator = 2 * max_outside - 1
drift_energy_denominator = drift_linear_denominator ** 2
root_mass_denominator = root_entries
root_max_denominator = root_mass_denominator * 2028
root_energy_denominator = root_mass_denominator ** 2 * 342732
require(drift_linear_denominator == 251, "251 drift factor")
require(drift_energy_denominator == 63001, "63001 drift energy")
require(root_max_denominator == 33226752, "root max denominator")
require(root_energy_denominator == 92001420705792, "root energy denominator")


# Four atoms are necessary.
atoms4 = ("k", "u", "v", "j")
tables4 = {
    "k": (F(1, 4), F(0)),
    "u": (F(0), F(1, 4)),
    "v": (F(1, 8), F(1, 8)),
    "j": (F(1, 4), F(0)),
}
total4 = vadd(*(tables4[a] for a in atoms4))
require(all(F(0) <= x <= F(1) for x in total4), "four-atom table capacity")
service_and_selected4 = []
working4 = []
for I in subsets(atoms4):
    service = "u" in I and "v" in I  # sole edge u->v
    selected = "k" in I
    drift = qpair(vadd(*(tables4[a] for a in I))) if I else (F(0), F(0))
    if service and selected:
        service_and_selected4.append((I, drift))
        if norm2(drift) > 0:
            working4.append(I)
require(service_and_selected4, "four-atom service candidates")
require(all(norm2(q) == 0 for I, q in service_and_selected4 if len(I) <= 3),
        "no at-most-three solution")
require(working4 == [frozenset(atoms4)], "unique four-atom solution")


# Aggregate drift is necessary: the minimal two-atom hostile.
atoms2 = ("k", "v")
tables2 = {
    "k": (F(1, 2), F(0)),
    "v": (F(0), F(1, 2)),
}
require(norm2(qpair(tables2["k"])) > 0, "selected two-atom drift")
require(qpair(vadd(tables2["k"], tables2["v"])) == (F(0), F(0)),
        "aggregate-zero hostile")
working2 = []
for I in subsets(atoms2):
    if "k" in I and "v" in I:  # sole edge k->v
        drift = qpair(vadd(*(tables2[a] for a in I)))
        if norm2(drift) > 0:
            working2.append(I)
require(not working2, "aggregate-zero hostile must have no solution")


# Optional observer lemma: one service cell plus one drift cell.
cell_q = ((F(0), F(0)), (F(3), F(-3)), (F(-1), F(1)))
cell_s = (F(2), F(0), F(0))
gate = (0, 1)
gate_q = vadd(*(cell_q[i] for i in gate))
gate_s = sum((cell_s[i] for i in gate), F(0))
require(norm2(gate_q) > 0 and gate_s > 0, "two-cell observer selector")
pointwise_zero = tuple(vadd(qpair(tables2["k"]), qpair(tables2["v"])) for _ in range(3))
require(all(q == (F(0), F(0)) for q in pointwise_zero),
        "pointwise root-constant hostile")


print("THM2459 FOUR-ATOM DRIFT/SERVICE AUDIT")
print(f"hilbert_identity_cases={identity_cases}")
print(f"sharp_nonnegative_realizations={sharp_nonnegative_cases}")
print(f"N={N} directed_entries={root_entries}")
print(f"drift_linear_denominator={drift_linear_denominator}")
print(f"drift_energy_denominator={drift_energy_denominator}")
print(f"root_mass_denominator={root_mass_denominator}")
print(f"root_max_denominator={root_max_denominator}")
print(f"root_energy_denominator={root_energy_denominator}")
print(f"four_atom_service_unions={len(service_and_selected4)}")
print(f"four_atom_working_unions={len(working4)} min_size={min(map(len, working4))}")
print(f"two_atom_working_unions={len(working2)}")
print(f"observer_cells={len(gate)} observer_service={gate_s} observer_drift_norm2={norm2(gate_q)}")
print("four_atom_hostile=PASS")
print("aggregate_sidecar_hostile=PASS")
print("optimized_require_path=PASS")
