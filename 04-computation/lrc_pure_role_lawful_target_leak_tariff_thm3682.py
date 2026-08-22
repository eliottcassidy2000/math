"""Exact controls for THM-3682's role-to-target leak tariff."""

from fractions import Fraction
import hashlib
import json


P = 13
GATES = 0


def mean(values):
    return sum(values, Fraction(0)) / len(values)


def centered_energy(values):
    """Normalized DFT energy off the zero character, by Parseval."""
    average = mean(values)
    return sum((value - average) ** 2 for value in values) / len(values)


def inner_centered(left, right):
    left_average = mean(left)
    right_average = mean(right)
    return sum(
        (a - left_average) * (b - right_average)
        for a, b in zip(left, right)
    ) / len(left)


def gate(condition, label):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


# D=(-1/14,1/14) on the circle.  Its overlap with D-s/13 has length
# max(0,1/7-min(s,13-s)/13).  These are exact away from null endpoints.
role = []
for shift in range(P):
    circular_distance = Fraction(min(shift, P - shift), P)
    role.append(max(Fraction(0), Fraction(1, 7) - circular_distance))

expected_role = [Fraction(0) for _ in range(P)]
expected_role[0] = Fraction(1, 7)
expected_role[1] = expected_role[-1] = Fraction(6, 91)
gate(role == expected_role, "danger overlap profile")

rho = Fraction(1, 7)
role_energy = centered_energy(role)
role_floor = Fraction(121, 2028) * rho * rho
gate(role_energy >= role_floor, "THM-2362 floor")

# Covary the weight by the same translation as the named danger factor.
# The lawful diagonal profile is then identically rho.  Consequently its
# leak is exactly the negative centered role profile and saturates the sharp
# cancellation threshold E_leak=E_role.
target = [rho for _ in range(P)]
leak = [c - r for c, r in zip(target, role)]
target_energy = centered_energy(target)
leak_energy = centered_energy(leak)
gate(target_energy == 0, "covariant target constancy")
gate(leak_energy == role_energy, "sharp leak threshold")
gate(
    target_energy
    == role_energy + leak_energy + 2 * inner_centered(role, leak),
    "polarization identity",
)
gate(inner_centered(role, leak) == -role_energy, "hostile antiparallelism")

# Deterministic rational controls for the reverse-triangle tariff.  Avoid
# square roots by checking the equivalent Cauchy/polarization inequality
# E_C >= E_R+E_L-2 sqrt(E_R E_L), squared when the RHS is positive.
controls = []
for seed in range(1, 65):
    r = [Fraction(((seed + 3) * (s + 2) * (s + seed)) % 37 - 18, 19) for s in range(P)]
    ell = [Fraction(((seed + 11) * (s + 5) * (2 * s + seed)) % 41 - 20, 23) for s in range(P)]
    c = [a + b for a, b in zip(r, ell)]
    er = centered_energy(r)
    el = centered_energy(ell)
    ec = centered_energy(c)
    cross = inner_centered(r, ell)
    gate(ec == er + el + 2 * cross, f"control polarization {seed}")
    gate(cross * cross <= er * el, f"control Cauchy {seed}")
    if ec < er + el:
        gate((er + el - ec) ** 2 <= 4 * er * el, f"control tariff {seed}")
    controls.append((str(er), str(el), str(ec)))

# A one-row variation inside a p by p successor table pays at least E_C/p
# in the globally normalized two-dimensional variance.  Equality occurs when
# the other p-1 rows are constant at the varying row's own mean.
row = role
row_average = mean(row)
table = [row] + [[row_average for _ in range(P)] for _ in range(P - 1)]
flat = [value for table_row in table for value in table_row]
table_variance = centered_energy(flat)
gate(table_variance == role_energy / P, "sharp row-to-table invoice")

# Exact p=13 tariff denominators after THM-3674.
drift_factor = Fraction(1, P * P * (P - 1))
deep_target_factor = Fraction(1, P * P * P * (P - 1))
gate(drift_factor == Fraction(1, 2028), "drift factor")
gate(deep_target_factor == Fraction(1, 26364), "deep target factor")

semantic = {
    "p": P,
    "role": [str(value) for value in role],
    "rho": str(rho),
    "role_energy": str(role_energy),
    "role_floor": str(role_floor),
    "target_energy": str(target_energy),
    "leak_energy": str(leak_energy),
    "row_to_table_factor": f"1/{P}",
    "drift_factor_from_row_energy": str(drift_factor),
    "deep_target_factor_from_row_energy": str(deep_target_factor),
    "control_count": len(controls),
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

print(f"danger_role={','.join(str(value) for value in role)}")
print(f"role_energy={role_energy};THM2362_floor={role_floor}")
print(f"covariant_target_energy={target_energy};leak_energy={leak_energy}")
print(f"row_to_table={table_variance};drift_factor={drift_factor};deep_target_factor={deep_target_factor}")
print(f"semantic_sha256={semantic_sha256}")
print(f"RESULT=PASS;gates={GATES}")
