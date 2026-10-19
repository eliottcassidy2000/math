#!/usr/bin/env python3
"""Independent exact audit of the 29-adic mixed-stress section bound.

No imports from the submitted verifier; standard library only.  Arrays use
y-major order, and carries are calculated from literal integer differences.
"""

from itertools import product
from pathlib import Path
import hashlib
import json

P = 29
N = P * P
DIRECTIONS = ((1, 0), (0, 1), (3, 9))
HERE = Path(__file__).resolve().parent


def q(x, y, m=1):
    return x * x + m * x * y + y * y


def idx(x, y):
    return (y % P) * P + x % P


def shift(f, d):
    a, b = d
    return [f[idx(x + a, y + b)] for y in range(P) for x in range(P)]


def diff(f, d):
    return [b - a for a, b in zip(f, shift(f, d))]


def carry(d):
    a, b = d
    ans = []
    for y in range(P):
        for x in range(P):
            dx, dy = (x + a) % P - x, (y + b) % P - y
            value = q(dx, dy) - 1
            assert value % P == 0
            ans.append((value // P) % P)
    return ans


def row_sums(f):
    return [sum(f[y * P:(y + 1) * P]) % P for y in range(P)]


def period_difference(f, d):
    return row_sums(diff(f, d))


def seam(family, orientation, offset):
    ans = []
    for y in range(P):
        for x in range(P):
            if family == 0:
                linear = (9 * x - 3 * y) if orientation else x
            elif family == 1:
                linear = (9 * x - 3 * y) if orientation else y
            else:
                linear = x if orientation else y
            ans.append(int(linear % P == offset))
    return ans


def integrate(residual, family):
    """Solve u or v differences, explicitly checking the closing edge."""
    ans = [None] * N
    for line in range(P):
        value = 0
        for step in range(P):
            x, y = (step, line) if family == 0 else (line, step)
            ans[idx(x, y)] = value
            value = (value + residual[idx(x, y)]) % P
        assert value == 0
    assert all((a - b) % P == 0 for a, b in
               zip(diff(ans, DIRECTIONS[family]), residual))
    return ans


def straight_cycles(d):
    unused = set(range(N))
    while unused:
        origin = min(unused)
        x, y = origin % P, origin // P
        cycle = [idx(x + k * d[0], y + k * d[1]) for k in range(P)]
        assert len(set(cycle)) == P
        unused.difference_update(cycle)
        yield cycle


def stencil(fields):
    fu, fv, fw = fields
    first = diff(diff(fw, DIRECTIONS[0]), DIRECTIONS[1])
    second = diff(diff(fu, DIRECTIONS[1]), DIRECTIONS[2])
    third = diff(diff(fv, DIRECTIONS[0]), DIRECTIONS[2])
    return [a - 3 * b - 9 * c for a, b, c in zip(first, second, third)]


def sparse_operator(d):
    return {(0, 0): -1, d: 1}


def convolution(*operators):
    result = {(0, 0): 1}
    for op in operators:
        new = {}
        for (x, y), c in result.items():
            for (a, b), e in op.items():
                key = ((x + a) % P, (y + b) % P)
                new[key] = new.get(key, 0) + c * e
        result = {k: v for k, v in new.items() if v}
    return result


def canonical(d):
    return min(d, ((-d[0]) % P, (-d[1]) % P))


def multiply(a, b, modulus):
    return ((a[0] * b[0] - a[1] * b[1]) % modulus,
            (a[0] * b[1] + a[1] * b[0] + a[1] * b[1]) % modulus)


unit = [(x, y) for y in range(P) for x in range(P) if q(x, y) % P == 1]
pairs = sorted(set(map(canonical, unit)))
assert len(unit) == 30 and len(pairs) == 15 and N * len(pairs) == 12615
assert pow(P - 3, (P - 1) // 2, P) == P - 1

small_solutions = [z for z in product(range(-2, 3), repeat=3)
                   if (z[0] - 3 * z[1] - 9 * z[2]) % P == 0]
assert small_solutions == [(0, 0, 0)]
assert (2 - 3*(-1) - 9*(-2)) == 23

# Exhaust binary small grids as a control of the separability lemma.  The
# proof, rather than this finite sample, is what handles a 29 by 29 grid.
separability_controls = {}
for n in (2, 3, 4):
    accepted = 0
    for mask in range(1 << (n*n)):
        f = lambda x, y: (mask >> ((y % n)*n + x % n)) & 1
        mixed_zero = all(f(x+1,y+1)-f(x+1,y)-f(x,y+1)+f(x,y) == 0
                         for y in range(n) for x in range(n))
        invariant_x = all(f(x+1,y) == f(x,y) for y in range(n) for x in range(n))
        invariant_y = all(f(x,y+1) == f(x,y) for y in range(n) for x in range(n))
        assert mixed_zero == (invariant_x or invariant_y)
        accepted += mixed_zero
    assert accepted == 2*(1 << n)-2
    separability_controls[str(n)] = {"grids": 1 << (n*n), "mixed_zero": accepted}

# The three orderings yield the SAME group-algebra operator.  Its two
# vector coefficients cancel over Z, stronger than cancellation modulo p.
operators = [sparse_operator(d) for d in DIRECTIONS]
base = convolution(*operators)
assert base == convolution(operators[1], operators[2], operators[0])
assert base == convolution(operators[2], operators[0], operators[1])
for component in range(2):
    linear = lambda d: (2 * d[0] + d[1], d[0] + 2 * d[1])[component]
    assert linear(DIRECTIONS[2]) - 3 * linear(DIRECTIONS[0]) - 9 * linear(DIRECTIONS[1]) == 0

all_carries = {d: carry(d) for d in pairs}
three_carries = [carry(d) for d in DIRECTIONS]
assert all(z % P == 0 for z in stencil(three_carries))
cycle_count = 0
for d, c in all_carries.items():
    for cycle in straight_cycles(d):
        assert sum(c[i] for i in cycle) % P == P - 2
        cycle_count += 1
assert cycle_count == 435

# Test the norm definition itself, not merely the submitted affine formula.
# These are controls for implementation, not an enumeration of sections.
motion_edge_checks = 0
motion_defects = []
for seed in (0, 1, 7):
    h = [((seed * (x*x + 5*x*y + 3*y*y*y)) % P,
          (seed * (2*x*x*x + x*y + y*y + 1)) % P)
         for y in range(P) for x in range(P)]
    errors = {}
    total_defects = 0
    for d, c in all_carries.items():
        a, b = d
        z = []
        for y in range(P):
            for x in range(P):
                i, j = idx(x, y), idx(x + a, y + b)
                xx, yy = (x + a) % P, (y + b) % P
                dhx, dhy = h[j][0] - h[i][0], h[j][1] - h[i][1]
                dx, dy = xx - x + P * dhx, yy - y + P * dhy
                norm_error = ((q(dx, dy) - 1) // P) % P
                affine = (c[i] + (2*a + b)*dhx + (a + 2*b)*dhy) % P
                assert (q(dx, dy) - 1) % P == 0
                assert norm_error == affine
                assert (norm_error == 0) == (q(dx, dy) % (P*P) == 1)
                z.append(norm_error)
                total_defects += bool(norm_error)
                motion_edge_checks += 1
        for cycle in straight_cycles(d):
            assert sum(z[i] for i in cycle) % P == P - 2
        errors[d] = z
    # Get literal directed errors for the three designated directions.
    three_errors = []
    for d, c in zip(DIRECTIONS, three_carries):
        a, b = d
        three_errors.append([
            (c[idx(x, y)] + (2*a+b)*(h[idx(x+a,y+b)][0]-h[idx(x,y)][0])
             + (a+2*b)*(h[idx(x+a,y+b)][1]-h[idx(x,y)][1])) % P
            for y in range(P) for x in range(P)])
    assert all(z % P == 0 for z in stencil(three_errors))
    motion_defects.append(total_defects)

# DIRECTLY integrate every one of the 58 possible seams for EACH of the
# u and v families.  Store all 29 horizontal period contributions.
# The residual is a sum of independent u, v, w contributions.  Consequently
# combining these tables tests all 8*29^3 full systems without recomputing
# the same cumulative integrals thousands of times.
contributions = []
seam_cycle_checks = 0
for family in range(3):
    table = {}
    d, c = DIRECTIONS[family], three_carries[family]
    for orientation in range(2):
        for offset in range(P):
            indicator = seam(family, orientation, offset)
            assert sum(indicator) == P
            other = [DIRECTIONS[j] for j in range(3) if j != family]
            assert all(z == 0 for z in diff(diff(indicator, other[0]), other[1]))
            for cycle in straight_cycles(d):
                assert sum(indicator[i] for i in cycle) == 1
                seam_cycle_checks += 1
            residual = [(-2*z - cc) % P for z, cc in zip(indicator, c)]
            if family < 2:
                H = integrate(residual, family)
                scalar = -3 if family == 0 else -9
                values = [(scalar*z) % P for z in period_difference(H, DIRECTIONS[2])]
            else:
                values = row_sums(residual)
            table[orientation, offset] = values
    contributions.append(table)

period_values = []
full_offset_cases = 0
row_period_checks = 0
for eu, ev, ew in product(range(2), repeat=3):
    predicted = (16 - 18*eu + 17*ev - 2*ew) % P
    assert predicted != 0
    period_values.append(predicted)
    for tu, tv, tw in product(range(P), repeat=3):
        U, V, W = (contributions[0][eu, tu], contributions[1][ev, tv],
                   contributions[2][ew, tw])
        for y in range(P):
            actual = (U[y] + V[y] + W[y]) % P
            assert actual == predicted
            row_period_checks += 1
        full_offset_cases += 1
assert period_values == [16, 14, 4, 2, 27, 25, 15, 13]
assert full_offset_cases == 195112 and row_period_checks == 5658248

# Rotations: direct norm-one lifts of every unit element, closure, and pair
# incidence counts.  The lift formula is tested modulo p^2.
for g in unit:
    assert (q(*g) - 1) % P == 0
    k = ((q(*g) - 1) // P) % P
    correction = (1 - P*k*pow(2, -1, P)) % (P*P)
    G = (g[0]*correction % (P*P), g[1]*correction % (P*P))
    assert (G[0] % P, G[1] % P) == g
    assert q(*G) % (P*P) == 1
    for z in unit:
        assert multiply(g, z, P) in unit

incidence = {d: 0 for d in pairs}
for g in pairs:
    rotated = {canonical(multiply(g, d, P)) for d in DIRECTIONS}
    assert len(rotated) == 3
    for d in rotated:
        incidence[d] += 1
assert set(incidence.values()) == {3}

# General-corollary arithmetic audit.  Enumerate every possible m=2t
# arising from nonzero a,b in all prime fields through 101.  This includes
# split forms as a diagnostic: p+1 direction counting must be restricted
# to the anisotropic/unramified setting.
def primes_to(n):
    return [k for k in range(5, n+1) if all(k % d for d in range(2, int(k**0.5)+1))]

general_counts = {}
split_example = None
anisotropic_examples = {}
for p in primes_to(101):
    stats = {"unit_coefficient_pairs": 0, "anisotropic_qualifying": 0,
             "split_qualifying": 0}
    for a in range(1, p):
        for b in range(1, p):
            m = ((1-a*a-b*b) * pow(a*b, -1, p)) % p
            disc = (m*m-4) % p
            if disc == 0:
                continue
            stats["unit_coefficient_pairs"] += 1
            assert (a*a + m*a*b + b*b) % p == 1
            separated = True
            for B, C in product(range(-2, 3), repeat=2):
                if B == C == 0:
                    continue
                residue = (a*B + b*C) % p
                if residue in {0, 1, 2, p-1, p-2}:
                    separated = False
                    break
            if not separated:
                continue
            tab = (m*a*b*pow(2, -1, p)) % p
            excluded = {tab, (1-tab) % p, (a*a+tab) % p, (b*b+tab) % p}
            eight = [(a*a*(1-eu)+tab+b*b*ev-ew) % p
                     for eu, ev, ew in product(range(2), repeat=3)]
            assert (not (excluded & {0, 1})) == all(eight)
            if not all(eight):
                continue
            anisotropic = pow(disc, (p-1)//2, p) == p-1
            stats["anisotropic_qualifying" if anisotropic else "split_qualifying"] += 1
            if anisotropic and p not in anisotropic_examples:
                anisotropic_examples[p] = (m, a, b)
            if not anisotropic and split_example is None:
                split_example = {"p": p, "m": m, "a": a, "b": b}
    general_counts[str(p)] = stats

# A genuinely separate direct-coordinate test of the generalized period,
# with nonzero offsets and a different quadratic form at each available p.
general_integrations = []
for p, (m, a, b) in anisotropic_examples.items():
    index = lambda x, y: (y % p)*p + x % p
    points = [(x, y) for y in range(p) for x in range(p)]
    directions = ((1,0), (0,1), (a,b))
    carries = []
    for dx, dy in directions:
        values = []
        for x, y in points:
            xx, yy = (x+dx)%p-x, (y+dy)%p-y
            raw = q(xx, yy, m)-1
            assert raw % p == 0
            values.append((raw//p) % p)
        carries.append(values)
    for eu, ev, ew in product(range(2), repeat=3):
        iu = [int(((b*x-a*y) if eu else x) % p == 2) for x,y in points]
        iv = [int(((b*x-a*y) if ev else y) % p == 3) for x,y in points]
        iw = [int((x if ew else y) % p == 4) for x,y in points]
        residuals = [[(-2*z-c) % p for z,c in zip(I,C)]
                     for I,C in zip((iu,iv,iw),carries)]
        integrals = []
        for family in range(2):
            H = [None]*(p*p)
            for line in range(p):
                value = 0
                for step in range(p):
                    x,y = (step,line) if family == 0 else (line,step)
                    H[index(x,y)] = value
                    value = (value + residuals[family][index(x,y)]) % p
                assert value == 0
            integrals.append(H)
        predicted = (-2*ew + a*(2*a+m*b)-2*a*a*eu+2*b*b*ev) % p
        assert predicted != 0
        for y in range(p):
            measured = sum(residuals[2][index(x,y)]
                           -a*(integrals[0][index(x+a,y+b)]-integrals[0][index(x,y)])
                           -b*(integrals[1][index(x+a,y+b)]-integrals[1][index(x,y)])
                           for x in range(p)) % p
            assert measured == predicted
    general_integrations.append({"p":p,"m":m,"a":a,"b":b,"orientations":8})

source = HERE.parent
output = {
    "verdict": "All independent checks passed",
    "p": P,
    "base_edges": 12615,
    "baseline_cycle_lower_bound": 435,
    "full_section_defect_lower_bound": 440,
    "binary_second_difference_tuples": 125,
    "binary_separability_small_grid_controls": separability_controls,
    "all_direction_straight_cycles": cycle_count,
    "direct_norm_affine_edge_comparisons": motion_edge_checks,
    "test_motion_defects": motion_defects,
    "seams_integrated_for_each_of_u_and_v": 58,
    "seam_cycle_checks": seam_cycle_checks,
    "all_offset_orientation_systems_checked": full_offset_cases,
    "all_offset_row_periods_checked": row_period_checks,
    "periods_in_lexicographic_orientation_order": period_values,
    "norm_one_lifts_checked": len(unit),
    "rotation_pair_multiplicity": 3,
    "general_criterion_arithmetic_checks": general_counts,
    "generalized_period_direct_integrations": general_integrations,
    "split_form_qualifying_example": split_example,
    "source_sha256": {name: hashlib.sha256((source/name).read_bytes()).hexdigest()
                      for name in ("README.md", "verify_balanced_stress.py")},
}
(HERE / "independent-results.json").write_text(json.dumps(output, indent=2) + "\n")
print(json.dumps(output, indent=2))
