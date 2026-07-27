#!/usr/bin/env python3
"""Exact referee for THM-2536.

Checks the two oriented selectors on the 23-ray anchored deep-comb cone,
their lossless start/end flow synthesis, target-line cancellation, the
constant-singleton zero-target hostile, the centered toothpick lift, and the
exact 169-cell positive control on the canonical typed non-cover word used by
the independent THM-2334 twist-variance computation.
"""

from bisect import bisect_right
from collections import Counter
from fractions import Fraction
from math import gcd


P = 13
Q = 7
checks = 0


def require(condition, label):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError("FAILED: " + label)


# ---------------------------------------------------------------------------
# The anchored 23-mask cone and its two honest selector marginals
# ---------------------------------------------------------------------------


def bit(mask, a):
    return (mask >> (a % P)) & 1


path_masks = [(1 << j, "singleton", j) for j in range(1, P)]
path_masks += [((1 << j) | (1 << (j + 1)), "pair", j) for j in range(1, P - 1)]

for mask, kind, j in path_masks:
    start = [bit(mask, a) * (1 - bit(mask, a - 1)) for a in range(P)]
    terminal = [bit(mask, a) * (1 - bit(mask, a + 1)) for a in range(P)]
    require(sum(start) == sum(terminal) == 1, "one start and one terminal")
    require(start[0] == terminal[0] == 0, "target anchor")
    require(start[j] == 1, "start coordinate")
    require(terminal[j if kind == "singleton" else j + 1] == 1, "terminal coordinate")

# Deliberately irregular integral cone coordinates.
alpha = [0] + [j * j + 2 * j + 1 for j in range(1, P)]
beta = [0] + [3 * j * j + 1 for j in range(1, P - 1)] + [0]
gamma_plus = [0] * P
gamma_minus = [0] * P
for j in range(1, P):
    gamma_plus[j] = alpha[j] + beta[j]
    gamma_minus[j] = alpha[j] + beta[j - 1]

require(sum(gamma_plus) == sum(gamma_minus), "start and terminal masses agree")
running = 0
beta_recovered = [0] * P
for j in range(1, P):
    running += gamma_plus[j] - gamma_minus[j]
    beta_recovered[j] = running
require(beta_recovered == beta, "anchored divergence recovers every pair mass")
alpha_recovered = [0] + [gamma_plus[j] - beta_recovered[j] for j in range(1, P)]
require(alpha_recovered == alpha, "start and pair masses recover every singleton mass")
delta_selector = [gamma_plus[j] - gamma_minus[j] for j in range(P)]
require(sum(delta_selector) == 0, "opposite-selector divergence has zero mass")
require(any(delta_selector) == any(beta), "anchored divergence detects pair mass")


def edge(u, v):
    matrix = [[0] * P for _ in range(P)]
    matrix[u][v] = 1
    matrix[v][u] = -1
    return matrix


def add_matrix(*matrices):
    return [[sum(matrix[u][v] for matrix in matrices) for v in range(P)]
            for u in range(P)]


def scale_matrix(matrix, scalar):
    return [[scalar * matrix[u][v] for v in range(P)] for u in range(P)]


def divergence(matrix):
    return [sum(row) for row in matrix]


for j in range(1, P - 1):
    anchored_path = add_matrix(edge(j, 0), edge(0, j + 1))
    direct_pair = edge(j, j + 1)
    triangle = add_matrix(anchored_path, scale_matrix(direct_pair, -1))
    require(divergence(anchored_path) == divergence(direct_pair), "path/direct divergence")
    require(divergence(triangle) == [0] * P, "target triangle is circulation")
    require(triangle[j][0] == triangle[0][j + 1] == triangle[j + 1][j] == 1,
            "oriented target triangle")


# ---------------------------------------------------------------------------
# Exact finite target transform and the diagonal-anchor line cancellation
# ---------------------------------------------------------------------------


MOD = 79
ZETA = 8
require(pow(ZETA, P, MOD) == 1 and ZETA != 1, "primitive thirteenth root in F_79")
require(all(pow(ZETA, k, MOD) != 1 for k in range(1, P)), "root order exactly thirteen")


def zpow(k):
    return pow(ZETA, k % P, MOD)


require(all(sum(delta_selector[a] * zpow(k * a) for a in range(P)) % MOD != 0
            for k in range(1, P)), "pair divergence retains all root colours")


# A nonuniform family of lawful anchored path-cone coordinates.  Its absolute
# tensor is Gamma(r,s,t)=gamma_(s,t)(r-t).  At fixed target (b,q), the colour-k
# entry has exponent k*a+b*s+q*t, because the endpoint colour is q-k.
family = {}
for s in range(P):
    for t in range(P):
        gp = [0] * P
        for j in range(1, P):
            aj = 1 + (j + 2 * s + 3 * t) % 11
            bj = 0 if j == P - 1 else (j * (s + 1) + t) % 7
            gp[j] = aj + bj
        family[(s, t)] = gp

for b in range(P):
    for q in range(P):
        line = 0
        star_line = 0
        star_colour_zero = 0
        for k in range(P):
            coeff = 0
            star_coeff = 0
            for s in range(P):
                for t in range(P):
                    gp = family[(s, t)]
                    rho = [-sum(gp)] + gp[1:]
                    for a in range(P):
                        coeff += gp[a] * zpow(k * a + b * s + q * t)
                        star_coeff += rho[a] * zpow(k * a + b * s + q * t)
            line = (line + coeff) % MOD
            star_line = (star_line + star_coeff) % MOD
            if k == 0:
                star_colour_zero = star_coeff % MOD
        require(line == 0, "anchored target-line cancellation")
        mass_hat = sum(sum(family[(s, t)]) * zpow(b * s + q * t)
                       for s in range(P) for t in range(P)) % MOD
        require(star_colour_zero == 0, "star boundary has zero root colour zero")
        require(star_line == (-P * mass_hat) % MOD,
                "unnormalized star line is minus thirteen times mass transform")

# Exact zero-target hostile: the same relative singleton in all 169 cells.
# It has every nonzero root colour but no nonzero target character.
for k in range(P):
    for b in range(P):
        for q in range(P):
            coeff = 0
            for s in range(P):
                for t in range(P):
                    coeff += zpow(k + b * s + q * t)
            coeff %= MOD
            if b == q == 0:
                require(coeff != 0, "singleton hostile retains root colour")
            else:
                require(coeff == 0, "singleton hostile has zero target colour")


# ---------------------------------------------------------------------------
# Selector/H-diagonal incomparability and the THM-2367 full-packet hostile
# ---------------------------------------------------------------------------


def path_coordinates(alpha_weights=None, beta_weights=None):
    alpha_weights = dict(alpha_weights or {})
    beta_weights = dict(beta_weights or {})
    gp = [Fraction(0) for _ in range(P)]
    diagonal = [Fraction(0) for _ in range(P)]
    for j in range(1, P):
        gp[j] = alpha_weights.get(j, 0) + beta_weights.get(j, 0)
        diagonal[j] = gp[j] + beta_weights.get(j - 1, 0)
    cell_mass = sum(alpha_weights.values(), Fraction(0)) + sum(beta_weights.values(), Fraction(0))
    require(cell_mass == sum(gp), "path selector partitions cell mass")
    return cell_mass, tuple(gp), tuple(diagonal)


single2 = path_coordinates({2: Fraction(1)}, {})
pair2 = path_coordinates({}, {2: Fraction(1)})
require(single2[:2] == pair2[:2] and single2[2] != pair2[2],
        "selector-flat but H-drifting control")

singles23 = path_coordinates({2: Fraction(1), 3: Fraction(1)}, {})
require(pair2[2] == singles23[2] and pair2[:2] != singles23[:2],
        "H-flat selector false-positive control")

left = path_coordinates({6: Fraction(1), 7: Fraction(1)}, {2: Fraction(1)})
right = path_coordinates({2: Fraction(1), 3: Fraction(1)}, {6: Fraction(1)})
require(left[0] == right[0] and left[2] == right[2] and left[1] != right[1],
        "mass-neutral H-flat selector variation")


def circular_distance(x):
    x %= 1
    return min(x, 1 - x)


def danger(y, root):
    return circular_distance(y - Fraction(root, P)) < Fraction(1, 14)


MESH = 2 * P * Q
safe_marker_counts = [0] * P
safe_diagonal_counts = [0] * P
singleton_counts = [0] * P
pair_counts = [0] * P
safe_cells = 0
for cell in range(MESH):
    y = Fraction(2 * cell + 1, 2 * MESH)
    bits = [int(danger(y, root)) for root in range(P)]
    support = [root for root, value in enumerate(bits) if value]
    require(len(support) in (1, 2), "deep mesh has one or two roots")
    starts = [root for root in range(P) if bits[root] and not bits[(root - 1) % P]]
    require(len(starts) == 1, "deep mesh has one positive start")
    start = starts[0]
    if len(support) == 1:
        singleton_counts[start] += 1
    else:
        pair_counts[start] += 1
    if not bits[0]:
        safe_cells += 1
        safe_marker_counts[start] += 1
        for root in support:
            safe_diagonal_counts[root] += 1

require(singleton_counts == [2] * P, "two singleton mesh cells per root")
require(pair_counts == [12] * P, "twelve pair mesh cells per edge")
require(safe_cells == 156, "target-safe deep mesh")
require(safe_marker_counts == [0] + [14] * 11 + [2], "full-packet selector profile")
require(safe_diagonal_counts == [0, 14] + [26] * 10 + [14], "full-packet H profile")

PEEL = Fraction(468, 12005)
full_packet_mu = PEEL * Fraction(safe_cells, MESH)
full_packet_gamma = tuple(PEEL * Fraction(value, MESH) for value in safe_marker_counts)
full_packet_h = tuple(PEEL * Fraction(value, MESH) for value in safe_diagonal_counts)
require(full_packet_mu == Fraction(2808, 84035), "THM-2367 selector mass")
require(full_packet_gamma[1:12] == (PEEL / 13,) * 11
        and full_packet_gamma[12] == PEEL / 91, "THM-2367 selector vector")
require(full_packet_h[1] == full_packet_h[12] == PEEL / 13
        and full_packet_h[2:12] == (PEEL / 7,) * 10, "THM-2367 H diagonal")


# ---------------------------------------------------------------------------
# Raw positive mass cannot be a toothpick output; centering is sufficient
# algebraically, but the lift uses a signed artificial ordered C7 pair.
# ---------------------------------------------------------------------------


mu = sum(gamma_plus)
centered13 = [13 * x - mu for x in gamma_plus]  # 13 times gamma - mu*1
star_boundary = [-mu] + gamma_plus[1:]
require(sum(centered13) == 0, "root centering")
require(sum(star_boundary) == 0, "anchored star boundary")
require(sum(gamma_plus) != 0, "raw selector cannot have zero output mass")
STAR_MOD = 131
STAR_ZETA = 39
require(pow(STAR_ZETA, P, STAR_MOD) == 1 and STAR_ZETA != 1,
        "primitive thirteenth root in F_131")
require(all(sum(star_boundary[a] * pow(STAR_ZETA, k * a, STAR_MOD)
                    for a in range(P)) % STAR_MOD != 0
            for k in range(1, P)), "anchored star boundary retains all root colours")

# Each individual occupied-to-target star edge has an exact two-tooth
# realization, at the root-dependent edge slope tau=-a.
for a in range(1, P):
    mass = gamma_plus[a]
    d_atom = [[0] * Q for _ in range(P)]
    d_atom[a][0] = mass
    d_atom[a][1] = -mass
    require(all(sum(row) == 0 for row in d_atom), "atomwise star row-zero")
    tau_edge = (-a) % P
    output_atom = [sum(d_atom[(v - tau_edge * r) % P][r] for r in range(Q))
                   for v in range(P)]
    expected_atom = [0] * P
    expected_atom[a] = mass
    expected_atom[0] = -mass
    require(output_atom == expected_atom, "atomwise two-tooth star boundary")
    require(d_atom[a][1] == -mass, "selected target tooth is signed nonzero")

for tau in range(1, P):
    primitive = [0] * P
    for j in range(1, P):
        v = (j * tau) % P
        primitive[v] = primitive[((j - 1) * tau) % P] + centered13[v]
    d = [[0] * Q for _ in range(P)]
    for h in range(P):
        d[h][0] = primitive[h]
        d[h][1] = -primitive[h]
        require(sum(d[h]) == 0, "two-column row-zero lift")
    output = [sum(d[(v - tau * r) % P][r] for r in range(Q)) for v in range(P)]
    require(output == centered13, "toothpick recovers centered selector")
    require(sum(output) == 0, "toothpick output zero mass")

    primitive_star = [0] * P
    for j in range(1, P):
        v = (j * tau) % P
        primitive_star[v] = primitive_star[((j - 1) * tau) % P] + star_boundary[v]
    d_star = [[0] * Q for _ in range(P)]
    for h in range(P):
        d_star[h][0] = primitive_star[h]
        d_star[h][1] = -primitive_star[h]
        require(sum(d_star[h]) == 0, "star-boundary row-zero lift")
    output_star = [sum(d_star[(v - tau * r) % P][r] for r in range(Q))
                   for v in range(P)]
    require(output_star == star_boundary, "toothpick recovers star boundary")

    primitive_delta = [0] * P
    for j in range(1, P):
        v = (j * tau) % P
        primitive_delta[v] = primitive_delta[((j - 1) * tau) % P] + delta_selector[v]
    d_delta = [[0] * Q for _ in range(P)]
    for h in range(P):
        d_delta[h][0] = primitive_delta[h]
        d_delta[h][1] = -primitive_delta[h]
        require(sum(d_delta[h]) == 0, "pair-divergence row-zero lift")
    output_delta = [sum(d_delta[(v - tau * r) % P][r] for r in range(Q))
                    for v in range(P)]
    require(output_delta == delta_selector, "toothpick recovers pair divergence")


# ---------------------------------------------------------------------------
# Canonical 169-cell positive control from the exact THM-2334 typed row
# ---------------------------------------------------------------------------


W = (1, 14, 27, 40, 53, 66, 13, 13**3, 2 * 13**5)
GUARD, OWNER, TA, TB = 0, 6, 7, 8
UNIT_IDX = (1, 2, 3, 4, 5)
RDIL = 169

LCM_W = 1
for value in W:
    LCM_W = LCM_W * value // gcd(LCM_W, value)
T_DEN = 182 * LCM_W
NN = RDIL * T_DEN


def in_comb(i, elli):
    """Intervals where ||W[i] x+elli/13||<1/14 on denominator T_DEN."""
    wsp = W[i]
    unit = T_DEN // (182 * wsp)
    lo = (-13 - 14 * elli) % 182
    out = []
    for n in range(wsp):
        left = (lo + 182 * n) * unit
        right = left + 26 * unit
        if right <= T_DEN:
            out.append((left, right))
        else:
            out.append((left, T_DEN))
            out.append((0, right - T_DEN))
    return sorted(out)


def subtract_comb(intervals, wsp, period_den, lo, hi):
    unit = T_DEN // (period_den * wsp)
    step = period_den * unit
    width = (hi - lo) * unit
    base = (lo % period_den) * unit
    out = []
    for left, right in intervals:
        first = left - width + 1
        k0 = -((base - first) // step)
        point = base + k0 * step
        cur = left
        while point < right:
            end = point + width
            if end > cur:
                if point > cur:
                    out.append((cur, point))
                cur = end
                if cur >= right:
                    break
            point += step
        if cur < right:
            out.append((cur, right))
    return out


def build_set(pattern, ell):
    inside = [i for i, mode in pattern.items() if mode == "in"]
    start = min(inside, key=lambda i: W[i])
    intervals = in_comb(start, ell[start])
    for i, mode in pattern.items():
        if mode == "gout":
            intervals = subtract_comb(intervals, W[i], 91, -13 - 7 * ell[i], 13 - 7 * ell[i])
    rest = sorted((W[i], i) for i, mode in pattern.items()
                  if i != start and mode in ("in", "out"))
    for _, i in rest:
        if pattern[i] == "out":
            intervals = subtract_comb(intervals, W[i], 182,
                                      -13 - 14 * ell[i], 13 - 14 * ell[i])
        else:
            intervals = subtract_comb(intervals, W[i], 182,
                                      13 - 14 * ell[i], 169 - 14 * ell[i])
    return intervals


PAT_E = {0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
         6: "in", 7: "out", 8: "out"}
PAT_QA = {0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
          6: "out", 7: "in", 8: "out"}
ZERO_ELL = (0,) * 9


def overlap_with_delayed_word(E, Qset, starts):
    """Length on denominator NN of E intersect (Qset o T^2)."""
    total = 0
    n_q = len(Qset)
    for left, right in E:
        lifted_left = RDIL * left
        start = lifted_left % T_DEN
        span = RDIL * (right - left)
        require(span < T_DEN, "one-wrap delayed-word sweep")
        end = start + span
        idx = bisect_right(starts, start) - 1
        offset = 0
        if idx < 0:
            idx = n_q - 1
            offset = -T_DEN
        while True:
            qa0, qb0 = Qset[idx]
            qa, qb = qa0 + offset, qb0 + offset
            if qa >= end:
                break
            if qb > start:
                lo = max(start, qa)
                hi = min(end, qb)
                if hi > lo:
                    total += hi - lo
            idx += 1
            if idx == n_q:
                idx = 0
                offset += T_DEN
    return total


def gf13_rank(vectors):
    matrix = [list(v) for v in vectors]
    rank = 0
    for col in range(9):
        pivot = next((r for r in range(rank, len(matrix)) if matrix[r][col] % P), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inv = pow(matrix[rank][col], P - 2, P)
        matrix[rank] = [(x * inv) % P for x in matrix[rank]]
        for r in range(len(matrix)):
            if r != rank and matrix[r][col] % P:
                factor = matrix[r][col]
                matrix[r] = [(matrix[r][c] - factor * matrix[rank][c]) % P
                             for c in range(9)]
        rank += 1
    return rank


def gf13_nullspace(rows):
    matrix = [list(r) for r in rows]
    pivots = []
    rank = 0
    for col in range(9):
        pivot = next((r for r in range(rank, len(matrix)) if matrix[r][col] % P), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inv = pow(matrix[rank][col], P - 2, P)
        matrix[rank] = [(x * inv) % P for x in matrix[rank]]
        for r in range(len(matrix)):
            if r != rank and matrix[r][col] % P:
                factor = matrix[r][col]
                matrix[r] = [(matrix[r][c] - factor * matrix[rank][c]) % P
                             for c in range(9)]
        pivots.append(col)
        rank += 1
    free = [c for c in range(9) if c not in pivots]
    basis = []
    for fc in free:
        v = [0] * 9
        v[fc] = 1
        for row, pc in enumerate(pivots):
            v[pc] = (-matrix[row][fc]) % P
        basis.append(tuple(v))
    return basis


# THM-2309 star+graft relation kernel and two quotient-dual generators.
u0 = 5
pivots = [OWNER, 0, 1, 2, 3, 4]
rows = []
for k in pivots:
    row = [0] * 9
    row[u0] = (row[u0] + W[k]) % P
    row[k] = (row[k] - W[u0]) % P
    rows.append(row)
ga = rows[pivots.index(1)]
ga[u0] = (ga[u0] + W[TA]) % P
ga[TA] = (ga[TA] - W[u0]) % P
gb = rows[pivots.index(2)]
gb[u0] = (gb[u0] + W[TB]) % P
gb[TB] = (gb[TB] - W[u0]) % P

wmod = tuple(x % P for x in W)
null_basis = gf13_nullspace(rows)
chosen = []
span = [wmod]
for vector in null_basis:
    if gf13_rank(span + [vector]) > gf13_rank(span):
        span.append(vector)
        chosen.append(vector)
require(len(chosen) == 2, "two target co-shift generators")
v1, v2 = chosen
reps = {(a, b): tuple((a * v1[i] + b * v2[i]) % P for i in range(9))
        for a in range(P) for b in range(P)}
require(len(set(reps.values())) == P * P, "169 lawful co-shifts")

QSET = build_set(PAT_QA, ZERO_ELL)
QSTARTS = [a for a, _ in QSET]
cell_mass = {}
for address, ell in reps.items():
    E = build_set(PAT_E, ell)
    cell_mass[address] = overlap_with_delayed_word(E, QSET, QSTARTS)

require(all(value > 0 for value in cell_mass.values()), "all canonical selector cells positive")
require(len(set(cell_mass.values())) == 75, "canonical selector mass has 75 values")
require(cell_mass[(0, 0)] == 60084076348296, "canonical base-cell mass")
require(cell_mass[(0, 1)] == 61887542465528, "first target-shift mass")

# A nonzero image under zeta_13 -> 8 in F_79 certifies a nonzero complex
# cyclotomic target coefficient.  Every nontrivial character survives here.
mass_spectrum = {}
for u in range(P):
    for v in range(P):
        mass_spectrum[(u, v)] = sum(
            cell_mass[(a, b)] * zpow(u * a + v * b)
            for a in range(P) for b in range(P)
        ) % MOD
require(all(mass_spectrum[q] != 0 for q in mass_spectrum if q != (0, 0)),
        "all 168 canonical target mass characters survive mod 79")


print("THM-2536 exact selector/target/toothpick referee")
print(f"path_masks={len(path_masks)} start_terminal_recovery=PASS pair_divergence=PASS")
print("target_triangle_holonomy_controls=11")
print("target_line_cancellation=PASS star_mass_line_identity=PASS singleton_zero_target_hostile=PASS")
print("selector_H_incomparability_controls=3")
print(f"thm2367_full_packet_mu={full_packet_mu} selector_and_H_target_invariant=PASS")
print("atomwise_root_to_target_two_tooth_boundaries=12")
print("uniform_centered_lifts=12 anchored_star_boundary_lifts=12 direct_pair_divergence_lifts=12")
print("raw_positive_output=IMPOSSIBLE")
print(f"canonical_cells={len(cell_mass)} positive={sum(v > 0 for v in cell_mass.values())}")
print(f"canonical_distinct_masses={len(set(cell_mass.values()))}")
print(f"canonical_overlap_min={min(cell_mass.values())} max={max(cell_mass.values())} denominator={NN}")
print(f"canonical_base_mass={Fraction(cell_mass[(0, 0)], NN)}")
print(f"canonical_first_shift_mass={Fraction(cell_mass[(0, 1)], NN)}")
print("canonical_nonzero_target_mass_characters_mod79=168/168 root=8")
print(f"checks={checks} all_passed=True")
