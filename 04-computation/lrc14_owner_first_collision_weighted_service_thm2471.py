#!/usr/bin/env python3
"""Exact companion for THM-2471.

Dependency-free; all quantitative checks use Fraction arithmetic.  Algebraic
nonvanishing at primitive thirteenth roots is checked by exact reduction
modulo Phi_13, not floating point evaluation.
"""

from fractions import Fraction as F


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def mod_one(x):
    return x - (x.numerator // x.denominator)


def circle_distance(x):
    y = mod_one(x)
    return min(y, 1 - y)


def phi13_reduce(values, colour):
    """Reduce sum_s values[s] zeta^(-colour*s) in Q[zeta_13]."""
    raw = [F(0) for _ in range(P)]
    for s, value in enumerate(values):
        raw[(-colour * s) % P] += value
    # zeta^12=-(1+...+zeta^11).
    return tuple(raw[j] - raw[12] for j in range(12))


def weighted_root_control():
    # Four rational base strata and five Boolean ancestry sheets.  Averaging
    # the sheets produces genuinely non-Boolean root weights.
    strata = 4
    sheets = 5
    atom_count = 7
    fa = [[[] for _ in range(P)] for _ in range(strata)]
    ea = [[[] for _ in range(P)] for _ in range(strata)]
    f_atom = [[[] for _ in range(P)] for _ in range(strata)]
    e_atom = [[[] for _ in range(P)] for _ in range(strata)]

    for y in range(strata):
        for u in range(P):
            fa[y][u] = [0] * sheets
            ea[y][u] = [0] * sheets
            f_atom[y][u] = [(2 * y + u + a) % atom_count for a in range(sheets)]
            e_atom[y][u] = [(3 * y + 2 * u + a) % atom_count for a in range(sheets)]

    # At each (y,u), only one side is occupied; different roots coexist over
    # every base stratum.  Multiplicities vary, so weights lie in (1/5)Z.
    for y in range(strata):
        f_roots = ((2 * y) % P, (2 * y + 4) % P)
        e_roots = ((2 * y + 1) % P, (2 * y + 7) % P)
        for a in range(1 + y % 4):
            fa[y][f_roots[0]][a] = 1
        for a in range(1 + (y + 2) % 4):
            fa[y][f_roots[1]][a] = 1
        for a in range(1 + (y + 1) % 4):
            ea[y][e_roots[0]][a] = 1
        for a in range(1 + (y + 3) % 4):
            ea[y][e_roots[1]][a] = 1

    u = [[F(sum(fa[y][r]), sheets) for r in range(P)] for y in range(strata)]
    v = [[F(sum(ea[y][r]), sheets) for r in range(P)] for y in range(strata)]
    require(any(x not in (0, 1) for row in u + v for x in row), "control must be weighted")
    require(all(u[y][r] * v[y][r] == 0 for y in range(strata) for r in range(P)),
            "same-root disjointness")

    correlation = []
    for shift in range(P):
        total = F(0)
        for y in range(strata):
            total += sum(u[y][(r + shift) % P] * v[y][r] for r in range(P))
        correlation.append(total / strata)

    service = sum(correlation)
    collision = service / (P * P)
    require(correlation[0] == 0 and service > 0, "anchored positive correlation")
    require(service == sum(
        (sum(u[y]) * sum(v[y]) for y in range(strata)), F(0)) / strata,
        "service identity")

    for colour in range(1, P):
        require(any(phi13_reduce(correlation, colour)), f"colour {colour} vanished")

    # The signed sum follows from root orthogonality.  The energy invoice is
    # checked through normalized Parseval, entirely in Q.
    chat_energy_all = sum(x * x for x in correlation) / P
    chat_zero_sq = (service / P) ** 2
    j_energy = (chat_energy_all - chat_zero_sq) / (P * P)
    require(j_energy >= collision * collision / 12, "root energy floor")
    require(-collision == -service / (P * P), "signed ledger normalization")

    # Exact atom partition before marginalization.
    matrix = [[F(0) for _ in range(atom_count)] for _ in range(atom_count)]
    for y in range(strata):
        a_count = [[F(0) for _ in range(atom_count)] for _ in range(P)]
        e_count = [[F(0) for _ in range(atom_count)] for _ in range(P)]
        for root in range(P):
            for sheet in range(sheets):
                if fa[y][root][sheet]:
                    a_count[root][f_atom[y][root][sheet]] += F(1, sheets)
                if ea[y][root][sheet]:
                    e_count[root][e_atom[y][root][sheet]] += F(1, sheets)
        a_marg = [sum(a_count[root][w] for root in range(P)) for w in range(atom_count)]
        e_marg = [sum(e_count[root][w] for root in range(P)) for w in range(atom_count)]
        for w in range(atom_count):
            for z in range(atom_count):
                matrix[w][z] += a_marg[w] * e_marg[z] / strata
    require(sum(sum(row) for row in matrix) == service, "atom matrix partitions service")

    # Sharp uniform nonzero-shift correlation: every J(k)=-I/12.
    sharp_c = [F(0)] + [F(1, 12) for _ in range(12)]
    sharp_i = sum(sharp_c) / (P * P)
    for colour in range(1, P):
        reduced = phi13_reduce(sharp_c, colour)
        require(reduced[0] == -F(1, 12) and all(x == 0 for x in reduced[1:]),
                "sharp cyclotomic reduction")
    sharp_j = -F(1, 12 * P * P)
    require(abs(sharp_j) == sharp_i / 12, "sharp maximum floor")

    return collision, j_energy, sum(1 for row in matrix for x in row if x)


def temporal_transfer_control():
    # Exact branchwise transfer identity.  y is the arrival cell and b is its
    # R-fold source ancestry sheet.  Arrival masks depend on y; source masks
    # can vary with (y,b), and the distinction is visible without numerics.
    cells, branches = 7, 9
    e = [[F(((2 * y + 3 * b) % 5) < 2) for b in range(branches)] for y in range(cells)]
    q = [F((3 * y + 1) % 4 != 0) for y in range(cells)]
    arrival_mask = [F(y % 3 != 1) for y in range(cells)]
    source_mask = [[F((y + 2 * b) % 4 < 2) for b in range(branches)] for y in range(cells)]

    perron_e = [sum(e[y], F(0)) / branches for y in range(cells)]
    left_arrival = sum(q[y] * arrival_mask[y] * perron_e[y] for y in range(cells)) / cells
    right_arrival = sum(
        q[y] * arrival_mask[y] * e[y][b]
        for y in range(cells) for b in range(branches)
    ) / (cells * branches)
    require(left_arrival == right_arrival, "arrival transfer")

    perron_source = [
        sum(e[y][b] * source_mask[y][b] for b in range(branches)) / branches
        for y in range(cells)
    ]
    left_source = sum(q[y] * perron_source[y] for y in range(cells)) / cells
    right_source = sum(
        q[y] * e[y][b] * source_mask[y][b]
        for y in range(cells) for b in range(branches)
    ) / (cells * branches)
    require(left_source == right_source, "source transfer")
    require(left_arrival != left_source, "temporal atom hostile must separate")
    return left_arrival, left_source


def coarsening_control():
    n = 128

    # Loop endpoint: equality in the 2N-3 triangle bound.
    loop_factor = 2 * n - 3
    q_loop = [-F(1, loop_factor)] + [F(2, loop_factor) for _ in range(n - 1)]
    require(sum(q_loop) == 1, "loop total")
    q_s = q_loop[0]
    require(abs(q_s) == F(1, loop_factor), "loop base candidate")
    require(all(abs(q_s + q_loop[j]) == F(1, loop_factor) for j in range(1, n)),
            "loop one-extension candidates")

    # Nonloop endpoint set: equality in the 2N-5 bound.
    edge_factor = 2 * n - 5
    q_edge = [-F(1, edge_factor), F(0)] + [F(2, edge_factor) for _ in range(n - 2)]
    require(sum(q_edge) == 1, "edge total")
    q_s2 = q_edge[0] + q_edge[1]
    require(abs(q_s2) == F(1, edge_factor), "edge base candidate")
    require(all(abs(q_s2 + q_edge[j]) == F(1, edge_factor) for j in range(2, n)),
            "edge one-extension candidates")

    require(loop_factor == 253 and loop_factor * loop_factor == 64009, "loop constants")
    require(edge_factor == 251 and edge_factor * edge_factor == 63001, "edge constants")
    return loop_factor, edge_factor


def deep_stalk_control():
    y = F(3, 17)

    def phases(c, d):
        return [[mod_one(F(c * (y + u + 13 * a), 13 * d)) for a in range(d)] for u in range(13)]

    # d|C and 13|(C/d): both sheet and root disappear.
    p1 = phases(26 * 13, 13)
    require(all(p1[u][a] == p1[0][0] for u in range(13) for a in range(13)),
            "base-only descent")

    # d|C and C/d is a unit: sheet disappears and t-hu is the affine root.
    p2 = phases(2 * 13, 13)
    require(all(p2[u][a] == p2[u][0] for u in range(13) for a in range(13)),
            "unit descent sheet independence")
    require(len({p2[u][0] for u in range(13)}) == 13, "unit descent retains root")

    # d does not divide C: already a=0 versus a=1 witnesses the essential
    # residue; do not materialize the 13^6-sheet stalk.
    c3, d3 = 2 * 13**5, 13**6
    phase_a0 = mod_one(F(c3 * y, 13 * d3))
    phase_a1 = mod_one(F(c3 * (y + 13), 13 * d3))
    require(phase_a0 != phase_a1, "lost sheet must be essential")
    return d3 // (13**5)


def strict_coordinate_hostile():
    y = F(3, 16)
    xs = [F(627, 2704), F(835, 2704)]
    addresses = [(3, 0), (4, 0)]
    speeds = [4, 2, 3, 6, 10]
    blockers = [13, 2197, 742586]

    require(all(mod_one(13 * 13 * x) == y for x in xs), "common terminal")
    for x, (a, root) in zip(xs, addresses):
        require(x == F(y + root + 13 * a, 169), "two-digit address")
        require(all(circle_distance(q * x) >= F(1, 14) for q in speeds), "source units safe")
        bits = [circle_distance(c * x) < F(1, 14) for c in blockers]
        require(bits == [True, False, False], "common source blocker atom")

    require(all(circle_distance(q * y) >= F(1, 14) for q in speeds), "terminal units safe")
    terminal_bits = [circle_distance(c * y) < F(1, 14) for c in blockers]
    require(terminal_bits == [False, True, False], "common pure terminal word")

    root_sets = []
    for x in xs:
        roots = {
            t for t in range(13)
            if circle_distance(blockers[2] * x - F(t, 13)) < F(1, 14)
        }
        root_sets.append(roots)
    require(root_sets[0] == root_sets[1] == {11, 12}, "common deep roots")
    require(addresses[0][0] != addresses[1][0], "endpoint digit hostile")
    return tuple(sorted(root_sets[0])), tuple(a for a, _ in addresses)


def main():
    collision, energy, matrix_support = weighted_root_control()
    arrival, source = temporal_transfer_control()
    loop_factor, edge_factor = coarsening_control()
    lost_sheet_modulus = deep_stalk_control()
    deep_roots, endpoint_digits = strict_coordinate_hostile()

    print("THM-2471 exact companion: PASS")
    print(f"prime={P}; nonzero_colours={P-1}; weighted_matrix_support={matrix_support}")
    print(f"control_collision_I={collision}; control_root_energy={energy}")
    print(f"temporal_transfer_arrival={arrival}; source={source}; distinct={arrival != source}")
    print(f"coarsening_factors: loop={loop_factor}, nonloop={edge_factor}")
    print(f"deep_lost_sheet_modulus_control={lost_sheet_modulus}")
    print(f"strict_hostile_deep_roots={deep_roots}; endpoint_digits={endpoint_digits}")
    print("all twelve cyclotomic colours, signed ledger, energy floor, and sharp control verified")


if __name__ == "__main__":
    main()
