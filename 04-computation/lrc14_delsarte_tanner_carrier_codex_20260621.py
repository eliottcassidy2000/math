#!/usr/bin/env python3
"""
LRC14 Delsarte/Tanner carrier audit -- codex 2026-06-21.

This script follows the Delsarte prompts by turning the THM-534 dual
polynomials into exact bipartite carriers:

  check nodes: factorial-moment rows r with y_r != 0
  variable nodes: missed-depth atoms t=0..6
  edge r--t: r <= t, weighted by y_r * binom(t,r)

It asks whether the exact Delsarte certificate has a Tanner/LDPC-like sparse
local realization, whether the sign data has a Doyle-Holt/half-arc flavor, and
which proof lens should be prioritized.  Everything numerical is exact over Q.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import comb


def fmt(x):
    if isinstance(x, F):
        return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"
    return str(x)


def lcm(a, b):
    from math import gcd
    return abs(a * b) // gcd(a, b)


def kraw(j, t, n=6):
    return sum((-1) ** i * comb(t, i) * comb(n - t, j - i) for i in range(j + 1))


KTAB = [[F(kraw(j, t)) for t in range(7)] for j in range(7)]


def kraw_expand(g):
    aug = [[KTAB[j][t] for j in range(7)] + [g[t]] for t in range(7)]
    for col in range(7):
        piv = next(r for r in range(col, 7) if aug[r][col] != 0)
        aug[col], aug[piv] = aug[piv], aug[col]
        pv = aug[col][col]
        aug[col] = [x / pv for x in aug[col]]
        for r in range(7):
            if r != col and aug[r][col] != 0:
                f = aug[r][col]
                aug[r] = [aug[r][i] - f * aug[col][i] for i in range(8)]
    return [aug[i][7] for i in range(7)]


def binom_basis_coeffs(g):
    """Newton/binomial-basis coefficients y_r with g(t)=sum_r y_r*C(t,r)."""
    y = []
    for r in range(7):
        y.append(sum((-1) ** (r - i) * F(comb(r, i)) * g[i] for i in range(r + 1)))
    return y


def scale_vector(v):
    q = 1
    for x in v:
        q = lcm(q, x.denominator)
    return q, [int(x * q) for x in v]


DUALS = {
    "K8": {
        "k_range": "k=8",
        "g": [F((t - 1) * (t - 2) * (t - 4) * (t - 5), 40) for t in range(7)],
    },
    "K9": {
        "k_range": "k=9,10",
        "g": [F(-(t - 2) * (t - 3) * (t - 6), 36) for t in range(7)],
    },
    "K11": {
        "k_range": "k=11,12,13",
        "g": [F((t - 3) * (t - 4), 12) for t in range(7)],
    },
}


def carrier_edges(y):
    rows = [r for r, yr in enumerate(y) if yr != 0]
    edges = []
    for r in rows:
        for t in range(r, 7):
            weight = y[r] * comb(t, r)
            if weight != 0:
                edges.append((r, t, weight))
    return rows, edges


def degree_hist(values):
    out = {}
    for v in values:
        out[v] = out.get(v, 0) + 1
    return dict(sorted(out.items()))


def four_cycle_count(rows, edges):
    adj_v = {t: set() for t in range(7)}
    for r, t, _ in edges:
        adj_v[t].add(r)
    count = 0
    witnesses = []
    for a, b in combinations(range(7), 2):
        common = sorted(adj_v[a] & adj_v[b])
        c = comb(len(common), 2)
        count += c
        if c and len(witnesses) < 5:
            witnesses.append((a, b, common, c))
    return count, witnesses


def automorphisms(rows, edges):
    """Automorphisms of the unweighted bipartite carrier preserving sides."""
    row_list = list(rows)
    edge_set = {(r, t) for r, t, _ in edges}
    row_deg = {r: sum(1 for rr, _, _ in edges if rr == r) for r in row_list}
    var_deg = {t: sum(1 for _, tt, _ in edges if tt == t) for t in range(7)}

    row_perm_candidates = []
    for perm in permutations(row_list):
        m = dict(zip(row_list, perm))
        if all(row_deg[r] == row_deg[m[r]] for r in row_list):
            row_perm_candidates.append(m)

    var_perm_candidates = []
    for perm in permutations(range(7)):
        m = dict(zip(range(7), perm))
        if all(var_deg[t] == var_deg[m[t]] for t in range(7)):
            var_perm_candidates.append(m)

    autos = []
    for rm in row_perm_candidates:
        for vm in var_perm_candidates:
            image = {(rm[r], vm[t]) for r, t in edge_set}
            if image == edge_set:
                autos.append((rm, vm))
    return autos


def edge_orbits(rows, edges):
    autos = automorphisms(rows, edges)
    edge_set = {(r, t) for r, t, _ in edges}
    weight = {(r, t): w for r, t, w in edges}
    unseen = set(edge_set)
    orbits = []
    while unseen:
        start = min(unseen)
        orb = set()
        frontier = {start}
        while frontier:
            e = frontier.pop()
            if e in orb:
                continue
            orb.add(e)
            r, t = e
            for rm, vm in autos:
                im = (rm[r], vm[t])
                if im not in orb:
                    frontier.add(im)
        unseen -= orb
        signs = sorted(set(1 if weight[e] > 0 else -1 for e in orb))
        orbits.append((sorted(orb), signs))
    return autos, orbits


def feasible_edge_flips(g, edges):
    target = [F(1)] + [F(0)] * 6
    feasible = []
    failing = []
    for r, t, w in edges:
        h = list(g)
        h[t] -= 2 * w
        ok = all(h[i] >= target[i] for i in range(7))
        (feasible if ok else failing).append((r, t, w, h[t]))
    return feasible, failing


def feasible_row_flips(g, y):
    target = [F(1)] + [F(0)] * 6
    out = []
    for r, yr in enumerate(y):
        if yr == 0:
            continue
        h = list(g)
        for t in range(r, 7):
            h[t] -= 2 * yr * comb(t, r)
        out.append((r, yr, all(h[i] >= target[i] for i in range(7)), h))
    return out


def has_mixed_sign_orbit(orbits):
    return any(len(signs) > 1 for _, signs in orbits)


def tournament_fingerprint(scores):
    names = list(scores)
    edges = {}
    score_hist = {name: 0 for name in names}
    flips = 0
    for a, b in combinations(names, 2):
        if scores[a] >= scores[b]:
            edges[(a, b)] = a
            score_hist[a] += 1
        else:
            edges[(a, b)] = b
            score_hist[b] += 1
            flips += 1
    cycles = 0
    for a, b, c in combinations(names, 3):
        wins = {
            (a, b): edges[tuple(sorted((a, b), key=names.index))],
            (a, c): edges[tuple(sorted((a, c), key=names.index))],
            (b, c): edges[tuple(sorted((b, c), key=names.index))],
        }
        out = {a: 0, b: 0, c: 0}
        for (u, v), winner in wins.items():
            out[winner] += 1
        if sorted(out.values()) == [1, 1, 1]:
            cycles += 1
    hp = 0
    for perm in permutations(names):
        ok = True
        for i in range(len(perm) - 1):
            a, b = perm[i], perm[i + 1]
            key = tuple(sorted((a, b), key=names.index))
            if edges[key] != a:
                ok = False
                break
        if ok:
            hp += 1
    return score_hist, cycles, hp, flips


def main():
    print("LRC14 Delsarte/Tanner carrier audit (exact fractions)")
    print("=" * 78)
    print("Vertex challenge:")
    print("  considered vertices: missed-depth atoms, factorial checks, signed edges,")
    print("  proof lenses, unit-distance edges, modular branch-address analogies.")
    print("  chosen carrier vertices below: atom depths t=0..6 and moment rows r.")
    print("  preserved predicate: the THM-534 Delsarte dominance g(t)>=1[t=0].")
    print("  destroyed information: actual generated speed words, relation lattice,")
    print("  and row-slice odometer geometry from HYP-2737.")

    global_sign_rigid = True
    global_sparse = True

    for label, data in DUALS.items():
        g = data["g"]
        y = binom_basis_coeffs(g)
        c = kraw_expand(g)
        y_scale, y_int = scale_vector(y)
        g_scale, g_int = scale_vector(g)
        c_scale, c_int = scale_vector(c)
        rows, edges = carrier_edges(y)
        autos, orbits = edge_orbits(rows, edges)
        mixed = has_mixed_sign_orbit(orbits)
        global_sign_rigid = global_sign_rigid and not mixed
        e_count = len(edges)
        density = F(e_count, len(rows) * 7)
        global_sparse = global_sparse and density <= F(1, 2)
        var_degrees = [sum(1 for _, t, _ in edges if t == v) for v in range(7)]
        row_degrees = [sum(1 for r, _, _ in edges if r == rr) for rr in rows]
        four_cycles, witnesses = four_cycle_count(rows, edges)
        feasible_flips, failing_flips = feasible_edge_flips(g, edges)
        row_flips = feasible_row_flips(g, y)
        sign_counts = {
            "+": sum(1 for _, _, w in edges if w > 0),
            "-": sum(1 for _, _, w in edges if w < 0),
        }
        tight = [t for t in range(7) if g[t] == (F(1) if t == 0 else F(0))]
        slack = [t for t in range(7) if g[t] > (F(1) if t == 0 else F(0))]
        print()
        print("-" * 78)
        print(f"{label} carrier ({data['k_range']})")
        print(f"  g scale {g_scale}: {g_int}")
        print(f"  y in binom basis, scale {y_scale}: {y_int}")
        print(f"  Krawtchouk coeffs, scale {c_scale}: {c_int}")
        print(f"  tight depths={tight}; positive-slack depths={slack}")
        print(f"  check rows={rows}; edges={e_count}; density={fmt(density)}; sign_counts={sign_counts}")
        print(f"  variable degree hist={degree_hist(var_degrees)}; row degree hist={degree_hist(row_degrees)}")
        print(f"  Tanner girth={'4' if four_cycles else 'acyclic/large'}; 4-cycles={four_cycles}")
        print(f"  sample 4-cycle witnesses={witnesses}")
        print(f"  automorphisms preserving sides={len(autos)}; edge orbits={len(orbits)}; mixed-sign orbits={mixed}")
        print("  first edge orbits:")
        for orbit, signs in orbits[:8]:
            print(f"    orbit={orbit} signs={signs}")
        print(f"  feasible single-edge sign flips={len(feasible_flips)}/{len(edges)}")
        if feasible_flips:
            sample = [(r, t, fmt(w), fmt(new)) for r, t, w, new in feasible_flips[:6]]
            print(f"    sample feasible flips (r,t,old_weight,new_g_t)={sample}")
        bad_sample = [(r, t, fmt(w), fmt(new)) for r, t, w, new in failing_flips[:6]]
        print(f"    sample failing flips (r,t,old_weight,new_g_t)={bad_sample}")
        print("  whole-row sign flips preserving dominance:")
        for r, yr, ok, h in row_flips:
            print(f"    row r={r}, y={fmt(yr)} -> feasible={ok}, new_g={[fmt(x) for x in h]}")

    print()
    print("=" * 78)
    print("Interpretation")
    print("  Tanner import test: negative.  The carrier is a nested Ferrers graph,")
    print("  dense at this scale, with many 4-cycles and non-biregular degree gradients.")
    print("  That is not an LDPC/Tanner expander mechanism; it is a global association")
    print("  scheme certificate written as local incidence.")
    print("  Doyle-Holt half-arc test: partly useful.  The undirected carrier has")
    print("  automorphisms, but no edge orbit mixes positive and negative signs.  The")
    print("  sign/orientation is rigid, like a half-arc split, and arbitrary arc flips")
    print("  usually break Delsarte dominance.")

    print()
    print("=" * 78)
    print("Tournament Analysis over proof lenses")
    print("  Pairwise observable: lexicographic proof value")
    print("    (preserves_bound, exact_audit, LRC_transfer, sign_structure, locality).")
    print("  Switch/gauge: keep generated atom words before scalarizing; Tanner locality")
    print("  is a tie-break, not the primary invariant.")
    scores = {
        "global_delsarte_lp": (5, 5, 5, 3, 1),
        "krawtchouk_parity_puncture": (4, 5, 4, 3, 1),
        "halfarc_sign_orbit": (3, 4, 3, 5, 2),
        "tanner_local_carrier": (2, 5, 2, 2, 1),
        "unit_distance_delsarte_transfer": (2, 3, 2, 2, 2),
        "modular_belyi_address": (1, 2, 1, 3, 1),
    }
    hist, cycles, hp, flips = tournament_fingerprint(scores)
    order = sorted(scores, key=lambda n: scores[n], reverse=True)
    print(f"  score order={order}")
    print(f"  score histogram={hist}")
    print(f"  directed_3cycles={cycles}; Hamiltonian_paths={hp}; edge_flips_vs_list_order={flips}")

    print()
    print("=" * 78)
    print("Next sharp target")
    print("  Do not try to prove LRC14 by importing a sparse Tanner-expander theorem.")
    print("  Use this audit as a guardrail: the Delsarte side is exact and globally")
    print("  positive, while HYP-2737's row-slice odometer is the active local proof")
    print("  obligation.  A productive hybrid target is a puncture/extend parity lemma:")
    print("  explain why the K8 even-only rows and K9/K11 mixed rows are compatible")
    print("  with the HYP-2737 single-row balance after quotienting by atom depth.")


if __name__ == "__main__":
    main()
