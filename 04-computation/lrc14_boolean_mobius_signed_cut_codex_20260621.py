#!/usr/bin/env python3
"""
LRC14 Boolean-Mobius signed low-depth cut -- codex 2026-06-21.

HYP-2744 showed that no single positive dihedral type atom certifies the
consecutive/AP row in the k=8 bounded bank.  This follow-up searches for the
small signed aggregate cut requested there.

The raw state is the 64-state Boolean atom law

    q[M] = meas{x : the exact missed inner-sector set is M}, M subset {1,...,6}.

We quotient by cyclic dihedral run type only after forming q[M].  The cut found
here lives on the three low-depth type sums

    T1    = q[type (1,(1))]
    T2sep = q[type (2,(1,1))]
    T2adj = q[type (2,(2))].

Positive coefficients make AP a strict minimizer of this low-depth miss mass.
Equivalently, negative coefficients give the small signed aggregate on which AP
is a strict maximizer.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd

import numpy as np
from scipy.optimize import linprog


LOW_FEATURES = ((1, (1,)), (2, (1, 1)), (2, (2,)))


def fmt(x):
    if isinstance(x, F):
        return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"
    return f"{x:.12g}"


def lcm(a, b):
    return abs(a * b) // gcd(a, b)


def mask_runs(mask):
    """Canonical cyclic run-length type under rotation and reflection."""
    bits = [(mask >> i) & 1 for i in range(6)]
    if sum(bits) == 0:
        return ()
    if all(bits):
        return (6,)
    candidates = []
    for seq in (bits, list(reversed(bits))):
        for shift in range(6):
            row = seq[shift:] + seq[:shift]
            if row[-1] == 0 and row[0] == 1:
                lens = []
                i = 0
                while i < 6:
                    if row[i]:
                        j = i
                        while j < 6 and row[j]:
                            j += 1
                        lens.append(j - i)
                        i = j
                    else:
                        i += 1
                candidates.append(tuple(lens))
    return min(candidates)


def exact_mask_atoms(E):
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    q = defaultdict(F)
    for a, b in zip(bps, bps[1:]):
        if a == b:
            continue
        mid = (a + b) / 2
        hit = {int(7 * e * mid) % 7 for e in E}
        mask = 0
        for s in range(1, 7):
            if s not in hit:
                mask |= 1 << (s - 1)
        q[mask] += b - a
    return dict(q)


def type_profile(q):
    out = defaultdict(F)
    for mask, val in q.items():
        out[(mask.bit_count(), mask_runs(mask))] += val
    return dict(out)


def equal_profiles(a, b):
    return all(a.get(k, F(0)) == b.get(k, F(0)) for k in set(a) | set(b))


def build_k8_bank(max_e=14):
    rows = []
    for rest in combinations(range(1, max_e + 1), 7):
        E = (0,) + rest
        rows.append((E, type_profile(exact_mask_atoms(E))))
    return rows


def vector(profile, features):
    return [profile.get(f, F(0)) for f in features]


def separator_lp(ap, rows, features, orbit):
    """Maximize margin gamma with L1-normalized free signed coefficients."""
    m = len(features)
    # Variables: w_0..w_{m-1}, gamma, u_0..u_{m-1}; |w_j| <= u_j, sum u <= 1.
    c = [0.0] * m + [-1.0] + [0.0] * m
    A_ub = []
    b_ub = []
    for idx, (_, prof) in enumerate(rows):
        if idx in orbit:
            continue
        diff = [float(ap.get(f, F(0)) - prof.get(f, F(0))) for f in features]
        A_ub.append([-d for d in diff] + [1.0] + [0.0] * m)
        b_ub.append(0.0)
    for j in range(m):
        row = [0.0] * (2 * m + 1)
        row[j] = 1.0
        row[m + 1 + j] = -1.0
        A_ub.append(row)
        b_ub.append(0.0)
        row = [0.0] * (2 * m + 1)
        row[j] = -1.0
        row[m + 1 + j] = -1.0
        A_ub.append(row)
        b_ub.append(0.0)
    A_ub.append([0.0] * (m + 1) + [1.0] * m)
    b_ub.append(1.0)
    bounds = [(None, None)] * m + [(0.0, None)] + [(0.0, None)] * m
    return linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), bounds=bounds, method="highs")


def positive_low_lp(ap, rows, orbit):
    """Optimize positive coefficients on LOW_FEATURES with sum one."""
    diffs = []
    kept = []
    for idx, (E, prof) in enumerate(rows):
        if idx in orbit:
            continue
        d = [prof.get(f, F(0)) - ap.get(f, F(0)) for f in LOW_FEATURES]
        diffs.append(d)
        kept.append((E, d))
    mat = np.array([[float(x) for x in d] for d in diffs])
    A_ub = [[-row[0], -row[1], -row[2], 1.0] for row in mat]
    res = linprog(
        [0.0, 0.0, 0.0, -1.0],
        A_ub=np.array(A_ub),
        b_ub=np.zeros(len(A_ub)),
        A_eq=np.array([[1.0, 1.0, 1.0, 0.0]]),
        b_eq=np.array([1.0]),
        bounds=[(0.0, None), (0.0, None), (0.0, None), (0.0, None)],
        method="highs",
    )
    return res, kept, mat


def solve_exact_active(active):
    """Solve coeff dot d_i = gamma for three active rows, sum coeff = 1."""
    if len(active) != 3:
        return None
    rows = [[active[i][1][j] for j in range(3)] + [F(-1)] + [F(0)] for i in range(3)]
    rows.append([F(1), F(1), F(1), F(0), F(1)])
    for col in range(4):
        pivot = None
        for r in range(col, 4):
            if rows[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            return None
        rows[col], rows[pivot] = rows[pivot], rows[col]
        pv = rows[col][col]
        rows[col] = [x / pv for x in rows[col]]
        for r in range(4):
            if r != col and rows[r][col] != 0:
                factor = rows[r][col]
                rows[r] = [rows[r][j] - factor * rows[col][j] for j in range(5)]
    return [rows[i][4] for i in range(4)]


def integer_coeffs(coeffs):
    den = 1
    for c in coeffs:
        den = lcm(den, c.denominator)
    return den, [int(c * den) for c in coeffs]


def validate_positive_cut(ap, rows, orbit, coeffs):
    vals = []
    for idx, (E, prof) in enumerate(rows):
        if idx in orbit:
            continue
        d = [prof.get(f, F(0)) - ap.get(f, F(0)) for f in LOW_FEATURES]
        vals.append((sum(c * x for c, x in zip(coeffs, d)), E, d))
    vals.sort(key=lambda row: row[0])
    return vals


def tournament_fingerprint(scores):
    names = list(scores)
    wins = {name: 0 for name in names}
    edge_winner = {}
    for i, j in combinations(range(len(names)), 2):
        a, b = names[i], names[j]
        winner = a if scores[a] >= scores[b] else b
        edge_winner[(a, b)] = winner
        wins[winner] += 1
    cycles = 0
    for a, b, c in combinations(names, 3):
        out = {a: 0, b: 0, c: 0}
        for u, v in ((a, b), (a, c), (b, c)):
            key = (u, v) if (u, v) in edge_winner else (v, u)
            out[edge_winner[key]] += 1
        if sorted(out.values()) == [1, 1, 1]:
            cycles += 1
    hp = 0
    for perm in permutations(names):
        if all(
            edge_winner.get((perm[i], perm[i + 1]), edge_winner.get((perm[i + 1], perm[i])))
            == perm[i]
            for i in range(len(perm) - 1)
        ):
            hp += 1
    return wins, cycles, hp


def print_cut(name, coeffs, vals, normalize=None):
    print(f"\n{name}")
    print(f"  coefficients on {LOW_FEATURES}: {[fmt(c) for c in coeffs]}")
    if normalize is not None:
        print(f"  L1/sum normalization: {fmt(normalize)}")
    print(f"  exact_min_margin = {fmt(vals[0][0])}")
    print(f"  normalized_margin = {fmt(vals[0][0] / sum(abs(c) for c in coeffs))}")
    print("  active/tight rows:")
    for val, E, diff in vals[:8]:
        print(f"    {fmt(val):>14}  E={E}  diffs={[fmt(x) for x in diff]}")


def main():
    print("LRC14 Boolean-Mobius signed low-depth cut (exact fractions)")
    print("=" * 80)
    rows = build_k8_bank()
    ap = type_profile(exact_mask_atoms(range(8)))
    features = sorted({f for _, prof in rows for f in prof} | set(ap))
    orbit = {i for i, (_, prof) in enumerate(rows) if equal_profiles(ap, prof)}
    print("Bank: E={0}+7-subsets of [1,14]")
    print(f"  rows={len(rows)}; non_AP_type_orbit_rows={len(rows)-len(orbit)}")
    print("  AP type orbit:")
    for i in sorted(orbit):
        print(f"    {rows[i][0]}")
    print(f"  dihedral run-type features={len(features)}: {features}")

    print("\n" + "=" * 80)
    print("LP separators in the 12-state dihedral quotient")
    all_lp = separator_lp(ap, rows, features, orbit)
    print(f"  all features: success={all_lp.success}; gamma={all_lp.x[len(features)] if all_lp.success else None}")
    if all_lp.success:
        nz = [(features[i], all_lp.x[i]) for i in range(len(features)) if abs(all_lp.x[i]) > 1e-8]
        print(f"    nonzero weights={nz}")
        q0_vals = []
        for i, (E, prof) in enumerate(rows):
            if i not in orbit:
                q0_vals.append((ap.get((0, ()), F(0)) - prof.get((0, ()), F(0)), E))
        q0_vals.sort()
        print(f"    q0-only exact min gap modulo AP orbit = {fmt(q0_vals[0][0])} at {q0_vals[0][1]}")

    no_q0 = [f for f in features if f != (0, ())]
    no_q0_lp = separator_lp(ap, rows, no_q0, orbit)
    print(f"\n  exclude q0: success={no_q0_lp.success}; gamma={no_q0_lp.x[len(no_q0)] if no_q0_lp.success else None}")
    if no_q0_lp.success:
        nz = [(no_q0[i], no_q0_lp.x[i]) for i in range(len(no_q0)) if abs(no_q0_lp.x[i]) > 1e-8]
        print(f"    nonzero signed weights={nz}")
        print("    sign reading: AP maximizes the negative of a positive low-depth miss functional.")

    print("\n" + "=" * 80)
    print("Three-atom positive low-depth functional")
    low_lp, kept, mat = positive_low_lp(ap, rows, orbit)
    print(f"  positive LP success={low_lp.success}; gamma={low_lp.x[3] if low_lp.success else None}")
    if not low_lp.success:
        raise SystemExit("low-depth LP failed")
    vals_float = mat.dot(low_lp.x[:3])
    active_idx = [i for i, v in enumerate(vals_float) if abs(v - vals_float.min()) < 1e-9]
    active = [kept[i] for i in active_idx]
    print(f"  float coeffs sum=1: {list(low_lp.x[:3])}; active_count={len(active)}")
    for E, d in active:
        print(f"    active E={E}; diffs={[fmt(x) for x in d]}")
    exact = solve_exact_active(active)
    if exact is None:
        raise SystemExit("could not solve exact active system")
    opt_coeffs, opt_gamma = exact[:3], exact[3]
    den, ints = integer_coeffs(opt_coeffs)
    opt_vals = validate_positive_cut(ap, rows, orbit, opt_coeffs)
    print("\n  exact optimum recovered from active rows:")
    print(f"    coeffs={list(map(fmt, opt_coeffs))}; gamma={fmt(opt_gamma)}")
    print(f"    integer coeffs={ints} over denominator {den}")
    print(f"    integer-margin={fmt(opt_gamma * den)}")
    print(f"    exact validation min={fmt(opt_vals[0][0])}; all_positive={opt_vals[0][0] > 0}")
    print("    This is the sharp L1-normalized cut on the three low-depth atoms in this bank.")

    compact = [F(21), F(57), F(2)]
    compact_vals = validate_positive_cut(ap, rows, orbit, compact)
    print_cut("  compact integer certificate", compact, compact_vals, normalize=sum(compact))
    print("  In signed-max form, use coefficients (-21,-57,-2); AP beats every non-orbit row")
    print(f"  by at least {fmt(compact_vals[0][0])} before normalization.")

    print("\n" + "=" * 80)
    print("64-state Boolean reading")
    print("  T1, T2sep, T2adj are sums of exact q[M] over 64 Boolean missed-sector atoms.")
    print("  The quotient preserves cyclic adjacency of the missed sectors and destroys:")
    print("    speed ownership, relation height, and the x-location of the atoms.")
    print("  Challenged assumption: vertices need not be runners or arcs; here the useful")
    print("    vertices are low-depth Boolean atom types, i.e. missed-sector shapes.")
    print("  What is proved here is finite bounded-bank separation modulo AP dilation,")
    print("    not a global LRC14 cap theorem.")

    print("\n" + "=" * 80)
    print("Tournament Analysis over proof lenses")
    scores = {
        "compact_low_depth_cut": (6, 5, 6, 4, 5),
        "sharp_three_atom_LP": (6, 5, 3, 5, 5),
        "q0_trivial_cut": (4, 4, 6, 2, 3),
        "relation_level2_pin": (5, 4, 4, 6, 5),
        "generic_CJJ_level2": (2, 2, 4, 2, 2),
        "full_64_boolean_law": (6, 6, 1, 6, 4),
    }
    wins, cycles, hp = tournament_fingerprint(scores)
    order = sorted(scores, key=lambda name: scores[name], reverse=True)
    print("  Pairwise observable: (separation, preserves_mask_shape, tractability,")
    print("    compatibility_with_incoming_hierarchy_work, formalizability).")
    print(f"  order={order}")
    print(f"  wins={wins}")
    print(f"  directed_3cycles={cycles}; Hamiltonian_paths={hp}")


if __name__ == "__main__":
    main()
