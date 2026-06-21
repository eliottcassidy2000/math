#!/usr/bin/env python3
"""
LRC14 Boolean-Mobius hierarchy scout -- codex 2026-06-21.

The Coregliano-Jeronimo-Jones hierarchy suggests a useful LRC analogue:
do not keep only the Delsarte size profile p_t of the missed-sector set.
Keep containment variables for actual missed-sector masks, then recover exact
atom probabilities by Möbius inversion on the Boolean lattice.

For the six inner sectors U={1,...,6}:

  q[M] = meas{x : the exact missed inner-sector set is M}
  a[A] = meas{x : A is contained in the missed set}
       = sum_{M superset A} q[M]
  q[M] = sum_{A superset M} (-1)^(|A|-|M|) a[A].

The size quotient sum_{|M|=t} q[M] is the THM-534/Delsarte depth law.  The
full Boolean lattice is complete for the sector atom law; intermediate
quotients such as dihedral run type are candidate LP cuts.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations


FULL = (1 << 6) - 1


def fmt(x):
    if isinstance(x, F):
        return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"
    return str(x)


def mask_bits(mask):
    return tuple(i + 1 for i in range(6) if (mask >> i) & 1)


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


def containment_from_atoms(q):
    a = {}
    for sub in range(1 << 6):
        a[sub] = sum(v for mask, v in q.items() if (mask & sub) == sub)
    return a


def mobius_atoms_from_containment(a):
    q = {}
    for mask in range(1 << 6):
        total = F(0)
        rest = FULL ^ mask
        sub = rest
        while True:
            sup = mask | sub
            total += ((-1) ** ((sup.bit_count() - mask.bit_count()))) * a[sup]
            if sub == 0:
                break
            sub = (sub - 1) & rest
        q[mask] = total
    return q


def size_profile(q):
    out = defaultdict(F)
    for mask, val in q.items():
        out[mask.bit_count()] += val
    return dict(out)


def type_profile(q):
    out = defaultdict(F)
    for mask, val in q.items():
        out[(mask.bit_count(), mask_runs(mask))] += val
    return dict(out)


def exact_equal_profiles(a, b):
    keys = set(a) | set(b)
    return all(a.get(k, F(0)) == b.get(k, F(0)) for k in keys)


def ap_profile_table():
    rows = []
    for k in range(8, 14):
        q = exact_mask_atoms(range(k))
        a = containment_from_atoms(q)
        inv = mobius_atoms_from_containment(a)
        ok = exact_equal_profiles(q, inv)
        rows.append((k, q, a, inv, ok, size_profile(q), type_profile(q)))
    return rows


def scan_k8_type_extrema(max_e=14):
    ap = type_profile(exact_mask_atoms(range(8)))
    features = sorted(ap)
    best = {f: (F(-1), None) for f in features}
    worst = {f: (F(2), None) for f in features}
    beat_max = defaultdict(int)
    beat_min = defaultdict(int)
    count = 0
    for rest in combinations(range(1, max_e + 1), 7):
        E = (0,) + rest
        prof = type_profile(exact_mask_atoms(E))
        count += 1
        for f in features:
            v = prof.get(f, F(0))
            if v > best[f][0]:
                best[f] = (v, E)
            if v < worst[f][0]:
                worst[f] = (v, E)
            if v > ap[f]:
                beat_max[f] += 1
            if v < ap[f]:
                beat_min[f] += 1
    return count, ap, features, best, worst, beat_max, beat_min


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
        for u, v in [(a, b), (a, c), (b, c)]:
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


def main():
    print("LRC14 Boolean-Mobius hierarchy scout (exact fractions)")
    print("=" * 78)
    print("Boolean containment/Mobius identity:")
    print("  a[A] = meas(A subset missed sectors)")
    print("  q[M] = sum_{A superset M} (-1)^(|A|-|M|) a[A]")
    print("  Size quotient sum_{|M|=t} q[M] is the THM-534 depth law.")
    print("  Full mask quotient has 64 states and is complete for sector atoms.")

    print()
    print("=" * 78)
    print("Consecutive/AP rows k=8..13: exact inversion and quotient sizes")
    for k, q, a, inv, ok, size, typ in ap_profile_table():
        support = sum(1 for v in q.values() if v)
        print()
        print(f"k={k}: mobius_inversion_ok={ok}; mask_support={support}; type_states={len(typ)}")
        print(f"  size profile={[fmt(size.get(t, F(0))) for t in range(7)]}")
        print("  dihedral run-type profile:")
        for key in sorted(typ):
            print(f"    {key}: {fmt(typ[key])}")

    print()
    print("=" * 78)
    print("k=8 bounded-bank scan: do individual dihedral-type atoms certify consec?")
    print("  Bank: E={0}+7-subsets of [1,14].  Observable: type atom mass.")
    count, ap, features, best, worst, beat_max, beat_min = scan_k8_type_extrema()
    print(f"  scanned_rows={count}")
    print("  feature | ap | is_AP_max | beat_AP_max | argmax | is_AP_min | beat_AP_min | argmin")
    for f in features:
        maxv, argmax = best[f]
        minv, argmin = worst[f]
        print(
            f"  {f}: ap={fmt(ap[f])}; "
            f"max={fmt(maxv)}; is_AP_max={beat_max[f] == 0}; beat_AP_max={beat_max[f]}; argmax={argmax}; "
            f"min={fmt(minv)}; is_AP_min={beat_min[f] == 0}; beat_AP_min={beat_min[f]}; argmin={argmin}"
        )

    print()
    print("=" * 78)
    print("Interpretation")
    print("  The full Boolean-Mobius lift is complete for the sector atom law, exactly")
    print("  mirroring the pseudoprobability/Mobius view of higher-order Delsarte LPs.")
    print("  But individual dihedral-type atoms do not give a clean monotone proof:")
    print("  in the k=8 bank, AP is maximal only for depth 0 and the deepest tail")
    print("  types (5) and (6).  Most type atoms have many AP beaters.")
    print("  This agrees with HYP-2738: the cap-closing certificate must be a signed")
    print("  aggregate cut, not a single positive monotone atom.")

    print()
    print("=" * 78)
    print("Tournament Analysis over hierarchy views")
    print("  Pairwise observable: (completeness, preserves_generated_masks, tractability, proof_signal).")
    scores = {
        "full_boolean_mobius": (6, 6, 2, 5),
        "dihedral_type_quotient": (4, 4, 4, 4),
        "size_delsarte_depth": (2, 2, 6, 4),
        "single_type_atom": (1, 3, 6, 1),
        "raw_tanner_carrier": (1, 1, 3, 2),
    }
    wins, cycles, hp = tournament_fingerprint(scores)
    order = sorted(scores, key=lambda name: scores[name], reverse=True)
    print(f"  order={order}")
    print(f"  wins={wins}")
    print(f"  directed_3cycles={cycles}; Hamiltonian_paths={hp}")
    print("  Next target: construct a small signed nonnegative-on-generated-laws cut")
    print("  in the type or full Boolean-Mobius basis, then certify it symbolically.")


if __name__ == "__main__":
    main()
