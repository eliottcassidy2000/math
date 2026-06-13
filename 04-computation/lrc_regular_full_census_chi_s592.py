#!/usr/bin/env python3
"""Full small regular-tournament chi census for the LRC tight-orbit question.

S592 follows the S581o handoff: compute the non-circulant regular classes on
7 and 9 vertices, then test whether the rotational/AP class is the unique
dichromatic-number-2 regular class.

Completion certificate.  The script counts labelled regular tournaments exactly
by degree-constrained row recursion.  It samples exact regular labelled
tournaments until the discovered unlabeled classes have orbit sizes summing to
that labelled count, using |orbit| = m! / |Aut(T)|.  Once the sum matches, the
unlabeled census is complete.

Tournament Analysis.  The TA vertices here are regular isomorphism classes, not
runners or arcs.  The pairwise observable is (chi(T), H(T)): lower chi is the
"closer to tight/AP" gauge, with higher Hamiltonian-path count as a deterministic
tie breaker and discovery index as the tie Hamiltonian path.  This quotient
preserves regular orbit invariants (chi, H, automorphism size) and destroys
labelled LRC geometry, observer/source marks, gap ownership, and wall positions.
Challenged assumption: "maximally cyclic" or vertex-transitive need not mean
Paley-like; the AP/interval orbit is regular but has strictly lower chi.
"""

from collections import Counter
from itertools import combinations, product
from math import factorial
from random import Random


RNG = Random(592)


def has_arc(adj, i, j):
    return (adj[i] >> j) & 1


def circulant(m, conn):
    conn = {c % m for c in conn}
    return tuple(
        sum(1 << j for j in range(m) if i != j and ((j - i) % m) in conn)
        for i in range(m)
    )


def qr_set(p):
    return sorted({(x * x) % p for x in range(1, p)})


def regular_labelled_count(m):
    """Count labelled regular tournaments with outdegree (m-1)/2."""
    d = (m - 1) // 2
    out = [0] * m
    count = 0

    def rec(i):
        nonlocal count
        if i == m - 1:
            if all(x == d for x in out):
                count += 1
            return

        later = list(range(i + 1, m))
        need = d - out[i]
        if need < 0 or need > len(later):
            return

        future_after_row = m - i - 2
        for outs_tuple in combinations(later, need):
            outs = set(outs_tuple)
            changed = []
            out[i] += need
            ok = True
            for j in later:
                if j not in outs:
                    out[j] += 1
                    changed.append(j)
                if out[j] > d or out[j] + future_after_row < d:
                    ok = False
                    break
            if ok:
                rec(i + 1)
            for j in changed:
                out[j] -= 1
            out[i] -= need

    rec(0)
    return count


def random_regular(m):
    """Construct one labelled regular tournament exactly, with randomized choices."""
    d = (m - 1) // 2
    out = [0] * m
    rows = [0] * m

    def rec(i):
        if i == m - 1:
            return all(x == d for x in out)

        later = list(range(i + 1, m))
        need = d - out[i]
        if need < 0 or need > len(later):
            return False

        choices = list(combinations(later, need))
        RNG.shuffle(choices)
        future_after_row = m - i - 2

        for outs_tuple in choices:
            outs = set(outs_tuple)
            changed = []
            old_row = rows[i]
            row_bits = 0
            out[i] += need
            ok = True
            for j in later:
                if j in outs:
                    row_bits |= 1 << j
                else:
                    rows[j] |= 1 << i
                    out[j] += 1
                    changed.append(j)
                if out[j] > d or out[j] + future_after_row < d:
                    ok = False
                    break
            rows[i] |= row_bits
            if ok and rec(i + 1):
                return True
            rows[i] = old_row
            for j in changed:
                rows[j] &= ~(1 << i)
                out[j] -= 1
            out[i] -= need
        return False

    if not rec(0):
        raise RuntimeError(f"failed to construct a regular tournament on {m} vertices")
    return tuple(rows)


def hamiltonian_paths(adj, m):
    full = (1 << m) - 1
    dp = [[0] * m for _ in range(1 << m)]
    for v in range(m):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << m):
        for v in range(m):
            val = dp[mask][v]
            if not val:
                continue
            outs = adj[v] & (full ^ mask)
            while outs:
                bit = outs & -outs
                u = bit.bit_length() - 1
                dp[mask | bit][u] += val
                outs -= bit
    return sum(dp[full])


def three_cycles(adj, m):
    count = 0
    for i, j, k in combinations(range(m), 3):
        s = has_arc(adj, i, j) + has_arc(adj, j, k) + has_arc(adj, k, i)
        if s == 0 or s == 3:
            count += 1
    return count


def acyclic_subset(adj, verts):
    for i, j, k in combinations(verts, 3):
        s = has_arc(adj, i, j) + has_arc(adj, j, k) + has_arc(adj, k, i)
        if s == 0 or s == 3:
            return False
    return True


def dichromatic_number(adj, m):
    """Minimum colors so every color class induces an acyclic subtournament."""
    classes = []

    def rec(v, k):
        if v == m:
            return True
        for color in range(k):
            classes[color].append(v)
            if acyclic_subset(adj, classes[color]) and rec(v + 1, k):
                return True
            classes[color].pop()
        return False

    for k in range(1, m + 1):
        classes = [[] for _ in range(k)]
        if rec(0, k):
            return k
    return m


def vertex_color(adj, m, v):
    """Cheap isomorphism invariant used to prune mapping search."""
    outs = [u for u in range(m) if has_arc(adj, v, u)]
    ins = [u for u in range(m) if u != v and not has_arc(adj, v, u)]

    def c3_sub(vertices):
        count = 0
        for a, b, c in combinations(vertices, 3):
            s = has_arc(adj, a, b) + has_arc(adj, b, c) + has_arc(adj, c, a)
            if s == 0 or s == 3:
                count += 1
        return count

    out_to_in = sum(has_arc(adj, o, i) for o in outs for i in ins)
    return (len(outs), c3_sub(outs), c3_sub(ins), out_to_in)


def isomorphisms(adj_a, adj_b, m, stop_after_one=False):
    colors_a = [vertex_color(adj_a, m, v) for v in range(m)]
    colors_b = [vertex_color(adj_b, m, v) for v in range(m)]
    if Counter(colors_a) != Counter(colors_b):
        return []

    order = sorted(range(m), key=lambda v: (colors_a[v], v))
    candidates = {
        v: [w for w in range(m) if colors_b[w] == colors_a[v]]
        for v in range(m)
    }
    used = [False] * m
    mapping = [-1] * m
    found = []

    def rec(pos):
        if pos == m:
            found.append(tuple(mapping))
            return stop_after_one

        best_v = None
        best_candidates = None
        for v in order:
            if mapping[v] != -1:
                continue
            feasible = []
            for w in candidates[v]:
                if used[w]:
                    continue
                ok = True
                for u in range(m):
                    wu = mapping[u]
                    if wu == -1 or u == v:
                        continue
                    if has_arc(adj_a, v, u) != has_arc(adj_b, w, wu):
                        ok = False
                        break
                    if has_arc(adj_a, u, v) != has_arc(adj_b, wu, w):
                        ok = False
                        break
                if ok:
                    feasible.append(w)
            if best_candidates is None or len(feasible) < len(best_candidates):
                best_v = v
                best_candidates = feasible
                if not feasible:
                    break
        if not best_candidates:
            return False

        for w in best_candidates:
            mapping[best_v] = w
            used[w] = True
            stop = rec(pos + 1)
            used[w] = False
            mapping[best_v] = -1
            if stop:
                return True
        return False

    rec(0)
    return found


def isomorphic(adj_a, adj_b, m):
    return bool(isomorphisms(adj_a, adj_b, m, stop_after_one=True))


def automorphism_data(adj, m):
    autos = isomorphisms(adj, adj, m, stop_after_one=False)
    seen = set()
    orbits = []
    for v in range(m):
        if v in seen:
            continue
        orbit = {p[v] for p in autos}
        seen |= orbit
        orbits.append(sorted(orbit))
    return len(autos), orbits


def scc_sizes(adj, m):
    graph = [[u for u in range(m) if has_arc(adj, v, u)] for v in range(m)]
    rev = [[] for _ in range(m)]
    for v in range(m):
        for u in graph[v]:
            rev[u].append(v)

    seen = [False] * m
    order = []

    def dfs(v):
        seen[v] = True
        for u in graph[v]:
            if not seen[u]:
                dfs(u)
        order.append(v)

    for v in range(m):
        if not seen[v]:
            dfs(v)

    seen = [False] * m
    sizes = []

    def rdfs(v):
        seen[v] = True
        size = 1
        for u in rev[v]:
            if not seen[u]:
                size += rdfs(u)
        return size

    for v in reversed(order):
        if not seen[v]:
            sizes.append(rdfs(v))
    return tuple(sorted(sizes, reverse=True))


def meta_tournament(rows):
    """Class-level TA: lower chi wins; ties prefer higher H; final tie by id."""
    m = len(rows)
    adj = [0] * m
    for i, j in combinations(range(m), 2):
        a = rows[i]
        b = rows[j]
        key_a = (a["chi"], -a["H"], a["id"])
        key_b = (b["chi"], -b["H"], b["id"])
        if key_a < key_b:
            adj[i] |= 1 << j
        else:
            adj[j] |= 1 << i
    return tuple(adj)


def relation_edge_flips(rows, adj):
    """Compare the chi-then-H gauge with a pure-H gauge."""
    flips = 0
    for i, j in combinations(range(len(rows)), 2):
        chi_edge = has_arc(adj, i, j)
        h_edge = rows[i]["H"] > rows[j]["H"]
        if rows[i]["H"] == rows[j]["H"]:
            h_edge = rows[i]["id"] < rows[j]["id"]
        if chi_edge != h_edge:
            flips += 1
    return flips


def discover_classes(m, labelled_count):
    reps = []

    def add_candidate(adj, seed_label):
        for rep in reps:
            if isomorphic(adj, rep["adj"], m):
                return False
        aut, orbits = automorphism_data(adj, m)
        rep = {
            "id": len(reps),
            "adj": adj,
            "aut": aut,
            "vt": len(orbits) == 1,
            "orbits": orbits,
            "orbit_size": factorial(m) // aut,
            "H": hamiltonian_paths(adj, m),
            "chi": dichromatic_number(adj, m),
            "c3": three_cycles(adj, m),
            "seed": seed_label,
        }
        reps.append(rep)
        return True

    half = (m - 1) // 2
    add_candidate(circulant(m, range(1, half + 1)), "AP/interval")
    if m in (3, 7, 11):
        add_candidate(circulant(m, qr_set(m)), "Paley/QR")

    for choice in product(*[(i, m - i) for i in range(1, half + 1)]):
        add_candidate(circulant(m, choice), "circulant")

    trials = 0
    while sum(rep["orbit_size"] for rep in reps) < labelled_count:
        trials += 1
        if trials > 100000:
            raise RuntimeError(f"failed to complete m={m} after {trials} trials")
        add_candidate(random_regular(m), "exact-random")
    return reps, trials


def label_rows(rows, m):
    ap = circulant(m, range(1, (m - 1) // 2 + 1))
    paley = circulant(m, qr_set(m)) if m in (3, 7, 11) else None
    for row in rows:
        labels = []
        if isomorphic(row["adj"], ap, m):
            labels.append("AP/interval")
        if paley is not None and isomorphic(row["adj"], paley, m):
            labels.append("Paley/QR")
        if not labels:
            labels.append(row["seed"])
        row["label"] = "+".join(dict.fromkeys(labels))


def print_census(m, expected_classes):
    labelled = regular_labelled_count(m)
    reps, trials = discover_classes(m, labelled)
    label_rows(reps, m)
    reps.sort(key=lambda r: (r["chi"], -r["vt"], -r["H"], r["aut"], r["id"]))
    for new_id, row in enumerate(reps):
        row["id"] = new_id

    orbit_sum = sum(rep["orbit_size"] for rep in reps)
    print(f"\n--- m={m} regular tournament full census ---")
    print(f"labelled_count={labelled}; classes_found={len(reps)}; "
          f"expected_classes={expected_classes}; orbit_sum={orbit_sum}; "
          f"completion={orbit_sum == labelled}; random_trials={trials}")
    print("id  label          aut  VT  chi  H     c3  orbit_size")
    for row in reps:
        print(f"{row['id']:2d}  {row['label'][:13]:13s} {row['aut']:3d}  "
              f"{str(row['vt'])[0]}   {row['chi']:3d}  {row['H']:4d}  "
              f"{row['c3']:3d}  {row['orbit_size']:10d}")

    chi2 = [row for row in reps if row["chi"] == 2]
    non_ap_chi2 = [row for row in chi2 if "AP/interval" not in row["label"]]
    print(f"chi=2 classes: {[row['label'] for row in chi2]}")
    print(f"non-AP chi=2 classes: {len(non_ap_chi2)}")

    adj = meta_tournament(reps)
    score_hist = Counter(sum(1 for u in range(len(reps)) if has_arc(adj, v, u))
                         for v in range(len(reps)))
    c3_meta = three_cycles(adj, len(reps))
    hp_meta = hamiltonian_paths(adj, len(reps))
    print("Tournament Analysis (class quotient; lower chi, then higher H):")
    print(f"  TA vertices={len(reps)}; score_hist={dict(sorted(score_hist.items()))}; "
          f"directed_3_cycles={c3_meta}; SCCs={scc_sizes(adj, len(reps))}; "
          f"Hamiltonian_paths={hp_meta}; edge_flips_vs_H_only={relation_edge_flips(reps, adj)}")

    assert len(reps) == expected_classes
    assert orbit_sum == labelled
    assert len(non_ap_chi2) == 0
    return reps


def main():
    print("Full regular-tournament chi census for S581o/HYP-2136")
    print("Completion uses exact labelled counts plus orbit-size sum m!/|Aut(T)|.")
    print_census(5, 1)
    print_census(7, 3)
    print_census(9, 15)

    print("\nConclusion:")
    print("  Through m=9, the AP/interval regular tournament is the unique chi=2 class.")
    print("  The m=7 non-circulant regular class exists, but chi=3; it is not a")
    print("  tight regular alternative with a new chi.  At m=9, all 14 non-AP")
    print("  regular classes have chi=3.  Thus chi does add beyond vertex-transitivity")
    print("  (AP vs Paley at m=7), and the full small census strengthens the LRC")
    print("  reading: regular tightness points to the minimally cyclic chi=2 AP orbit.")


if __name__ == "__main__":
    main()
