#!/usr/bin/env python3
"""S673: Borel antidiagonalization, constructivity, tangible incompleteness.

The finite labs here are intentionally small.

1. Antidiagonal channel audit:
   Given an m x m binary matrix, the constructive diagonal witness is the
   anti-diagonal word a_i = 1 - M_ii.  We ask which quotient channels determine
   the witness.  Row/column weights leak; the embedded diagonal address is pure.

2. Regressive shift-collision toy:
   Friedman-style Borel/subtle-cardinal diagonalization uses regressive maps and
   shifted tuples.  On finite intervals, the predecessor boundary gives a
   canonical avoiding branch.  If that boundary address is forbidden, the branch
   collapses immediately.  This is a toy, not the theorem.

3. Method tournament:
   Compare transfer lanes by constructible witness, embedded address, Borel-code
   naturality, recursive strength, repo transfer, and overclaim control.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations, product


def banner(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def bits_of_int(x: int, width: int) -> tuple[int, ...]:
    return tuple((x >> i) & 1 for i in range(width))


def matrix_from_mask(mask: int, m: int) -> tuple[tuple[int, ...], ...]:
    bits = bits_of_int(mask, m * m)
    return tuple(tuple(bits[r * m + c] for c in range(m)) for r in range(m))


def antiword(mat: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    return tuple(1 - mat[i][i] for i in range(len(mat)))


def channel_audit(items: list[tuple[object, object]]) -> dict[str, object]:
    groups: defaultdict[object, set[object]] = defaultdict(set)
    sizes: Counter[object] = Counter()
    for key, value in items:
        groups[key].add(value)
        sizes[key] += 1
    mixed = sum(1 for values in groups.values() if len(values) > 1)
    return {
        "groups": len(groups),
        "mixed": mixed,
        "max_bucket": max(sizes.values()) if sizes else 0,
        "pure": mixed == 0,
    }


def antidiagonal_channel_rows(max_m: int = 4) -> list[dict[str, object]]:
    rows = []
    for m in range(1, max_m + 1):
        matrices = [matrix_from_mask(mask, m) for mask in range(1 << (m * m))]
        channels = {
            "row_weight_multiset": lambda mat: tuple(sorted(sum(row) for row in mat)),
            "row_weight_sequence": lambda mat: tuple(sum(row) for row in mat),
            "column_weight_sequence": lambda mat: tuple(
                sum(mat[r][c] for r in range(m)) for c in range(m)
            ),
            "row_and_column_weights": lambda mat: (
                tuple(sum(row) for row in mat),
                tuple(sum(mat[r][c] for r in range(m)) for c in range(m)),
            ),
            "diagonal_vector": lambda mat: tuple(mat[i][i] for i in range(m)),
            "diagonal_plus_row_weights": lambda mat: (
                tuple(mat[i][i] for i in range(m)),
                tuple(sum(row) for row in mat),
            ),
            "full_matrix": lambda mat: mat,
        }
        for name, fn in channels.items():
            audit = channel_audit([(fn(mat), antiword(mat)) for mat in matrices])
            rows.append({"m": m, "channel": name, "matrices": len(matrices), **audit})
    return rows


def first_mixed_example(m: int, channel_name: str) -> dict[str, object] | None:
    channels = {
        "row_weight_sequence": lambda mat: tuple(sum(row) for row in mat),
        "column_weight_sequence": lambda mat: tuple(
            sum(mat[r][c] for r in range(m)) for c in range(m)
        ),
        "row_and_column_weights": lambda mat: (
            tuple(sum(row) for row in mat),
            tuple(sum(mat[r][c] for r in range(m)) for c in range(m)),
        ),
    }
    fn = channels[channel_name]
    bucket: defaultdict[object, list[tuple[tuple[tuple[int, ...], ...], tuple[int, ...]]]] = defaultdict(list)
    for mask in range(1 << (m * m)):
        mat = matrix_from_mask(mask, m)
        bucket[fn(mat)].append((mat, antiword(mat)))
    for key, values in bucket.items():
        seen: dict[tuple[int, ...], tuple[tuple[int, ...], ...]] = {}
        for mat, aw in values:
            if aw not in seen:
                seen[aw] = mat
            if len(seen) >= 2:
                first = list(seen.items())[:2]
                return {
                    "m": m,
                    "channel": channel_name,
                    "key": key,
                    "antiwords": [aw for aw, _ in first],
                    "matrices": [mat for _, mat in first],
                }
    return None


def increasing_tuples(n: int, k: int) -> tuple[tuple[int, ...], ...]:
    return tuple(t for t in combinations(range(n), k) if t[0] > 0)


def regressive_domains(n: int, k: int, forbid_predecessor: bool = False) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    out = []
    for t in increasing_tuples(n, k):
        domain = tuple(range(t[0]))
        if forbid_predecessor and t[0] > 0:
            domain = tuple(x for x in domain if x != t[0] - 1)
        out.append((t, domain))
    return out


def shifted_constraints(n: int, k: int) -> tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]:
    cons = []
    for block in combinations(range(n), k + 1):
        left = block[:-1]
        right = block[1:]
        if left[0] > 0 and right[0] > 0:
            cons.append((left, right))
    return tuple(cons)


def count_shift_avoiders(n: int, k: int, forbid_predecessor: bool = False, node_limit: int = 2_000_000) -> dict[str, object]:
    vars_with_domains = regressive_domains(n, k, forbid_predecessor)
    variables = [t for t, _ in vars_with_domains]
    domains = {t: d for t, d in vars_with_domains}
    if any(len(d) == 0 for d in domains.values()):
        return {
            "N": n,
            "k": k,
            "forbid_predecessor": forbid_predecessor,
            "variables": len(variables),
            "constraints": len(shifted_constraints(n, k)),
            "avoiders": 0,
            "nodes": 1,
            "cutoff": False,
            "empty_domain": True,
        }
    idx = {t: i for i, t in enumerate(variables)}
    constraints = [
        (idx[a], idx[b])
        for a, b in shifted_constraints(n, k)
        if a in idx and b in idx
    ]
    incident: list[list[int]] = [[] for _ in variables]
    for a, b in constraints:
        incident[a].append(b)
        incident[b].append(a)
    order = sorted(range(len(variables)), key=lambda i: -len(incident[i]))
    assigned: list[int | None] = [None] * len(variables)
    count = 0
    nodes = 0
    cutoff = False

    def rec(pos: int) -> None:
        nonlocal count, nodes, cutoff
        nodes += 1
        if nodes > node_limit:
            cutoff = True
            return
        if pos == len(order):
            count += 1
            return
        i = order[pos]
        for value in domains[variables[i]]:
            if all(assigned[j] != value for j in incident[i]):
                assigned[i] = value
                rec(pos + 1)
                assigned[i] = None
                if cutoff:
                    return

    rec(0)
    return {
        "N": n,
        "k": k,
        "forbid_predecessor": forbid_predecessor,
        "variables": len(variables),
        "constraints": len(constraints),
        "avoiders": count,
        "nodes": nodes,
        "cutoff": cutoff,
        "empty_domain": False,
    }


def predecessor_assignment(n: int, k: int) -> dict[tuple[int, ...], int]:
    return {t: t[0] - 1 for t in increasing_tuples(n, k)}


def has_shift_collision(assignment: dict[tuple[int, ...], int], n: int, k: int) -> bool:
    for left, right in shifted_constraints(n, k):
        if left in assignment and right in assignment and assignment[left] == assignment[right]:
            return True
    return False


def regressive_rows() -> list[dict[str, object]]:
    rows = []
    for k in (1, 2, 3):
        for n in range(k + 1, k + 7):
            row = count_shift_avoiders(n, k)
            pred = predecessor_assignment(n, k)
            row["predecessor_avoids"] = not has_shift_collision(pred, n, k)
            rows.append(row)
    gated = []
    for k in (1, 2, 3):
        for n in range(k + 1, k + 5):
            gated.append(count_shift_avoiders(n, k, forbid_predecessor=True))
    return rows + gated


def constructivist_status_table() -> list[dict[str, object]]:
    return [
        {
            "principle": "Cantor antidiagonal",
            "classical_existence": 2,
            "constructive_witness": 2,
            "needs_embedded_address": 2,
            "recursion_strength": 1,
            "repo_transfer": 2,
            "risk_control": 2,
            "note": "Explicit witness if diagonal bits are named.",
        },
        {
            "principle": "Borel diagonalization",
            "classical_existence": 2,
            "constructive_witness": 1,
            "needs_embedded_address": 2,
            "recursion_strength": 2,
            "repo_transfer": 2,
            "risk_control": 1,
            "note": "Needs Borel code plus invariant embedding/selection address.",
        },
        {
            "principle": "Constructive mathematics",
            "classical_existence": 1,
            "constructive_witness": 2,
            "needs_embedded_address": 2,
            "recursion_strength": 1,
            "repo_transfer": 2,
            "risk_control": 2,
            "note": "Existence claims become witness-producing programs.",
        },
        {
            "principle": "Tangible incompleteness",
            "classical_existence": 2,
            "constructive_witness": 1,
            "needs_embedded_address": 2,
            "recursion_strength": 2,
            "repo_transfer": 2,
            "risk_control": 2,
            "note": "Concrete finite statements expose proof-strength cost.",
        },
        {
            "principle": "Embedded maximality",
            "classical_existence": 1,
            "constructive_witness": 2,
            "needs_embedded_address": 2,
            "recursion_strength": 2,
            "repo_transfer": 2,
            "risk_control": 2,
            "note": "A maximum is stable only relative to allowed extensions.",
        },
        {
            "principle": "LRC14 owner/carry rank",
            "classical_existence": 1,
            "constructive_witness": 1,
            "needs_embedded_address": 2,
            "recursion_strength": 2,
            "repo_transfer": 2,
            "risk_control": 1,
            "note": "Need rank drop after owner-private deletion.",
        },
        {
            "principle": "A000568 endpoint recursion",
            "classical_existence": 2,
            "constructive_witness": 2,
            "needs_embedded_address": 2,
            "recursion_strength": 2,
            "repo_transfer": 2,
            "risk_control": 2,
            "note": "Half-filter state should survive n->n+1.",
        },
        {
            "principle": "Unit-distance spine recursion",
            "classical_existence": 1,
            "constructive_witness": 1,
            "needs_embedded_address": 2,
            "recursion_strength": 2,
            "repo_transfer": 2,
            "risk_control": 1,
            "note": "Spine owner must beat bulk tail escape.",
        },
        {
            "principle": "Raw quotient/cardinality",
            "classical_existence": 2,
            "constructive_witness": 0,
            "needs_embedded_address": 0,
            "recursion_strength": 0,
            "repo_transfer": 0,
            "risk_control": 1,
            "note": "Counts prove too little without a selected address.",
        },
    ]


def route_beats(a: dict[str, object], b: dict[str, object]) -> bool:
    criteria = [
        "constructive_witness",
        "needs_embedded_address",
        "recursion_strength",
        "repo_transfer",
        "risk_control",
        "classical_existence",
    ]
    aw = sum(int(a[c]) > int(b[c]) for c in criteria)
    bw = sum(int(b[c]) > int(a[c]) for c in criteria)
    if aw != bw:
        return aw > bw
    return str(a["principle"]) < str(b["principle"])


def tournament_fingerprint(vertices: list[dict[str, object]]) -> dict[str, object]:
    n = len(vertices)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if route_beats(vertices[i], vertices[j]):
            adj[i][j] = True
        else:
            adj[j][i] = True
    scores = [sum(row) for row in adj]
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[b][a] and adj[c][b] and adj[a][c]
        ):
            c3 += 1

    graph = [[j for j in range(n) if adj[i][j]] for i in range(n)]
    rev = [[i for i in range(n) if adj[i][j]] for j in range(n)]
    seen = [False] * n
    order = []

    def dfs(v: int) -> None:
        seen[v] = True
        for w in graph[v]:
            if not seen[w]:
                dfs(w)
        order.append(v)

    for v in range(n):
        if not seen[v]:
            dfs(v)

    seen = [False] * n
    sccs = []

    def rdfs(v: int, acc: list[int]) -> None:
        seen[v] = True
        acc.append(v)
        for w in rev[v]:
            if not seen[w]:
                rdfs(w, acc)

    for v in reversed(order):
        if not seen[v]:
            acc: list[int] = []
            rdfs(v, acc)
            sccs.append(tuple(vertices[i]["principle"] for i in acc))

    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for used in range(1 << n):
        for v in range(n):
            if not dp[used][v]:
                continue
            for w in range(n):
                if ((used >> w) & 1) == 0 and adj[v][w]:
                    dp[used | (1 << w)][w] += dp[used][v]

    return {
        "score_histogram": dict(sorted(Counter(scores).items())),
        "directed_3cycles": c3,
        "scc_sizes": [len(s) for s in sccs],
        "sccs": sccs,
        "hamiltonian_paths": sum(dp[-1]),
        "ranking": sorted(
            [(scores[i], vertices[i]["principle"]) for i in range(n)],
            reverse=True,
        ),
    }


def main() -> None:
    banner("S673 Borel antidiagonalization and constructive embedded maximality")
    print("HYP-2248 / T743.")
    print(
        "Working synthesis: diagonalization is constructive only after the "
        "diagonal address is embedded; incompleteness becomes tangible when a "
        "finite-looking extension law asks for a uniform bound beyond the local "
        "quotient."
    )

    banner("Antidiagonal witness channel audit")
    for row in antidiagonal_channel_rows():
        print(row)
    print("Mixed examples for m=4:")
    for channel in ("row_weight_sequence", "column_weight_sequence", "row_and_column_weights"):
        print(first_mixed_example(4, channel))
    print(
        "Conclusion: row/column summaries can be exact scalar data and still "
        "fail to construct the antiword.  The diagonal vector is the embedded "
        "address."
    )

    banner("Regressive shift-collision toy")
    for row in regressive_rows():
        print(row)
    print(
        "Finite intervals have a constructive predecessor escape "
        "f(t)=min(t)-1, so the shifted diagonal collision is avoided in every "
        "tested interval.  Forbidding that predecessor boundary collapses the "
        "toy immediately.  This is the finite shadow of why endpoint/no-endpoint "
        "embedding matters in k-critical linear-order diagonalization."
    )

    banner("Constructive status table and Tournament Analysis")
    vertices = constructivist_status_table()
    for vertex in vertices:
        print(vertex)
    print(tournament_fingerprint(vertices))

    banner("Repo transfer notes")
    print(
        "HYP-2242: embedded maximality names the ambient extension before "
        "calling something maximal."
    )
    print(
        "HYP-2243: outer-extension usability says a quotient is useful only if "
        "its retained address keeps target fibers pure under growth."
    )
    print(
        "HYP-2247/T742: recursion is the transition law; constructive witness "
        "extraction asks whether that transition is programmatic."
    )
    print(
        "LRC14 next move: define diagonal/cut addresses for the C=27 owner-carry "
        "fiber, then prove the local rank drop constructively rather than by "
        "post-hoc classification."
    )
    print(
        "A000568 next move: pair the half-filter deck with a constructive child "
        "selection rule; without selection, purity can be true but unusable."
    )
    print(
        "Unit distance next move: treat a unit-spine owner as the diagonal "
        "address; bulk edge counts alone are row/column weights."
    )


if __name__ == "__main__":
    main()
