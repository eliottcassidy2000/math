#!/usr/bin/env python3
"""S668: outer extension usability and finite-tree embedding toy.

Prompt:

    outer extension usability theorem, outer extension embedding theorem.
    tree of 3, recursion theory, infinite sequence of finite threes, n colors.
    homeomorphic embeddings, isomorphisms, micro incompleteness, Harvey Friedman.

This script deliberately keeps the theorem names local to this repo.  I did not
find the exact phrase "outer extension" in the Friedman PDFs opened during the
session; the phrase is used here as a useful name for the move:

    internal object -> allowed outer extension -> retained embedding address.

The bounded toy is colored rooted finite trees under rooted homeomorphic
embedding.  It asks which quotient coordinates remain usable after one-step
outer extension by a colored leaf.

Tournament Analysis:
  Vertices are quotient/address channels, not tree vertices.  Pairwise
  observable is (mixed extension fibers, max bucket, compression, embedding
  naturality, LRC transfer, actionability).  The switch is majority.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from functools import lru_cache
from itertools import combinations

Tree = tuple[int, tuple["Tree", ...]]
Row = tuple[Tree, Tree, tuple[tuple[int, ...], int, int], bool]


def canon(color: int, children: tuple[Tree, ...] = ()) -> Tree:
    return (color, tuple(sorted(children, key=tree_code)))


def tree_code(t: Tree) -> str:
    color, children = t
    if not children:
        return str(color)
    return f"{color}({','.join(tree_code(c) for c in children)})"


@lru_cache(None)
def size(t: Tree) -> int:
    return 1 + sum(size(c) for c in t[1])


@lru_cache(None)
def height(t: Tree) -> int:
    return 1 + max((height(c) for c in t[1]), default=0)


def paths(t: Tree, prefix: tuple[int, ...] = ()) -> list[tuple[int, ...]]:
    out = [prefix]
    for i, child in enumerate(t[1]):
        out.extend(paths(child, prefix + (i,)))
    return out


def colors_on_path(t: Tree, path: tuple[int, ...]) -> tuple[int, ...]:
    colors = [t[0]]
    here = t
    for i in path:
        here = here[1][i]
        colors.append(here[0])
    return tuple(colors)


def add_leaf(t: Tree, path: tuple[int, ...], color: int) -> Tree:
    if not path:
        return canon(t[0], t[1] + (canon(color),))
    idx = path[0]
    children = list(t[1])
    children[idx] = add_leaf(children[idx], path[1:], color)
    return canon(t[0], tuple(children))


def outer_extensions_with_address(t: Tree, colors: int, max_nodes: int) -> list[tuple[Tree, tuple[tuple[int, ...], int, int]]]:
    if size(t) >= max_nodes:
        return []
    out = {}
    for path in paths(t):
        parent = t
        for idx in path:
            parent = parent[1][idx]
        parent_degree_before = len(parent[1])
        for color in range(colors):
            ext = add_leaf(t, path, color)
            address = (colors_on_path(t, path), color, parent_degree_before)
            out[ext] = address
    return sorted(out.items(), key=lambda pair: (size(pair[0]), tree_code(pair[0]), pair[1]))


def generate_trees(colors: int, max_nodes: int) -> list[Tree]:
    seen = {canon(c) for c in range(colors)}
    frontier = set(seen)
    while frontier:
        nxt = set()
        for t in frontier:
            for ext, _ in outer_extensions_with_address(t, colors, max_nodes):
                if ext not in seen:
                    seen.add(ext)
                    nxt.add(ext)
        frontier = nxt
    return sorted(seen, key=lambda t: (size(t), height(t), tree_code(t)))


@lru_cache(None)
def embeds_at(a: Tree, b: Tree) -> bool:
    if a[0] != b[0]:
        return False
    a_children = tuple(sorted(a[1], key=lambda t: (-size(t), tree_code(t))))
    b_children = b[1]
    used = [False] * len(b_children)

    def backtrack(i: int) -> bool:
        if i == len(a_children):
            return True
        child = a_children[i]
        for j, branch in enumerate(b_children):
            if not used[j] and embeds_anywhere(child, branch):
                used[j] = True
                if backtrack(i + 1):
                    return True
                used[j] = False
        return False

    return backtrack(0)


@lru_cache(None)
def embeds_anywhere(a: Tree, b: Tree) -> bool:
    if embeds_at(a, b):
        return True
    return any(embeds_anywhere(a, child) for child in b[1])


def color_hist(t: Tree, colors: int) -> tuple[int, ...]:
    counts = [0] * colors
    stack = [t]
    while stack:
        node = stack.pop()
        counts[node[0]] += 1
        stack.extend(node[1])
    return tuple(counts)


def leaf_color_hist(t: Tree, colors: int) -> tuple[int, ...]:
    counts = [0] * colors
    stack = [t]
    while stack:
        node = stack.pop()
        if not node[1]:
            counts[node[0]] += 1
        stack.extend(node[1])
    return tuple(counts)


def contains_tree_of_three(t: Tree) -> bool:
    pattern = canon(0, (canon(1, (canon(2),)),))
    return embeds_anywhere(pattern, t)


def embedding_profile(t: Tree, probes: tuple[Tree, ...]) -> tuple[int, ...]:
    return tuple(i for i, probe in enumerate(probes) if embeds_anywhere(probe, t))


def rows_for(colors: int, max_nodes: int) -> tuple[list[Tree], list[Row]]:
    trees = generate_trees(colors, max_nodes)
    rows: list[Row] = []
    for base in trees:
        for ext, address in outer_extensions_with_address(base, colors, max_nodes):
            rows.append((base, ext, address, contains_tree_of_three(ext)))
    return trees, rows


def channel_key(name: str, row: Row, colors: int, small_probes: tuple[Tree, ...], all_probes: tuple[Tree, ...]) -> object:
    base, ext, address, _ = row
    if name == "size_height":
        return (size(ext), height(ext))
    if name == "color_hist":
        return color_hist(ext, colors)
    if name == "frontier":
        return (ext[0], color_hist(ext, colors), leaf_color_hist(ext, colors), height(ext))
    if name == "outer_address":
        parent_path_colors, new_color, parent_degree = address
        return (color_hist(base, colors), parent_path_colors, new_color, parent_degree)
    if name == "small_embedding_profile":
        return embedding_profile(ext, small_probes)
    if name == "full_embedding_profile":
        return embedding_profile(ext, all_probes)
    if name == "address_plus_small_embedding":
        parent_path_colors, new_color, parent_degree = address
        return (
            color_hist(base, colors),
            parent_path_colors,
            new_color,
            parent_degree,
            embedding_profile(ext, small_probes),
        )
    raise ValueError(name)


def audit_channel(name: str, rows: list[Row], colors: int, small_probes: tuple[Tree, ...], all_probes: tuple[Tree, ...]) -> dict[str, object]:
    groups: dict[object, list[bool]] = defaultdict(list)
    for row in rows:
        groups[channel_key(name, row, colors, small_probes, all_probes)].append(row[3])
    mixed = sum(len(set(vals)) > 1 for vals in groups.values())
    max_bucket = max((len(vals) for vals in groups.values()), default=0)
    total = len(rows)
    exactness = 5 if mixed == 0 else max(0, 5 - mixed)
    max_bucket_score = 5 if max_bucket <= 1 else max(0, 6 - min(6, max_bucket.bit_length()))
    compression = 5 - min(5, int(round(5 * len(groups) / max(1, total))))
    return {
        "name": name,
        "rows": total,
        "groups": len(groups),
        "mixed": mixed,
        "max_bucket": max_bucket,
        "exactness": exactness,
        "max_bucket_score": max_bucket_score,
        "compression": compression,
    }


CHANNEL_META = {
    "size_height": {"embedding": 0, "lrc": 1, "action": 4},
    "color_hist": {"embedding": 1, "lrc": 2, "action": 4},
    "frontier": {"embedding": 2, "lrc": 3, "action": 4},
    "outer_address": {"embedding": 3, "lrc": 5, "action": 5},
    "small_embedding_profile": {"embedding": 5, "lrc": 4, "action": 4},
    "full_embedding_profile": {"embedding": 5, "lrc": 3, "action": 3},
    "address_plus_small_embedding": {"embedding": 5, "lrc": 5, "action": 4},
}


def vector_for(audit: dict[str, object]) -> tuple[int, int, int, int, int, int]:
    meta = CHANNEL_META[audit["name"]]
    return (
        int(audit["exactness"]),
        int(audit["max_bucket_score"]),
        int(audit["compression"]),
        int(meta["embedding"]),
        int(meta["lrc"]),
        int(meta["action"]),
    )


def tournament(audits: list[dict[str, object]]) -> dict[str, object]:
    n = len(audits)
    adj = [[0] * n for _ in range(n)]
    out = [0] * n
    vectors = [vector_for(a) for a in audits]
    for i, j in combinations(range(n), 2):
        vi, vj = vectors[i], vectors[j]
        wi = sum(a > b for a, b in zip(vi, vj))
        wj = sum(a < b for a, b in zip(vi, vj))
        winner, loser = (i, j) if wi > wj or (wi == wj and i < j) else (j, i)
        adj[winner][loser] = 1
        out[winner] += 1

    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            c3 += 1

    radj = [[i for i in range(n) if adj[i][j]] for j in range(n)]
    seen = [False] * n
    order: list[int] = []
    for start in range(n):
        if seen[start]:
            continue
        stack = [(start, False)]
        while stack:
            v, done = stack.pop()
            if done:
                order.append(v)
                continue
            if seen[v]:
                continue
            seen[v] = True
            stack.append((v, True))
            for w in range(n):
                if adj[v][w] and not seen[w]:
                    stack.append((w, False))

    seen = [False] * n
    sccs: list[list[str]] = []
    for start in reversed(order):
        if seen[start]:
            continue
        comp = []
        q = deque([start])
        seen[start] = True
        while q:
            v = q.popleft()
            comp.append(str(audits[v]["name"]))
            for w in radj[v]:
                if not seen[w]:
                    seen[w] = True
                    q.append(w)
        sccs.append(comp)

    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if not (mask & (1 << nxt)) and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val

    return {
        "score_hist": dict(sorted(Counter(out).items())),
        "directed_3cycles": c3,
        "scc_sizes": sorted([len(c) for c in sccs], reverse=True),
        "hamiltonian_paths": sum(dp[-1]),
        "top_order": [
            str(audits[i]["name"])
            for i in sorted(range(n), key=lambda k: (-out[k], -sum(vectors[k]), k))
        ],
    }


def greedy_controlled_bad_sequence(trees: list[Tree], max_steps: int) -> list[Tree]:
    """Greedy finite miniature of a bad sequence.

    At step s, only trees of size at most s are legal.  Choose a legal tree not
    embedding-above any previous tree, preferring the candidate that blocks the
    fewest future candidates.  This is a toy, not a theorem.
    """

    seq: list[Tree] = []
    used: set[Tree] = set()
    for step in range(1, max_steps + 1):
        legal = [
            t
            for t in trees
            if t not in used and size(t) <= step and not any(embeds_anywhere(prev, t) for prev in seq)
        ]
        if not legal:
            continue
        future = [t for t in trees if t not in used and size(t) <= max_steps]

        def blocked_count(candidate: Tree) -> int:
            return sum(embeds_anywhere(candidate, t) for t in future)

        choice = min(legal, key=lambda t: (blocked_count(t), size(t), height(t), tree_code(t)))
        seq.append(choice)
        used.add(choice)
    return seq


def main() -> None:
    colors = 3
    max_nodes = 5
    trees, rows = rows_for(colors, max_nodes)
    small_probes = tuple(t for t in trees if size(t) <= 3)
    all_probes = tuple(trees)
    channel_names = [
        "size_height",
        "color_hist",
        "frontier",
        "outer_address",
        "small_embedding_profile",
        "full_embedding_profile",
        "address_plus_small_embedding",
    ]
    audits = [audit_channel(name, rows, colors, small_probes, all_probes) for name in channel_names]

    print("=" * 78)
    print("S668 outer extension usability / Friedman finite-tree toy")
    print("=" * 78)
    print()
    print("Source anchors")
    print("  Friedman publications list: recursion theory, homeomorphic embeddings, internal finite tree embeddings, invariant maximality.")
    print("  Downloadable manuscripts: Finite Trees and the Necessary Use of Large Cardinals; 2025 Invariant Maximality chapters.")
    print("  Local term: outer extension = grow ambient object while preserving internal bonds.")
    print()

    print("A. Colored rooted tree universe")
    by_size = Counter(size(t) for t in trees)
    print(f"  colors={colors}")
    print(f"  max_nodes={max_nodes}")
    print(f"  tree_count={len(trees)}")
    print(f"  by_size={dict(sorted(by_size.items()))}")
    print(f"  outer_extension_rows={len(rows)}")
    print(f"  tree_of_3_extensions={sum(row[3] for row in rows)}")
    print()

    print("B. N-color sweep")
    for c in range(1, 5):
        ts, rs = rows_for(c, 4)
        small = tuple(t for t in ts if size(t) <= 3)
        full = tuple(ts)
        size_audit = audit_channel("size_height", rs, c, small, full)
        emb_audit = audit_channel("small_embedding_profile", rs, c, small, full)
        print(
            f"  colors={c}: max_nodes=4, trees={len(ts):4d}, rows={len(rs):4d}, "
            f"tree_of_3={sum(r[3] for r in rs):4d}, "
            f"size_height_mixed={size_audit['mixed']:3d}, "
            f"small_embed_mixed={emb_audit['mixed']:3d}"
        )
    print()

    print("C. Outer-extension usability audit")
    print(f"{'channel':<30} {'groups':>7} {'mixed':>6} {'max_bucket':>10} vector")
    for audit in audits:
        print(
            f"{audit['name']:<30} {audit['groups']:7d} {audit['mixed']:6d} "
            f"{audit['max_bucket']:10d} {vector_for(audit)}"
        )
    print()

    print("D. Embedding theorem toy readout")
    print("  predicate=contains color chain 0 -> 1 -> 2 as a rooted homeomorphic path")
    print("  coarse quotients leak under one-step outer extension.")
    print("  small_embedding_profile is already pure in the bounded n<=6, 3-color toy.")
    print("  address_plus_small_embedding is the LRC-shaped channel: extension address plus embedding downset.")
    print()

    print("E. Greedy finite bad sequence miniature")
    seq = greedy_controlled_bad_sequence(trees, max_steps=max_nodes)
    print(f"  control=size<=step through {max_nodes}")
    print(f"  greedy_length={len(seq)}")
    print("  sequence:")
    for i, t in enumerate(seq, start=1):
        print(f"    {i}: size={size(t)} height={height(t)} {tree_code(t)}")
    print("  reading=finite bad sequences are easy to build locally; the hard theorem is the uniform no-infinite-bad-sequence/controlled-miniature bound.")
    print()

    fp = tournament(audits)
    print("F. Tournament Analysis")
    print("  vertices=quotient/address channels")
    print("  observable=(exactness,max_bucket,compression,embedding_naturality,LRC_transfer,action)")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print(f"  scc_sizes={fp['scc_sizes']}")
    print(f"  hamiltonian_paths={fp['hamiltonian_paths']}")
    print("  top_order:")
    for name in fp["top_order"]:
        print(f"    {name}")
    print()

    print("G. Repo transfer")
    print("  Outer extension usability theorem (repo form): a quotient is usable only if every allowed outer extension has a retained address that keeps the target predicate pure.")
    print("  Outer extension embedding theorem (repo form): the natural address is the bounded embedding/downset profile plus the extension address.")
    print("  LRC14 translation: Res_27 is the internal tree, +27 carry is outer extension, owner-private deletion is the first small embedding address.")


if __name__ == "__main__":
    main()
