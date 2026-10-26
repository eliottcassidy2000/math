#!/usr/bin/env python3
"""Exact tournament isoclasses by canonical vertex augmentation, not Burnside.

Upper-pair code: bit (i,j), in lexicographic i<j order, is 1 iff i->j.
Every gate uses an explicit exception and remains active under python -O.
"""
from dataclasses import dataclass
from itertools import combinations
import argparse
import json


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def rows(n, code):
    need(n >= 0 and 0 <= code < (1 << (n * (n - 1) // 2)), "invalid code")
    a = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        if (code >> k) & 1:
            a[i] |= 1 << j
        else:
            a[j] |= 1 << i
    return tuple(a)


def encode(a, order=None):
    if order is None:
        order = tuple(range(len(a)))
    return sum(((a[u] >> v) & 1) << k
               for k, (u, v) in enumerate(combinations(order, 2)))


def relabel(n, code, order):
    """order[new vertex] = old vertex; returns the transported code."""
    need(sorted(order) == list(range(n)), "invalid relabeling")
    return encode(rows(n, code), order)


@dataclass(frozen=True)
class Canonical:
    code: int
    orders: tuple
    leaves: int

    def automorphisms(self):
        """Automorphisms on the input vertex labels, as image permutations."""
        base = self.orders[0]
        ans = []
        for order in self.orders:
            p = [0] * len(base)
            for u, v in zip(base, order):
                p[u] = v
            ans.append(tuple(p))
        return tuple(ans)

    def canonical_automorphisms(self):
        """Automorphisms on the returned canonical-code labels."""
        base = self.orders[0]
        inverse = {v: i for i, v in enumerate(base)}
        return tuple(tuple(inverse[order[i]] for i in range(len(base)))
                     for order in self.orders)


def canonical(n, code, colors=None):
    """Exact equivariant individualization/refinement; all minima retained.

    A color is an ordered integer tag. Colors are never discarded by the
    search. The return code plus colors transported by orders[0] is the
    colored canonical key. The uncolored case is used for A000568.
    """
    a = rows(n, code)
    if colors is None:
        colors = (0,) * n
    need(len(colors) == n, "color length")
    cells = tuple(tuple(v for v in range(n) if colors[v] == c)
                  for c in sorted(set(colors)))

    def refine(partition):
        while True:
            masks = tuple(sum(1 << v for v in cell) for cell in partition)
            new = []
            for cell in partition:
                pieces = {}
                for v in cell:
                    signature = tuple((a[v] & mask).bit_count() for mask in masks)
                    pieces.setdefault(signature, []).append(v)
                new.extend(tuple(pieces[s]) for s in sorted(pieces))
            new = tuple(new)
            if len(new) == len(partition):
                return new
            partition = new

    best = None
    orders = []
    leaves = 0

    def visit(partition):
        nonlocal best, orders, leaves
        partition = refine(partition)
        for i, cell in enumerate(partition):
            if len(cell) > 1:
                for v in cell:
                    rest = tuple(w for w in cell if w != v)
                    visit(partition[:i] + ((v,), rest) + partition[i + 1:])
                return
        order = tuple(cell[0] for cell in partition)
        value = encode(a, order)
        leaves += 1
        if best is None or value < best:
            best, orders = value, [order]
        elif value == best:
            orders.append(order)

    visit(cells)
    return Canonical(best, tuple(orders), leaves)


def subset_image(mask, permutation):
    return sum(1 << permutation[i] for i in range(len(permutation)) if (mask >> i) & 1)


def subset_orbits(n, group):
    """Least masks, one per orbit of the supplied complete group action."""
    unseen = set(range(1 << n))
    while unseen:
        mask = min(unseen)
        orbit = {subset_image(mask, p) for p in group}
        need(mask in orbit and orbit <= unseen, "group does not partition subset action")
        unseen.difference_update(orbit)
        yield mask, len(orbit)


def extend(n, code, mask):
    """Adjoin vertex n; mask bit i=1 means old vertex i -> new vertex."""
    need(0 <= mask < (1 << n), "extension mask")
    a = list(rows(n, code)) + [0]
    for i in range(n):
        if (mask >> i) & 1:
            a[i] |= 1 << n
        else:
            a[n] |= 1 << i
    return encode(a)


@dataclass(frozen=True)
class Record:
    code: int
    automorphisms: tuple
    parent: int | None
    extension_mask: int | None


def enumerate_classes(max_n):
    """Return (levels, receipts); levels[n] maps canonical code to Record.

    Exactly one accepted extension for every isomorphism class. Repeated
    acceptance is a hard error, not silently repaired by a set insertion.
    """
    need(max_n >= 0, "negative order")
    levels = [{0: Record(0, ((),), None, None)}]
    receipts = []
    for n in range(1, max_n + 1):
        children = {}
        candidates = leaves = orbit_mass = 0
        for parent in levels[-1].values():
            for mask, mass in subset_orbits(n - 1, parent.automorphisms):
                candidates += 1
                orbit_mass += mass
                child = extend(n - 1, parent.code, mask)
                cf = canonical(n, child)
                leaves += cf.leaves
                # The last canonical vertex defines one intrinsic deletion
                # orbit. Any minimizing labeling gives the same orbit.
                deletion_orbit = {p[-1] for p in cf.orders}
                if n - 1 not in deletion_orbit:
                    continue
                need(cf.code not in children, "canonical augmentation accepted a duplicate")
                children[cf.code] = Record(cf.code, cf.canonical_automorphisms(), parent.code, mask)
        levels.append(dict(sorted(children.items())))
        receipts.append({"n": n, "classes": len(children),
                         "rooted_extensions": candidates, "canonical_search_leaves": leaves,
                         "labelled_parent_masks": orbit_mass})
    return levels, receipts


def hamilton_count(n, code):
    a = rows(n, code)
    if n == 0:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            value = dp[mask][v]
            if not value:
                continue
            todo = a[v] & ~mask
            while todo:
                bit = todo & -todo
                todo -= bit
                dp[mask | bit][bit.bit_length() - 1] += value
    return sum(dp[-1])


def hamilton_path(n, code):
    """Deterministic insertion proof of Redei's existence theorem."""
    a = rows(n, code)
    path = []
    for v in range(n):
        j = 0
        while j < len(path) and not ((a[v] >> path[j]) & 1):
            j += 1
        path.insert(j, v)
    need(all((a[u] >> v) & 1 for u, v in zip(path, path[1:])), "insertion path")
    return tuple(path)


def hamilton_paths(n, code):
    """Every literal Hamiltonian word, including the empty word at n=0."""
    a = rows(n, code)
    full = (1 << n) - 1

    def visit(path, used):
        if used == full:
            yield path
            return
        todo = (a[path[-1]] if path else full) & ~used
        while todo:
            bit = todo & -todo
            todo -= bit
            yield from visit(path + (bit.bit_length() - 1,), used | bit)

    yield from visit((), 0)


def chart(n, code, path=None):
    """Hamiltonian spine order plus C(n-1,2) free bits. Return (word,path)."""
    if path is None:
        path = hamilton_path(n, code)
    need(sorted(path) == list(range(n)), "chart path permutation")
    a = rows(n, code)
    need(all((a[u] >> v) & 1 for u, v in zip(path, path[1:])), "chart nonpath")
    pairs = ((i, j) for i in range(n) for j in range(i + 2, n))
    word = sum(((a[path[i]] >> path[j]) & 1) << k for k, (i, j) in enumerate(pairs))
    return word, path


def unchart(n, word):
    width = max(0, (n - 1) * (n - 2) // 2) if n else 0
    need(0 <= word < 1 << width, "chart word")
    a = [0] * n
    for i in range(n - 1):
        a[i] |= 1 << (i + 1)
    for k, (i, j) in enumerate((i, j) for i in range(n) for j in range(i + 2, n)):
        if (word >> k) & 1:
            a[i] |= 1 << j
        else:
            a[j] |= 1 << i
    return encode(a)


def chart_complement(n, word):
    """Flip every nonspine arc, retaining the specified Hamiltonian spine."""
    unchart(n, word)  # input validation
    width = max(0, (n - 1) * (n - 2) // 2) if n else 0
    return word ^ ((1 << width) - 1)


def chart_reverse_converse(n, word):
    """Reverse the spine and take the converse: reflected free-pair bits."""
    unchart(n, word)  # input validation
    pairs = tuple((i, j) for i in range(n) for j in range(i + 2, n))
    where = {pair: k for k, pair in enumerate(pairs)}
    return sum(((word >> k) & 1) << where[n - 1 - j, n - 1 - i]
               for k, (i, j) in enumerate(pairs))


def reframe_ports(path, ports):
    """Transport old vertex labels to new positions in the retained path."""
    inverse = {v: i for i, v in enumerate(path)}
    need(all(v in inverse for v in ports), "port outside path")
    return tuple(inverse[v] for v in ports)


def vertex_orbits(n, group):
    unseen = set(range(n))
    out = []
    while unseen:
        v = min(unseen)
        orbit = {p[v] for p in group}
        need(orbit <= unseen, "vertex group action")
        unseen.difference_update(orbit)
        out.append(tuple(sorted(orbit)))
    return tuple(out)


def clone_core(n, code):
    """Literal ordered homogeneous-pair components, quotient code, weights.

    Components are in minimum-original-label order; vertices within each
    block are in their unique transitive order. No quotient labels erased.
    """
    a = rows(n, code)
    neighbors = [set() for _ in range(n)]
    full = (1 << n) - 1
    for u, v in combinations(range(n), 2):
        exterior = full ^ (1 << u) ^ (1 << v)
        if not ((a[u] ^ a[v]) & exterior):
            neighbors[u].add(v)
            neighbors[v].add(u)
    unseen = set(range(n))
    blocks = []
    while unseen:
        start = min(unseen)
        found, queue = {start}, [start]
        while queue:
            u = queue.pop()
            for v in neighbors[u] - found:
                found.add(v)
                queue.append(v)
        unseen.difference_update(found)
        blockmask = sum(1 << v for v in found)
        block = tuple(sorted(found, key=lambda v: -(a[v] & blockmask).bit_count()))
        need(all((a[u] >> v) & 1 for u, v in combinations(block, 2)), "clone block is not transitive")
        blocks.append(block)
    qrows = [0] * len(blocks)
    for i, block in enumerate(blocks):
        for j, other in enumerate(blocks):
            if i == j:
                continue
            bit = (a[block[0]] >> other[0]) & 1
            need(all(((a[u] >> v) & 1) == bit for u in block for v in other), "clone block not a module")
            qrows[i] |= bit << j
    return tuple(blocks), encode(qrows), tuple(map(len, blocks))


def inflate(q, weights):
    """Inflate a literal weighted quotient into transitive blocks."""
    need(all(w >= 1 for w in weights), "nonpositive block size")
    a = rows(len(weights), q)
    blocks, nextv = [], 0
    for w in weights:
        blocks.append(tuple(range(nextv, nextv + w)))
        nextv += w
    expanded = [0] * nextv
    for i, block in enumerate(blocks):
        for u, v in combinations(block, 2):
            expanded[u] |= 1 << v
        for j, other in enumerate(blocks):
            if (a[i] >> j) & 1:
                for u in block:
                    for v in other:
                        expanded[u] |= 1 << v
    return encode(expanded)


def weighted_key(q, weights):
    cf = canonical(len(weights), q, tuple(weights))
    return cf.code, tuple(weights[v] for v in cf.orders[0])


def port_coordinates(blocks, ports):
    """Ordered port occurrences, including coincidences and internal ranks."""
    position = {v: (i, j) for i, block in enumerate(blocks) for j, v in enumerate(block)}
    need(all(v in position for v in ports), "unknown port")
    return tuple(position[v] for v in ports)


def marked_key(n, code, ports):
    """Canonical key with an ordered port list; repeated endpoints retained."""
    need(all(0 <= v < n for v in ports), "unknown marked vertex")
    colors = tuple(tuple(int(v == port) for port in ports) for v in range(n))
    cf = canonical(n, code, colors)
    return cf.code, tuple(colors[v] for v in cf.orders[0])


def compositions(total, length):
    if length == 0:
        if total == 0:
            yield ()
        return
    if length == 1:
        if total >= 1:
            yield (total,)
        return
    for first in range(1, total - length + 2):
        for rest in compositions(total - first, length - 1):
            yield (first,) + rest


def weight_image(weights, permutation):
    result = [0] * len(weights)
    for v, w in enumerate(weights):
        result[permutation[v]] = w
    return tuple(result)


def weighted_core_allocations(levels, n):
    """Generate each canonical pair-free core with each weight orbit once.

    This is a second constructive route. Core catalogs are required, and the
    order-n all-one allocation still contains every order-n pair-free core.
    """
    need(0 <= n < len(levels), "core catalog horizon")
    if n == 0:
        yield 0, 0, ()
        return
    for q in range(1, n + 1):
        for rec in levels[q].values():
            blocks, _, _ = clone_core(q, rec.code)
            if len(blocks) != q:
                continue
            for weights in compositions(n, q):
                if weights == min(weight_image(weights, g) for g in rec.automorphisms):
                    yield q, rec.code, weights


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=8)
    parser.add_argument("--catalog", help="optional JSON of codes, groups, parents and chart words")
    args = parser.parse_args()
    levels, receipts = enumerate_classes(args.max_n)
    if args.catalog:
        catalog = []
        for n, level in enumerate(levels):
            catalog.append([{"code": r.code, "automorphisms": r.automorphisms,
                             "parent": r.parent, "extension_mask": r.extension_mask,
                             "chart_word": chart(n, r.code)[0]} for r in level.values()])
        with open(args.catalog, "w", encoding="utf-8", newline="\n") as out:
            json.dump(catalog, out, separators=(",", ":"))
            out.write("\n")
    print(json.dumps({"status": "PASS", "method": "canonical vertex augmentation",
                      "counts": [len(level) for level in levels], "levels": receipts}, indent=2))


if __name__ == "__main__":
    main()
