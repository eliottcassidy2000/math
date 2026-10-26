"""Exact family enumeration with a canonical deletion owner and metagraph export.

No previous tournament enumerator is imported. Standard library only.
Finite output bounds are explicit; the algorithm and proof work at every order.
"""
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations, product
from math import comb, factorial
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import time


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def arc(rows, u, v):
    return (rows[u] >> v) & 1


def decode(n, code):
    out = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            out[i if (code >> k) & 1 else j] |= 1 << (j if (code >> k) & 1 else i)
            k += 1
    return tuple(out)


def encode(rows, order=None):
    order = tuple(range(len(rows))) if order is None else order
    result = k = 0
    for a, i in enumerate(order):
        for j in order[a + 1:]:
            result |= arc(rows, i, j) << k
            k += 1
    return result


def canonical(rows):
    """Min invariant-refinement leaf code; return every minimizing labeling.

    Each labeling maps a new position to its old vertex. All minimizers, rather
    than the first, are essential for the invariant deletion-owner orbit.
    """
    groups = defaultdict(list)
    for v, row in enumerate(rows):
        groups[row.bit_count()].append(v)
    def refine(cells):
        while True:
            masks = tuple(sum(1 << v for v in cell) for cell in cells)
            updated = []
            for cell in cells:
                parts = defaultdict(list)
                for v in cell:
                    parts[tuple((rows[v] & mask).bit_count() for mask in masks)].append(v)
                updated.extend(tuple(parts[key]) for key in sorted(parts))
            updated = tuple(updated)
            if len(updated) == len(cells):
                return updated
            cells = updated

    def leaves(cells):
        cells = refine(cells)
        index = next((i for i, cell in enumerate(cells) if len(cell) > 1), None)
        if index is None:
            yield tuple(cell[0] for cell in cells)
            return
        for v in cells[index]:
            rest = tuple(w for w in cells[index] if w != v)
            yield from leaves(cells[:index] + ((v,), rest) + cells[index+1:])

    best, owners = None, []
    initial = tuple(tuple(groups[d]) for d in sorted(groups))
    for order in leaves(initial):
        code = encode(rows, order)
        if best is None or code < best:
            best, owners = code, [order]
        elif code == best:
            owners.append(order)
    return best, tuple(owners)


def automorphisms(rows):
    """Return automorphisms in canonical coordinates, including for noncanonical input."""
    code, labels = canonical(rows)
    first = labels[0]
    inverse = {v: i for i, v in enumerate(first)}
    # Input core rows are canonical in production. For arbitrary input, these
    # are automorphisms of decode(n, code), not of the original labeling.
    return tuple(tuple(inverse[p[i]] for i in range(len(rows))) for p in labels)


def homogeneous(rows, u, v):
    mask = ((1 << len(rows)) - 1) ^ (1 << u) ^ (1 << v)
    return ((rows[u] ^ rows[v]) & mask) == 0


def chains(rows):
    n = len(rows)
    after, before = {}, {}
    for u in range(n):
        for v in range(u + 1, n):
            if homogeneous(rows, u, v):
                a, b = (u, v) if arc(rows, u, v) else (v, u)
                need(a not in after and b not in before, "homogeneous pair graph branches")
                after[a], before[b] = b, a
    answer = []
    for v in range(n):
        if v not in before:
            block = [v]
            while block[-1] in after:
                block.append(after[block[-1]])
            answer.append(tuple(block))
    need(sum(map(len, answer)) == n, "homogeneous pair graph has a cycle")
    return tuple(answer)


def compositions(n, q):
    if q == 0:
        if n == 0:
            yield ()
    elif q == 1:
        if n >= 1:
            yield (n,)
    else:
        for x in range(1, n - q + 2):
            for rest in compositions(n - x, q - 1):
                yield (x,) + rest


def expand(core, weights):
    labels = tuple(i for i, w in enumerate(weights) for _ in range(w))
    rows = tuple(sum(((i < j) if labels[i] == labels[j]
                      else arc(core, labels[i], labels[j])) << j
                     for j in range(len(labels)) if i != j)
                 for i in range(len(labels)))
    offset = 0
    blocks = []
    for w in weights:
        blocks.append(tuple(range(offset, offset + w)))
        offset += w
    return rows, tuple(blocks)


def extend(rows, mask):
    """Bit v is the arc new vertex -> v."""
    n = len(rows)
    return tuple(row | (((mask >> v) & 1 ^ 1) << n)
                 for v, row in enumerate(rows)) + (mask,)


def permute_bits(mask, p):
    return sum(((mask >> p[i]) & 1) << i for i in range(len(p)))


@dataclass(frozen=True)
class Core:
    code: int
    rows: tuple
    group: tuple


@dataclass(frozen=True)
class FamilyMember:
    core: Core
    weights: tuple
    stabilizer: tuple
    rows: tuple
    blocks: tuple


def members(n, bank):
    for q in sorted(bank):
        if q > n or (q == 0 and n != 0):
            continue
        for core in bank[q]:
            for w in compositions(n, q):
                if any(tuple(w[p[i]] for i in range(q)) < w for p in core.group):
                    continue
                stabilizer = tuple(p for p in core.group
                                   if all(w[p[i]] == w[i] for i in range(q)))
                rows, blocks = expand(core.rows, w)
                yield FamilyMember(core, w, stabilizer, rows, blocks)


def build(max_order):
    empty = Core(0, (), ((),))
    bank = {0: [empty]}
    census = []
    parents = list(members(0, bank))
    receipts = {0: [{"core_order": 0, "core": 0, "weights": [], "code": 0, "aut": 1}]}
    for n in range(1, max_order + 1):
        found = {}
        stats = Counter()
        births = []
        for parent in parents:
            q = len(parent.weights)
            for phase in range(1 << q):
                stats["raw_phases"] += 1
                if any(permute_bits(phase, p) < phase for p in parent.stabilizer):
                    continue
                stats["phase_orbits"] += 1
                mask = sum((((phase >> i) & 1) ^ (j & 1)) << v
                           for i, block in enumerate(parent.blocks)
                           for j, v in enumerate(block))
                child = extend(parent.rows, mask)
                if any(homogeneous(child, n - 1, v) for v in range(n - 1)):
                    continue
                stats["accepted_rooted_orbits"] += 1
                # The canonical vertex zero has minimum outdegree.
                if child[-1].bit_count() != min(r.bit_count() for r in child):
                    continue
                code, labels = canonical(child)
                if not any(p[0] == n - 1 for p in labels):
                    continue
                need(code not in found, "owner construction generated a duplicate")
                rows = decode(n, code)
                core = Core(code, rows, automorphisms(rows))
                need(len(chains(rows)) == n, "accepted child is not pair-free")
                found[code] = core
                births.append({"core": code, "parent_core_order": q,
                               "parent_core": parent.core.code,
                               "parent_weights": list(parent.weights), "phase": phase})
        bank[n] = [found[k] for k in sorted(found)]
        parents = list(members(n, bank))
        records, codes = [], set()
        strata = Counter()
        labelled_strata = Counter()
        for member in parents:
            code, labels = canonical(member.rows)
            need(code not in codes, "weighted families overlap")
            codes.add(code)
            need(len(labels) == len(member.stabilizer), "inflation automorphism mismatch")
            need(tuple(sorted(map(len, chains(member.rows)))) == tuple(sorted(member.weights)),
                 "expansion lost its canonical chains")
            q = len(member.weights)
            strata[q] += 1
            labelled_strata[q] += factorial(n) // len(member.stabilizer)
            records.append({"core_order": q, "core": member.core.code,
                            "weights": list(member.weights), "code": code,
                            "aut": len(member.stabilizer)})
        need(sum(labelled_strata.values()) == 1 << comb(n, 2), "labelled orbit mass")
        need(labelled_strata[n] == labelled_core(n), "independent labelled core mass")
        for q, value in labelled_strata.items():
            expected = labelled_core(q) * factorial(n) // factorial(q) * comb(n - 1, q - 1)
            need(value == expected, "labelled core stratum mass")
        receipts[n] = sorted(records, key=lambda r: r["code"])
        census.append({"order": n, "classes": len(parents), "cores": len(found),
                       "by_core_order": dict(sorted(strata.items())), **dict(stats),
                       "births": sorted(births, key=lambda r: r["core"])})
        print(json.dumps({k: v for k, v in census[-1].items() if k != "births"}), flush=True)
    return bank, receipts, census


def labelled_core(n):
    if n == 0:
        return 1
    return sum((-1) ** (n-k) * comb(n-1, k-1) * factorial(n) // factorial(k)
               * (1 << comb(k, 2)) for k in range(1, n+1))


def hamilton_paths(rows):
    n = len(rows)
    def visit(path, remaining):
        if not remaining:
            yield path
            return
        candidates = remaining if not path else remaining & rows[path[-1]]
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            yield from visit(path + (bit.bit_length() - 1,), remaining ^ bit)
    yield from visit((), (1 << n) - 1)


def free_coordinates(n):
    return tuple((i, j) for i in range(n) for j in range(i+2, n))


def path_mask(rows, path):
    return sum(arc(rows, path[i], path[j]) << k
               for k, (i, j) in enumerate(free_coordinates(len(rows))))


def fixed_path(n, mask):
    rows = [0] * n
    for i in range(n-1):
        rows[i] |= 1 << (i+1)
    for k, (i, j) in enumerate(free_coordinates(n)):
        a, b = (i, j) if (mask >> k) & 1 else (j, i)
        rows[a] |= 1 << b
    return tuple(rows)


def converse(rows):
    full = (1 << len(rows)) - 1
    return tuple(full ^ row ^ (1 << v) for v, row in enumerate(rows))


def metagraph(n, records):
    """Recover every tiling and literal complement line from exact families."""
    m = len(free_coordinates(n))
    owner, opposite, fibres = {}, {}, {}
    total_paths = 0
    for record in records:
        code = record["code"]
        rows = decode(n, code)
        paths = list(hamilton_paths(rows))
        masks = {path_mask(rows, p) for p in paths}
        need(len(paths) == record["aut"] * len(masks), "Hamilton path / Aut fibre law")
        need(len(masks) & 1, "unmerged fibre parity")
        need(not (owner.keys() & masks), "different classes share a tiling")
        owner.update({mask: code for mask in masks})
        fibres[code] = len(masks)
        opposite[code] = canonical(converse(rows))[0]
        total_paths += len(paths)
    need(set(owner) == set(range(1 << m)), "family fibres do not cover fixed-path cube")
    # Independently decode the ambient fixed-path cube, with no Hamilton search.
    for mask, code in owner.items():
        need(canonical(fixed_path(n, mask))[0] == code, "direct path-cube representative mismatch")
    coords = free_coordinates(n)
    coord_index = {pair: k for k, pair in enumerate(coords)}
    reflection = tuple(coord_index[(n-1-j, n-1-i)] for i, j in coords)
    def reflect(mask):
        return permute_bits(mask, reflection)
    node = lambda code: min(code, opposite[code])
    edges = Counter()
    blue_tilings = 0
    full = (1 << m) - 1
    for t in range(1 << m):
        need(node(owner[t]) == node(owner[reflect(t)]), "reflection changes merged node")
        blue_tilings += reflect(t) == t
        u = t ^ full
        if t > u:
            continue
        a, b = sorted((node(owner[t]), node(owner[u])))
        edges[(a, b, reflect(t) == t)] += 1
    for (a, b, blue), multiplicity in edges.items():
        if not blue:
            need(multiplicity % 2 == 0, "black line reflection pairing")
    sc = sum(code == op for code, op in opposite.items())
    need(len(set(node(c) for c in opposite)) == (len(records) + sc) // 2, "converse node count")
    return {"order": n, "classes": len(records), "self_converse": sc,
            "merged_nodes": (len(records) + sc) // 2,
            "tilings": len(owner), "literal_lines": sum(edges.values()),
            "blue_tilings": blue_tilings, "hamilton_paths_over_representatives": total_paths,
            "fibres": {str(k): v for k, v in sorted(fibres.items())},
            "edges": [[a, b, "blue" if blue else "black", value]
                      for (a, b, blue), value in sorted(edges.items())]}


def independent_small(records, bound=5):
    totals = {}
    for n in range(bound+1):
        # All labelings, with no degree restriction, are an independent normal form.
        unrestricted = set()
        for mask in range(1 << comb(n, 2)):
            rows = decode(n, mask)
            unrestricted.add(min(encode(rows, p) for p in permutations(range(n))))
        generated = {min(encode(decode(n, r["code"]), p)
                         for p in permutations(range(n))) for r in records[n]}
        need(generated == unrestricted, "independent all-labelled all-permutation census")
        totals[n] = len(unrestricted)
    return totals


def hostile_controls():
    # Frozen exact witnesses inherited from weighted_core_transport/core_count.md.
    witnesses = [(0, 33, 11, 19, 7, 29), (10, 20, 33, 38, 13, 19)]
    data = []
    for rows in witnesses:
        code, labels = canonical(rows)
        core = Core(code, decode(6, code), automorphisms(decode(6, code)))
        need(len(core.group) == 3 and len(chains(core.rows)) == 6, "group-order hostile scope")
        response = [sum(1 for _ in members(n, {6: [core]})) for n in range(6, 13)]
        data.append(response)
    need(data == [[1, 4, 11, 26, 54, 102, 180], [1, 2, 7, 20, 42, 84, 156]],
         "same group order does not retain weight response")
    # A marked root is not an unrooted object: the two order-three classes have
    # respectively three and one vertex orbits.
    rooted = []
    for rows in [(6, 4, 0), (2, 4, 1)]:
        code, _ = canonical(rows)
        group = automorphisms(decode(3, code))
        rooted.append(len({tuple(sorted({p[v] for p in group})) for v in range(3)}))
    need(rooted == [3, 1], "root-forgetting multiplicity hostile")
    return {"same_group_order_weight_responses": data, "root_orbit_counts_TT3_C3": rooted}


def main():
    source_digest = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-order", type=int, default=8)
    parser.add_argument("--metagraph-order", type=int, default=7)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("normal.json"))
    parser.add_argument("--bank", type=Path, default=Path(__file__).with_name("representatives.json"))
    args = parser.parse_args()
    need(args.max_order >= 0 and 3 <= args.metagraph_order <= args.max_order, "invalid bounds")
    start = time.time()
    bank, records, census = build(args.max_order)
    known = [1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536]
    need(all(len(records[n]) == known[n] for n in range(min(args.max_order, 9)+1)), "external A000568 control")
    report = {"status": "FINITE-EXACT; implementation pending independent audit",
              "max_order": args.max_order, "source_sha256": source_digest,
              "census": [{k: v for k, v in row.items() if k != "births"} for row in census],
              "independent_all_labelled": independent_small(records, min(5, args.max_order)),
              "hostiles": hostile_controls(),
              "metagraph": metagraph(args.metagraph_order, records[args.metagraph_order])}
    need(hashlib.sha256(Path(__file__).read_bytes()).hexdigest() == source_digest,
         "producer source changed during run; refusing to publish mixed-provenance output")
    bank_bytes = (json.dumps({"encoding": "upper triangle i<j, increasing i then j, bit1 means i->j; invariant individualization/refinement canonical representatives",
                             "members": records, "births": {row["order"]: row["births"] for row in census}},
                            indent=2) + "\n").encode()
    if args.bank.suffix == ".gz":
        args.bank.write_bytes(gzip.compress(bank_bytes, mtime=0))
    else:
        args.bank.write_bytes(bank_bytes)
    report["bank_sha256"] = hashlib.sha256(args.bank.read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"result": "PASS", "elapsed_seconds": round(time.time()-start, 3),
                      "classes": len(records[args.max_order]), "bank": str(args.bank), "receipt": str(args.output)}))


if __name__ == "__main__":
    main()
