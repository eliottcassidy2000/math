"""Exact small black distances between rigid SCC words, inherited full tables."""
from collections import defaultdict, deque
from itertools import permutations, product
from pathlib import Path
import csv
import json

ROOT = Path(__file__).resolve().parents[3]
TABLES = ROOT / "04-computation/experiments"


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def decode(n, code):
    out = [0] * n
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bit, (i, j) in enumerate(reversed(pairs)):
        a, b = (i, j) if (code >> bit) & 1 else (j, i)
        out[a] |= 1 << b
    return tuple(out)


def canon(out):
    groups = defaultdict(list)
    for v, row in enumerate(out):
        groups[row.bit_count()].append(v)
    best = None
    for parts in product(*(permutations(groups[s]) for s in sorted(groups))):
        perm = tuple(x for part in parts for x in part)
        code = 0
        for i in range(len(out)):
            for j in range(i + 1, len(out)):
                code = (code << 1) | ((out[perm[i]] >> perm[j]) & 1)
        best = code if best is None else min(best, code)
    return best


def word_tournament(word):
    n = sum(word)
    out = [sum(1 << j for j in range(i + 1, n)) for i in range(n)]
    start = 0
    for length in word:
        if length == 3:
            out[start] ^= 1 << (start + 2)
            out[start + 2] ^= 1 << start
        start += length
    return tuple(out)


def words(n):
    if n == 0:
        yield ()
    for k in (1, 3):
        if n >= k:
            for rest in words(n-k):
                yield (k,) + rest


def paths(out):
    n = len(out)
    def rec(path, unseen):
        if not unseen:
            yield tuple(path)
        else:
            candidates = unseen & out[path[-1]] if path else unseen
            while candidates:
                bit = candidates & -candidates
                candidates ^= bit
                yield from rec(path + [bit.bit_length()-1], unseen ^ bit)
    yield from rec([], (1 << n)-1)


def blue(out, path):
    reverse = {v: path[-i-1] for i, v in enumerate(path)}
    return all(((out[u] >> v) & 1) == ((out[reverse[v]] >> reverse[u]) & 1)
               for u in range(len(out)) for v in range(len(out)) if u != v)


def reverse_path(out, path):
    changed = list(out)
    for u, v in zip(path, path[1:]):
        need((out[u] >> v) & 1, "invalid Hamiltonian path")
        changed[u] ^= 1 << v
        changed[v] ^= 1 << u
    return tuple(changed)


def distances(graph, root):
    parent = {root: None}
    distance = {root: 0}
    queue = deque([root])
    while queue:
        u = queue.popleft()
        for v in graph[u]:
            if v not in parent:
                parent[v] = u
                distance[v] = distance[u]+1
                queue.append(v)
    return distance, parent


def main():
    reports = []
    for n in range(4, 8):
        with (TABLES / f"tiling_line_metagraph_census_20261020_classes_n{n}.tsv").open() as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        by_id = {int(r["class"]): r for r in rows}
        merged = {i: min(i, int(r["sigmaClass"])) for i, r in by_id.items()}
        address = {int(r["canon"]): merged[i] for i, r in by_id.items()}
        graph = defaultdict(set)
        with (TABLES / f"tiling_line_metagraph_census_20261020_lines_n{n}.tsv").open() as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if int(r["black"]):
                    u, v = merged[int(r["classA"])], merged[int(r["classB"])]
                    graph[u].add(v)
                    graph[v].add(u)
        rigid = {}
        for w in words(n):
            code = canon(word_tournament(w))
            need(code in address, "word absent from complete universe")
            key = address[code]
            need(int(by_id[key]["tc"]) == 1, "rigid ordinary fibre not one")
            if w <= w[::-1] and w != w[::-1]:
                rigid[w] = key
                need(key in graph and len(graph[key]-{key}) == 1,
                     "nonsymmetric rigid word is not a black support leaf")
        pair_data = []
        for w, u in rigid.items():
            dist, parent = distances(graph, u)
            for z, v in rigid.items():
                if w >= z:
                    continue
                need(v in dist, "rigid basins disconnected")
                route = [v]
                while route[-1] != u:
                    route.append(parent[route[-1]])
                route.reverse()
                witnesses = []
                literal = word_tournament(w)
                literal_witnesses = []
                for a, b in zip(route, route[1:]):
                    out = decode(n, int(by_id[a]["canon"]))
                    found = None
                    for path in paths(out):
                        if blue(out, path):
                            continue
                        target = reverse_path(out, path)
                        target_id = address[canon(target)]
                        if target_id == b:
                            delta = [target[i].bit_count()-out[i].bit_count()
                                     for i in range(n)]
                            need(delta == [int(i == path[-1])-int(i == path[0])
                                           for i in range(n)], "endpoint charge")
                            found = list(path)
                            break
                    need(found is not None, "table edge lacks an actual black path")
                    witnesses.append(found)
                    for path in paths(literal):
                        if not blue(literal, path):
                            target = reverse_path(literal, path)
                            if address[canon(target)] == b:
                                literal_witnesses.append(list(path))
                                literal = target
                                break
                    else:
                        raise RuntimeError("no literal lift")
                common = set.intersection(*[
                    {tuple(sorted(e)) for e in zip(p, p[1:])}
                    for p in literal_witnesses])
                pair_data.append({"source":w, "target":z, "distance":dist[v],
                                  "route_class_ids":route, "black_paths":witnesses,
                                  "literal_paths":literal_witnesses,
                                  "literal_final_out":literal,
                                  "common_path_edges":sorted(common)})
        report = {"n":n, "black_support_nodes":len(graph),
                  "active_rigid_words":len(rigid),
                  "rigid_word_pairs":pair_data}
        reports.append(report)
        print(json.dumps(report), flush=True)
    (Path(__file__).with_suffix(".json")).write_text(
        json.dumps(reports, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
