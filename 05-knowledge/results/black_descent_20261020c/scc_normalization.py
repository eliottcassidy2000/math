"""Exact audits of SCC minima, descent basins, and a uniform black macro."""
from collections import Counter, defaultdict, deque
import csv
import json
from math import comb
from pathlib import Path
from scc_probe import decode, paths, blue, reverse_path, need

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def strong(out):
    n = len(out)
    for backward in (False, True):
        seen = 1
        while True:
            reach = seen
            for u in range(n):
                if seen >> u & 1:
                    reach |= (sum(1 << v for v in range(n)
                                  if out[v] >> u & 1) if backward else out[u])
            if reach == seen:
                break
            seen = reach
        if seen != (1 << n)-1:
            return False
    return True


def triangles(out):
    return comb(len(out), 3)-sum(comb(row.bit_count(), 2) for row in out)


def macro(prefix, prefix_path):
    """Source A -> x -> C3; target A -> C3 -> c, preserving A literally."""
    m = len(prefix)
    x, a, b, c = range(m, m+4)
    suffix = sum(1 << v for v in (x, a, b, c))
    source = tuple(row | suffix for row in prefix) + (
        (1 << a) | (1 << b) | (1 << c), 1 << b, 1 << c, 1 << a)
    target = tuple(row | suffix for row in prefix) + (
        (1 << a) | (1 << c), (1 << b) | (1 << c),
        (1 << x) | (1 << c), 0)
    A = tuple(prefix_path)
    ps = (A+(x, c, a, b), A[::-1]+(c, x, b, a),
          (x, a, b, c)+A, (c, b, a, x)+A[::-1])
    actual = source
    colours = []
    for p in ps:
        need(sorted(p) == list(range(m+4)), "macro path misses a vertex")
        colours.append("blue" if blue(actual, p) else "black")
        actual = reverse_path(actual, p)
    need(actual == target, "macro changes the wrong arcs")
    return colours


def macro_audit():
    receipt = []
    # The residual B is arbitrary, including nonrigid and self-converse types.
    for k in range(6):
        total = 0
        for code in range(1 << comb(k, 2)):
            B = decode(k, code)
            A = ((1 << (k+2))-2, (1 << (k+2))-4) + tuple(row << 2 for row in B)
            for p in paths(B):
                need(macro(A, (0, 1)+tuple(v+2 for v in p)) == ["black"]*4,
                     "two-singleton guard fails")
                total += 1
        receipt.append({"residual_order":k, "labelled_tournaments":1 << comb(k,2),
                        "prefix_path_certificates":total})
    need(macro((0,), (0,)) == ["black", "black", "black", "blue"],
         "one-singleton hostile changed")
    return receipt


def census_audit():
    receipt = []
    for n in range(3, 8):
        stem = ROOT/"04-computation/experiments"/"tiling_line_metagraph_census_20261020"
        with Path(f"{stem}_classes_n{n}.tsv").open() as f:
            data = list(csv.DictReader(f, delimiter="\t"))
        by_id = {int(r["class"]):r for r in data}
        merged = {i:min(i, int(r["sigmaClass"])) for i,r in by_id.items()}
        reps = sorted(set(merged.values()))
        out = {i:decode(n, int(by_id[i]["canon"])) for i in reps}
        is_strong = {i:strong(out[i]) for i in reps}
        potential = {i:triangles(out[i]) for i in reps}
        all_neighbours = defaultdict(set)
        black_neighbours = defaultdict(set)
        with Path(f"{stem}_lines_n{n}.tsv").open() as f:
            for r in csv.DictReader(f, delimiter="\t"):
                u,v = merged[int(r["classA"])],merged[int(r["classB"])]
                all_neighbours[u].add(v)
                all_neighbours[v].add(u)
                if int(r["black"]):
                    black_neighbours[u].add(v)
                    black_neighbours[v].add(u)
        full_minima = {u for u in reps if all(potential[v] > potential[u]
                                            for v in all_neighbours[u])}
        need(full_minima == {u for u in reps if not is_strong[u]},
             "full strict minima differ from SCC periphery")
        active = set(black_neighbours)
        terminals = {u for u in active if not is_strong[u]}
        next_step = {}
        for u in active-terminals:
            lower = [v for v in black_neighbours[u] if potential[v] < potential[u]]
            need(lower, "strong active node has no strict black descent")
            next_step[u] = min(lower, key=lambda v:(potential[v],v))
        basin = {}
        for u in sorted(active, key=lambda v:potential[v]):
            basin[u] = u if u in terminals else basin[next_step[u]]
        quotient = defaultdict(set)
        for u in active:
            for v in black_neighbours[u]:
                quotient[basin[u]].add(basin[v])
        unseen = set(terminals)
        components = []
        while unseen:
            found = {min(unseen)}
            queue = deque(found)
            while queue:
                u = queue.popleft()
                for v in quotient[u]-found:
                    found.add(v)
                    queue.append(v)
            unseen -= found
            components.append(len(found))
        # Type-count recurrence, separate from the fixed-path tiling recurrence.
        a = len(data)
        b = sum(int(r["SC"]) for r in data)
        s = sum(strong(decode(n,int(r["canon"]))) for r in data)
        r = sum(int(row["SC"]) and strong(decode(n,int(row["canon"])))
                for row in data)
        p = sum(int(row["tc"]) == int(row["bm"]) for row in data)
        expected = (a-s+b-r)//2-p+int(n in (3,5))
        need(expected == len(terminals), "terminal type-count formula")
        receipt.append({"n":n, "ordinary_types":a, "self_converse_types":b,
                        "strong_types":s, "strong_self_converse_types":r,
                        "pure_blue_types":p, "full_local_minima":len(full_minima),
                        "active_black_nodes":len(active), "terminal_nodes":len(terminals),
                        "terminal_quotient_component_sizes":components,
                        "basin_sizes":sorted(Counter(basin.values()).values())})
    return receipt


def main():
    result = {"uniform_macro":macro_audit(), "complete_census":census_audit()}
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
