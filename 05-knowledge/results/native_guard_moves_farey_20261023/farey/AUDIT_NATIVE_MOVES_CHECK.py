"""Independent literal-path census and threshold BFS for the native move receipt."""
import itertools
import json
from pathlib import Path

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def components(paths, threshold):
    remaining = set(range(len(paths)))
    groups = []
    while remaining:
        start = min(remaining)
        remaining.remove(start)
        todo, found = [start], {start}
        while todo:
            u = todo.pop()
            neighbors = [v for v in remaining
                         if sum(x != y for x, y in zip(paths[u], paths[v])) <= threshold]
            for v in neighbors:
                remaining.remove(v)
                found.add(v)
                todo.append(v)
        groups.append(found)
    return groups


def main():
    packet = Path(__file__).resolve().parents[1]
    receipt = json.loads((packet / "native_moves/receipt.json").read_text())
    fibre_count, census = 0, []
    for source in receipt["sources"]:
        pairs = list(itertools.combinations(range(8), 2))
        code = source["canon"]
        adjacency = [[False]*8 for _ in range(8)]
        for k, (u, v) in enumerate(pairs):
            adjacency[u][v] = bool(code & (1 << (27-k)))
            adjacency[v][u] = not adjacency[u][v]
        literal = {p for p in itertools.permutations(range(8))
                   if all(adjacency[u][v] for u, v in zip(p, p[1:]))}
        rows = source["path_rows"]
        require(literal == {tuple(row["path"]) for row in rows},
                "independent complete Hamiltonian path census")
        for fibre in source["endpoint_phase_fibres"] + source["nontrivial_indexed_arc_fibres"]:
            fibre_count += 1
            key, ids = fibre["pins_and_phase"], fibre["path_indices"]
            paths = [rows[i]["path"] for i in ids]
            expected = []
            for i, row in enumerate(rows):
                path = row["path"]
                if [path[0], path[-1], row["phase_mod8"]] != key[:3]:
                    continue
                if len(key) > 3 and path[key[3]:key[3]+2] != key[4:6]:
                    continue
                expected.append(i)
            require(ids == expected, "complete accepted fibre")
            budget = next(b for b in range(9) if len(components(paths, b)) == 1)
            require(budget == fibre["sharp_clock_budget"], "BFS exact threshold versus Kruskal")
            below = components(paths, budget-1)
            require(sorted(len(c) for c in below) == fibre["component_sizes_at_budget_minus_one"],
                    "BFS lower partition sizes")
            require(sorted(sorted(ids[i] for i in c) for c in below)
                    == fibre["components_at_budget_minus_one"], "BFS lower partition owners")
        census.append([source["merged_source"], len(rows)])
    print(json.dumps({"status": "PASS", "checks": CHECKS,
                      "fibres": fibre_count, "path_census": census}, indent=2))


if __name__ == "__main__":
    main()
