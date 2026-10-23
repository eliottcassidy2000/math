#!/usr/bin/env python3
"""Exact G25 block transfer, standard-library only; no Git/network operations."""
from collections import Counter
from itertools import product
from math import isqrt
from pathlib import Path
import hashlib
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def graph(n):
    return [(u, v) for u in range(1, n + 1) for v in range(u + 1, n + 1)
            if isqrt(u + v) ** 2 == u + v]


def components(vertices, edges):
    adj = {v: [] for v in vertices}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    remaining = set(vertices)
    result = []
    while remaining:
        todo = [min(remaining)]
        part = set()
        while todo:
            v = todo.pop()
            if v in part:
                continue
            part.add(v)
            todo.extend(adj[v])
        remaining -= part
        result.append(sorted(part))
    return result


def subsets(values):
    for mask in range(1 << len(values)):
        yield [value for i, value in enumerate(values) if mask & (1 << i)]


def profiles(vertices, internal, ports):
    """Enumerate spanning linear forests with named boundary and endpoint tags.

    Each pair in ports is (inside vertex, outside/core vertex). A component
    without a selected boundary incidence is rejected because the core is
    nonempty. Endpoint identities are retained, not merely their number.
    """
    result = []
    for forest in subsets(internal):
        deg = Counter(v for edge in forest for v in edge)
        if any(deg[v] > 2 for v in vertices):
            continue
        parts = components(vertices, forest)
        if len(forest) != len(vertices) - len(parts):
            continue
        for chosen in subsets(list(range(len(ports)))):
            total = deg + Counter(ports[i][0] for i in chosen)
            if any(total[v] not in (1, 2) for v in vertices):
                continue
            ends = [v for v in vertices if total[v] == 1]
            if len(ends) > 2:
                continue
            terminal_pairs = []
            for part in parts:
                q = [i for i in chosen if ports[i][0] in part]
                e = [v for v in ends if v in part]
                if not q:
                    break
                require(len(q) + len(e) == 2, "Component terminal identity")
                terminal_pairs.append({"ports": q, "endpoints": e})
            else:
                result.append({"internal_edges": forest, "selected_ports": chosen,
                               "endpoints": ends, "pairs": terminal_pairs})
    return result


def lifted_edges(blocks, states, core_edges):
    result = list(core_edges)
    for name, state in zip(blocks, states):
        result.extend(state["internal_edges"])
        result.extend(blocks[name]["ports"][i] for i in state["selected_ports"])
    return [tuple(sorted(edge)) for edge in result]


def direct_paths(n, edges, start):
    adj = {v: [] for v in range(1, n + 1)}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    result = []
    def visit(path, seen):
        if len(path) == n:
            result.append(path.copy())
            return
        for nxt in adj[path[-1]]:
            if nxt not in seen:
                seen.add(nxt)
                path.append(nxt)
                visit(path, seen)
                path.pop()
                seen.remove(nxt)
    visit([start], {start})
    return sorted(result)


def main():
    edges = graph(25)
    sets = {"A": [2, 7, 9, 16, 18], "B": [4, 5, 11, 12],
            "C": [1, 3, 6, 8, 10, 17]}
    occupied = set(v for values in sets.values() for v in values)
    core = sorted(set(range(1, 26)) - occupied)
    core_edges = [edge for edge in edges if all(v in core for v in edge)]
    blocks = {}
    for name, vertices in sets.items():
        internal = [edge for edge in edges if all(v in vertices for v in edge)]
        ports = sorted(edge if edge[0] in vertices else edge[::-1]
                       for edge in edges if sum(v in vertices for v in edge) == 1)
        require(all(v in core for u, v in ports), "No interblock edges")
        states = profiles(vertices, internal, ports)
        keys = [json.dumps([s["selected_ports"], s["endpoints"], s["pairs"]],
                           sort_keys=True) for s in states]
        require(len(keys) == len(set(keys)), "Distinct signatures for these trees")
        breakdown = Counter((len(s["endpoints"]), len(s["selected_ports"]),
                             len(s["pairs"])) for s in states)
        anonymous = Counter(tuple(sorted(tuple(pair["ports"]) for pair in s["pairs"]))
                            for s in states)
        blocks[name] = {"vertices": vertices, "internal_edges": internal,
                        "ports": ports, "profile_count": len(states),
                        "anonymous_endpoint_profile_count": len(anonymous),
                        "anonymous_multiplicity_histogram": dict(sorted(Counter(
                            anonymous.values()).items())),
                        "breakdown_endpoint_port_component": [
                            {"parameters": key, "count": value}
                            for key, value in sorted(breakdown.items())],
                        "states": states}

    valid = []
    degree_feasible = 0
    budget_triples = 0
    degree_feasible_disconnected_examples = []
    for ids in product(*(range(block["profile_count"]) for block in blocks.values())):
        states = [block["states"][i] for block, i in zip(blocks.values(), ids)]
        local_ends = [v for state in states for v in state["endpoints"]]
        if len(local_ends) > 2:
            continue
        budget_triples += 1
        virtual = []
        for name, state in zip(blocks, states):
            ports = blocks[name]["ports"]
            for pair in state["pairs"]:
                terminals = [ports[i][1] for i in pair["ports"]]
                terminals.extend(25 + v for v in pair["endpoints"])
                require(len(terminals) == 2, "Two virtual endpoints")
                virtual.append(tuple(terminals))
        virtual_vertices = core + [25 + v for v in local_ends]
        for chosen_core in subsets(core_edges):
            reduced_edges = chosen_core + virtual
            degrees = Counter(v for edge in reduced_edges for v in edge)
            if any(degrees[v] not in (1, 2) for v in core):
                continue
            core_ends = [v for v in core if degrees[v] == 1]
            if len(local_ends) + len(core_ends) != 2:
                continue
            degree_feasible += 1
            if len(components(virtual_vertices, reduced_edges)) != 1:
                if len(degree_feasible_disconnected_examples) < 3:
                    degree_feasible_disconnected_examples.append({
                        "profile_ids": ids, "core_edges": chosen_core,
                        "reduced_components": components(virtual_vertices, reduced_edges)})
                continue
            full_edges = lifted_edges(blocks, states, chosen_core)
            full_paths = direct_paths(25, full_edges, 18)
            require(len(full_paths) == 1, "Unique lifted Hamiltonian path")
            valid.append({"profile_ids": ids, "core_edges": chosen_core,
                          "path": full_paths[0], "other_endpoint": full_paths[0][-1]})

    direct = direct_paths(25, edges, 18)
    require(sorted(item["path"] for item in valid) == direct,
            "Transfer table equals independent DFS census")
    require(len(direct) == 10, "G25 count")
    require([blocks[name]["profile_count"] for name in blocks] == [5, 101, 149],
            "Local profile census")
    require(len(edges) == 32 and len(components(range(1, 26), edges)) == 1,
            "Graph size and connectivity")
    # Contract the entire induced trees, retaining parallel edges. The bare
    # quotient fails the necessary path cut bound despite G25 being traceable.
    owner = {v: 26 + i for i, values in enumerate(sets.values()) for v in values}
    quotient = [tuple(sorted((owner.get(u, u), owner.get(v, v)))) for u, v in edges
                if owner.get(u, u) != owner.get(v, v)]
    quotient_vertices = core + [26, 27, 28]
    require(len(quotient) == 20 and len(set(quotient)) == 18,
            "Parallel edges are part of the cycle-rank sidecar")
    require(len(quotient) - len(quotient_vertices) + 1 == 8,
            "Tree contraction preserves cycle rank as a multigraph")
    require(len(components(core, core_edges)) == 6 > 3 + 1,
            "Bare quotient violates Hamiltonian path cut condition")
    # Minimal nonpath-tree hostile to edge-count-only compression.
    p4 = [(1, 2), (2, 3), (3, 4)]
    star = [(1, 2), (1, 3), (1, 4)]
    require(bool(direct_paths(4, p4, 1)), "P4 positive control")
    require(not any(direct_paths(4, star, v) for v in range(1, 5)), "K1,3 hostile")
    # Four-port pairing is necessary even if selected ports/degrees agree.
    attachment = [(2, 5), (3, 5)]
    require(bool(direct_paths(5, [(1, 2), (3, 4)] + attachment, 1)),
            "Pairing positive control")
    require(not any(direct_paths(5, [(1, 4), (2, 3)] + attachment, v)
                    for v in range(1, 6)), "Pairing hostile")
    # Degree-two endpoints are not forced to use both their edges.
    degree_two_endpoint = next(path for path in direct if path[-1] == 8)
    require(sum(8 in edge for edge in edges) == 2, "Endpoint hostile degree")
    # Every reflection formula is checked on every length-two G25 walk.
    adjacency = {v: [u for u in range(1, 26) if tuple(sorted((u, v))) in edges]
                 for v in range(1, 26)}
    reflection_checks = 0
    for x in adjacency:
        for y in adjacency[x]:
            for z in adjacency[y]:
                a, b = x + y, y + z
                require(z == x + b - a, "Two-reflection translation")
                reflection_checks += 1
    receipt = {"status": "FINITE-EXACT", "n": 25, "edges": edges,
               "edge_count": len(edges), "cycle_rank": len(edges) - 25 + 1,
               "core": core, "core_edges": core_edges, "blocks": blocks,
               "bare_tree_quotient": {"vertices": quotient_vertices, "edges": quotient,
                   "multigraph_edge_count": 20, "simple_edge_count": 18,
                   "core_components_after_deleting_3_block_vertices": 6,
                   "hamiltonian_path_exists": False},
               "endpoint_budget_profile_triples": budget_triples,
               "degree_feasible_gluings": degree_feasible,
               "disconnected_degree_feasible_gluings": degree_feasible - len(valid),
               "disconnected_examples": degree_feasible_disconnected_examples,
               "hamiltonian_path_count_mod_reversal": len(direct),
               "solutions": sorted(valid, key=lambda item: item["other_endpoint"]),
               "degree_two_endpoint8_witness": degree_two_endpoint,
               "reflection_checks": reflection_checks,
               "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    out = Path(__file__).with_name("certificate.json")
    out.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({k: v for k, v in receipt.items() if k not in
                      {"edges", "blocks", "disconnected_examples", "solutions"}}, indent=2))
    print("Local endpoint-marked signature counts:",
          {name: block["profile_count"] for name, block in blocks.items()})
    print("Solution profile IDs:", [(s["other_endpoint"], s["profile_ids"])
                                      for s in receipt["solutions"]])
    print("PASS: exact block-transfer census agrees with independent DFS")


if __name__ == "__main__":
    main()
