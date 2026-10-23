"""Exact binary-packet assembly, implication-block moves, and native all-size seeds.

No imports from another checker. Normal and -O runs compare a frozen receipt.
"""
from collections import deque
from itertools import combinations, product, permutations
from pathlib import Path
import hashlib
import json
import random
import sys

CHECKS = 0


def require(condition, why):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(why)


def accepted(bits, forbidden):
    return all(not (bits[i] == a and bits[j] == b) for i, a, j, b in forbidden)


def implication_model(n, forbidden, weights):
    """Literal (i,b) means packet i is in state b; forbidden pairs give 2-CNF."""
    edges = [set() for _ in range(2*n)]
    for i, a, j, b in forbidden:
        edges[2*i+a].add(2*j+1-b)
        edges[2*j+b].add(2*i+1-a)
    reach = []
    for root in range(2*n):
        seen, todo = {root}, [root]
        while todo:
            for v in edges[todo.pop()]:
                if v not in seen:
                    seen.add(v)
                    todo.append(v)
        reach.append(seen)
    classes, owner = [], {}
    for u in range(2*n):
        if u in owner:
            continue
        members = tuple(v for v in range(2*n) if v in reach[u] and u in reach[v])
        for v in members:
            owner[v] = len(classes)
        classes.append(members)
    if any(owner[2*i] == owner[2*i+1] for i in range(n)):
        return {"satisfiable": False}
    opposite = {c: owner[members[0] ^ 1] for c, members in enumerate(classes)}
    blocks, seen = [], set()
    for c, members in enumerate(classes):
        if c in seen:
            continue
        d = opposite[c]
        seen.update((c, d))
        variables = sorted({v//2 for v in members})
        free = (classes[d][0] not in reach[members[0]]
                and members[0] not in reach[classes[d][0]])
        blocks.append({"classes": (c, d), "variables": variables,
                       "weight": sum(weights[i] for i in variables), "free": free})
    return {"satisfiable": True, "classes": classes, "owner": owner,
            "opposite": opposite, "reach": reach, "blocks": blocks,
            "budget": max((b["weight"] for b in blocks if b["free"]), default=0)}


def block_path(start, target, model, forbidden):
    """Remove an implication-minimal true SCC that must become false."""
    current, answer = tuple(start), [tuple(start)]
    target_true = {model["owner"][2*i+b] for i, b in enumerate(target)}
    while current != tuple(target):
        true = {model["owner"][2*i+b] for i, b in enumerate(current)}
        difference = true-target_true
        roots = {c: model["classes"][c][0] for c in difference}
        minimal = [c for c in sorted(difference)
                   if all(roots[c] not in model["reach"][roots[d]]
                          for d in difference if d != c)]
        require(bool(minimal), "nonempty difference has an implication-minimal SCC")
        chosen = minimal[0]
        block = next(b for b in model["blocks"] if chosen in b["classes"])
        next_state = list(current)
        for i in block["variables"]:
            next_state[i] ^= 1
        next_state = tuple(next_state)
        require(accepted(next_state, forbidden), "entire SCC-pair flip stays accepted")
        require(all(next_state[i] == target[i] for i in block["variables"]),
                "chosen block reaches its final value without recycling")
        answer.append(next_state)
        current = next_state
    return answer


def literal_budget(states, weights):
    if len(states) <= 1:
        return 0
    distances = {(i, j): sum(w for a, b, w in zip(states[i], states[j], weights) if a != b)
                 for i, j in combinations(range(len(states)), 2)}
    for bound in sorted(set(distances.values())):
        seen, todo = {0}, [0]
        while todo:
            u = todo.pop()
            for v in range(len(states)):
                if v not in seen and distances[tuple(sorted((u, v)))] <= bound:
                    seen.add(v)
                    todo.append(v)
        if len(seen) == len(states):
            return bound
    raise RuntimeError("complete accepted graph was disconnected")


def check_formula(n, forbidden, weights):
    states = [x for x in product(range(2), repeat=n) if accepted(x, forbidden)]
    model = implication_model(n, forbidden, weights)
    require(model["satisfiable"] == bool(states), "SCC obstruction versus literal satisfiability")
    if not states:
        return False
    require(model["budget"] == literal_budget(states, weights), "sharp SCC weight versus independent threshold BFS")
    for block in model["blocks"]:
        i = block["variables"][0]
        require(block["free"] == (len({s[i] for s in states}) == 2), "exact frozen/free SCC criterion")
    pairs = list(product(states, repeat=2))
    if len(pairs) > 64:
        pairs = pairs[::max(1, len(pairs)//64)]
    for start, end in pairs:
        path = block_path(start, end, model, forbidden)
        total = sum(sum(w for a, b, w in zip(x, y, weights) if a != b)
                    for x, y in zip(path, path[1:]))
        endpoint = sum(w for a, b, w in zip(start, end, weights) if a != b)
        require(total == endpoint, "monotone coupled path is additive-cost geodesic")
        require(all(sum(w for a, b, w in zip(x, y, weights) if a != b) <= model["budget"]
                    for x, y in zip(path, path[1:])), "all proposed move costs obey sharp budget")
    return True


def compile_packets(packets, used=()):
    """Each record keeps actual owner values and a set of literal port identities."""
    forbidden = []
    used = set(used)
    for i, packet in enumerate(packets):
        for b in range(2):
            if b >= len(packet) or (b < len(packet) and used.intersection(packet[b]["ports"])):
                forbidden.append((i, b, i, b))
    for i, j in combinations(range(len(packets)), 2):
        for a, left in enumerate(packets[i]):
            for b, right in enumerate(packets[j]):
                shared = set(left["owners"]) & set(right["owners"])
                mismatch = any(left["owners"][v] != right["owners"][v] for v in shared)
                collision = set(left["ports"]) & set(right["ports"])
                if mismatch or collision:
                    forbidden.append((i, a, j, b))
    return forbidden


def direct_packet_acceptance(packets, bits, used):
    owners, inventory = {}, set(used)
    for packet, bit in zip(packets, bits):
        if bit >= len(packet):
            return False
        record = packet[bit]
        for variable, value in record["owners"].items():
            if variable in owners and owners[variable] != value:
                return False
            owners[variable] = value
        if inventory.intersection(record["ports"]):
            return False
        inventory.update(record["ports"])
    return True


def abstract_controls():
    formulas = satisfiable = 0
    clauses2 = [(i, a, i, a) for i in range(2) for a in range(2)]
    clauses2 += [(0, a, 1, b) for a, b in product(range(2), repeat=2)]
    for mask in range(1 << len(clauses2)):
        f = [c for j, c in enumerate(clauses2) if mask >> j & 1]
        satisfiable += check_formula(2, f, [1, 3])
        formulas += 1
    clauses3 = [(i, a, i, a) for i in range(3) for a in range(2)]
    clauses3 += [(i, a, j, b) for i, j in combinations(range(3), 2)
                 for a, b in product(range(2), repeat=2)]
    for size in range(5):
        for f in combinations(clauses3, size):
            satisfiable += check_formula(3, f, [2, 1, 3])
            formulas += 1
    rng = random.Random(20261023)
    for n in (4, 5, 6):
        clauses = [(i, a, i, a) for i in range(n) for a in range(2)]
        clauses += [(i, a, j, b) for i, j in combinations(range(n), 2)
                    for a, b in product(range(2), repeat=2)]
        for _ in range(160):
            f = rng.sample(clauses, rng.randrange(1, min(len(clauses), 15)))
            satisfiable += check_formula(n, f, [1+rng.randrange(4) for _ in range(n)])
            formulas += 1
    for _ in range(400):
        n = rng.randrange(1, 7)
        packets = []
        for i in range(n):
            observed = [v for v in range(3) if rng.randrange(2)]
            packets.append([{"owners": {v: rng.randrange(3) for v in observed},
                             "ports": rng.sample(range(9), rng.randrange(1, 4))}
                            for _ in range(1 if rng.randrange(8) == 0 else 2)])
        used = rng.sample(range(9), rng.randrange(3))
        f = compile_packets(packets, used)
        for bits in product(range(2), repeat=n):
            require(accepted(bits, f) == direct_packet_acceptance(packets, bits, used),
                    "lossless shared-coordinate and resource compiler")
    # Each source admits an injection, but their required shared owner values conflict.
    def rec(value, ports):
        return {"owners": {"shared_owner": value}, "ports": ports}
    packets = [[rec(0, ["p"]), rec(1, ["q"])],
               [{"owners": {}, "ports": ["p"]}],
               [rec(0, ["r"]), rec(1, ["s"])],
               [{"owners": {}, "ports": ["s"]}]]
    require(any(direct_packet_acceptance(packets[:2], b, ()) for b in product(range(2), repeat=2)),
            "source A separately has an injective matching")
    require(any(direct_packet_acceptance(packets[2:], b, ()) for b in product(range(2), repeat=2)),
            "source B separately has an injective matching")
    require(not any(direct_packet_acceptance(packets, b, ()) for b in product(range(2), repeat=4)),
            "source matchings fail after shared owner identification")
    return {"formulas": formulas, "satisfiable_formulas": satisfiable,
            "random_full_record_packets": 400,
            "hostile": "A requires shared owner1 while B requires shared owner0"}


CORE_CODE = 29558785
P0 = (3, 0, 6, 7, 1, 5, 4, 2)
P1 = (3, 1, 6, 7, 0, 5, 4, 2)
G = (4, 0, 2, 3, 1, 5, 6, 7)


def core_matrix():
    A = [[0]*8 for _ in range(8)]
    for bit, (i, j) in enumerate(combinations(range(8), 2)):
        A[i][j] = (CORE_CODE >> (27-bit)) & 1
        A[j][i] = 1-A[i][j]
    return A


def partner(A, path):
    retained = set(zip(path, path[1:]))
    return [[0 if i == j else (A[i][j] if (i, j) in retained or (j, i) in retained else A[j][i])
             for j in range(len(A))] for i in range(len(A))]


def strongly_connected(A):
    for root in range(len(A)):
        seen, todo = {root}, [root]
        while todo:
            u = todo.pop()
            for v in range(len(A)):
                if A[u][v] and v not in seen:
                    seen.add(v)
                    todo.append(v)
        if len(seen) != len(A):
            return False
    return True


def normalized_word(A, path):
    chords = [(i, j) for i in range(len(A)) for j in range(i+2, len(A))]
    return sum(A[path[i]][path[j]] << b for b, (i, j) in enumerate(chords))


def native_controls():
    A = core_matrix()
    paths = [p for p in permutations(range(8)) if all(A[i][j] for i, j in zip(p, p[1:]))]
    require(len(paths) == 409, "independent core full Hamiltonian census")
    require(sum(all(A[i][j] == A[g[i]][g[j]] for i, j in combinations(range(8), 2))
                for g in permutations(range(8))) == 1, "trivial actual core automorphism kernel")
    cycle = (0, 4, 1, 5, 2, 6, 7, 3)
    require(all(A[i][j] for i, j in zip(cycle, cycle[1:]+cycle[:1])), "explicit strong core cycle")
    target_cycle = (0, 2, 3, 4, 6, 7, 1, 5)
    K0 = partner(A, P0)
    require(all(K0[i][j] for i, j in zip(target_cycle, target_cycle[1:]+target_cycle[:1])),
            "explicit strong complement-target core cycle")
    eligible = [p for p in paths if all(p[i] == P0[i] for i in range(8) if i not in (1, 4))]
    require(eligible == [P0, P1], "complete two-state eligibility comes from literal pins")
    hosts = 0
    summaries = []
    for k in range(1, 9):
        n, L = 2*k+7, 2*k-1
        packet_states = []
        profiles = []
        source_words = [[], []]
        full = (1 << ((n-1)*(n-2)//2))-1
        for r in range(k):
            pre, post = tuple(range(8, 8+r)), tuple(range(8+r, n))
            order = pre+tuple(range(8))+post
            rank = {v: j for j, v in enumerate(order)}
            B = [[0 if i == j else A[i][j] if i < 8 and j < 8 else int(rank[i] < rank[j])
                  for j in range(n)] for i in range(n)]
            states = [pre+p+post for p in (P0, P1)]
            g = G+tuple(range(8, n))
            K = [partner(B, p) for p in states]
            require(all(K[0][i][j] == K[1][g[i]][g[j]] for i, j in combinations(range(n), 2)),
                    "literal endpoint-fixed target isomorphism after padding")
            require(not strongly_connected(B), "padded source is decomposable")
            for bit, (p, Kp) in enumerate(zip(states, K)):
                require(all(B[i][j] for i, j in zip(p, p[1:])), "actual padded path")
                converse = [[B[j][i] for j in range(n)] for i in range(n)]
                require(partner(converse, p[::-1]) == [list(row) for row in zip(*Kp)],
                        "literal converse-partner square")
                require(strongly_connected(Kp), "padded inverse target is strong")
                word, converse_word = normalized_word(B, p), normalized_word(converse, p[::-1])
                source_words[bit].extend((word, converse_word))
                require(normalized_word(Kp, p) == (full ^ word), "inverse endpoint word is literal complement")
                require(p[r+2:r+4] == (6, 7) and p[r+5:r+7] == (5, 4),
                        "retained actual shared arcs at shifted clocks")
            require(sum(x != y for x, y in zip(*states)) == 2, "two changed actual owner clocks per source")
            # All other clocks are literally pinned: the two remaining owners have exactly two orders.
            only = []
            for a, b in ((0, 1), (1, 0)):
                q = list(states[0]);q[r+1], q[r+4] = a, b
                if all(B[i][j] for i, j in zip(q, q[1:])):
                    only.append(tuple(q))
            require(only == states, "exact padded two-hole eligibility")
            profiles.append((r, 8, L-r))
            packet_states.append(states)
            hosts += 1
        require(len(set(profiles)) == k, "ordered SCC profiles separate sources")
        require(not set(profiles).intersection((s, m, r) for r, m, s in profiles),
                "sources remain distinct even after converse")
        legal = []
        for bits in product(range(2), repeat=k):
            values = [packet_states[i][bits[i]][i+1] for i in range(k)]
            if len(set(values)) == 1:
                legal.append(bits)
        require(legal == [(0,)*k, (1,)*k], "connected shared owners leave exactly two native global states")
        require(sum(a != b for a, b in zip(*legal)) == k, "sharp native simultaneous slot budget")
        reservoir = set()
        for words in source_words:
            inverse_words = [full ^ word for word in words]
            require(len(set(words)) == len(set(inverse_words)) == 2*k, "source and inverse target injection")
            require(set(words).isdisjoint(inverse_words), "source and inverse endpoint inventories disjoint")
            reservoir.update(words)
            reservoir.update(inverse_words)
        require(len(reservoir) == 8*k, "complete two-state converse-closed endpoint reservoir is injective")
        summaries.append({"sources": k, "order": n, "legal_packets": 2,
                          "sharp_source_slot_budget": k, "actual_clock_occurrences": 2*k,
                          "distinct_changed_shared_variables": 2, "reserved_endpoint_ports": 8*k})
    # Three native binary slots can satisfy every separate edge and still have no common assembly.
    parity_constraints = [(0, 1, 0), (1, 2, 0), (2, 0, 1)]
    require(all(any((x[i] ^ x[j]) == d for x in product(range(2), repeat=3))
                for i, j, d in parity_constraints), "each native equality constraint separately feasible")
    require(not any(all((x[i] ^ x[j]) == d for i, j, d in parity_constraints)
                    for x in product(range(2), repeat=3)), "native shared-clock odd holonomy obstructs global assembly")
    return {"core_code": CORE_CODE, "core_paths": 409, "core_automorphisms": 1,
            "strong_cycle": cycle, "strong_target_cycle": target_cycle,
            "hosts": hosts, "families": summaries}


def main():
    here = Path(__file__).resolve().parent
    results = here.parents[1]
    source_paths = ["session44_guarded_native_moves_20261023/coupled-moves/THEOREM.md",
                    "session44_guarded_native_moves_20261023/stable-padding/THEOREM.md"]
    hashes = {p: hashlib.sha256((results/p).read_bytes()).hexdigest() for p in source_paths}
    result = {"status": "PROVED and FINITE-EXACT on declared carriers",
              "source_sha256": hashes, "abstract": abstract_controls(), "native": native_controls()}
    result["checks"] = CHECKS
    output = here/"receipt.json"
    rendered = json.dumps(result, indent=2)+"\n"
    if "--write" in sys.argv:
        output.write_text(rendered, encoding="utf-8")
    else:
        require(json.loads(output.read_text(encoding="utf-8")) == json.loads(rendered), "frozen assembly receipt")
    print(json.dumps({"status": "PASS", "checks": CHECKS,
                      "formulas": result["abstract"]["formulas"], "native_hosts": result["native"]["hosts"]}))


if __name__ == "__main__":
    main()
