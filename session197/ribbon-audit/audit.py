"""Independent audit: concrete graph merges, literal covers, bounded closures.

The source is copied read-only by convention into source/. No source outputs
are rewritten. This audit intentionally does not call transfer.transitions or
transfer.advance while building its independent automata and closures.
"""
import collections
import itertools
import json
import pathlib
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "source"))
import transfer as src


def graph_step(state, width, letter):
    # Lift each abstract component to an actual directed graph. Distinct
    # forgotten endpoints receive distinct negative vertex IDs. A forgotten
    # singleton may equally be lifted as a longer path: only its component
    # identity and free endpoints affect future attachments.
    vertex = set()
    edges = []
    counter = -1
    for head, tail in state:
        if head < 0:
            head = counter
            counter -= 1
        if tail < 0:
            tail = counter
            counter -= 1
        vertex.update((head, tail))
        if head != tail:
            edges.append((head, tail))
    starts = vertex - {b for a, b in edges}
    ends = vertex - {a for a, b in edges}
    incoming = [None] + [u for u in sorted(ends) if u < 0 or not (letter & (1 << u))]
    outgoing = [None] + [u for u in sorted(starts) if u >= 0 and letter & (1 << u)]
    result = collections.Counter()
    for left, right in itertools.product(incoming, outgoing):
        es = list(edges)
        if left is not None:
            es.append((left, width))
        if right is not None:
            es.append((width, right))
        vs = vertex | {width}
        successors = dict(es)
        assert len(successors) == len(es)
        predecessors = {v: u for u, v in es}
        assert len(predecessors) == len(es)
        heads = vs - predecessors.keys()
        covered = set()
        components = []
        for h in heads:
            u = h
            while u not in covered:
                covered.add(u)
                if u not in successors:
                    break
                u = successors[u]
            components.append((max(-1, h - 1), max(-1, u - 1)))
        if covered != vs:  # Any vertices missed are in a directed cycle.
            continue
        if sum(h < 0 for h, t in components) >= 2:
            continue
        result[tuple(sorted(components))] += 1
    return result


def independent_automaton(w):
    todo = [()]
    found = {()}
    transitions = {}
    for st in todo:
        for a in range(1 << w):
            row = graph_step(st, w, a)
            transitions[st, a] = row
            for nx in row:
                if nx not in found:
                    todo.append(nx)
                    found.add(nx)
    # Deliberately choose a different coordinate order from source automaton.
    states = sorted(found, reverse=True)
    index = {st: i for i, st in enumerate(states)}
    matrices = [
        [[(index[nx], c) for nx, c in transitions[st, a].items()] for st in states]
        for a in range(1 << w)
    ]
    return states, matrices


def closure(w, M):
    states, matrices = independent_automaton(w)
    size = len(states)
    final = [i for i, s in enumerate(states) if len(s) == 1]
    init = tuple(int(s == ()) for s in states)
    seen = {init: ()}
    queue = collections.deque([init])
    counts = {}
    processed = 0
    rejected = 0
    t = time.time()
    while queue:
        vector = queue.popleft()
        processed += 1
        word = seen[vector]
        val = sum(vector[i] for i in final)
        if val:
            counts.setdefault(val, word)
        support = [(i, n) for i, n in enumerate(vector) if n]
        for a, matrix in enumerate(matrices):
            next_vector = [0] * size
            for i, n in support:
                for j, c in matrix[i]:
                    next_vector[j] += n * c
            next_vector = tuple(min(n, M + 1) for n in next_vector)
            if sum(next_vector[i] for i in final) > M:
                rejected += 1
                continue
            if next_vector not in seen:
                seen[next_vector] = word + (a,)
                queue.append(next_vector)
    return {
        "states": len(seen), "processed": processed,
        "rejected_transitions": rejected, "complete": not queue,
        "counts": {str(h): list(wd) for h, wd in sorted(counts.items())},
        "odd_gaps": [h for h in range(1, M + 1, 2) if h not in counts],
        "even_outputs": [h for h in counts if not h % 2],
        "seconds": time.time() - t,
    }


def literal_covers(n):
    covers = set()
    for p in itertools.permutations(range(n)):
        for cutset in range(1 << (n - 1)):
            pieces = []
            start = 0
            for i in range(n - 1):
                if cutset & (1 << i):
                    pieces.append(p[start:i + 1])
                    start = i + 1
            pieces.append(p[start:])
            covers.add(tuple(sorted(pieces)))
    return covers


def literal_vector(graph, w, covers):
    n = len(graph)
    vector = collections.Counter()
    for pieces in covers:
        if not all(graph[u][v] for piece in pieces for u, v in zip(piece, piece[1:])):
            continue
        state = tuple(sorted((max(-1, piece[0] - n + w), max(-1, piece[-1] - n + w)) for piece in pieces))
        if sum(s < 0 for s, e in state) <= 1:
            vector[state] += 1
    return vector


def literal_count(graph):
    return sum(all(graph[a][b] for a, b in zip(p, p[1:])) for p in itertools.permutations(range(len(graph))))


def sccs(graph):
    n = len(graph)
    reach = [set([i]) | {j for j in range(n) if graph[i][j]} for i in range(n)]
    for k in range(n):
        for i in range(n):
            if k in reach[i]:
                reach[i].update(reach[k])
    todo = set(range(n))
    ans = []
    while todo:
        a = min(todo)
        component = sorted(b for b in todo if b in reach[a] and a in reach[b])
        ans.append(component)
        todo.difference_update(component)
    return ans


def main():
    report = {}
    for w in range(1, 5):
        states, matrices = independent_automaton(w)
        src_states, src_mat = src.automaton(w)
        assert set(states) == set(src_states)
        checked = 0
        for i, st in enumerate(src_states):
            for letter in range(1 << w):
                original = collections.Counter({src_states[j]: c for j, c in src_mat[i][letter]})
                assert graph_step(st, w, letter) == original, (w, st, letter)
                checked += 1
        report[f"independent_span{w}_transition_rows"] = {"states": len(states), "rows": checked}
    print("Every transition independently reconstructed through width 4", flush=True)

    covers = literal_covers(5)
    report["literal_directed_linear_forests_order5"] = len(covers)
    for w in range(1, 5):
        states, mat = src.automaton(w)
        checked = 0
        for wd in itertools.product(*(range(0, 1 << w, 1 << max(0, w - j)) for j in range(5))):
            vector = (1,) + (0,) * (len(states) - 1)
            for letter in wd:
                vector = src.advance(vector, mat, letter)
            abstract = collections.Counter({st: x for st, x in zip(states, vector) if x})
            literal = literal_vector(src.graph(wd, w), w, covers)
            assert abstract == literal, (w, wd, abstract - literal, literal - abstract)
            checked += 1
        report[f"literal_all_components_span{w}_order5"] = checked
    print("Literal path-cover vectors agree on every ordered tournament of order 5, all widths", flush=True)

    for w, M in [(1, 101), (2, 47), (2, 1001), (3, 101)]:
        result = closure(w, M)
        assert result["complete"] and result["processed"] == result["states"]
        if w > 1:
            prior = json.loads((ROOT / "source" / f"closure-{w}-{M}.json").read_text())
            assert result["states"] == prior["states"]
            assert set(result["counts"]) == set(prior["counts"])
            for value, word in prior["counts"].items():
                assert len(word) == len(result["counts"][value])
                # Distinct independent Held-Karp via predecessor removal.
                a = src.graph(tuple(map(int, word)), w)
                from functools import lru_cache
                @lru_cache(None)
                def ends(subset, v):
                    if subset == (1 << v):
                        return 1
                    rest = subset ^ (1 << v)
                    return sum(ends(rest, u) for u in range(len(a)) if rest & (1 << u) and a[u][v])
                assert sum(ends((1 << len(a)) - 1, v) for v in range(len(a))) == int(value)
        (ROOT / f"independent-closure-{w}-{M}.json").write_text(json.dumps(result, indent=2) + "\n")
        report[f"independent_closure_span{w}_cap{M}"] = {k: v for k, v in result.items() if k != "counts"}
        print(w, M, report[f"independent_closure_span{w}_cap{M}"], flush=True)

    witness = tuple(map(int, "0001321"))
    graph = src.graph(witness, 3)
    assert literal_count(graph) == 47
    assert sccs(graph) == [list(range(7))]
    report["literal_47_witness"] = {"word": list(witness), "count": 47, "sccs": sccs(graph)}

    # Find and record the smallest natural-order noninterval SCC as a warning
    # against replacing graph reachability by uncovered natural-order cuts.
    found = None
    for n in range(2, 7):
        for word in itertools.product(*(range(0, 4, 1 << max(0, 2 - j)) for j in range(n))):
            components = sccs(src.graph(word, 2))
            if any(c != list(range(min(c), max(c) + 1)) for c in components):
                found = {"word": list(word), "sccs_1based": [[i + 1 for i in c] for c in components]}
                break
        if found:
            break
    report["noninterval_scc_counterexample"] = found
    (ROOT / "AUDIT.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
