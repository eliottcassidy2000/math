"""Explicit terminal certificates complete one frozen finite Collatz atlas.

Production functions only verify supplied words and retain observed edges.
The separately labelled discover_selected function is the sole orbit discovery
routine. Frozen terminal words below were discovered only at the 37 selected
missing-edge sinks; they are data, never assumed ROOT claims.
"""
from collections import defaultdict, deque
from dataclasses import dataclass, field
from hashlib import sha256
import json

import collatz_uncovered_join_routes_20261007 as routes


CHECKS = 0
LIMIT = 32767
REQUESTED = tuple(range(1, LIMIT+1, 2))


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def odd(n):
    need(type(n) is int and n > 0 and n % 2 == 1, 'exact positive odd source')


def valuation_word(word):
    need(type(word) is tuple and all(type(a) is int and a > 0 for a in word),
         'tuple of exact positive valuation letters')


def literal_step(n):
    odd(n)
    value, a = 3*n+1, 0
    while value % 2 == 0:
        value //= 2
        a += 1
    return value, a


def replay(source, word):
    odd(source)
    valuation_word(word)
    states, current = [source], source
    for a in word:
        need(current != 1, 'strict first-hit word has no ROOT padding')
        current, actual = literal_step(current)
        need(actual == a, 'independent while-even replay authenticates valuation')
        states.append(current)
    return tuple(states)


@dataclass
class Graph:
    edges: dict = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)
    submissions: int = 0

    def insert(self, source, word, origin):
        """Atomic insertion: a corrupt word cannot partially mutate the graph."""
        path = replay(source, word)
        for n, target, a in zip(path, path[1:], word):
            need(n not in self.edges or self.edges[n] == (target, a),
                 'deterministic edge compatibility')
        for i, (n, target, a) in enumerate(zip(path, path[1:], word)):
            self.edges[n] = (target, a)
            self.provenance.setdefault(n, set()).add((*origin, i))
            self.submissions += 1
        return path[-1]

    def copy(self):
        return Graph(dict(self.edges), {n: set(p) for n, p in self.provenance.items()},
                     self.submissions)

    def vertices(self):
        return {1, *self.edges, *(target for target, _ in self.edges.values())}


def build_observations():
    """Frozen B8 plus two inverse rows and paid final-frontier inspections.

    Exactly one selector call per requested source. No orbit discovery, seeded
    ROOT words, learned rows, or unrelated extra policy is consulted here.
    """
    graph, dependencies, kinds = Graph(), {}, defaultdict(int)
    for n in REQUESTED:
        selected = routes.prior.select(n, 8, 'adaptive')
        need(graph.insert(n, selected.word, ('B8-source', n)) == selected.endpoint,
             'selected endpoint is authenticated')
        kinds[selected.status] += 1
        if selected.kind == 'common_future':
            for child in selected.dependencies:
                if child == selected.endpoint or child == 1:
                    continue
                target, a = literal_step(child)
                need(target == selected.endpoint, 'advertised sibling joins endpoint')
                graph.insert(child, (a,), ('B8-child', n, child))
        added = routes.new_receipts(n, True, selected) if selected.status == 'PENDING' else ()
        for rec in added:
            left = graph.insert(rec.source, rec.source_word, ('new-source', n, rec.child))
            right = graph.insert(rec.child, rec.child_word, ('new-child', n, rec.child))
            need(left == right == rec.endpoint and 0 < rec.child < n,
                 'new receipt joins exact strictly smaller source')
        children = selected.dependencies if selected.status == 'REDUCED' else tuple(z.child for z in added)
        need(all(type(c) is int and c % 2 and 0 < c < n for c in children),
             'the frozen dependency relation is a strict-size DAG')
        dependencies[n] = tuple(sorted(set(children)))
    return graph, dependencies, dict(kinds)


def rooted_ranks(edges):
    """Least rooted closure on a finite relation; this does not validate edges."""
    reverse = defaultdict(list)
    for n, (target, _) in edges.items():
        reverse[target].append(n)
    ranks, queue = {1: 0}, deque([1])
    while queue:
        target = queue.popleft()
        for n in reverse[target]:
            if n not in ranks:
                ranks[n] = ranks[target]+1
                queue.append(n)
    return ranks


def walk(edges, source):
    odd(source)
    path, word, seen = [source], [], set()
    while path[-1] != 1 and path[-1] in edges and path[-1] not in seen:
        n = path[-1]
        seen.add(n)
        target, a = edges[n]
        path.append(target)
        word.append(a)
    return tuple(path), tuple(word), path[-1] in seen


def frontier(graph, requested=REQUESTED):
    demand, cycles = defaultdict(list), []
    ranks = rooted_ranks(graph.edges)
    for n in requested:
        odd(n)
        if n in ranks:
            continue
        path, _, cycle = walk(graph.edges, n)
        if cycle:
            cycles.append(n)
        else:
            need(path[-1] != 1 and path[-1] not in graph.edges, 'missing outgoing edge')
            demand[path[-1]].append(n)
    return dict(demand), tuple(cycles)


def complete(graph, supplied_terminal_words):
    """Verify explicit point receipts, then compile all ROOT-connected words.

    No orbit discovery occurs. Pending vertices remain absent from the result.
    The initial graph is reauthenticated, so forged Graph values are rejected.
    """
    need(isinstance(graph, Graph) and type(supplied_terminal_words) is dict,
         'graph and explicit finite terminal dictionary')
    result = Graph()
    for n, (target, a) in sorted(graph.edges.items()):
        odd(target)
        need(result.insert(n, (a,), ('retained-edge', n)) == target,
             'retained graph target is exact')
        result.provenance[n].update(graph.provenance.get(n, ()))
    for terminal, word in sorted(supplied_terminal_words.items()):
        path = replay(terminal, word)
        need(path[-1] == 1, 'a supplied terminal is proved ROOT, not merely assumed')
        result.insert(terminal, word, ('supplied-terminal', terminal))
    ranks = rooted_ranks(result.edges)
    certificates = {1: ()}
    for n in sorted(ranks, key=lambda z: (ranks[z], z)):
        if n == 1:
            continue
        target, a = result.edges[n]
        need(target in certificates, 'smaller actual ROOT rank')
        certificates[n] = (a,)+certificates[target]
    return result, certificates


def discover_selected(source, cap=10000):
    """EXPERIMENT ONLY: bounded literal search for one explicitly selected sink."""
    odd(source)
    need(type(cap) is int and cap >= 0, 'exact finite discovery cap')
    word, current = [], source
    for _ in range(cap):
        if current == 1:
            return tuple(word)
        current, a = literal_step(current)
        word.append(a)
    return tuple(word) if current == 1 else None


# FROZEN_TERMINALS_BEGIN
TERMINAL_WORDS = {49897: (2, 1, 1, 1, 2, 2, 1, 1, 2, 8, 3, 1, 1, 1, 5, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2, 4, 3, 1, 1,
         5, 4),
 50153: (2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 3, 3, 1, 7, 2, 1, 2, 1, 1, 4, 1, 2, 1, 3, 3, 1, 10),
 57187: (1, 6, 5, 2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2, 4, 3, 1, 1, 5, 4),
 58159: (1, 1, 1, 2, 2, 3, 3, 1, 1, 2, 1, 2, 1, 1, 2, 1, 3, 1, 1, 1, 3, 2, 1, 4, 1, 1, 1, 2, 2, 2, 1, 1, 1, 4, 3, 3,
         1, 3, 6, 2, 1, 3, 1, 2, 3, 4),
 70417: (2, 3, 3, 1, 4, 2, 1, 1, 1, 1, 2, 2, 2, 6, 3, 1, 2, 1, 4, 1, 3, 1, 2, 3, 4),
 91931: (1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 2, 2, 4, 1, 1, 2, 1, 1, 1, 1, 3, 7, 1, 1, 1, 2, 2, 2, 1, 4, 5, 2, 1,
         3, 2, 1, 3, 1, 1, 3, 4, 1, 3, 1, 2, 3, 4),
 98851: (1, 5, 3, 3, 1, 5, 1, 3, 1, 2, 1, 1, 1, 2, 3, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2, 4, 3, 1, 1,
         5, 4),
 109147: (1, 2, 1, 1, 2, 1, 2, 2, 1, 3, 2, 1, 1, 1, 1, 2, 2, 1, 2, 1, 2, 3, 2, 1, 1, 1, 5, 2, 1, 1, 1, 3, 2, 1, 2,
          4, 2, 2, 1, 1, 1, 3, 1, 7, 2, 1, 4, 1, 3, 1, 2, 3, 4),
 114043: (1, 2, 1, 2, 2, 3, 1, 3, 1, 3, 3, 2, 3, 1, 2, 5, 1, 1, 1, 4, 1, 2, 1, 1, 2, 1, 1, 1, 2, 3, 1, 1, 2, 1, 2,
          1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2, 4, 3, 1, 1, 5, 4),
 131041: (2, 2, 1, 1, 1, 2, 1, 1, 4, 1, 2, 1, 2, 2, 1, 1, 1, 1, 3, 1, 2, 1, 1, 1, 1, 7, 1, 1, 3, 2, 1, 1, 1, 3, 1,
          1, 1, 1, 1, 3, 2, 1, 2, 3, 1, 6, 2, 2, 3, 1, 2, 4, 1, 2, 2, 4, 1, 1, 2, 3, 4),
 131387: (1, 2, 1, 6, 2, 1, 3, 2, 2, 1, 5, 1, 2, 1, 1, 2, 6, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 2, 1, 1, 1, 2, 3, 1, 1,
          2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2, 4, 3, 1, 1, 5, 4),
 134951: (1, 1, 2, 1, 3, 1, 6, 1, 1, 1, 2, 2, 3, 3, 4, 1, 4, 1, 1, 1, 1, 1, 1, 2, 4, 3, 3, 3, 1, 2, 3, 4),
 146777: (2, 1, 4, 2, 5, 2, 2, 2, 1, 1, 6, 2, 2, 1, 1, 3, 1, 1, 1, 2, 2, 1, 2, 1, 1, 2, 1, 1, 1, 2, 3, 1, 1, 2, 1,
          2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2, 4, 3, 1, 1, 5, 4),
 164359: (1, 1, 2, 3, 2, 1, 1, 3, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 2, 5, 1, 1, 3, 3, 2, 1, 6, 1, 4, 1, 2, 2, 1, 1,
          2, 2, 3, 1, 1, 1, 1, 4, 1, 2, 1, 1, 2, 1, 1, 1, 2, 3, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2,
          4, 3, 1, 1, 5, 4),
 214111: (1, 1, 1, 1, 4, 1, 1, 2, 1, 1, 3, 1, 6, 1, 2, 1, 2, 2, 1, 3, 2, 7, 1, 1, 1, 1, 1, 2, 4, 3, 3, 3, 1, 2, 3,
          4),
 220249: (2, 1, 4, 1, 2, 1, 2, 2, 4, 3, 3, 1, 2, 3, 1, 2, 2, 1, 1, 1, 1, 2, 1, 2, 2, 5, 3, 1, 1, 1, 2, 2, 1, 2, 1,
          1, 2, 1, 1, 1, 2, 3, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2, 4, 3, 1, 1, 5, 4),
 245935: (1, 1, 1, 2, 3, 1, 1, 2, 1, 1, 4, 1, 3, 3, 2, 5, 2, 2, 1, 2, 2, 1, 3, 1, 3, 1, 3, 2, 3, 1, 2, 1, 4, 1, 3,
          1, 2, 3, 4),
 272857: (2, 1, 6, 5, 4, 3, 1, 1, 1, 1, 4, 1, 2, 1, 1, 2, 1, 1, 1, 2, 3, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1,
          4, 2, 2, 4, 3, 1, 1, 5, 4),
 284305: (2, 3, 2, 1, 4, 1, 7, 2, 3, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2, 4, 3, 1, 1, 5, 4),
 297067: (1, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 4, 4, 2, 1, 3, 1, 1, 5, 2, 2, 6, 2, 1, 1, 2, 2, 3, 1, 1, 1, 1, 4, 1, 2,
          1, 1, 2, 1, 1, 1, 2, 3, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2, 4, 3, 1, 1, 5, 4),
 380743: (1, 1, 2, 2, 1, 1, 4, 3, 1, 2, 4, 2, 3, 1, 2, 2, 2, 1, 1, 1, 2, 2, 3, 1, 1, 3, 2, 1, 2, 3, 2, 2, 1, 3, 3,
          2, 2, 4, 1, 1, 2, 3, 4),
 386839: (1, 1, 5, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 2, 1, 3, 3, 3, 1, 1, 1, 3, 4, 2, 1, 4, 1, 3, 1, 1, 4, 1, 1, 1, 1,
          2, 1, 2, 2, 5, 3, 1, 1, 1, 2, 2, 1, 2, 1, 1, 2, 1, 1, 1, 2, 3, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1,
          4, 2, 2, 4, 3, 1, 1, 5, 4),
 402733: (3, 2, 5, 1, 2, 1, 2, 2, 2, 2, 1, 1, 2, 1, 6, 2, 2, 1, 2, 1, 2, 2, 2, 1, 1, 1, 6, 1, 2, 2, 4, 1, 1, 2, 3,
          4),
 408185: (2, 1, 2, 1, 1, 4, 2, 3, 2, 2, 1, 2, 4, 2, 5, 10),
 446107: (1, 2, 1, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 3, 1, 1, 1, 4, 3, 5, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 3, 4, 3, 1,
          2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 2, 1, 3, 2, 3, 1, 2, 1, 2, 3, 2, 1, 2, 5, 4, 4, 1, 1, 2, 3, 4),
 454531: (1, 4, 3, 1, 4, 1, 2, 2, 2, 10, 2, 1, 3, 1, 2, 3, 4),
 471359: (1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 3, 2, 1, 2, 1, 4, 2, 4, 1, 4, 1, 5, 5, 2, 1, 6, 2, 1, 3, 1, 2, 3, 4),
 552703: (1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 7, 4, 2, 9, 1, 2, 3, 2, 1, 1, 3, 2, 3, 1, 2, 1, 4, 1, 3, 1, 2,
          3, 4),
 765449: (2, 1, 1, 2, 4, 3, 6, 2, 1, 2, 1, 1, 2, 4, 9, 4),
 805949: (3, 1, 1, 8, 2, 2, 1, 4, 3, 2, 1, 6, 2, 1, 3, 1, 2, 3, 4),
 841265: (2, 4, 1, 1, 1, 1, 2, 2, 2, 1, 2, 1, 2, 1, 1, 2, 1, 3, 2, 1, 2, 5, 1, 2, 4, 2, 1, 1, 4, 1, 4, 1, 4, 1, 1,
          1, 1, 1, 1, 2, 4, 3, 3, 3, 1, 2, 3, 4),
 1072397: (3, 4, 1, 2, 3, 1, 4, 1, 3, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 6, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2, 4, 3,
           1, 1, 5, 4),
 1378889: (2, 1, 1, 3, 1, 4, 4, 1, 1, 1, 1, 1, 2, 1, 2, 1, 3, 7, 1, 1, 1, 2, 2, 7, 3, 2, 2, 4, 1, 1, 2, 3, 4),
 1421995: (1, 2, 2, 2, 2, 3, 3, 1, 2, 1, 3, 1, 2, 1, 1, 1, 2, 1, 1, 5, 3, 2, 1, 5, 1, 3, 1, 5, 2, 8),
 2717873: (2, 4, 2, 1, 2, 1, 3, 6, 1, 1, 3, 1, 1, 2, 3, 3, 1, 1, 2, 1, 2, 1, 5, 1, 1, 2, 2, 2, 2, 4, 3, 1, 1, 5, 4),
 9275309: (3, 2, 2, 1, 1, 1, 1, 2, 1, 1, 5, 2, 2, 2, 3, 3, 1, 5, 2, 3, 3, 1, 6, 2, 2, 1, 1, 3, 1, 1, 1, 2, 2, 1, 2,
           1, 1, 2, 1, 1, 1, 2, 3, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2, 4, 3, 1, 1, 5, 4),
 14005993: (2, 1, 1, 1, 2, 2, 2, 5, 4, 3, 1, 1, 1, 1, 2, 4, 3, 2, 3, 1, 2, 2, 1, 1, 2, 1, 1, 1, 2, 2, 2, 4, 1, 1, 1,
            1, 1, 1, 1, 2, 5, 1, 1, 3, 1, 1, 2, 1, 3, 1, 1, 7, 1, 1, 1, 4, 1, 2, 1, 1, 2, 1, 1, 1, 2, 3, 1, 1, 2, 1,
            2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2, 4, 3, 1, 1, 5, 4)}
# FROZEN_TERMINALS_END


def digest(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def rejected(call):
    try:
        call()
    except (ValueError, TypeError):
        return
    raise ValueError('hostile unexpectedly accepted')


def main():
    global CHECKS
    CHECKS = 0
    graph, dependencies, kinds = build_observations()
    demand, cycles = frontier(graph)
    initial_ranks = rooted_ranks(graph.edges)
    need(not cycles and len(demand) == 37 and len(graph.edges) == 25663,
         'declared finite observation universe')
    need(len(graph.vertices()) == 25701 and sum(n in initial_ranks for n in REQUESTED) == 16258,
         'declared initial rooted closure')
    leaves = tuple(n for n in REQUESTED if n != 1 and not dependencies[n])
    need(len(leaves) == 212, 'no-route leaves in strict-size dependency DAG')
    dag_grounded = {1}
    for n in REQUESTED:
        if any(c in dag_grounded for c in dependencies[n]):
            dag_grounded.add(n)
    need(len(dag_grounded) == 10474, 'ROOT-only strict-child fixed point')
    assumed = {1, *leaves}
    for n in REQUESTED:
        if any(c in assumed for c in dependencies[n]):
            assumed.add(n)
    need(assumed == set(REQUESTED), 'finite conditional closure from all no-route leaves')
    need(set(TERMINAL_WORDS) == set(demand), 'frozen supplied terminals are exactly missing sinks')

    # Inherited uncapped reset switch: all its edges already occur elsewhere
    # in this particular retained graph, even at the 22 newly paid requests.
    reset_graph, reset_count, reset_no_route = graph.copy(), 0, 0
    for n in REQUESTED[1:]:
        k = ((n+1) & -(n+1)).bit_length()-1
        if k < 2:
            continue
        z = 3**k*((n+1)//2**k)-1
        r = (z & -z).bit_length()-1
        if r < 2:
            continue
        child = (n-1)//2
        first = (1,)*(k-1)+(r+1,)
        second = () if child == 1 else (1,)*(k-2)+(2, r-1)
        left = reset_graph.insert(n, first, ('reset-source', n))
        right = reset_graph.insert(child, second, ('reset-child', n))
        need(left == right, 'uncapped inherited reset has exact common future')
        reset_count += 1
        reset_no_route += not dependencies[n]
    need(reset_count == 4096 and reset_no_route == 22 and reset_graph.edges == graph.edges,
         'stronger local reset dispatcher adds no literal evidence to this graph')

    # Independent explicit discovery is only a reproducibility control on these
    # 37 preselected terminals; the production compiler consumes the frozen data.
    discoveries = {n: discover_selected(n) for n in sorted(demand)}
    need(discoveries == TERMINAL_WORDS, '37 declared bounded discovery controls match frozen receipts')
    for terminal, word in TERMINAL_WORDS.items():
        states = replay(terminal, word)
        need(set(states) & set(demand) == {terminal}, 'pairwise incomparable terminal futures')

    completed, certs = complete(graph, TERMINAL_WORDS)
    need(all(n in certs for n in REQUESTED), 'all requested inputs have compiled ROOT receipts')
    need(len(completed.edges) == 26098 and len(certs) == 26099,
         'complete finite proof DAG size')
    need(completed.vertices() == set(certs), 'all retained proof vertices are grounded')
    independent_union = {}
    for n in REQUESTED:
        states = replay(n, certs[n])
        need(states[-1] == 1, 'every requested receipt independently reaches first ROOT')
        for x, y, a in zip(states, states[1:], certs[n]):
            need(x not in independent_union or independent_union[x] == (y, a), 'route determinism')
            independent_union[x] = (y, a)
    need(independent_union == completed.edges, 'final graph is exactly canonical requested-route union')
    walk_ranks = {}
    for n in graph.vertices():
        states, word, cycle = walk(graph.edges, n)
        if not cycle and states[-1] == 1:
            walk_ranks[n] = len(word)
    need(walk_ranks == initial_ranks, 'independent path walks equal reverse least closure')

    # Two independent finite minimum controls: reverse deletion, and a demand
    # greedy ordering. No further orbit discoveries occur in either control.
    deletion = []
    for removed in sorted(demand):
        remaining = {n: w for n, w in TERMINAL_WORDS.items() if n != removed}
        test = graph.copy()
        for n, w in remaining.items():
            test.insert(n, w, ('deletion', removed, n))
        remaining_ranks = rooted_ranks(test.edges)
        missing = tuple(n for n in REQUESTED if n not in remaining_ranks)
        need(missing == tuple(demand[removed]), 'deleting one receipt loses exactly its original demands')
        deletion.append((removed, len(missing)))
    greedy, trace = graph.copy(), []
    while True:
        remaining, cyclic = frontier(greedy)
        need(not cyclic, 'greedy control has no ungrounded cycle')
        if not remaining:
            break
        terminal = min(remaining, key=lambda n: (-len(remaining[n]), n))
        old = len(greedy.edges)
        greedy.insert(terminal, TERMINAL_WORDS[terminal], ('greedy-terminal', terminal))
        trace.append((terminal, len(remaining[terminal]), len(greedy.edges)-old))
    need(len(trace) == 37 and greedy.edges == completed.edges, 'greedy and simultaneous certificate imports agree')

    # All-height *descent* lifting; positive lifts are not asserted ROOT here.
    lifts = 0
    for terminal, word in TERMINAL_WORDS.items():
        p, q, b = routes.carrier(word)
        need(p*terminal+b == q and q > p, 'ROOT certificate implies contracting affine slope')
        for t in (0, 1, 2, 17):
            n, endpoint = terminal+2*q*t, 1+2*p*t
            need(replay(n, word)[-1] == endpoint < n, 'native-cylinder all-height descent control')
            lifts += 1

    for bad in (True, 1.0, 0, 2):
        rejected(lambda bad=bad: replay(bad, ()))
    rejected(lambda: complete(graph, {3: ()}))
    rejected(lambda: complete(graph, {1: (2,)}))
    rejected(lambda: complete(Graph({3: (1, 1)}), {}))
    rejected(lambda: complete(Graph({5: (True, 4)}), {}))
    rejected(lambda: complete(Graph({3: (5.0, 1)}), {}))
    probe = Graph()
    rejected(lambda: probe.insert(3, (1, 5), ('corrupt',)))
    need(not probe.edges, 'failed word insertion leaves no partial observation')
    need(rooted_ranks({3: (5, 1), 5: (3, 1)}) == {1: 0},
         'generic closure cannot self-ground an unrooted cycle')
    need(discover_selected(27, 0) is None, 'bounded discovery expiry is pending, not nonconvergence')
    need(discover_selected(1, 0) == (), 'ROOT has its unique empty word')

    representatives = {min(v): sink for sink, v in demand.items()}
    print('TERMINAL BASIS: finite ROOT1-only observation closure; no universal coverage claim')
    print('Universe: odd1..32767; frozen adaptive budget8; two inverse rows; all paid final-frontier siblings')
    print('Dependency DAG:', json.dumps(dict(statuses=kinds, no_route_leaves=len(leaves),
          no_route_leaf_hash=digest(leaves), conditional_leaf_closure=len(assumed),
          root_only_strict_child_closure=len(dag_grounded)), sort_keys=True))
    selected_max = max(p[-1]+1 for marks in graph.provenance.values() for p in marks if p[0] == 'B8-source')
    print('Initial observation graph:', json.dumps(dict(edges=len(graph.edges), vertices=len(graph.vertices()),
          submitted_edges=graph.submissions, grounded_requested=sum(n in initial_ranks for n in REQUESTED),
          unresolved_requested=sum(map(len, demand.values())), missing_sinks=len(demand),
          max_selected_prefix=selected_max, ungrounded_cycles=len(cycles),
          edges_hash=digest(sorted((n,*z) for n,z in graph.edges.items()))), sort_keys=True))
    print('Inherited uncapped reset:', reset_count, 'guarded inputs;', reset_no_route,
          'new local paid requests; exactly0 new observed edges')
    print('TERMINAL_RECEIPTS_JSON', json.dumps({str(n): list(w) for n,w in sorted(TERMINAL_WORDS.items())}, sort_keys=True))
    print('COMPONENTS_JSON', json.dumps({str(n): demand[n] for n in sorted(demand)}, sort_keys=True))
    print('Smallest requested representatives:', tuple(sorted(representatives)))
    print('Greedy (sink, requested demand, new distinct edges):', trace)
    print('37/37 receipt-deletion controls lose exactly their component; terminal futures pairwise incomparable')
    print('Final proof memory:', json.dumps(dict(requested=len(REQUESTED), certificates=len(certs),
          edges=len(completed.edges), new_edges=len(completed.edges)-len(graph.edges),
          selected_terminal_queries=len(discoveries), discovered_terminal_steps=sum(map(len, discoveries.values())),
          largest_terminal_word=max(map(len, discoveries.values())),
          separate_requested_word_edges=sum(len(certs[n]) for n in REQUESTED),
          largest_requested_word=max(len(certs[n]) for n in REQUESTED),
          all_receipts_hash=digest(sorted((n,certs[n]) for n in REQUESTED)),
          all_edges_hash=digest(sorted((n,*z) for n,z in completed.edges.items()))), sort_keys=True))
    print('All-height descent lifts:', lifts, 'finite controls of the symbolic cylinder theorem; endpoints remain explicit obligations')
    print('Checks:', CHECKS)


if __name__ == '__main__':
    main()
