"""Monotone reuse of incomplete exact Collatz observations.

Only ROOT is initially certified. Missing outgoing edges are queried explicitly,
with a declared cap. Every inserted edge receives an independent while-even
check; graph closure never invokes an orbit oracle or certifies an unrooted cycle.
"""
from collections import Counter, defaultdict, deque
from dataclasses import dataclass, field
from hashlib import sha256
import json

import adaptive_boundary_selector_20261004 as selector
import inverse_ray_ternary_addresses_20261004 as codec


def need(ok, message):
    if not ok:
        raise ValueError(message)


def literal_step(n):
    need(type(n) is int and n > 0 and n % 2, 'positive odd exact source')
    value, exponent = 3*n+1, 0
    while value % 2 == 0:
        value //= 2
        exponent += 1
    return value, exponent


@dataclass
class ObservationGraph:
    edges: dict = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)
    submitted_edges: int = 0

    def add_edge(self, source, target, exponent, origin):
        need(type(source) is int and source > 1 and source % 2,
             'nonroot positive odd observation source')
        need(type(target) is int and target > 0 and target % 2,
             'positive odd observation target')
        need(type(exponent) is int and exponent > 0, 'positive exact valuation')
        need(literal_step(source) == (target, exponent), 'independent literal edge validation')
        need(source not in self.edges or self.edges[source] == (target, exponent),
             'deterministic odd map has one observed outgoing edge')
        new = source not in self.edges
        self.edges[source] = (target, exponent)
        self.provenance.setdefault(source, set()).add(origin)
        self.submitted_edges += 1
        return new

    def add_word(self, source, word, origin, root_tail=False):
        need(len(word) <= 128, 'declared observed-word cap128')
        current, new = source, 0
        for position, exponent in enumerate(word):
            target, actual = literal_step(current)
            need(actual == exponent, 'supplied word uses actual valuations')
            if current == 1:
                need(root_tail and exponent == 2, 'only explicit reset proof may have formal root tail')
            else:
                new += self.add_edge(current, target, exponent, (*origin, position))
            current = target
        return current, new

    def copy(self):
        return ObservationGraph(dict(self.edges), {n: set(v) for n, v in self.provenance.items()},
                                self.submitted_edges)

    def vertices(self):
        return {1, *self.edges, *(edge[0] for edge in self.edges.values())}


def least_closure(edges):
    """Generic finite rooted reachability; no validation or new observations."""
    reverse = defaultdict(list)
    for source, (target, _) in sorted(edges.items()):
        reverse[target].append(source)
    ranks, queue = {1: 0}, deque([1])
    while queue:
        target = queue.popleft()
        for source in reverse[target]:
            if source not in ranks:
                ranks[source] = ranks[target]+1
                queue.append(source)
    return ranks


def route_walk(edges, source):
    """Read the finite observed prefix until root, a missing edge, or a cycle."""
    path, word, positions = [source], [], {}
    current = source
    while current != 1 and current in edges and current not in positions:
        positions[current] = len(word)
        current, exponent = edges[current]
        path.append(current)
        word.append(exponent)
    return tuple(path), tuple(word), current in positions


def missing_frontier(graph, requested, ranks=None):
    ranks = least_closure(graph.edges) if ranks is None else ranks
    demand, cyclic = defaultdict(list), []
    for source in requested:
        if source in ranks:
            continue
        path, _, cycle = route_walk(graph.edges, source)
        if cycle:
            cyclic.append(source)
        else:
            need(path[-1] not in graph.edges and path[-1] != 1, 'genuinely missing outgoing edge')
            demand[path[-1]].append(source)
    return dict(demand), tuple(cyclic)


def add_selection_result(graph, result, policy):
    source = result.source
    endpoint, new = graph.add_word(source, result.word, (policy, source, 'actual'))
    need(endpoint == result.endpoint, 'stored selected endpoint')
    if result.kind == 'common_future':
        for child in result.dependencies:
            if child == result.endpoint:
                continue
            target, exponent = literal_step(child)
            need(target == result.endpoint, 'advertised sibling meets the actual endpoint')
            if child != 1:
                new += graph.add_edge(child, target, exponent, (policy, source, 'sibling'))
            else:
                need(target == 1, 'root child needs no self-edge')
    return new


def add_selected(graph, source, policy):
    return add_selection_result(graph, selector.select(source, 8, policy, word_cap=128), policy)


def recorded_seed_trace():
    """Observe the frozen original one-pass algorithm without changing choices."""
    original_select = selector.select
    graph, calls = ObservationGraph(), []

    def recorded(*args, **kwargs):
        result = original_select(*args, **kwargs)
        calls.append(result.source)
        add_selection_result(graph, result, 'original-seed-trace')
        return result

    selector.select = recorded
    try:
        certified, pending, unresolved = selector.seed_closure(1023, 8, 'adaptive', suffix_memory=True)
    finally:
        selector.select = original_select
    return graph, tuple(calls), certified, pending, unresolved


def add_reset(graph, source):
    if source == 1:
        return 0, False
    run = selector.v2(source+1)-1
    if run < 1 or run+1 > 128:
        return 0, False
    checkpoint = (3**run*(source+1)//2**run)-1
    endpoint, exponent = literal_step(checkpoint)
    if exponent < 3:
        return 0, False
    first = (1,)*run+(exponent,)
    second = (1,)*(run-1)+(2, exponent-2)
    child = (source-1)//2
    left, count1 = graph.add_word(source, first, ('reset', source, 'actual'), root_tail=True)
    right, count2 = graph.add_word(child, second, ('reset', source, 'child'), root_tail=True)
    need(0 < child < source and child % 2 and left == right == endpoint, 'reset common-future witness')
    return count1+count2, True


def graph_certificates(graph, ranks):
    """Build ASTs after least closure, using only the retained certified edges."""
    certs = {1: codec.ROOT}
    for source in sorted(ranks, key=lambda n: (ranks[n], n)):
        if source == 1:
            continue
        target, exponent = graph.edges[source]
        need(target in certs and ranks[target]+1 == ranks[source], 'proved smaller graph rank')
        least = codec.kappa(target % 9, source % 3)
        need(exponent >= least and (exponent-least) % 6 == 0, 'canonical inverse address')
        certs[source] = codec.extend(certs[target], source % 3, (exponent-least)//6)
        need(codec.expand(certs[source]) == source, 'AST retains exact source label')
    return certs


def describe(name, graph, requested):
    ranks = least_closure(graph.edges)
    frontier, cyclic = missing_frontier(graph, requested, ranks)
    missing = [n for n in requested if n not in ranks]
    print(name, json.dumps(dict(edges=len(graph.edges), vertices=len(graph.vertices()),
                               submitted_edges=graph.submitted_edges,
                               all_certified_vertices=len(ranks), requested_certified=len(requested)-len(missing),
                               missing=missing, frontier={str(n): v for n, v in sorted(frontier.items())},
                               ungrounded_cycles=list(cyclic)), sort_keys=True))
    return ranks


def frontier_expansion(graph, requested, cap=128):
    need(type(cap) is int and cap >= 0, 'nonnegative exact new-edge query cap')
    graph = graph.copy()
    trace = []
    for _ in range(cap):
        before = least_closure(graph.edges)
        demand, cyclic = missing_frontier(graph, requested, before)
        if not demand:
            break
        chosen = min(demand, key=lambda n: (-len(demand[n]), n))
        target, exponent = selector.step(chosen)
        need(graph.add_edge(chosen, target, exponent, ('frontier-query', len(trace))),
             'a query must add a previously unknown edge')
        after = least_closure(graph.edges)
        need(set(before) <= set(after), 'closure is monotone under observations')
        trace.append(dict(source=chosen, target=target, exponent=exponent,
                          demand=len(demand[chosen]), requested=demand[chosen],
                          newly_certified=[n for n in requested if n in after and n not in before]))
    return graph, tuple(trace)


def separate_completion_queries(initial, sources, cap=128):
    """Each input gets the same initial cache, with no sharing of new edges."""
    grounded = least_closure(initial.edges)
    counts, all_new = {}, set()
    for source in sources:
        graph, current, queries, seen = initial.copy(), source, 0, set()
        while current not in grounded:
            need(current not in seen, 'declared separate controls have no ungrounded cycle')
            seen.add(current)
            if current not in graph.edges:
                need(queries < cap, 'separate-query comparison cap')
                target, exponent = selector.step(current)
                graph.add_edge(current, target, exponent, ('separate-query', source, queries))
                all_new.add(current)
                queries += 1
            current = graph.edges[current][0]
        counts[source] = queries
    return counts, all_new


def main():
    requested = tuple(range(1, 1024, 2))
    print('ADAPTIVE OBSERVATION UNION: only ROOT1 is initially certified')
    print('Universe:512 odd sources1..1023; selector budget8; per observed word cap128')
    graph = ObservationGraph()
    for source in requested:
        add_selected(graph, source, 'adaptive')
    initial = graph.copy()
    ranks = describe('adaptive observations + advertised sibling edges', graph, requested)
    need(len(graph.edges) == 802 and sum(n in ranks for n in requested) == 509,
         'declared base census')
    walk_ranks = {}
    for source in initial.vertices():
        path, word, cyclic = route_walk(initial.edges, source)
        if not cyclic and path[-1] == 1:
            walk_ranks[source] = len(word)
    need(walk_ranks == ranks, 'independent per-vertex walk and reverse fixed-point agreement')
    snapshot = dict(requested=[1, 1023, 2], macro_budget=8, word_cap=128,
                    edges=[[source, target, exponent]
                           for source, (target, exponent) in sorted(initial.edges.items())],
                    frontier={str(n): v for n, v in sorted(missing_frontier(initial, requested)[0].items())})
    encoded_snapshot = json.dumps(snapshot, sort_keys=True, separators=(',', ':'))
    print('Initial graph snapshot SHA256:', sha256(encoded_snapshot.encode()).hexdigest())
    print('INITIAL_GRAPH_JSON', encoded_snapshot)
    recorded, calls, old_certified, old_pending, old_unresolved = recorded_seed_trace()
    need(recorded.edges == initial.edges, 'same checked transition set as frozen one-pass trace')
    need(len(old_certified) == 354 and len(calls) == 353 and recorded.submitted_edges == 1388,
         'declared identical-observation control')
    print('Same-observation original trace:', json.dumps(dict(selector_calls=len(calls),
          nonroot_selector_calls=sum(n != 1 for n in calls), submitted_edges=recorded.submitted_edges,
          distinct_edges=len(recorded.edges), original_certified=len(old_certified),
          original_pending=len(old_pending), original_unavailable_child=len(old_unresolved),
          union_certified=sum(n in least_closure(recorded.edges) for n in requested)), sort_keys=True))
    for policy in ('inherited', 'literal'):
        before = len(graph.edges)
        for source in requested:
            add_selected(graph, source, policy)
        describe('union '+policy+'; new distinct edges '+str(len(graph.edges)-before), graph, requested)
    reset_count, reset_edges = 0, 0
    for source in requested:
        added, accepted = add_reset(graph, source)
        reset_count += accepted
        reset_edges += added
    describe('union reset witnesses; accepted '+str(reset_count)+'; new edges '+str(reset_edges), graph, requested)
    learned = selector.discover_descent(27, cap=128)
    need(learned is not None, 'declared inherited learned-cylinder control')
    learned_additions = 0
    for source in requested:
        if source % learned.modulus == learned.residue:
            _, new = graph.add_word(source, learned.word, ('learned27', source))
            learned_additions += new
    describe('union learned27; new edges '+str(learned_additions), graph, requested)
    need(graph.edges == initial.edges, 'all optional witnesses are redundant in this finite graph')

    expanded, trace = frontier_expansion(initial, requested, cap=128)
    final_ranks = describe('after demand-scheduled missing-edge queries', expanded, requested)
    print('Explicit new queries:', len(trace))
    for record in trace:
        print('  query', json.dumps(record, sort_keys=True))
    required_edges = set()
    for source in requested:
        path, _, cyclic = route_walk(expanded.edges, source)
        need(not cyclic and path[-1] == 1, 'declared requested route is complete')
        required_edges.update(path[:-1])
    need(required_edges == set(expanded.edges), 'final graph equals requested canonical route union')
    need(set(initial.edges) <= required_edges and len(required_edges-set(initial.edges)) == len(trace) == 32,
         'exact minimum32 missing edges for literal-edge completion of this graph')
    missing_sources = tuple(n for n in requested if n not in ranks)
    separate_counts, separate_distinct = separate_completion_queries(initial, missing_sources)
    need(separate_distinct == required_edges-set(initial.edges), 'same required missing edge set')
    print('Separate missing-input completion, same initial cache each time:',
          json.dumps(separate_counts, sort_keys=True), '; total new-edge queries', sum(separate_counts.values()),
          '; distinct required', len(separate_distinct))
    certificates = graph_certificates(expanded, final_ranks)
    for source in requested:
        if source in certificates:
            need(certificates[source] == codec.encode_source(source),
                 'post-closure independent first-hit encoder audit')
    print('Post-closure audited requested ASTs:', sum(n in certificates for n in requested),
          '; all graph certificates:', len(certificates),
          '; separate requested odd-edge count:', sum(final_ranks.get(n, 0) for n in requested))

    # Generic reachability control only: this formal cycle is NOT Collatz data.
    need(least_closure({7: (11, 1), 11: (7, 1)}) == {1: 0},
         'ungrounded formal cycles cannot certify themselves')
    hostile = ObservationGraph()
    for args in ((1, 1, 2, ('bad',)), (3, 7, 1, ('bad',)), (True, 1, 2, ('bad',))):
        try:
            hostile.add_edge(*args)
        except ValueError:
            pass
        else:
            raise ValueError('invalid or root edge accepted')
    capped, empty = frontier_expansion(initial, requested, cap=0)
    need(not empty and capped.edges == initial.edges, 'zero cap adds no observation')
    reset_hostile = ObservationGraph()
    need(add_reset(reset_hostile, 7) == (0, False) and not reset_hostile.edges,
         'reset exponent2 does not create the forbidden exponent0 edge')
    # Insertion order affects neither the fixed point nor the exact edge labels.
    reverse = ObservationGraph()
    for source, (target, exponent) in reversed(list(expanded.edges.items())):
        reverse.add_edge(source, target, exponent, ('reordered', source))
    need(least_closure(reverse.edges) == final_ranks, 'order-independent least closure')
    print('Hostiles passed:ungrounded formal cycle, root self-edge, false edge, bool source, zero cap, reset2, reordered insertion')
    print('Every extra edge is an explicitly counted observation; no new oracle seeds or universal coverage claim.')


if __name__ == '__main__':
    main()
