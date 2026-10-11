"""Last-consumer scheduling for an authenticated Collatz receipt DAG.

Only supplied finite words are replayed. No orbit discovery, no new source
coverage, and no claim about Python byte memory or multiplication runtime.
The resource model counts live immutable word payloads, excludes the retained
recipe and output sink, forbids recomputation and in-place result construction.
"""
from collections import Counter
from dataclasses import dataclass, replace
from hashlib import sha256
from itertools import permutations

import paper_local_global_grounding_20261007 as proof
from collatz_receipt_flow_geometry_20261009 import ROOT27

CHECKS = 0


def require(ok, why):
    if not ok:
        raise ValueError(why)


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


@dataclass(frozen=True)
class Node:
    key: str
    parent: str | None
    receipt: proof.Receipt | None


@dataclass(frozen=True)
class Event:
    kind: str
    key: str
    slot: int


def source(node):
    return 1 if node.key == 'ROOT' else node.receipt.left


def audit_program(program, order, requested):
    require(type(program) is tuple and type(order) is tuple and type(requested) is tuple,
            'exact tuple program/order/requested')
    nodes = {}
    for node in program:
        require(type(node) is Node and type(node.key) is str and node.key,
                'exact Node and nonempty string key')
        require(node.key not in nodes, 'duplicate key')
        if node.key == 'ROOT':
            require(node.parent is None and node.receipt is None, 'unique empty ROOT recipe')
        else:
            require(type(node.parent) is str and node.parent and node.parent != node.key,
                    'distinct named parent')
            proof.verify_receipt(node.receipt)
            require(node.receipt.left != 1, 'only ROOT recipe may claim source one')
        nodes[node.key] = node
    require('ROOT' in nodes, 'supplied ROOT1 required')
    require(all(type(k) is str for k in order + requested), 'exact string names')
    require(len(set(order)) == len(order) and set(order) == set(nodes), 'complete unique order')
    require(len(set(requested)) == len(requested) and set(requested) <= set(nodes),
            'distinct listed outputs')
    positions = {k: i for i, k in enumerate(order)}
    for node in program:
        if node.key != 'ROOT':
            require(node.parent in nodes and positions[node.parent] < positions[node.key],
                    'every dependency earlier: circular assumptions rejected')
            require(node.receipt.right == source(nodes[node.parent]), 'immutable source binding')
    return nodes


def use_counts(nodes, requested):
    counts = Counter({k: 0 for k in nodes})
    for node in nodes.values():
        if node.parent is not None:
            counts[node.parent] += 1
    counts.update(requested)
    return counts


def plan(program, order, requested):
    """Greedy slot reuse for this order; never reorders a valuation word."""
    nodes = audit_program(program, order, requested)
    remaining = use_counts(nodes, requested)
    occupied, events = {}, []
    for key in order:
        slot = 0
        while slot in occupied.values():
            slot += 1
        occupied[key] = slot
        events.append(Event('build', key, slot))
        parent = nodes[key].parent
        if parent is not None:
            remaining[parent] -= 1
        if key in requested:
            events.append(Event('emit', key, slot))
            remaining[key] -= 1
        for dead in tuple(occupied):
            if remaining[dead] == 0:
                events.append(Event('drop', dead, occupied.pop(dead)))
    return tuple(events)


def splice(node, parent_word):
    """Transport a supplied parent's first-hit word across an actual seam."""
    rec = node.receipt
    require(proof.replay(rec.right, parent_word) == 1, 'parent is not a ROOT word')
    k = len(rec.right_word)
    require(parent_word[:k] == rec.right_word, 'right receipt must be the actual parent prefix')
    result = rec.left_word + parent_word[k:]
    require(proof.replay(rec.left, result) == 1, 'source first-hit export failed')
    return result


def word_digest(word):
    return sha256(','.join(map(str, word)).encode('ascii')).hexdigest()


def execute(program, order, requested, events, sink=None):
    """Reauthenticate retained recipes and independently check the event plan.

    sink receives (key, source, first_hit_word) at emission; externally retained
    outputs are outside the workspace metric. Metadata alone is not a proof.
    """
    nodes = audit_program(program, order, requested)
    require(type(events) is tuple and (sink is None or callable(sink)), 'tuple events/callable sink')
    remaining = use_counts(nodes, requested)
    slots, located, built, emitted = {}, {}, [], set()
    peak_slots = peak_cells = allocations = 0
    summaries = []
    for event in events:
        require(type(event) is Event and type(event.kind) is str and type(event.key) is str,
                'exact typed Event')
        require(type(event.slot) is int and event.slot >= 0, 'exact nonnegative physical slot')
        require(event.key in nodes and event.kind in ('build', 'emit', 'drop'), 'known event')
        key, slot = event.key, event.slot
        if event.kind == 'build':
            require(len(built) < len(order) and key == order[len(built)], 'declared order exactly once')
            require(slot not in slots, 'cannot overwrite a live payload')
            node = nodes[key]
            if node.parent is None:
                word = ()
            else:
                require(node.parent in located, 'consumer requires a live parent payload')
                parent_key, parent_word = slots[located[node.parent]]
                require(parent_key == node.parent, 'physical label/source must travel together')
                word = splice(node, parent_word)
            # Allocate before consuming/releasing the parent: non-in-place model.
            slots[slot] = (key, word)
            located[key] = slot
            built.append(key)
            allocations = max(allocations, slot + 1)
            peak_slots = max(peak_slots, len(slots))
            peak_cells = max(peak_cells, sum(1 + len(w) for _, w in slots.values()))
            if node.parent is not None:
                remaining[node.parent] -= 1
        elif event.kind == 'emit':
            require(key in requested and key not in emitted, 'requested output exactly once')
            require(located.get(key) == slot and slot in slots, 'emit correct live payload')
            payload_key, word = slots[slot]
            require(payload_key == key and proof.replay(source(nodes[key]), word) == 1,
                    'independent output authentication')
            if sink is not None:
                sink(key, source(nodes[key]), word)
            summaries.append((key, source(nodes[key]), len(word), sum(word), word_digest(word)))
            emitted.add(key)
            remaining[key] -= 1
        else:
            require(located.get(key) == slot and slot in slots and slots[slot][0] == key,
                    'release correct live payload')
            require(remaining[key] == 0, 'release only after last declared consumer/output')
            del slots[slot]
            del located[key]
    require(tuple(built) == order and emitted == set(requested), 'all obligations discharged')
    require(not slots and not any(remaining.values()), 'complete workspace release')
    return {'peak_slots': peak_slots, 'peak_word_cells': peak_cells,
            'allocated_slots': allocations, 'outputs': tuple(summaries)}


def interval_bound(program, order, requested):
    """Independent fixed-order lower bound from closed live intervals."""
    nodes = audit_program(program, order, requested)
    births = {key: i for i, key in enumerate(order)}
    deaths = dict(births)
    for node in program:
        if node.parent is not None:
            deaths[node.parent] = max(deaths[node.parent], births[node.key])
    widths = tuple(sum(births[k] <= i <= deaths[k] for k in nodes) for i in range(len(order)))
    return max(widths), tuple((k, births[k], deaths[k]) for k in order)


def retain_all_program(program, order, requested):
    """Control schedule: emit immediately, retain every constructed payload."""
    audit_program(program, order, requested)
    events = []
    for i, key in enumerate(order):
        events.append(Event('build', key, i))
        if key in requested:
            events.append(Event('emit', key, i))
    events.extend(Event('drop', key, i) for i, key in enumerate(order))
    return tuple(events)


def example():
    R = proof.Receipt
    return (Node('ROOT', None, None),
            Node('five', 'ROOT', R(5, 1, (4,), ())),
            Node('three', 'five', R(3, 5, (1,), ())),
            Node('seven', 'three', R(7, 3, (1, 1, 2, 3), (1,))),
            Node('nine', 'seven', R(9, 7, (2,), ())),
            Node('eleven', 'seven', R(11, 7, (1, 2, 3), (1, 1, 2, 3))),
            Node('twenty_seven', 'ROOT', R(27, 1, ROOT27, ())),
            Node('forty_one', 'twenty_seven', R(41, 27, (2,), (1, 2))),
            Node('thirty_one', 'twenty_seven', R(31, 27, (1,), (1, 2, 1))))


def rejects(fn, why):
    try:
        fn()
    except ValueError:
        check(True, why)
        return
    check(False, why)


def carrier(word):
    p, q, b = 1, 1, 0
    for a in word:
        p, q, b = 3*p, q*2**a, 3*b+q
    return p, q, b


def main():
    program = example()
    order = tuple(node.key for node in program)
    requested = order
    events = plan(program, order, requested)
    emitted = {}
    def sink(key, n, word):
        check(proof.replay(n, word) == 1, 'sink strict replay')
        emitted[key] = word
    report = execute(program, order, requested, events, sink)
    naive = execute(program, order, requested, retain_all_program(program, order, requested))
    check(report['outputs'] == naive['outputs'], 'same exact outputs')
    check(report['peak_slots'] == interval_bound(program, order, requested)[0], 'interval optimum')
    _, live_intervals = interval_bound(program, order, requested)
    weighted_bound = max(sum(1 + len(emitted[key]) for key, birth, death in live_intervals
                             if birth <= i <= death) for i in range(len(order)))
    check(report['peak_word_cells'] == weighted_bound, 'weighted interval optimum')
    check(naive['peak_word_cells'] == sum(1 + len(w) for w in emitted.values()), 'retained payload bill')
    check(report['allocated_slots'] == report['peak_slots'], 'dense greedy physical slots')
    check(report['peak_slots'] < naive['peak_slots'], 'declared model workspace saving')
    check(emitted['seven'] == (1, 1, 2, 3, 4), 'common-future transport')
    check(emitted['thirty_one'] == ROOT27[2:], 'supplied suffix retained')
    # Verify every topological order of the seven-node sub-DAG, not every DAG.
    small = tuple(n for n in program if n.key not in ('nine', 'eleven'))
    first = small[0].key
    valid_orders, peaks = 0, Counter()
    for tail in permutations(tuple(n.key for n in small[1:])):
        candidate = (first,) + tail
        positions = {k: i for i, k in enumerate(candidate)}
        if any(n.parent is not None and positions[n.parent] >= positions[n.key] for n in small):
            continue
        got = execute(small, candidate, candidate, plan(small, candidate, candidate))
        lower, intervals = interval_bound(small, candidate, candidate)
        check(got['peak_slots'] == got['allocated_slots'] == lower, 'all-order interval bound')
        check({(n, rank, cost, digest) for _, n, rank, cost, digest in got['outputs']}
              == {(n, rank, cost, digest) for key, n, rank, cost, digest in report['outputs']
                  if key not in ('nine', 'eleven')}, 'schedule does not alter proofs')
        valid_orders += 1
        peaks[lower] += 1
    # Finite all-height one-step ROOT schema, not fresh orbit seeds.
    rays = (Node('ROOT', None, None),) + tuple(
        Node('ray'+str(k), 'ROOT', proof.Receipt((4**k-1)//3, 1, (2*k,), ()))
        for k in range(2, 66))
    ray_order = tuple(n.key for n in rays)
    ray_report = execute(rays, ray_order, ray_order, plan(rays, ray_order, ray_order))
    check(ray_report['peak_slots'] == 2, 'shared ROOT payload remains live through final consumer')
    check(len(ray_report['outputs']) == 65, 'explicit ray test universe')
    check(carrier((1, 2)) == (9, 8, 5) and carrier((2, 1)) == (9, 8, 7),
          'same local profile does not retain ordered semantics')
    check(proof.replay(11, (1, 2)) == 13, 'positive chronological word')
    rejects(lambda: proof.replay(11, (2, 1)), 'invalid chronological reordering')
    rejects(lambda: execute(program, order, requested,
                            (events[0], Event('drop', 'ROOT', 0)) + events[1:]),
            'premature release with another consumer')
    second_build = next(i for i, e in enumerate(events) if e.kind == 'build' and e.key == 'five')
    bad = events[:second_build] + (replace(events[second_build], slot=0),) + events[second_build+1:]
    rejects(lambda: execute(program, order, requested, bad), 'live slot overwritten')
    rejects(lambda: execute(program, order, requested, (replace(events[0], slot=True),)+events[1:]),
            'bool slot alias')
    rejects(lambda: execute(program, order, requested, (replace(events[0], slot=0.0),)+events[1:]),
            'float slot alias')
    rejects(lambda: execute(program, order, requested, tuple(e for e in events if e.kind != 'emit')),
            'dropped outputs remain obligations')
    rejects(lambda: execute(program, order, requested, events[:-1]), 'retained final live payload')
    changed = (program[0], replace(program[1], parent='three')) + program[2:]
    rejects(lambda: plan(changed, order, requested), 'cycle cannot authenticate itself')
    changed = program[:3] + (replace(program[3], parent='five'),) + program[4:]
    rejects(lambda: plan(changed, order, requested), 'parent source mismatch')
    for malformed in (True, 5.0):
        changed = (program[0], replace(program[1], receipt=proof.Receipt(malformed, 1, (4,), ())))
        rejects(lambda p=changed: plan(p, ('ROOT', 'five'), ('five',)), 'exact source types')
    padded = (program[0], Node('bad', 'ROOT', proof.Receipt(5, 1, (4, 2), (2,))))
    rejects(lambda: plan(padded, ('ROOT', 'bad'), ('bad',)), 'ROOT padding rejected')
    wrong = (program[0], Node('bad', 'ROOT', proof.Receipt(5, 1, (3,), ())))
    rejects(lambda: plan(wrong, ('ROOT', 'bad'), ('bad',)), 'false valuation rejected')
    print('Multiplication-inspired retained-recipe / last-consumer certificate routing')
    print('Pinned external main: 12a4934d15006e8c10f9b2264f13f8f341be802b (2026-10-10T14:13:23Z)')
    print('9-node supplied receipt DAG: retain-all slots', naive['peak_slots'],
          '; reused slots', report['peak_slots'], '; word cells',
          naive['peak_word_cells'], '->', report['peak_word_cells'])
    print('Same exact requested output metadata:', report['outputs'] == naive['outputs'])
    print('7-node all-topological-order controls:', valid_orders, '; peak histogram', sorted(peaks.items()))
    print('Root-ray schema controls k=2..65: 65 nodes including ROOT; peak slots', ray_report['peak_slots'])
    print('All chronological words and full immutable recipes retained; no new source coverage.')
    print('Cost model: payload slots/cells only; excludes recipe, verifier transients and emitted output storage.')
    print('Checks:', CHECKS)


if __name__ == '__main__':
    main()
