"""Finite controls for the all-clock S182 +4 family's native implementation.

Prints one deterministic JSON receipt. The compiler and the independent
elementary interpreter are separate implementations. No closure search.
"""
from pathlib import Path
import hashlib
import importlib.util
import json

from followup import direct_compile, count, nine_fields, pair_fields, phase

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent.parent
OLD = RESULTS / 'zero_ready_component_returns_20261103'
ELEMENTARY = OLD / 'component/audit_native.py'
spec = importlib.util.spec_from_file_location('uniform_elementary_audit', ELEMENTARY)
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)


def strong(rows):
    for start in range(len(rows)):
        reached = 1 << start
        while True:
            more = reached
            for u, row in enumerate(rows):
                if reached >> u & 1:
                    more |= row
            if more == reached:
                break
            reached = more
        if reached != (1 << len(rows)) - 1:
            return False
    return True


def graph(prefix, suffix, ports, names):
    """Direct named M(TT_t;TT_m,W), independent of compiler operations."""
    order = prefix + suffix
    position = {v: j for j, v in enumerate(order)}
    def arrow(u, v):
        if u == v:
            return False
        if v == 39:
            return u in ports
        if u == 39:
            return v not in ports
        return position[u] < position[v]
    rows = [sum(1 << j for j, v in enumerate(names) if arrow(u, v)) for u in names]
    A.valid(rows)
    return rows


def subset_h(rows):
    """Independent endpoint subset DP, not using a feedback representation."""
    n = len(rows)
    dp = [{} for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v, number in dp[mask].items():
            possible = rows[v] & (full ^ mask)
            while possible:
                bit = possible & -possible
                possible -= bit
                w = bit.bit_length() - 1
                child = dp[mask | bit]
                child[w] = child.get(w, 0) + number
    return sum(dp[full].values())


def state_record(names, rows, expected_H, expected_defect, include_pair):
    rows = tuple(rows)
    r = names.index(39)
    feedback = (r,)
    A.need(strong(rows), 'all completed adjacency stages are strong')
    nine = nine_fields(rows, feedback, r)
    h, s, q, F, e, P2, S2, E3, Q2 = nine
    eight = [h, e, S2, F, Q2, q, P2,
             count(rows, feedback, 3, starts=((0, r),))]
    A.need(h == expected_H == subset_h(rows), 'formula/gap/subset full H agreement')
    A.need(F - h == expected_defect, 'actual transported root defect')
    A.need((S2 - e) % 2 == 0, 'integral current rooted golden')
    golden = [h, (S2 - e) // 2, (S2 + e) // 2]
    result = dict(names=names[:], rows=list(rows), H=h, root=39,
                  root_inports=[v for i, v in enumerate(names) if rows[i] >> r & 1],
                  root_fields=[e, s, q, F, S2], root_eight=eight,
                  root_defect=F-h, golden=golden, charts=phase(golden),
                  strong=True, independent_subset_H=h)
    if include_pair:
        f = pair_fields(rows, feedback, r, names.index(1))
        A.need(f[0] == h and f[2:7] == result['root_fields'], 'named root-pair binding')
        result['root_old_port_pair'] = dict(ports=[39, 1], H=f[0], P2=f[1],
                                            port_data=[f[2:7], f[7:12]], ordered_W2=f[12:])
    return result


def one_clock(t):
    A.need(t >= 1, 'nonempty transitive prefix')
    prefix = list(range(100, 100+t))
    names = [0, 1, 2] + prefix + [39]
    original_names = names[:]
    rows = graph(prefix, [0, 1, 2], {0, 1, 2}, names)
    source_rows = rows[:]
    stages = [state_record(names, rows, 7 * (1 << t) - 5, 6, True)]
    births = []
    for newborn, beaten, suffix, ports, coefficient, intercept, defect in (
        (4, [39], [0, 1, 2, 4], {0, 1, 2, 4}, 15, -13, 14),
        (3, [4], [0, 1, 2, 3, 4], {0, 1, 2, 4}, 23, -19, 20),
    ):
        cut = sum(1 << names.index(v) for v in beaten)
        before = rows[:]
        rows = A.add_vertex(rows, cut)
        names.append(newborn)
        A.need([row & ((1 << len(before)) - 1) for row in rows[:-1]] == before,
               'permanent birth preserves every old arc')
        A.need(rows == graph(prefix, suffix, ports, names), 'literal named birth attachment')
        births.append(dict(name=newborn, out_mask=cut, beaten_old_names=beaten,
                           names=names[:], rows=rows[:]))
        stages.append(state_record(names, rows, coefficient * (1 << t) + intercept,
                                   defect, False))

    before = rows[:]
    U, V, W, p, q, case, inverse, pair_word = direct_compile(tuple(rows), names.index(1), names.index(39))
    helper = -1
    local = names + [helper]
    words = [[local[v] for v in word] for word in (p, q)]
    A.need(words == [[1, 39] + prefix + [0, 2, 3, 4, helper],
                     [1, helper, 39] + prefix + [0, 2, 3, 4]],
           'literal all-clock guard words')
    A.need(U[-1] == 1 << names.index(1), 'helper beats exactly old port1')
    state = A.add_vertex(rows, U[-1])
    A.need(state == list(U), 'independent helper birth equals compiler input')
    native_rows = []
    for word, expected in zip(words, (V, W)):
        state = A.native(state, local, word)
        A.need(state == list(expected), 'independent literal native guard replay')
        native_rows.append(state[:])
    rows = [row & ((1 << len(names)) - 1) for row in state[:-1]]
    A.need(A.differences(before, rows, names) == [[1, 39]], 'only the specified root arc changes')
    A.need(rows == graph(prefix, [0, 1, 2, 3, 4], {0, 2, 4}, names), 'exact inherited target family')
    A.need(names[:len(original_names)] == original_names and set(names) - set(original_names) == {3, 4},
           'all old physical names retained and exactly two permanent births')
    induced_old = [row & ((1 << len(original_names)) - 1) for row in rows[:len(original_names)]]
    A.need(A.differences(source_rows, induced_old, original_names) == [[1, 39]],
           'all source incidences except the declared root arc survive at the endpoint')
    target = state_record(names, rows, 21 * (1 << t) - 11, 12, True)
    A.need(target['H'] == 3 * stages[0]['H'] + 4, 'actual +4 return')
    stages.append(target)
    return dict(clock=t, prefix=prefix, stage_labels=['source', 'birth_terminal_port4',
                'birth_nonport3_before4', 'reverse_old_port1_root39'], stages=stages,
                permanent_births=births,
                native_program=dict(changed_pair=[1, 39], helper=helper,
                    helper_out_mask=U[-1], pair_guard_names=[names[v] for v in pair_word],
                    paths=words, native_rows=native_rows, case=case, inverse=inverse,
                    after_delete_rows=rows, temporary_births=1, helper_deletions=1,
                    native_complements=2, peak_order=len(names)+1),
                old_names=original_names, deleted_old_names=[])


def main():
    control = [6, 4, 0]
    A.need(A.native(control, [0, 1, 2], [0, 1, 2]) == [2, 4, 1], 'positive native control')
    A.need(A.native([2, 4, 1], [0, 1, 2], [0, 1, 2]) == control, 'native involution control')
    for bad in ([2, 1, 0], [0, 1, 1]):
        rejected = False
        try:
            A.native(control, [0, 1, 2], bad)
        except RuntimeError:
            rejected = True
        A.need(rejected, 'illegal directed and repeated-name guards rejected')
    records = [one_clock(t) for t in (1, 2, 3, 4, 8)]
    pins = [Path(__file__), HERE/'UNIFORM.md', HERE/'followup.py', HERE/'FOLLOWUP.md',
            ELEMENTARY, OLD/'interface/gap.py', OLD/'interface/check.py',
            OLD/'inverse/native.py', OLD/'continuation/native_tree.py',
            RESULTS/'adaptive_metagraph_reformulation_20261103/metagraph_probe.py',
            RESULTS/'adaptive_metagraph_reformulation_20261103/METAGRAPH.md',
            RESULTS/'session182_physical_returns_20261104/induction/PROOF.md',
            RESULTS/'session182_physical_returns_20261104/induction_audit/AUDIT.md']
    print(json.dumps(dict(status='PASS',
        scope='PROVED uniform native implementation of the inherited shared-clock +4 family; FINITE-EXACT controls at clocks1,2,3,4,8. No returning family, +2 map, or unrestricted spectrum coverage.',
        records=records, root_eight_fields=['H', 'E1', 'S2', 'E2', 'deletion_P2', 'deletion_H', 'P2', 'S3'],
        totals=dict(clocks=5, completed_strong_stages=20, independent_subset_counts=20,
                    named_endpoint_pairs=10, native_complements=10, temporary_births_and_deletions=5,
                    permanent_births=10, positive_hostile_controls=4),
        source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins}),
        indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
