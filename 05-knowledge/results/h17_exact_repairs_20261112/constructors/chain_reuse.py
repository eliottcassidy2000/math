"""Universal coefficient certificate and literal retained-chain renewals."""
from collections import defaultdict
from math import gcd
from pathlib import Path
import json
from wrapper_repair import birth, require, state, strong, wrapper

HERE = Path(__file__).resolve().parent
ZERO = (0, 0, 0, 0, 0)


def add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def kernel(carrier, root, newborn, rename):
    """Exact unordered endpoint-cover birth kernel, with forgotten mark -1."""
    result = defaultdict(lambda: ZERO)
    for signature, weight in carrier.items():
        incoming = [None] + [i for i, (_, end) in enumerate(signature) if end == root]
        outgoing = [None] + [i for i, (start, _) in enumerate(signature) if start != root]
        for left in incoming:
            for right in outgoing:
                if left is not None and left == right:
                    continue
                start = newborn if left is None else signature[left][0]
                end = newborn if right is None else signature[right][1]
                components = [p for i, p in enumerate(signature) if i not in (left, right)]
                components.append((start, end))
                renamed = tuple(sorted((rename.get(a, a), rename.get(b, b))
                                       for a, b in components))
                if sum(b == -1 for _, b in renamed) > 1:
                    continue
                result[renamed] = add(result[renamed], weight)
    return dict(result)


def observation(carrier, root):
    answer = [ZERO] * 4
    for sig, weight in carrier.items():
        flags = [len(sig) == 1,
                 len(sig) == 1 and sig[0][0] == root,
                 len(sig) == 2 and (root, root) in sig,
                 len(sig) <= 2 and any(b == root for _, b in sig)]
        for i, flag in enumerate(flags):
            if flag:
                answer[i] = add(answer[i], weight)
    return answer


def universal_certificate():
    # Formal variable order h,s,q,F,e; e is the old root-ending path count.
    # These five signature weights are the exact universal one-port cover
    # decomposition in ROOTED_RESETS R1, not fitted from source examples.
    carrier = {((-1, -1),): (1, -1, 0, 0, -1),
               ((-1, 0),): (0, 0, 0, 0, 1),
               ((0, -1),): (0, 1, 0, 0, 0),
               ((-1, -1), (0, 0)): (0, 0, 1, 0, 0),
               ((-1, -1), (-1, 0)): (0, 0, -1, 1, -1)}
    # Add w retaining (oldroot,w); add y at oldroot, forgetting oldroot;
    # add z at y, forgetting z.  Final slots are (y,w), and we observe w.
    stages = []
    for root, new, rename in [(0, 1, {}), (0, 2, {0: -1, 2: 0}),
                              (0, 2, {2: -1})]:
        carrier = kernel(carrier, root, new, rename)
        stages.append(len(carrier))
    matrix = observation(carrier, 1)
    expected = [(1, -1, -2, 5, 0), (0, 0, -1, 1, 0),
                (1, 0, 0, 2, 0), (4, -1, -1, 10, 0)]
    require(matrix == expected, 'universal coefficient identity failed')
    return {'formal_variables': ['h', 's', 'q', 'F', 'e'],
            'three_birth_chain_mark_matrix': matrix,
            'carrier_sizes_after_births': stages,
            'status': 'Universal coefficient identity from the endpoint-cover bijection'}


def apply_chain(rows, root, k):
    require(k >= 2, 'closed chart requires k>=2')
    source = state(rows, root)
    old_root = root
    for _ in range(k):
        rows, root = birth(rows, old_root)
    mark = len(rows) - 1
    rows, y = birth(rows, old_root, True)
    rows, _ = birth(rows, y)
    observed = state(rows, mark)
    h, s, q, F = source
    a, d, v = h - s, F - q, 2**(k - 1)
    new_s = v * d
    new_q = a + 3 * q + (3 * v - 1) * d
    expected = (new_q + 3 * new_s, new_s, new_q, 4 * new_q + 2 * new_s)
    require(observed == expected and strong(rows), 'literal retained-chain law failed')
    indegree = sum(r >> mark & 1 for r in rows)
    require(indegree == 3, 'wrong physical chain-root indegree')
    return rows, mark, {'k': k, 'rows': rows, 'source_root': old_root,
                        'retained_ports': [y, mark], 'root': mark,
                        'source_state': source, 'state_h_s_q_F': observed,
                        'root_indegree': indegree}


def inverse(q, s):
    require(q > 0 and q % 2 == 1 and s > 0 and s % 2 == 0, 'outside positive chart')
    v = s & -s
    k = v.bit_length()
    old_q = q - 3 * s
    numerator = s // v - 3 * old_q
    if old_q <= 0 or numerator <= 0 or numerator % 4:
        return None
    return k, old_q, numerator // 2


def replay_fixed_slots(rows, slots, active, k):
    initial_slots = tuple(slots)
    slots = list(slots)
    word = [1] * k + [0, 2] if active == 0 else [3] * k + [4, 5]
    for command in word:
        root_slot, discard = divmod(command, 3)
        newborn = len(rows)
        rows, _ = birth(rows, slots[root_slot])
        if discard < 2:
            slots[discard] = newborn
    trace = {'input_slots': initial_slots, 'input_active_slot': active,
             'word': word, 'output_slots': slots, 'output_active_slot': 1 - active}
    return rows, slots[1 - active], trace


def connector_certificate():
    records = []
    for b in range(3):
        for t in range(3):
            if b == t == 0:
                continue
            x, y = b, t
            guard_word = []
            if (x - 2 * y) % 3 == 0:
                x = (2 * x - y) % 3
                guard_word.append('R')
            require((x - 2 * y) % 3 != 0, 'C3 guard bridge failed')
            x, y = (61 * x - 2 * y) % 3, (9 * x) % 3
            guard_word.append('C3')
            require(x != 0 and y == 0, 'usual controller guard not reached')
            x, y = b, t
            chain_word = []
            if (x - y) % 3 == 0:
                require((x - 2 * y) % 3 != 0, 'chain bridge chose illegal C3')
                x, y = (61 * x - 2 * y) % 3, (9 * x) % 3
                chain_word.append('C3')
            require((x - y) % 3 != 0, 'primitive chain guard not reached')
            chain_word.append('L_k')
            records.append({'residue_b_t': [b, t], 'usual_guard_word': guard_word,
                            'primitive_chain_word': chain_word})
    return records


def main():
    certificate = universal_certificate()
    source, _ = wrapper(4, 1)
    rows, root, entry = apply_chain(source, 0, 2)
    replay, replay_root, trace = replay_fixed_slots(source, (0, 5), 0, 2)
    require(replay == rows and replay_root == root, 'entry fixed-slot replay failed')
    entry['fixed_slot_trace'] = trace
    require(entry['state_h_s_q_F'] == (249, 38, 135, 616), 'wrong H17 entry')
    require(inverse(135, 38) is None, 'entry should be outside the positive chart image')
    renewals = []
    for k in (2, 3):
        target, target_root, receipt = apply_chain(rows, root, k)
        replay, replay_root, trace = replay_fixed_slots(rows, entry['retained_ports'], 1, k)
        require(replay == target and replay_root == target_root, 'renewal fixed-slot replay failed')
        receipt['fixed_slot_trace'] = trace
        h, s, q, F = receipt['state_h_s_q_F']
        require(gcd(q, s) == 1 and inverse(q, s) == (k, 135, 38),
                'dyadic content or exact inverse failure')
        renewals.append(receipt)
    h9_repaired = (22, 8, 2, 21, 6)
    h9_parent, h9_root = birth(h9_repaired, 0, True)
    require(state(h9_parent, h9_root) == (15, 5, 9, 19), 'H9 preparatory N failed')
    h9_rows, h9_mark, h9_entry = apply_chain(h9_parent, h9_root, 2)
    require(h9_entry['state_h_s_q_F'] == (147, 20, 87, 388), 'H9 chain entry failed')
    replay, replay_root, trace = replay_fixed_slots(h9_parent, (0, 5), 1, 2)
    require(replay == h9_rows and replay_root == h9_mark, 'H9 fixed-slot replay failed')
    h9_entry['fixed_slot_trace'] = trace
    # Many iterations use the proved exact two-coordinate law, not large DP.
    q, s = 135, 38
    word = [2, 3, 2, 5, 4, 2, 3, 2]
    address = []
    for k in word:
        v = 2**(k - 1)
        next_s = v * (3 * q + 2 * s)
        next_q = q + 3 * next_s
        require(inverse(next_q, next_s) == (k, q, s), 'exact address decoding failed')
        require(gcd(next_q, next_s) == 1, 'renewal content changed')
        q, s = next_q, next_s
        address.append({'k': k, 'q': q, 's': s, 'H': q + 3 * s})
    out = {'coefficient_certificate': certificate, 'literal_H17_entry': entry,
           'literal_H9_N_entry': h9_entry,
           'all_nonzero_mod3_connector_cases': connector_certificate(),
           'literal_renewals': renewals, 'formula_only_long_address': address,
           'status': 'Universal law proved in CHAIN_REUSE.md; bounded literal DP controls'}
    (HERE / 'chain_reuse_receipt.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'passed': True, 'entry_H': 249, 'H9_N_entry_H': 147,
                      'literal_renewal_H': [r['state_h_s_q_F'][0] for r in renewals],
                      'universal_matrix': certificate['three_birth_chain_mark_matrix']}))


if __name__ == '__main__':
    main()
