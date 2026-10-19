#!/usr/bin/env python3
"""Finite closure certificate and independent direct-triangle checks.

No assertion statements: every verification remains active under python -O.
Uses Python standard library only. No conjectural decision horizon is used.
"""
from collections import deque
from itertools import product
from math import comb
import json

ALPHABET = ((0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2))
Z, A, B = (0, 0), (2, 0), (2, 2)


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def advance(state, pair):
    """Three-state DFA for avoiding A Z* B and B Z* A; -1 rejects."""
    if pair not in ALPHABET:
        return -1
    if pair == Z:
        return state
    if pair == A:
        return -1 if state == 2 else 1
    if pair == B:
        return -1 if state == 1 else 2
    return 0


def pair_transition(q, pair):
    x = abs(q - pair[0])
    y = abs(x - pair[1])
    return (x, y), y


def belongs(vector):
    need(len(vector) % 2 == 0, 'The invariant checker needs whole pairs.')
    state = 0
    for j in range(0, len(vector), 2):
        state = advance(state, tuple(vector[j:j+2]))
        if state == -1:
            return False
    return True


def T(vector, bit):
    output = []
    q = bit
    for x in vector:
        q = abs(q - x)
        output.append(q)
    return tuple(output)


def certify_closure():
    """Explore the finite product: input DFA, output DFA, transducer state.

    Every allowed input of arbitrary finite pair length is represented by a
    path from (0,0,b). Proving no reachable edge rejects the output proves
    invariant closure at all lengths, rather than checking sample lengths.
    """
    initial = {(0, 0, 0), (0, 0, 1)}
    seen = set(initial)
    queue = deque(sorted(initial))
    edges = []
    while queue:
        state = queue.popleft()
        source_monitor, output_monitor, q = state
        for pair in ALPHABET:
            next_source = advance(source_monitor, pair)
            if next_source < 0:
                continue
            output_pair, next_q = pair_transition(q, pair)
            next_output = advance(output_monitor, output_pair)
            need(next_output >= 0,
                 f'Invariant closure fails: {state}, {pair} -> {output_pair}')
            target = next_source, next_output, next_q
            edges.append([list(state), list(pair), list(output_pair), list(target)])
            if target not in seen:
                seen.add(target)
                queue.append(target)
    return {'reachable_product_states': [list(s) for s in sorted(seen)],
            'reachable_product_state_count': len(seen),
            'transition_count': len(edges), 'transitions': edges}


def check_initial_phase():
    count = 0
    vector = (3,) * 128
    for z in range(64):
        observed = T(vector, 1)
        expected = tuple(1 + (comb(z+j+1, j) & 1) for j in range(128))
        need(observed == expected, f'First-one formula fails at z={z}')
        need(belongs(observed), f'Initial invariant fails at z={z}')
        permitted = {(1, 1), (2, 1)} if z % 2 == 0 else {(1, 1), (2, 2)}
        need(all(observed[j:j+2] in permitted for j in range(0,128,2)),
             f'Initial pair parity classification fails at z={z}')
        vector = T(vector, 0)
        count += 1
    return {'zero_prefixes_tested': count, 'age_coordinates_per_prefix': 128}


def direct_edge(word):
    row = tuple(word)
    answer = []
    while row:
        answer.append(row[0])
        row = tuple(abs(a-b) for a,b in zip(row,row[1:]))
    return tuple(answer)


def pascal_front(word):
    reverse = tuple(reversed(word))
    return tuple(sum(comb(s,j)*reverse[j] for j in range(s+1)) & 1
                 for s in range(len(word)))


def check_direct_triangles():
    count = 0
    # This is an independent orientation/index check, not the proof of closure.
    for L in range(11):
        for word in product((0,1), repeat=L):
            vector = (3,)*32
            alpha = pascal_front(word)
            for bit in alpha:
                vector = T(vector, bit)
            edge = direct_edge(word + (3,) + (0,)*31)
            need(vector == edge[L:], f'Rotated/direct mismatch: {word}')
            if 1 in alpha:
                need(belongs(vector), f'Observed orbit violates invariant: {word}')
            first = next((j for j,value in enumerate(vector) if value >= 2), None)
            need(first is None or first % 2 == 0,
                 f'Odd first defect in direct triangle: {word}, age={first}')
            count += 1
    return {'left_words_checked': count, 'maximum_front_length': 10,
            'age_coordinates_checked': 32}


def check_hostiles():
    first = (2,0,2,2)
    second = (2,2,2,0)
    need(T(first,0) == (2,2,0,2), 'First hostile calculation changed.')
    need(T(second,0) == first, 'Reverse hostile calculation changed.')
    need(not belongs(first) and not belongs(second), 'Hostile accepted by invariant.')
    # For unrestricted boundary q=2, even the legal one-pair word 22 fails.
    need(pair_transition(2,B)[0] == (0,2), 'Boundary hostile changed.')
    return {'without_A_Zstar_B': {'input': first, 'T0': T(first,0)},
            'without_B_Zstar_A': {'input': second, 'T0': T(second,0),
                                   'T0_squared': T(T(second,0),0)},
            'nonbinary_boundary': {'input': B, 'boundary': 2, 'output': (0,2)}}


def main():
    table = [{'input': pair,
              'outputs_by_incoming_state_0_1_2':
                 [pair_transition(q,pair)[0] for q in range(3)]}
             for pair in ALPHABET]
    certificate = {'status': 'ALL-LENGTH FINITE INVARIANT CERTIFICATE',
                   'pair_transducer_table': table,
                   'closure': certify_closure(),
                   'initial_formula_checks': check_initial_phase(),
                   'independent_direct_triangle_checks': check_direct_triangles(),
                   'hostile_controls': check_hostiles()}
    print(json.dumps(certificate, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
