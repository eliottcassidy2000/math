#!/usr/bin/env python3
"""Independent bounded audit of ROOTED_RESETS; never imports its producer code."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product as words
import importlib.util
import json

HERE = Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE/(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


G = load('observer_induction_probe')
B = load('observer_induction_basis')


def positive(profile):
    h, s, q, F = profile
    return h-s, q-s, s, F-h-q+s


def scalar_rows(M):
    return [[int(x) if x.denominator == 1 else str(x) for x in row] for row in M]


def main():
    source = HERE.parent/'adaptive_state_induction_20261103'
    old = json.loads((source/'live_port.json').read_text())
    common = json.loads((source/'two_port_envelope.json').read_text())
    receipt = json.loads((HERE/'rooted_resets.json').read_text())
    queries = old['rank']['suffixes']
    A = [[[Q(x) for x in row] for row in M] for M in old['rank']['query_matrices']]
    F = [list(map(Q, row)) for row in receipt['sealing_factor']]
    # Four actual rooted sources span the profile coordinates; direct graph
    # queries recover the entire factor without using its five-cover code.
    sources = [((0,), 0), ((2, 0), 1), ((2, 4, 1), 0), ((6, 4, 0), 2)]
    U, Y = [], []
    for T, root in sources:
        U.append(positive(G.profile(T, root)))
        child = G.birth(T, root)
        Y.append([G.H(G.run(child, (root, len(T)), word)[0]) for word in queries])
    recovered = B.product(B.inverse(U), Y)
    G.need(recovered == F, 'factor recovered from 68 independent actual graph counts')
    # Independently construct port swap on the universal 25 endpoint carrier.
    def column(word):
        v = common['observer']
        for a in reversed(word):
            w = [0]*25
            for i, j, count in common['matrices'][a]:
                w[i] += count*v[j]
            v = w
        return v
    C = B.transpose([column([3*(1 << (a//3))+a%3 for a in word]) for word in queries])
    states = [tuple(map(tuple, state)) for state in common['states']]
    indices = {state: i for i, state in enumerate(states)}
    perm = [indices[tuple(sorted((1-s if s >= 0 else s, 1-e if e >= 0 else e)
                                 for s, e in state))] for state in states]
    _, rows = B.reduce(B.transpose(C))
    swap = B.product(B.inverse([C[i] for i in rows]), [C[perm[i]] for i in rows])
    G.need(B.product(C, swap) == [C[i] for i in perm], 'full carrier port permutation')
    cols = list(range(4))
    minor_inverse = B.inverse([[row[j] for j in cols] for row in F])
    macros = {}
    for length in range(3):
        for prefix in words((0, 2, 4, 5), repeat=length):
            for last in (1, 3):
                word = prefix+(last,)
                FA = F
                for a in word:
                    FA = B.product(FA, A[a])
                if last == 3:
                    FA = B.product(FA, swap)
                M = B.product([[row[j] for j in cols] for row in FA], minor_inverse)
                G.need(B.product(M, F) == FA, 'macro equality on all 17 queries')
                G.need(all(x.denominator == 1 and x >= 0 for row in M for x in row), 'positive integral macro')
                macros[word] = M
    G.need(len({tuple(map(tuple, M)) for M in macros.values()}) == 28, '28 distinct matrices')
    for item in receipt['grammar']['alphabet']:
        G.need(macros[tuple(item['word'])] == item['matrix'], 'every supplied macro matrix')
    # A separate physical-slot controller replays every length<=3 macro,
    # then two composite words exercise repeated normalization.
    cases = [((3,), (word,)) for word in macros]
    cases += [((1,), ((0, 4, 1), (4, 0, 3))), ((3,), ((4, 0, 3), (0, 4, 1)))]
    replays = []
    for entry, sequence in cases:
        T, slots, controller = (2, 4, 1), [0, 1], [0, 1]
        raw = []
        expected = positive(G.profile(T, 0 if entry[-1] == 1 else 1))
        for word in (entry,)+sequence:
            for letter in word:
                root, discard = divmod(letter, 3)
                actual_root = controller[root]
                actual_discard = controller[discard] if discard < 2 else 2
                before, root_name = T, slots[actual_root]
                T = G.birth(T, root_name)
                raw.append(3*actual_root+actual_discard)
                if actual_discard < 2:
                    slots[actual_discard] = len(before)
            if word[-1] == 3:
                controller.reverse()
            if word != entry or len(raw) > len(entry):
                expected = tuple(B.product([expected], macros[word])[0])
        G.need(slots[controller[1]] == len(T)-1, 'last interface child is newborn')
        G.need(slots[controller[0]] == root_name, 'retained root is parent root')
        G.need(positive(G.profile(before, root_name)) == expected, 'literal parent profile after omitting final birth')
        replays.append({'entry': entry, 'macros': sequence, 'raw': raw, 'parent_H': G.H(before)})
    # Scope hostile: universal divisor monotonicity alone does not imply
    # strict count growth outside the strong source language.
    weak = (30, 28, 8, 16, 4)  # TT2 -> C3, the two TT vertices are live.
    weak_queries = [G.H(G.run(weak, (0, 1), word)[0]) for word in queries]
    from math import gcd
    G.need(gcd(*weak_queries) == 3 and G.H(G.run(weak, (0, 1), (2,))[0]) == 3,
           'weak divisor-three source can preserve prime three')
    print(json.dumps({
        'scope': 'independent graph/profile and coefficient audit; no bounded target saturation',
        'actual_factor_counts': 68, 'profile_source_det': str(B.determinant(U)),
        'all_macro_words': len(macros), 'distinct_macros': 28,
        'literal_compiler_replays': replays,
        'weak_prime_scope_hostile': {'rows': weak, 'ports': [0, 1], 'D': 3, 'child_command': 2, 'child_H': 3},
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
