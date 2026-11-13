#!/usr/bin/env python3
"""Exact, bounded controls for SHARED_DICTIONARY.md; standard library only."""
from collections import Counter
from itertools import product
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


F = {-1: 1, 1: 1, 2: 1, 3: 2}
for index in range(4, 24):
    F[index] = F[index - 1] + F[index - 2]


def legal_support_bank(top, deleted):
    """Independent exhaustive slot-subset enumeration, not the fibre formula."""
    slots = [-1] + list(range(1, top + 1))
    bank = [[] for _ in range(F[top] + 1)]
    for mask in range(1 << len(slots)):
        chosen_positions = [j for j in range(len(slots)) if mask >> j & 1]
        chosen = [slots[j] for j in chosen_positions]
        if deleted:
            legal = all(b - a != 1 for a, b in zip(chosen_positions, chosen_positions[1:]))
        else:
            legal = all(b - a != 1 for a, b in zip(chosen, chosen[1:]))
        if not legal:
            continue
        value = sum(F[j] for j in chosen)
        if value <= F[top]:
            bank[value].append(tuple(reversed(chosen)))
    return bank


def tagged_blocks(top):
    result = {}
    for n in range(3, top + 1):
        if n % 2:
            result[n] = sum((result[i] for i in range(n - 1, 3, -2)), [])
            result[n] += [('C', n), ('A', 0)]
        else:
            result[n] = sum((result[i] for i in range(n - 1, 2, -2)), [])
            result[n] += [('B', n)]
        require(len(result[n]) == F[n], 'tagged length')
        if n > 3:
            require(result[n][:F[n - 1]] == result[n - 1], 'tagged nesting')
    return result


def all_dictionary_controls(top, deleted):
    bank = legal_support_bank(top, deleted)
    tags = tagged_blocks(top)[top]
    proper = {}
    for n in range(3, top + 1):
        proper[n] = [s for s in bank[F[n]] if max(s) < n]
        require(len(proper[n]) == (1 if deleted and n % 2 else 2), 'proper count')
        if n >= 4:
            require(all(s[0] == n - 1 for s in proper[n]), 'largest proper block')
    free = [n for n in range(3, top + 1) if len(proper[n]) == 2]
    prefixes = [set() for _ in range(F[top] + 1)]
    dictionaries = 0
    row_checks = 0
    for bits in product((0, 1), repeat=len(free)):
        choices = dict(zip(free, bits))
        blocks = {-1: 'A', 1: 'B', 2: 'C'}
        for n in range(3, top + 1):
            blocks[n] = ''.join(blocks[i] for i in proper[n][choices.get(n, 0)])
            if n > 3:
                require(blocks[n].startswith(blocks[n - 1]), 'dictionary nesting')
        stripe = blocks[top]
        owner_values = {}
        for position, (canonical, owner) in enumerate(tags):
            if owner == 0:
                require(stripe[position] == 'A', 'forced A')
            else:
                require(stripe[position] in ({'C', 'B'} if owner % 2 else {'B', 'A'}), 'owner alphabet')
                if deleted and owner % 2:
                    require(stripe[position] == 'C', 'deleted odd owner')
                previous = owner_values.setdefault(owner, stripe[position])
                require(previous == stripe[position], 'owner synchronization')
        for N in range(1, len(stripe) + 1):
            languages = {''.join(blocks[i] for i in support) for support in bank[N]}
            require(stripe[:N] in languages, 'prefix has no legal outer support')
            prefixes[N].add(stripe[:N])
            row_checks += 1
        dictionaries += 1
    for N in range(1, F[top] + 1):
        visible = sum((F[n] - (n % 2) <= N) for n in free)
        require(len(prefixes[N]) == 2 ** visible, 'exact prefix count')
    return {'grammar': 'D' if deleted else 'R', 'top_index': top,
            'horizon': F[top], 'dictionaries': dictionaries, 'row_checks': row_checks,
            'legal_supports': sum(map(len, bank)), 'prefix_count_checks': F[top]}


def pair_polynomial(tags, N):
    constant = 0
    linear = Counter()
    quadratic = Counter()
    for a in range(1, (N + 1) // 2):
        b = N - a
        ca, ia = tags[a - 1]
        cb, ib = tags[b - 1]
        if 'C' in (ca, cb):
            continue
        if ca == cb == 'A':
            constant += 1
        elif ca == cb == 'B':
            if ia == ib:
                linear[ia] += 1
            else:
                quadratic[tuple(sorted((ia, ib)))] += 1
        else:
            linear[ia if ca == 'B' else ib] += 1
    return constant, linear, quadratic


def pair_controls():
    tags = tagged_blocks(16)[16]
    first_gap = None
    scalar_checks = 0
    gap_values = None
    for N in range(3, 129):
        owners = sorted({i for c, i in tags[:N - 1] if c == 'B'})
        constant, linear, quadratic = pair_polynomial(tags, N)
        values = set()
        for bits in product((0, 1), repeat=len(owners)):
            select = dict(zip(owners, bits))
            predicted = constant + sum(v * select[i] for i, v in linear.items())
            predicted += sum(v * select[i] * select[j] for (i, j), v in quadratic.items())
            A = lambda pos: tags[pos - 1][0] == 'A' or (tags[pos - 1][0] == 'B' and select[tags[pos - 1][1]])
            actual = sum(bool(A(a) and A(N - a)) for a in range(1, (N + 1) // 2))
            require(predicted == actual, 'owner polynomial differs from literal pair count')
            values.add(actual)
            scalar_checks += 1
        if values != set(range(min(values), max(values) + 1)) and first_gap is None:
            first_gap, gap_values = N, sorted(values)
    require(first_gap == 16 and gap_values == [0, 2], 'first selector gap')
    c, linear, quadratic = pair_polynomial(tags, 16)
    require(c == 0 and linear == {4: 2} and not quadratic, 'target16 polynomial')
    c18, linear18, quadratic18 = pair_polynomial(tags, 18)
    require(c18 == 1 and linear18 == {4: 3, 6: 1} and not quadratic18, 'content-one target18 polynomial')
    values18 = sorted({c18 + 3 * a + b for a, b in product((0, 1), repeat=2)})
    require(values18 == [1, 2, 4, 5], 'content-one gap')
    independent = set()
    Bpos = [n for n in range(1, 16) if tags[n - 1][0] == 'B']
    for bits in product((0, 1), repeat=len(Bpos)):
        Apos = {n for n in range(1, 16) if tags[n - 1][0] == 'A'}
        Apos.update(n for n, bit in zip(Bpos, bits) if bit)
        independent.add(sum(a in Apos and 16 - a in Apos for a in range(1, 8)))
    require(independent == {0, 1, 2}, 'independent repair')
    require(tags[0] == tags[34] == ('C', 3), 'historical35 same owner')
    return {'targets_checked': 126, 'owner_assignments_checked': scalar_checks,
            'first_noninterval_target': first_gap, 'shared_values': gap_values,
            'independent_values': sorted(independent),
            'content_one_target18_values': values18,
            'historical35_owner': 3}


def main():
    receipt = {'status': 'PASS', 'dictionary_controls': [all_dictionary_controls(11, False),
               all_dictionary_controls(12, True)], 'pair_controls': pair_controls(),
               'scope': 'Exact finite controls for the proved shared-dictionary classification; no HP spectrum claim.'}
    print(json.dumps(receipt, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
