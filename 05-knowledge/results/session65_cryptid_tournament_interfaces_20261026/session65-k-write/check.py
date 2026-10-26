#!/usr/bin/env python3
"""Independent finite controls for a fixed-label native-K tape-write obstruction.

No producer imports, graph libraries, random samples, or optimization-sensitive
assertions. The all-n / all-finite-chain result is proved in REPORT.md.
"""

import itertools
import json
import math


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sign(exponent):
    return -1 if exponent % 2 else 1


def pairs(n):
    return itertools.combinations(range(n), 2)


def from_mask(n, mask):
    matrix = [[0] * n for _ in range(n)]
    for bit, (a, b) in enumerate(pairs(n)):
        matrix[a][b] = (mask >> bit) & 1
        matrix[b][a] = 1 - matrix[a][b]
    return tuple(tuple(row) for row in matrix)


def chi(matrix):
    return sign(sum(matrix[a][b] for a, b in pairs(len(matrix))))


def changed_pairs(left, right):
    require(len(left) == len(right), "different orders")
    return tuple((a, b) for a, b in pairs(len(left))
                 if left[a][b] != right[a][b])


def directed(matrix, word):
    return all(matrix[a][b] for a, b in zip(word, word[1:]))


def native_k(matrix, word):
    n = len(matrix)
    require(sorted(word) == list(range(n)), "word is not Hamiltonian")
    require(directed(matrix, word), "word is not directed in current source")
    path_pairs = {frozenset((a, b)) for a, b in zip(word, word[1:])}
    result = [[0] * n for _ in range(n)]
    for a, b in pairs(n):
        result[a][b] = (matrix[a][b] if frozenset((a, b)) in path_pairs
                        else matrix[b][a])
        result[b][a] = 1 - result[a][b]
    return tuple(tuple(row) for row in result)


def audit_move(matrix, word):
    n = len(matrix)
    result = native_k(matrix, word)
    difference = changed_pairs(matrix, result)
    path_pairs = {tuple(sorted((a, b))) for a, b in zip(word, word[1:])}
    require(set(difference) == set(pairs(n)) - path_pairs,
            "native move does not flip exactly the non-path pairs")
    require(len(difference) == (n - 1) * (n - 2) // 2, "flip count")
    require(len(set(pairs(n)) - set(difference)) == n - 1, "retained count")
    require(chi(result) == sign((n - 1) * (n - 2) // 2) * chi(matrix),
            "native K parity law")
    require(directed(result, word), "path lost after native move")
    require(native_k(result, word) == matrix, "native K not involutive")
    return result


def permutation_sign(permutation):
    return sign(sum(permutation[a] > permutation[b]
                    for a, b in pairs(len(permutation))))


def relabel(matrix, old_to_new):
    n = len(matrix)
    result = [[0] * n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            result[old_to_new[a]][old_to_new[b]] = matrix[a][b]
    return tuple(tuple(row) for row in result)


def blue(matrix, word):
    reflection = dict(zip(word, reversed(word)))
    return all(matrix[reflection[a]][reflection[b]] == matrix[b][a]
               for a in range(len(matrix)) for b in range(len(matrix)))


def all_hamiltonian_words(matrix):
    """Independent adjacency-guided exhaustive enumerator for larger tapes."""
    n = len(matrix)

    def extend(prefix, unused):
        if not unused:
            yield prefix
            return
        for vertex in sorted(unused):
            if not prefix or matrix[prefix[-1]][vertex]:
                yield from extend(prefix + (vertex,), unused - {vertex})

    yield from extend((), set(range(n)))


def tape(word):
    n = 3 * len(word)
    matrix = [[int(a < b) for b in range(n)] for a in range(n)]
    for cell, bit in enumerate(word):
        require(bit in (0, 1), "nonbinary tape")
        if bit:
            a, b = 3 * cell, 3 * cell + 2
            matrix[a][b], matrix[b][a] = 0, 1
    return tuple(tuple(row) for row in matrix)


def exhaust_small():
    rows = []
    for n in range(1, 6):
        N = n * (n - 1) // 2
        C = (n - 1) * (n - 2) // 2
        words = list(itertools.permutations(range(n)))
        legal_moves = relabellings = odd_relabs = 0
        for mask in range(1 << N):
            matrix = from_mask(n, mask)
            for word in words:
                renamed = relabel(matrix, word)
                factor = permutation_sign(word)
                require(chi(renamed) == factor * chi(matrix), "relabel sign law")
                relabellings += 1
                odd_relabs += int(factor == -1)
                if directed(matrix, word):
                    audit_move(matrix, word)
                    legal_moves += 1
        expected = math.factorial(n) * (1 << C)
        require(legal_moves == expected, "exhaustive legal-move count")
        rows.append(dict(n=n, labelled_tournaments=1 << N,
                         legal_K_moves=legal_moves, retained_arcs=n - 1,
                         reversed_arcs=C, chi_multiplier=sign(C),
                         relabellings=relabellings, odd_relabellings=odd_relabs))
    return rows


def tape_hostiles():
    rows = []
    for m in (2, 3):
        n = 3 * m
        write_count = move_count = odd_relabs = padding_writes = 0
        padding_cells = next(r for r in (1, 2) if (m + r) % 4 in (0, 1))
        for word in itertools.product((0, 1), repeat=m):
            matrix = tape(word)
            paths = list(all_hamiltonian_words(matrix))
            require(len(paths) == 3 ** sum(word), "ordered-join path count")
            targets = []
            for cell in range(m):
                written_word = list(word)
                written_word[cell] ^= 1
                target = tape(written_word)
                require(changed_pairs(matrix, target) == ((3 * cell, 3 * cell + 2),),
                        "tape write not exactly its designated arc flip")
                require(chi(target) == -chi(matrix), "single-write parity")
                targets.append(target)
                write_count += 1

                padded = tape(word + (0,) * padding_cells)
                padded_target = tape(tuple(written_word) + (0,) * padding_cells)
                Npad = len(padded)
                Cpad = (Npad - 1) * (Npad - 2) // 2
                require(Cpad % 2 == 1, "padding did not remove conserved parity")
                require(len(changed_pairs(padded, padded_target)) == 1,
                        "common padding altered literal write size")
                require(chi(padded_target) == -chi(padded), "padded write parity")
                output = audit_move(padded, tuple(range(Npad)))
                require(chi(output) == chi(padded_target), "odd K not parity-compatible")
                require(output != padded_target, "padding falsely gave a one-K write")
                padding_writes += 1

            for path in paths:
                output = audit_move(matrix, path)
                require(chi(output) == chi(matrix), "blocked tape order not conserved")
                require(all(output != target for target in targets), "one-K hostile failed")
                move_count += 1

            permutation = (1, 0) + tuple(range(2, n))
            renamed = relabel(matrix, permutation)
            require(chi(renamed) == -chi(matrix), "odd relabel hostile failed")
            require(relabel(renamed, permutation) == matrix, "not the same isomorphism class")
            odd_relabs += 1

        rows.append(dict(m=m, n=n, words=1 << m, directed_one_bit_writes=write_count,
                         legal_K_moves=move_count, odd_relabel_hostiles=odd_relabs,
                         padding_zero_cells=padding_cells, padded_order=12,
                         padded_write_controls=padding_writes,
                         padding_claim="only removes conserved-parity obstruction"))
    return rows


def legal_positive_control():
    source, target = tape((0,)), tape((1,))
    path = (0, 1, 2)
    require(audit_move(source, path) == target, "n=3 legal all-K write")
    require(blue(source, path), "n=3 control unexpectedly black")
    require(len(list(all_hamiltonian_words(source))) == 1, "transitive path count")
    return dict(n=3, source="G0", target="G1", word=list(path),
                legal=True, blue=True, black_witness=False,
                role="all-K legality control only, not a black realization")


def main():
    small = exhaust_small()
    hostiles = tape_hostiles()
    result = dict(verdict="PASS", implementation="independent adjacency matrices",
                  exhaustive_small_orders=small,
                  total_small_legal_K_moves=sum(row["legal_K_moves"] for row in small),
                  total_small_relabellings=sum(row["relabellings"] for row in small),
                  tape_write_hostiles=hostiles,
                  odd_C_control=legal_positive_control(),
                  sequence_scope="all finite fixed-label chains, by REPORT.md proof",
                  not_claimed=["unlabelled-class obstruction", "relabelling-enabled obstruction",
                               "padding or larger-gadget obstruction", "black realization at odd C",
                               "universal Turing-machine simulation impossibility"])
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
