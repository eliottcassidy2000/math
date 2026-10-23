"""Exact labelled curl fibres and authenticated Hamiltonian port transport.

Independent of repository enumeration producers; Python standard library only.
Run: python -B .../check.py --out .../receipt.json
No assert statements: checks survive python -O.
"""
from collections import Counter, defaultdict
from itertools import combinations, permutations, product
import argparse
import hashlib
import json
from pathlib import Path

CHECKS = 0


def require(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def edges(n):
    return list(combinations(range(n), 2))


def matrix(n, mask):
    a = [[0] * n for _ in range(n)]
    for b, (i, j) in enumerate(edges(n)):
        a[i][j] = 1 if mask >> b & 1 else -1
        a[j][i] = -a[i][j]
    return a


def key(a):
    return tuple(a[i][j] for i, j in edges(len(a)))


def curl(a, i, j, k):
    return a[i][j] + a[j][k] + a[k][i]


def rooted(a):
    return tuple(curl(a, 0, i, j) for i, j in combinations(range(1, len(a)), 2))


def score(a, r=0):
    return sum(v == 1 for v in a[r])


def components(a):
    """SCCs in source-to-sink order, using exact transitive closure."""
    n = len(a)
    reach = [(1 << i) | sum(1 << j for j in range(n) if a[i][j] == 1) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    left = set(range(n)); blocks = []
    while left:
        i = min(left)
        block = tuple(j for j in sorted(left) if reach[i] >> j & 1 and reach[j] >> i & 1)
        left.difference_update(block); blocks.append(block)
    order = {b: -sum(a[b[0]][c[0]] == 1 for c in blocks if c != b) for b in blocks}
    blocks.sort(key=order.__getitem__)
    return blocks


def decode(n, kappa, root_score=None):
    """All realizations; with a score, the result has size zero or one.

    The sign tournament gives a totally ordered SCC list. All feasible root
    stars are prefix cuts. Input rejection is part of the representation.
    """
    if type(n) is not int or n < 1:
        return []
    if root_score is not None and (type(root_score) is not int or not 0 <= root_score < n):
        return []
    pairs = list(combinations(range(1, n), 2))
    if len(kappa) != len(pairs) or any(type(c) is not int or c not in (-3, -1, 1, 3) for c in kappa):
        return []
    if n == 1:
        return [[[0]]] if root_score in (None, 0) else []
    cdict = dict(zip(pairs, kappa))
    s = [[0] * (n - 1) for _ in range(n - 1)]
    for (i, j), c in cdict.items():
        s[i - 1][j - 1] = 1 if c > 0 else -1
        s[j - 1][i - 1] = -s[i - 1][j - 1]
    blocks = components(s)
    ans = []; plus = set()
    for cut in range(len(blocks) + 1):
        if root_score is None or len(plus) == root_score:
            x = {i: (1 if i - 1 in plus else -1) for i in range(1, n)}
            a = [[0] * n for _ in range(n)]
            good = True
            for i in range(1, n):
                a[0][i] = x[i]; a[i][0] = -x[i]
            for (i, j), c in cdict.items():
                a[i][j] = c - x[i] + x[j]
                a[j][i] = -a[i][j]
                if abs(a[i][j]) != 1:
                    good = False
            if good:
                ans.append(a)
        if cut < len(blocks):
            plus.update(blocks[cut])
    return ans


def rotated(a, blocks, cut):
    n = len(a); b = [row[:] for row in a]
    prefix = {v for block in blocks[:cut] for v in block}
    for i, j in edges(n):
        if (i in prefix) != (j in prefix):
            b[i][j] = -b[i][j]; b[j][i] = -b[j][i]
    return b


def path_partner(a):
    n = len(a); b = [row[:] for row in a]
    for i, j in edges(n):
        if j != i + 1:
            b[i][j] = -b[i][j]; b[j][i] = -b[j][i]
    return b


def blue(a):
    n = len(a)
    return all(a[i][j] == -a[n - 1 - i][n - 1 - j] for i, j in edges(n))


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    rows = []
    for n in range(2, 7):
        fibres = defaultdict(list)
        for mask in range(1 << len(edges(n))):
            a = matrix(n, mask); fibres[rooted(a)].append(a)
        sizes = Counter()
        for kappa, members in fibres.items():
            a = members[0]; blocks = components(a)
            actual = {key(b) for b in members}
            expected = {key(rotated(a, blocks, cut)) for cut in range(len(blocks))}
            require(actual == expected, ('SCC rotation fibre', n, kappa))
            decoded = {key(b) for b in decode(n, kappa)}
            require(decoded == actual, ('decoder fibre', n, kappa))
            require(len({score(b) for b in members}) == len(members), ('root score collision', n, kappa))
            for b in members:
                d = decode(n, kappa, score(b))
                require(len(d) == 1 and key(d[0]) == key(b), ('pointed reconstruction', n, kappa))
                require(len(components(b)) == len(members), ('fibre cardinality', n, kappa))
            sizes[len(members)] += 1
        rows.append({'n': n, 'labelled_tournaments': 1 << len(edges(n)),
                     'root_curl_fibres': len(fibres), 'fibre_histogram': dict(sorted(sizes.items()))})
    ports = []
    for n in range(3, 8):
        free = [(i, j) for i, j in edges(n) if j != i + 1]
        signatures = set(); blue_count = 0
        for mask in range(1 << len(free)):
            a = [[0] * n for _ in range(n)]
            for i in range(n - 1):
                a[i][i + 1] = 1; a[i + 1][i] = -1
            for bit, (i, j) in enumerate(free):
                a[i][j] = 1 if mask >> bit & 1 else -1; a[j][i] = -a[i][j]
            kappa = rooted(a)
            require(kappa not in signatures, ('fixed-path curl collision', n, mask))
            signatures.add(kappa)
            b = path_partner(a)
            pcurl = tuple((1 if i == 1 else 0) + (1 if j == i + 1 else 0)
                          for i, j in combinations(range(1, n), 2))
            target = tuple(-c + 2 * p for c, p in zip(kappa, pcurl))
            target_score = n - score(a)  # root is the start, so div(p)=1
            require(target == rooted(b), ('affine curl transport', n, mask))
            require(target_score == score(b), ('score transport', n, mask))
            d = decode(n, target, target_score)
            require(len(d) == 1 and key(d[0]) == key(b), ('authenticated target', n, mask))
            require(key(path_partner(b)) == key(a), ('involution', n, mask))
            require((a[n - 1][0] == 1) != (b[n - 1][0] == 1), ('closed/open', n, mask))
            require(blue(a) == blue(b), ('colour', n, mask))
            blue_count += blue(a)
        m = len(free); f = (n - 1) // 2
        require(blue_count == 1 << ((m + f) // 2), ('blue orbit count', n))
        ports.append({'n': n, 'ports': 1 << m, 'closed_open_pairs': 1 << (m - 1), 'blue_ports': blue_count})
    # Hostile: unsigned cycle support cannot distinguish the two C3 orientations.
    c3 = [[0, 1, -1], [-1, 0, 1], [1, -1, 0]]
    opposite = [[-v for v in row] for row in c3]
    require(score(c3) == score(opposite) == 1, 'cycle hostile score')
    require(abs(rooted(c3)[0]) == abs(rooted(opposite)[0]) == 3, 'cycle hostile support')
    require(rooted(c3) != rooted(opposite), 'cycle hostile sign')
    # Hostile: allowed local curl alphabet does not ensure global realizability.
    invalid = (3, 3, 3)
    require(decode(4, invalid) == [], 'inconsistent three positive curl maxima')
    # Hostile: scalar score is indispensable without a directed path marking.
    transitive = [[0 if i == j else (1 if i < j else -1) for j in range(4)] for i in range(4)]
    fibre = decode(4, rooted(transitive))
    require(len(fibre) == 4 and sorted(score(a) for a in fibre) == [0, 1, 2, 3], 'missing root score')
    require(decode(3, (1.0,), 2) == [], 'floating curl rejected')
    require(decode(3, (True,), 2) == [], 'boolean curl rejected')
    require(decode(3, (1,), 2.0) == [], 'floating score rejected')
    require(decode(3, (1,), True) == [], 'boolean score rejected')
    # Literal five-vertex parity witness, independently counting permutations.
    witness = [[0 if i == j else (1 if i < j else -1) for j in range(5)] for i in range(5)]
    witness[2][0] = 1; witness[0][2] = -1
    witness_rotations = [rotated(witness, [(0, 1, 2), (3,), (4,)], c) for c in range(3)]
    parity_rows = []
    for a in witness_rotations:
        pp = list(permutations(range(5)))
        h = sum(all(a[p[i]][p[i + 1]] == 1 for i in range(4)) for p in pp)
        aut = sum(all(a[i][j] == a[p[i]][p[j]] for i, j in edges(5)) for p in pp)
        anti = sum(all(a[i][j] == -a[p[i]][p[j]] for i, j in edges(5)) for p in pp)
        require((h, aut) == (3, 3), 'parity witness H and Aut')
        parity_rows.append({'H': h, 'Aut': aut, 'anti': anti, 'root_score': score(a)})
    require([r['anti'] for r in parity_rows] == [0, 0, 3], 'SC cut changes')
    # Dihedral law on realized strong types: 0=K1, 1=the chiral matrix(6,18), 2=its converse.
    chiral = matrix(6, 18)
    require(len(components(chiral)) == 1, 'chiral letter strong')
    require(not any(all(chiral[i][j] == -chiral[p[i]][p[j]] for i, j in edges(6))
                    for p in permutations(range(6))), 'chiral letter not self converse')
    dual = {0: 0, 1: 2, 2: 1}; necklace_cases = 0
    for k in range(1, 9):
        for w in product(range(3), repeat=k):
            rotations = list(dict.fromkeys(w[t:] + w[:t] for t in range(k)))
            if w != min(rotations):
                continue
            necklace_cases += 1; d = len(rotations)
            reverse_dual = tuple(dual[v] for v in reversed(w))
            actual_sc = sum(v == tuple(dual[x] for x in reversed(v)) for v in rotations)
            if reverse_dual in rotations:
                a = rotations.index(reverse_dual)
                predicted_sc = sum((2 * t - a) % d == 0 for t in range(d))
                require(actual_sc == predicted_sc, ('dihedral fixed cuts', w))
                require(actual_sc == (1 if d % 2 else (2 if a % 2 == 0 else 0)), ('odd even buckets', w))
            else:
                require(actual_sc == 0, ('paired necklaces no SC', w))
    result = {'scope': 'FINITE-EXACT: all labelled tournaments n=2..6; all fixed-path ports n=3..7',
              'labelled': rows, 'ports': ports, 'hostiles': {'unsigned_cycle_support': 'C3 and its converse, same root score',
              'unrealizable_root_curls': [3, 3, 3], 'missing_score': 'four transitive cyclic rotations'},
              'five_vertex_cut_witness': parity_rows, 'dihedral_necklaces_through_length_eight': necklace_cases,
              'chiral_strong_letter': {'n': 6, 'mask': 18},
              'checks': CHECKS, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    payload = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.out:
        args.out.write_text(payload, encoding='utf-8')
    print(payload)


if __name__ == '__main__':
    main()
