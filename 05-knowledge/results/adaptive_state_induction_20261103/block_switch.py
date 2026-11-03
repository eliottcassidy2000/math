#!/usr/bin/env python3
"""Exact delayed marking of a transitive twin block; finite source-image probes."""
from pathlib import Path
import importlib.util
from collections import deque
from fractions import Fraction
import json

provider = Path(__file__).resolve().parents[1]/'adaptive_metagraph_reformulation_20261103'/'adaptive_probe.py'
spec = importlib.util.spec_from_file_location('prior_adaptive', provider)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def block(z, A, B):
    h, s, q, F = z
    a, d, k = h-s, F-q, A+B+1
    hp = a+q+(2**k-1)*d
    sp = 2**B*d if A else a+(2**B-1)*d
    qp = h if k == 1 else a+q+(2**(k-1)-1)*d
    Fp = 2**A*a+(2**(B+1)*3**A-2**A)*d+2**(A+1)*q+(s-q if B == 0 else 0)
    return hp, sp, qp, Fp


def positive(z):
    h, s, q, F = z
    return h-s, q-s, s, F-h-q+s


def profile(v):
    a, t, s, delta = v
    return a+s, s, t+s, a+t+s+delta


def det(M):
    rows = [list(map(Fraction, r)) for r in M]
    answer = Fraction(1)
    for j in range(len(rows)):
        k = next((i for i in range(j, len(rows)) if rows[i][j]), None)
        if k is None:
            return Fraction(0)
        if k != j:
            rows[k], rows[j] = rows[j], rows[k]
            answer = -answer
        pivot = rows[j][j]
        answer *= pivot
        for i in range(j+1, len(rows)):
            c = rows[i][j]/pivot
            rows[i] = [x-c*y for x, y in zip(rows[i], rows[j])]
    return answer


def direct_block(T, v, A, B):
    n = len(T)
    for _ in range(A+B+1):
        T = old.add(T, v)
    return T, n+B


def checks():
    total = 0
    for layer in old.inventories(5)[1:]:
        for T in sorted(layer):
            for v, z in enumerate(old.profiles(T)):
                need(min(positive(z)) >= 0, 'physical positive coordinates')
                for k in range(1, 5):
                    for B in range(k):
                        A = k-B-1
                        U, w = direct_block(T, v, A, B)
                        want = old.profiles(U)[w]
                        got = block(z, A, B)
                        need(got == want, 'actual block profile')
                        rz = z
                        for _ in range(B):
                            rz = old.keep(rz)
                        need(block(rz, A, 0) == got, 'finite-alphabet decomposition')
                        if A:
                            need(got[0]-got[2] == 2**A*got[1], 'marked block chart')
                        else:
                            need(positive(got)[3] == 0, 'newborn chart')
                        total += 1
                old.table.cache_clear()
                old.profiles.cache_clear()
    return total


def decode(seed_root, word):
    T, root = old.add((2, 4, 1), 0), seed_root
    roots = [0]
    for letter in word:
        if letter == 'R':
            T = old.add(T, root)
            roots.append(root)
        else:
            A = int(letter[1:])
            n = len(T)
            for _ in range(A+1):
                T = old.add(T, root)
                roots.append(root)
            root = n
    return T, root, roots


def saturation(bound, cap):
    T = old.add((2, 4, 1), 0)
    seen = {z: (v, ()) for v, z in enumerate(old.profiles(T))
            if cap is None or 3-T[v].bit_count() <= cap}
    todo = deque(seen)
    hits = {3: {'roots': [], 'rows': (2, 4, 1)}}
    tested = 0
    while todo:
        z = todo.popleft()
        source, word = seen[z]
        if z[0] not in hits:
            U, root, roots = decode(source, word)
            need(old.profiles(U)[root] == z, 'named program full-state replay')
            hits[z[0]] = {'seed_root': source, 'program': word, 'roots': roots}
            old.table.cache_clear()
            old.profiles.cache_clear()
        options = [('R', old.keep(z))]
        A = 0
        while cap is None or A < cap:
            zz = block(z, A, 0)
            if zz[0] > bound:
                break
            options.append(('J'+str(A), zz))
            A += 1
        for letter, zz in options:
            tested += 1
            need(zz[0] > z[0], 'strict source growth')
            if zz[0] <= bound and zz not in seen:
                seen[zz] = (source, word+(letter,))
                todo.append(zz)
    return {'cutoff': bound, 'degree_cap': cap, 'states': len(seen),
            'tested_transitions': tested,
            'missing_allowed': [h for h in range(3, bound+1, 2) if h not in hits and h not in (7, 21)],
            'witnesses': {str(h): hits[h] for h in sorted(hits)}}


def main():
    matrices = {}
    for letter, operation in [('R', old.keep)] + [('J'+str(A), lambda z, A=A: block(z, A, 0)) for A in range(3)]:
        M = [positive(operation(profile(tuple(int(i == j) for i in range(4))))) for j in range(4)]
        need(all(x >= 0 for row in M for x in row), 'positive four-coordinate matrix')
        matrices[letter] = M
    z = (3, 1, 1, 3)
    prefixes = [z, old.keep(z), old.new(z), block(z, 1, 0)]
    hankel = [[v[0], old.keep(v)[0], old.keep(old.keep(v))[0], old.keep(old.new(v))[0]] for v in prefixes]
    need(det(hankel) == -32, 'exact four-dimensional lower bound')
    result = {'scope': 'arbitrary actual rooted sources; fresh-block grammar, not arbitrary old-root switching',
              'profile_checks': checks(), 'positive_coordinates': ['h-s', 'q-s', 's', 'F-h-q+s'],
              'matrices_row_vectors': matrices,
              'hankel': {'prefixes': ['empty', 'R', 'N', 'J1'], 'suffixes': ['empty', 'R', 'RR', 'NR'],
                         'matrix': hankel, 'determinant': -32},
              'coverage': [saturation(199, 2), saturation(199, 3), saturation(199, None)]}
    for h, word in [(39, (0, 3, 0, 0)), (67, (0, 1, 3, 4))]:
        T = old.history(word)
        need(old.table(T)[0][-1] == h, 'actual excluded-grammar target witness')
        result.setdefault('older_port_controls', []).append({'H': h, 'roots': word, 'rows': T})
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
