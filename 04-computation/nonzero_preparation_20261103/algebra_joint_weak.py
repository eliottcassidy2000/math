#!/usr/bin/env python3
"""Independent all-tournament guard and literal nonreachability injection.

No imports of producer code. All path objects are enumerated by permutations.
"""
from itertools import combinations, permutations
from pathlib import Path
import argparse
import hashlib
import json

GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def graph(n, code):
    rows = [0] * n
    for i, (u, v) in enumerate(combinations(range(n), 2)):
        if code >> i & 1:
            rows[u] |= 1 << v
        else:
            rows[v] |= 1 << u
    return tuple(rows)


def multiply(A, B):
    C = [0] * len(A)
    for mask in range(len(A)):
        s = mask
        while True:
            C[mask] += A[s] * B[mask ^ s]
            if not s:
                break
            s = (s - 1) & mask
    return C


def full(A, B):
    mask = len(A) - 1
    return sum(A[s] * B[mask ^ s] for s in range(mask + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, default=Path(__file__).with_suffix('.json'))
    args = ap.parse_args()
    cases = nonreach_cases = cover_images = weak_positive = 0
    census = []
    for n in range(1, 6):
        at_n = 0
        for code in range(1 << (n * (n - 1) // 2)):
            T = graph(n, code)
            size, mask = 1 << n, (1 << n) - 1
            P = [0] * size
            P[0] = 1
            A = [[0] * size for _ in T]
            E = [[0] * size for _ in T]
            Q = [[[0] * size for _ in T] for _ in T]
            literals = [[] for _ in range(size)]
            for k in range(1, n + 1):
                for path in permutations(range(n), k):
                    if all(T[u] >> v & 1 for u, v in zip(path, path[1:])):
                        s = sum(1 << v for v in path)
                        literals[s].append(path)
                        P[s] += 1
                        A[path[0]][s] += 1
                        E[path[-1]][s] += 1
                        Q[path[0]][path[-1]][s] += 1
            reach = [r | (1 << i) for i, r in enumerate(T)]
            for k in range(n):
                for i in range(n):
                    if reach[i] >> k & 1:
                        reach[i] |= reach[k]
            strong = n >= 3 and all(r == mask for r in reach)
            P2 = multiply(P, P)
            for a in range(n):
                B = {b for b in range(n) if T[a] >> b & 1}
                EB = [sum(E[b][s] for b in B) for s in range(size)]
                sigma = full(P, A[a]) - P[mask]
                need(sigma == full(EB, A[a]), 'exact scalar failed-join identity')
                for w in range(n):
                    F = full(E[w], P)
                    G = full(multiply(E[w], P), A[a]) + full(P2, Q[a][w])
                    tau = G - F
                    right = a != w and B == {w} and E[w][mask ^ (1 << a)] > 0
                    need((sigma > 0 and tau == 0) == right, 'all-tournament positive joint iff')
                    if not (reach[a] >> w & 1):
                        nonreach_cases += 1
                        images = set()
                        for s in range(size):
                            for failed_prefix in literals[s]:
                                if failed_prefix[-1] not in B:
                                    continue
                                for Apath in literals[mask ^ s]:
                                    if Apath[0] != a:
                                        continue
                                    need(w in failed_prefix and w != failed_prefix[-1], 'unreachable root belongs inside failed prefix')
                                    cut = failed_prefix.index(w) + 1
                                    image = (failed_prefix[:cut], failed_prefix[cut:], Apath)
                                    need(image[0][-1] == w and image[1][-1] in B and all(image), 'literal L-surplus target')
                                    need(image not in images, 'injective typed cut')
                                    images.add(image)
                        need(len(images) == sigma and tau >= sigma, 'nonreachability injection full counts')
                        cover_images += len(images)
                    weak_positive += (not strong and sigma > 0 and tau == 0)
                    cases += 1
                    at_n += 1
        census.append({'n': n, 'marked_cases': at_n})
    need(cases == 26705, 'complete declared universe')
    result = {'scope': 'all labelled tournaments through5, all ordered marks; exact nonreachability injection and full positive guard',
              'gates': GATES, 'marked_cases': cases, 'nonreachable_marked_cases': nonreach_cases,
              'literal_injection_images': cover_images, 'weak_positive_guards': weak_positive, 'census': census}
    raw = (json.dumps(result, indent=2, sort_keys=True)+'\n').encode()
    args.output.write_bytes(raw)
    print(json.dumps({'gates':GATES,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))


if __name__ == '__main__':
    main()
