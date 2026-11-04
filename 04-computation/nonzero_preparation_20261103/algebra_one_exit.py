#!/usr/bin/env python3
"""Exact one-exit identities, named persistence certificates, and clone control.

Standalone: no imported producer or input receipt. The finite universes here
control the analytic all-order statements in algebra/ONE_EXIT_AND_PERSISTENCE.md.
"""
from array import array
from itertools import combinations, permutations
from math import factorial
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


def tournament(n, code):
    rows = [0] * n
    for k, (a, b) in enumerate(combinations(range(n), 2)):
        if code >> k & 1:
            rows[a] |= 1 << b
        else:
            rows[b] |= 1 << a
    return tuple(rows)


def birth(rows, incoming):
    n = len(rows)
    return tuple(r | (((incoming >> i) & 1) << n) for i, r in enumerate(rows)) + (((1 << n) - 1) ^ incoming,)


def strong(rows):
    reach = [r | (1 << i) for i, r in enumerate(rows)]
    for k in range(len(rows)):
        for i in range(len(rows)):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return all(r == (1 << len(rows)) - 1 for r in reach)


def induced(rows, keep):
    return tuple(sum(((rows[a] >> b) & 1) << j for j, b in enumerate(keep)) for a in keep)


def enumerate_paths(rows):
    """Independent literal permutations: all nonempty simple directed paths."""
    n, size = len(rows), 1 << len(rows)
    P = [0] * size
    P[0] = 1
    starts = [[0] * size for _ in rows]
    ends = [[0] * size for _ in rows]
    joint = [[[0] * size for _ in rows] for _ in rows]
    for length in range(1, n + 1):
        for path in permutations(range(n), length):
            if all(rows[a] >> b & 1 for a, b in zip(path, path[1:])):
                mask = sum(1 << v for v in path)
                P[mask] += 1
                starts[path[0]][mask] += 1
                ends[path[-1]][mask] += 1
                joint[path[0]][path[-1]][mask] += 1
    return P, starts, ends, joint


def conv(A, B):
    C = [0] * len(A)
    for mask in range(len(A)):
        sub = mask
        while True:
            C[mask] += A[sub] * B[mask ^ sub]
            if not sub:
                break
            sub = (sub - 1) & mask
    return C


def fullprod(A, B, full):
    return sum(A[s] * B[full ^ s] for s in range(full + 1))


def counts(rows, root=None):
    """Separate forward endpoint DP for larger literal named controls."""
    n = len(rows)
    need(n <= 18 and factorial(n) < 2 ** 64, 'uint64 exact universe')
    size, full = 1 << n, (1 << n) - 1
    dp = array('Q', [0]) * (size * n)
    P = array('Q', [0]) * size
    P[0] = 1
    for v in range(n):
        dp[(1 << v) * n + v] = 1
    for mask in range(1, size):
        vertices, total = mask, 0
        while vertices:
            bit = vertices & -vertices
            vertices -= bit
            v = bit.bit_length() - 1
            value = dp[mask * n + v]
            total += value
            if value:
                choices = rows[v] & (full ^ mask)
                while choices:
                    bit2 = choices & -choices
                    choices -= bit2
                    dp[(mask | bit2) * n + bit2.bit_length() - 1] += value
        P[mask] = total
    F = None if root is None else sum(int(dp[s * n + root]) * int(P[full ^ s]) for s in range(1, size))
    return int(P[full]), F, [int(dp[full * n + v]) for v in range(n)]


def source(t):
    k, D = t - 2, (2, 4, 57, 3, 11, 27)
    suffix = tuple(row | 64 for row in D) + (0,)
    root, ports, suffix_mask = k + 7, 80, 127 << k
    rows = [sum(1 << j for j in range(i + 1, k)) | suffix_mask for i in range(k)]
    rows += [(row << k) | ((1 << root) if ports >> i & 1 else 0) for i, row in enumerate(suffix)]
    rows += [((1 << k) - 1) | (suffix_mask ^ (ports << k))]
    return tuple(rows), root, k + 6, (k, k + 1, k + 2)


def c2(rows, aux):
    n = len(rows)
    X = birth(rows, ((1 << n) - 1) ^ (1 << aux))
    T = birth(X, (1 << n) - 1)
    need(strong(X) and strong(T), 'literal C2 strong at both births')
    need(induced(T, tuple(range(n))) == rows, 'C2 retains all old arcs')
    return X, T, n, n + 1


def triangle(rows, tri):
    a, b, c = tri
    return bool(rows[a] >> b & 1 and rows[b] >> c & 1 and rows[c] >> a & 1)


def certify(rows, v, tri, witnesses):
    """Dual of inherited I2: actual old names, not scalar readiness."""
    out = tuple(i for i in range(len(rows)) if rows[v] >> i & 1)
    if len(out) >= 3:
        return {'v': v, 'kind': 1, 'exits': out}
    need(len(out) in (1, 2), 'certificate positive outdegree')
    need(not (set(tri) & ({v} | set(out))) and triangle(rows, tri), 'literal residual triangle')
    if len(out) == 2:
        return {'v': v, 'kind': 2, 'exits': out, 'triangle': tri}
    keep, expected = witnesses[v]
    need(v not in keep and out[0] in keep, 'endpoint witness names')
    H, F, ends = counts(induced(rows, keep))
    actual = ends[keep.index(out[0])]
    need(actual == expected and actual >= 3, 'endpoint witness exact count')
    return {'v': v, 'kind': 3, 'exit': out[0], 'triangle': tri,
            'witness': keep, 'end_count': actual}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, default=Path(__file__).with_suffix('.json'))
    args = ap.parse_args()
    census, examples = [], {}
    for n in range(2, 6):
        cases = 0
        for code in range(1 << (n * (n - 1) // 2)):
            rows = tournament(n, code)
            P, A, E, Q = enumerate_paths(rows)
            full = (1 << n) - 1
            P2 = conv(P, P)
            for a in range(n):
                if rows[a].bit_count() != 1:
                    continue
                b = rows[a].bit_length() - 1
                U = full ^ (1 << a)
                S = fullprod(P, A[a], full)
                sigma = S - P[full]
                need(sigma == E[b][U], 'one-exit sigma')
                for w in range(n):
                    if w == a:
                        continue
                    cases += 1
                    F = fullprod(E[w], P, full)
                    G = fullprod(conv(E[w], P), A[a], full) + fullprod(P2, Q[a][w], full)
                    joint = conv(E[w], E[b])[U]
                    need(G - F == joint, 'one-exit joint tau')
                    need(joint == 0 if b == w else joint >= sigma, 'named-exit cut injection')
                    need(not (sigma > 0 and G == F) or b == w, 'positive preparation exit forced')
                    kind = ('root_exit_' + str(sigma)) if b == w and sigma in (1, 2) and strong(rows) else None
                    if b != w and sigma > 0 and strong(rows):
                        kind = 'wrong_exit_positive_joint'
                    if kind and kind not in examples:
                        examples[kind] = {'rows': rows, 'root': w, 'aux': a, 'exit': b, 'H': P[full], 'F': F, 'sigma': sigma, 'tau': G - F}
        need(cases == n * (n - 1) ** 2 * 2 ** ((n - 1) * (n - 2) // 2), 'complete one-exit marked universe')
        census.append({'order': n, 'graphs': 2 ** (n * (n - 1) // 2), 'marked_cases': cases})
    need(all(k in examples for k in ('root_exit_1', 'root_exit_2', 'wrong_exit_positive_joint')), 'positive and hostile literal controls')

    family, clone = [], []
    for t in (3, 4, 5):
        B, root, u, tri = source(t)
        h, f, ends = counts(B, root)
        need((h, f, ends[u]) == (63 * 2 ** t - 29, 63 * 2 ** t + 17, 28 * 2 ** (t - 2) + 15), 'source full endpoint formula')
        X, T, aux, y = c2(B, u)
        witnesses = {aux: (tuple(range(len(B))), ends[u]), y: (tuple(range(len(X))), h - ends[u])}
        stages = []
        oldH = h
        for j in range(3):
            H, F, end = counts(T, root)
            need((H, F) == (3 ** j * (189 * 2 ** t - 83), 3 ** j * (189 * 2 ** t + 51)), 'all-clock zero formulas')
            need(end[aux] == oldH, 'transported auxiliary end count')
            certs = [certify(T, v, tri, witnesses) for v in range(len(T))]
            stages.append({'level': j, 'rows': T, 'H': H, 'F': F, 'aux': aux, 'certificates': certs})
            if j < 2:
                old, old_aux = T, aux
                first, T, aux, y = c2(old, old_aux)
                witnesses[aux] = (tuple(range(len(old))), end[old_aux])
                witnesses[y] = (tuple(range(len(first))), H - end[old_aux])
                oldH = H
        family.append({'clock': t, 'root': root, 'source': B, 'stages': stages})
        # The clone blocks have the exact named order y_1<...<y_a.
        T, a = c2(B, u)[1], 1
        x = len(B)
        seq = []
        for a in range(1, 6):
            H, F, _ = counts(T, root)
            need((H, F) == (h + 2 ** a * (h + 2), (1 + 2 ** a) * f), 'all-clock clone formula')
            need(strong(T), 'clone family strong')
            clone_witnesses = {x: (tuple(range(len(B))), ends[u]),
                               len(T) - 1: (tuple(range(len(B) + 1)), h - ends[u])}
            certs = [certify(T, v, tri, clone_witnesses) for v in range(len(T))]
            seq.append({'a': a, 'rows': T, 'H': H, 'F': F, 'certificates': certs})
            if a < 5:
                old = T
                T = birth(T, ((1 << len(T)) - 1) ^ (1 << x))
                need(induced(T, tuple(range(len(old)))) == old and strong(T), 'clone birth retains old names and is strong')
        need((seq[2]['H'], seq[2]['F']) == (3 * seq[0]['H'] + 4, 3 * seq[0]['F']), 'literal non-C2 clone digit four')
        clone.append({'clock': t, 'sequence': seq})

    # All 2^11 one-birth cuts of the first prepared output: certificates
    # stay literal, even if the target is weak; no Hamilton-response inference.
    initial = tuple(family[0]['stages'][0]['rows'])
    checked = 0
    for cut in range(1 << len(initial)):
        U = birth(initial, cut)
        need(induced(U, tuple(range(len(initial)))) == initial, 'complete extension cut preserves certificate graph')
        for cert in family[0]['stages'][0]['certificates']:
            v = cert['v']
            exits = {i for i in range(len(U)) if U[v] >> i & 1}
            if len(exits) >= 3:
                need(True, 'extended exit-degree certificate')
            elif len(exits) == 2:
                need(not (set(cert['triangle']) & (exits | {v})) and triangle(U, cert['triangle']),
                     'extended two-exit residual triangle')
            else:
                need(cert['kind'] == 3 and exits == {cert['exit']} and cert['end_count'] >= 3,
                     'extended unique-exit witness certificate')
                need(induced(U, tuple(cert['witness'])) == induced(initial, tuple(cert['witness'])),
                     'extended endpoint witness unchanged')
        checked += 1
    result = {'status': 'FINITE-EXACT controls; all-order proof is separate', 'gates': GATES,
              'one_exit_census': census, 'examples': examples, 'family': family,
              'clone_controls': clone, 'all_one_birth_cuts': checked,
              'scope': 'Old saved vertices permanently nonready; fresh exterior auxiliaries and old-edge changes are not excluded.'}
    text = json.dumps(result, indent=2, sort_keys=True) + '\n'
    args.output.write_text(text)
    print(json.dumps({'gates': GATES, 'cases': sum(c['marked_cases'] for c in census),
                      'bytes': len(text.encode()), 'sha256': hashlib.sha256(text.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
