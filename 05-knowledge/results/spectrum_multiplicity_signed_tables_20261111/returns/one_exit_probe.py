"""Exact fixed-gadget gate; inherited ordinary tournament bank through order seven."""
from pathlib import Path
import hashlib
import json
import subprocess

from native_return import construct, converse, ends, strong

HERE = Path(__file__).resolve().parent
REV = '63b86513f8be'
BANK = '05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/extended_census.json'


def need(condition, label):
    if not condition:
        raise RuntimeError(label)


def moments(G):
    n = len(G)
    full = (1 << n) - 1
    e, s = ends(G), ends(converse(G))
    h = [sum(row) for row in e]
    h[0] = 1
    u = [sum(s[a][j] * h[full ^ a] for a in range(1, full + 1)) for j in range(n)]
    k = [[sum(e[a][i] * s[full ^ a][j] for a in range(1, full))
          for j in range(n)] for i in range(n)]
    return h[full], e[full], u, k


def selected(data, A):
    H, e, u, k = data
    ids = [i for i in range(len(e)) if A >> i & 1]
    return H, sum(u[i] for i in ids), sum(e[i] for i in ids), sum(k[i][j] for i in ids for j in ids)


def attach(X, x, G, A):
    n = len(X)
    T = list(X) + [row << n for row in G]
    for j in range(len(G)):
        for i in range(n):
            if i == x and A >> j & 1:
                T[i] |= 1 << (n + j)
            else:
                T[n + j] |= 1 << i
    return T


def source_moments(X, x):
    full = (1 << len(X)) - 1
    e, s = ends(X), ends(converse(X))
    h = [sum(row) for row in e]
    h[0] = 1
    F = sum(e[a][x] * h[full ^ a] for a in range(1, full + 1))
    return h[full], F, s[full][x], h[full ^ (1 << x)]


def delete(T, r):
    ids = [i for i in range(len(T)) if i != r]
    return [sum(1 << j for j, v in enumerate(ids) if T[u] >> v & 1) for u in ids], ids


def main():
    raw = subprocess.check_output(['git', 'show', REV + ':' + BANK], cwd=HERE)
    inherited = json.loads(raw)
    bank = [(0,), (2, 0)]
    for r in inherited['records']:
        T = tuple(r['adjacency_outmasks'])
        bank.append(T)
        if r['ordinary_multiplicity'] == 2:
            bank.append(tuple(converse(T)))
    counts = {n: sum(len(T) == n for T in bank) for n in range(1, 8)}
    need(counts == {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456}, 'ordinary bank universe')
    compatible = []
    cut_count = 0
    for G in bank:
        data = moments(G)
        for A in range(1 << len(G)):
            cut_count += 1
            H, U, E, W = selected(data, A)
            need(U >= E + W, 'selected capacity inequality')
            if H > E + 1 and (H - E - 1) % 45 == 0:
                t = (H - E - 1) // 45
                if W == 99 * t - E:
                    compatible.append({'lambda': U + 10*t + 1 + E, 't': t,
                                       'G': list(G), 'A': A, 'HUEW': [H, U, E, W]})
    need(sorted(r['lambda'] for r in compatible) == [375, 377], 'complete small compatible-gadget list')

    # Independent literal path counts exercise arbitrary old-root incidences.
    small_controls = 0
    for X in bank[:4]:
        for x in range(len(X)):
            h, F, s, d = source_moments(X, x)
            for G in bank[:4]:
                data = moments(G)
                for A in range(1 << len(G)):
                    H, U, E, W = selected(data, A)
                    actual = sum(ends(attach(X, x, G, A))[-1])
                    need(actual == h*H + F*U - s*E - d*W, 'literal one-exit identity')
                    small_controls += 1

    hostile_G = (6, 4, 0)
    hostile_moments = selected(moments(hostile_G), 3)
    hostile_actual = sum(ends(attach((0,), 0, hostile_G, 3))[-1])
    need(hostile_moments == (1, 6, 0, 2) and hostile_actual == 5, 'singleton hostile')

    X = construct(1)[0]
    x = len(X) - 1
    literal = []
    for row in compatible:
        compressed = attach((0,), 0, row['G'], row['A'])
        cm = moments(compressed)
        compressed_state = [cm[0], cm[1][0], cm[2][0]]
        need(compressed_state == [row['lambda'] - 64, 46, row['lambda'] + 35], 'rooted compression')
        T = attach(X, x, row['G'], row['A'])
        H = sum(ends(T)[-1])
        need(H == row['lambda'] * 666 - 1 and strong(T), 'literal strong Q1 construction')
        root_defects = []
        for r in range(len(T)):
            predecessors = [i for i, out in enumerate(T) if out >> r & 1]
            if len(predecessors) == 1:
                p = predecessors[0]
                P, ids = delete(T, r)
                defect = ends(converse(P))[-1][ids.index(p)]
                root_defects.append({'root': r, 'predecessor': p, 'V_minus_H': defect})
        literal.append({'lambda': row['lambda'], 'H': H, 'strong': True,
                        'compressed_root_H_E_U': compressed_state,
                        'rows': T, 'indegree_one_root_defects': root_defects})

    tt5 = tuple(sum(1 << j for j in range(i + 1, 5)) for i in range(5))
    naive = selected(moments(tt5), 11)
    naive_actual = sum(ends(attach(X, x, tt5, 11))[-1])
    need(naive == (1, 26, 0, 8) and naive_actual == 16101, 'nonuniform transitive-cut hostile')

    receipt = {
        'scope': 'Complete inherited ordinary-isomorphism bank of orders 1..7, every selected subset. Bank completeness inherited; endpoint and two-cover moments computed exactly. Not an all-order gadget exclusion for lambda243.',
        'bank_commit': REV, 'bank_path': BANK, 'bank_sha256': hashlib.sha256(raw).hexdigest(),
        'ordinary_counts_by_order': counts, 'all_selected_cuts': cut_count,
        'positive_t_uniform_gadgets': compatible, 'lambda243_hits': 0,
        'small_literal_formula_controls': small_controls,
        'singleton_hostile': {'HUEW': hostile_moments, 'actual_H': hostile_actual,
                              'incorrect_omitted_penalty_value': 7},
        'literal_Q1_extensions': literal,
        'naive_lambda27_tt5': {'A': 11, 'HUEW': naive, 'actual_H': naive_actual,
                              'desired_H': 17981},
    }
    (HERE / 'one_exit_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': 'pass', 'bank_count': len(bank), 'cuts': cut_count,
                      'positive_lambdas': [r['lambda'] for r in compatible],
                      'small_literal_controls': small_controls}))


if __name__ == '__main__':
    main()
