"""Exact unique-start star faces and full newborn target response.

The one-ear polynomial is inherited; this audits the ready-face restriction
and the additional ending-cover field needed for the second birth.
"""
from pathlib import Path
from itertools import combinations
import argparse
import hashlib
import importlib.util
import json
import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('spectrum_for_star', HERE.parent / 'containment/count_spectrum.py')
counter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(counter)
need = counter.need


def induced(rows, keep):
    return tuple(sum(1 << j for j, v in enumerate(keep) if rows[u] >> v & 1) for u in keep)


def decode(n, code):
    rows = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        a, b = (i, j) if code >> k & 1 else (j, i)
        rows[a] |= 1 << b
    return tuple(rows)


def kernel(rows):
    n = len(rows)
    need(n > 0, 'nonempty deletion graph required')
    full = (1 << n) - 1
    end = counter.backward(np.array(rows, dtype=np.int64))
    start = counter.backward(np.array([full ^ (1 << i) ^ row for i, row in enumerate(rows)], dtype=np.int64))
    h = end.sum(axis=1)
    h[0] = 1
    F = [sum(int(end[mask, i]) * int(h[full ^ mask]) for mask in range(1, full + 1)) for i in range(n)]
    Q = [[sum(int(end[mask, i]) * int(start[full ^ mask, j]) for mask in range(1, full))
          for j in range(n)] for i in range(n)]
    return dict(H=int(h[full]), E=list(map(int, end[full])), S=list(map(int, start[full])),
                deletion_H=[int(h[full ^ (1 << i)]) for i in range(n)], F=F, Q=Q)


def evaluate(k, mask):
    n = len(k['S'])
    out = [j for j in range(n) if mask >> j & 1]
    inc = [i for i in range(n) if not (mask >> i & 1)]
    s = sum(k['S'][j] for j in out)
    h = sum(k['E'][i] for i in inc) + s + sum(k['Q'][i][j] for i in inc for j in out)
    F = k['H'] + sum(k['F'][i] for i in inc)
    target = k['H'] + sum(k['E'][i] + k['F'][i] for i in inc) + sum(k['Q'][i][j] for i in inc for j in out)
    need(target == h + F - s, 'newborn cancellation')
    return h, s, k['H'], F, target


def star(rows, mask):
    n = len(rows)
    return tuple(row | (0 if mask >> i & 1 else (1 << n)) for i, row in enumerate(rows)) + (mask,)


def state(rows, root):
    full = (1 << len(rows)) - 1
    e = counter.backward(np.array(rows, dtype=np.int64))
    h = e.sum(axis=1)
    h[0] = 1
    s = counter.backward(np.array([full ^ (1 << i) ^ row for i, row in enumerate(rows)], dtype=np.int64))
    F = sum(int(e[mask, root]) * int(h[full ^ mask]) for mask in range(1, full + 1))
    return [int(h[full]), int(s[full, root]), int(h[full ^ (1 << root)]), F]


def faces(k):
    zero = [v for v, value in enumerate(k['S']) if value == 0]
    for w, value in enumerate(k['S']):
        if value != 1:
            continue
        yield w, zero, [(1 << w) | sum((1 << v) for i, v in enumerate(zero) if bits >> i & 1)
                        for bits in range(1 << len(zero))]


def all_port_repairs(rows, target):
    records = []
    total = 0
    for v in range(len(rows)):
        keep = [i for i in range(len(rows)) if i != v]
        U = induced(rows, keep)
        k = kernel(U)
        for w, zero, masks in faces(k):
            values = [evaluate(k, mask)[-1] for mask in masks]
            minimum = min(values)
            minimizing = [mask for mask, value in zip(masks, values) if value == minimum]
            need(all(a & b in (a, b) for a in minimizing for b in minimizing), 'minimizers must be a chain')
            for mask, value in zip(masks, values):
                total += 1
                P = star(U, mask)
                need(state(P, len(U))[1] == 1, 'literal unique-start star')
                child = counter.birth(P, len(U))
                child_state = state(child, len(P))
                need(child_state[0] == value and child_state[3] - value == 1, 'literal full ready child')
                if value == target:
                    changed = [keep[i] for i in range(len(U)) if ((rows[v] >> keep[i]) & 1) != ((mask >> i) & 1)]
                    records.append(dict(selected_old_port=v, changed_old_neighbours=changed,
                                        parent_rows=P, parent_state=state(P, len(U)),
                                        child_rows=child, child_state=child_state,
                                        deletion_rows=U, kernel=k, unique_entry=keep[w]))
    return dict(source_rows=rows, target=target, total_unique_start_stars=total, hits=records)


def run(output):
    parents = stars = chains = 0
    for n in range(1, 5):
        for code in range(1 << (n * (n - 1) // 2)):
            U = decode(n, code)
            k = kernel(U)
            parents += 1
            need(all(k['Q'][i][j] + k['Q'][j][i] >= 2 for i in range(n) for j in range(i)), 'strict curvature')
            feasible = {mask for _, _, masks in faces(k) for mask in masks}
            direct_feasible = set()
            for mask in range(1 << n):
                stars += 1
                P = star(U, mask)
                expected = evaluate(k, mask)
                need(state(P, n) == list(expected[:4]), 'literal parent state')
                child = counter.birth(P, n)
                c = state(child, n + 1)
                need(c[0] == expected[-1] and c[3] - c[0] == expected[1], 'literal newborn count and defect')
                if expected[1] == 1:
                    direct_feasible.add(mask)
            need(feasible == direct_feasible, 'complete unique-start faces')
            for w, zero, masks in faces(k):
                values = [evaluate(k, mask)[-1] for mask in masks]
                minimizers = [mask for mask, value in zip(masks, values) if value == min(values)]
                need(len(minimizers) <= len(zero) + 1, 'minimum-chain size')
                need(all(a & b in (a, b) for a in minimizers for b in minimizers), 'minimum chain')
                largest = masks[-1]
                need(minimizers == [largest], 'unique minimum at all zero-start vertices outward')
                for mask, value in zip(masks, values):
                    incoming_tail = sum(not (mask >> i & 1) for i in zero)
                    need(value - min(values) >= 2 * incoming_tail * (k['H'] + k['deletion_H'][w]),
                         'two triangle families lower-bound exact target cost')
                chains += 1
    A2 = (518,388,456,1011,487,455,387,257,1,502)
    P153 = (130,388,457,1011,487,455,387,256,1,503)
    T1691 = counter.birth(P153, 0)
    first = all_port_repairs(A2, 1691)
    second = all_port_repairs(T1691, 3311)
    need(first['total_unique_start_stars'] == 12 and len(first['hits']) == 1, 'complete A2 star repair')
    need(first['hits'][0]['selected_old_port'] == 0 and first['hits'][0]['changed_old_neighbours'] == [2,7,9], 'selective star witness')
    need(second['total_unique_start_stars'] == 6 and not second['hits'], 'first child star continuation hostile')
    root = HERE.parents[3]
    source_paths = ['05-knowledge/results/session124_tournament_cut_response_20261031/README.md',
                    '05-knowledge/results/session136_full_tournament_spectrum_20261102/ear_geometry/positive_cut_geometry.md']
    pins = {p:hashlib.sha256((root / p).read_bytes()).hexdigest() for p in source_paths}
    report = dict(status='FINITE-EXACT controls for the proved ready-star kernel',
                  universe='All labelled nonempty deletion tournaments through order4, every star, every newborn; all unique-start star choices at every port of A2 and its selected1691 child',
                  inherited_filters='No strongness filter; weak candidates retained',
                  deletion_tournaments=parents, stars=stars, minimum_faces=chains,
                  source_pins=pins, first=first, second=second)
    (HERE / output).write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('status','deletion_tournaments','stars','minimum_faces')}))
    print('A2 repairs', len(first['hits']), 'next repairs', len(second['hits']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='receipt.json')
    run(parser.parse_args().output)
