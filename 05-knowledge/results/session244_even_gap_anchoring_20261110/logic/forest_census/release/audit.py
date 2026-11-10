#!/usr/bin/env python3
"""Independent Python census of labelled graphs through six and sampled witnesses.

Run: python3 -I -B audit.py [CSV_DIRECTORY]
This verifier never edits the CSV directory. All checks are explicit exceptions.
"""
import csv
import hashlib
import itertools
import json
import math
import sys
from functools import lru_cache
from pathlib import Path

HEADER = ['H', 'chronological_code', 'min_forest_depth', 'min_depth_code']


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def decode(n, code):
    need(0 <= code < math.factorial(n), 'code range')
    choices = [0] * n
    for size in range(n, 1, -1):
        choices[size - 1] = code % size
        code //= size
    need(code == 0, 'code remainder')
    a, depths = [0] * n, [0] * n
    for v in range(1, n):
        choice = choices[v]
        depths[v] = depths[choice - 1] + 1 if choice else 0
        for u in range(v):
            if choice == u + 1:
                a[v] |= 1 << u
            else:
                a[u] |= 1 << v
    return tuple(a), max(depths)


def hp(a):
    @lru_cache(None)
    def f(mask, last):
        if mask == 1 << last:
            return 1
        rest = mask ^ (1 << last)
        return sum(f(rest, u) for u in range(len(a))
                   if rest >> u & 1 and a[u] >> last & 1)
    return sum(f((1 << len(a)) - 1, v) for v in range(len(a)))


def degone(a):
    mask = (1 << len(a)) - 1
    while mask:
        candidates = [v for v in range(len(a)) if mask >> v & 1
                      and (a[v] & mask).bit_count() <= 1]
        if not candidates:
            return False
        mask ^= 1 << candidates[0]
    return True


def labelled(n, code):
    a = [0] * n
    for i, j in itertools.combinations(range(n), 2):
        if code & 1:
            a[i] |= 1 << j
        else:
            a[j] |= 1 << i
        code >>= 1
    return tuple(a)


def audit(root):
    paths = sorted(root.glob('n*_witnesses.csv'),
                   key=lambda p: int(p.name.split('_')[0][1:]))
    need(paths, 'missing CSV files')
    entries, report = {}, []
    summaries = {r['n']: r for r in map(json.loads,
                 (root / 'SUMMARY.jsonl').read_text().splitlines())}
    need(set(summaries) == set(range(1, len(paths) + 1)), 'summary orders')
    for path in paths:
        n = int(path.name.split('_')[0][1:])
        need(1 <= n <= 12, 'order bound')
        with path.open(newline='') as stream:
            reader = csv.reader(stream)
            need(next(reader) == HEADER, 'CSV header')
            data = []
            for row in reader:
                need(len(row) == 4 and all(x.isascii() and x.isdigit() for x in row), 'CSV format')
                data.append(tuple(map(int, row)))
        need(data, 'empty CSV')
        vals = [r[0] for r in data]
        need(vals == sorted(set(vals)), 'duplicate or unordered values')
        need(all(h % 2 and h not in (7, 21) and h <= math.factorial(n)
                 and 0 <= code < math.factorial(n)
                 and 0 <= mincode < math.factorial(n) and 0 <= depth < n
                 for h, code, depth, mincode in data), 'CSV value/witness bounds')
        rows = {r[0]: r for r in data}
        entries[n] = rows
        if n <= 7:
            selected = vals
        else:
            indexes = {0, len(vals) - 1, *range(min(24, len(vals)))}
            indexes.update(i * (len(vals) - 1) // 100 for i in range(101))
            selected = [vals[i] for i in sorted(indexes)]
        for h in selected:
            _, code, recorded_depth, mincode = rows[h]
            a, depth = decode(n, code)
            need(degone(a) and hp(a) == h,
                 f'first witness n={n} H={h} code={code}')
            b, depth_min = decode(n, mincode)
            need(degone(b) and hp(b) == h and depth_min == recorded_depth <= depth,
                 f'minimum-depth witness n={n} H={h} code={mincode}')
        missing = [h for h in range(1, max(vals) + 3, 2)
                   if h not in rows and h not in (7, 21)]
        first_miss = missing[0] if missing else max(vals) + 2
        while first_miss in (7, 21):
            first_miss += 2
        depth_counts = {str(d): sum(r[2] == d for r in data) for d in range(n)}
        first_deep = next((r[0] for r in data if r[2] > 2), 0)
        actual = {'n': n, 'histories': math.factorial(n), 'distinct_values': len(data),
                  'max_H': max(vals), 'first_allowed_miss': first_miss,
                  'depth2_values': sum(r[2] <= 2 for r in data),
                  'depth3_values': sum(r[2] <= 3 for r in data),
                  'first_needing_depth_gt2': first_deep}
        need(actual == summaries[n], f'summary mismatch n={n}')
        report.append({**actual, 'first_20_allowed_misses': missing[:20],
                       'rows_sampled': len(selected), 'codes_recounted': 2 * len(selected),
                       'csv_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                       'depth_counts': depth_counts})
    need(set(entries) == set(range(1, max(entries) + 1)), 'nonconsecutive orders')
    labelled_reports = []
    for n in range(1, min(6, max(entries)) + 1):
        values, count = set(), 0
        for code in range(1 << (n * (n - 1) // 2)):
            a = labelled(n, code)
            if degone(a):
                values.add(hp(a))
                count += 1
        need(values == set(entries[n]), f'complete labelled class spectrum order {n}')
        labelled_reports.append({'n': n, 'labelled_deg_one': count,
                                 'spectrum_size': len(values)})
    for n in range(2, max(entries) + 1):
        need(set(entries[n - 1]) <= set(entries[n]), 'sink-padding spectrum nesting')
    return {'status': 'PASS', 'order_results': report,
            'independent_labelled_census': labelled_reports,
            'scope': 'Complete labelled graph-class spectra independently checked through order 6; deterministic literal witness samples above order 6. Summary minimal depths rely on the exhaustive history census, not these positive recounts.'}


if __name__ == '__main__':
    need(len(sys.argv) <= 2, 'usage: audit.py [CSV_DIRECTORY]')
    root = Path(sys.argv[1]) if len(sys.argv) == 2 else Path(__file__).resolve().parent
    print(json.dumps(audit(root), indent=2))
