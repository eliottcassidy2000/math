#!/usr/bin/env python3
"""Exact coprime branch compatibility; inherited complete strong H<375 bank.

All assertions are explicit runtime checks and survive python -O. No claim of
re-enumerating the inherited bank is made. Canonicalisation uses every actual
Hamilton word, so neither an external graph package nor an orbit heuristic is
needed for the selected small-count sources.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BANK = ROOT / '04-computation/experiments/h_spectrum_multiplicities_20261021/audit_S50/inventory_375.json'

def require(q, msg):
    if not q:
        raise RuntimeError(msg)

def decode(n, code):
    rows = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code >> k & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
            k += 1
    return tuple(rows)

def induced(rows, vertices):
    return tuple(sum(((rows[v] >> w) & 1) << j for j, w in enumerate(vertices)) for v in vertices)

def hdp(rows):
    n = len(rows)
    if not n:
        return 1
    d = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        d[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for i, x in enumerate(d[mask]):
            if x:
                out = rows[i] & ~mask
                while out:
                    bit = out & -out
                    out -= bit
                    d[mask | bit][bit.bit_length() - 1] += x
    return sum(d[-1])

def words(rows):
    n = len(rows)
    def walk(word, unused):
        if not unused:
            yield word
            return
        out = unused if not word else unused & rows[word[-1]]
        while out:
            bit = out & -out
            out -= bit
            yield from walk(word + (bit.bit_length() - 1,), unused ^ bit)
    yield from walk((), (1 << n) - 1)

def code_at(rows, p):
    return sum(((rows[p[i]] >> p[j]) & 1) << k for k, (i, j) in enumerate((i, j) for i in range(len(p)) for j in range(i + 1, len(p))))

@lru_cache(None)
def canonical(rows):
    return min((code_at(rows, p), p) for p in words(rows))

def components(rows):
    n = len(rows)
    reach = [r | (1 << i) for i, r in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    left = set(range(n)); ans = []
    while left:
        v = min(left)
        part = tuple(w for w in sorted(left) if reach[v] >> w & 1 and reach[w] >> v & 1)
        ans.append(part)
        left.difference_update(part)
    parts = tuple(ans)
    ans.sort(key=lambda p: -sum(bool(rows[p[0]] & (1 << q[0])) for q in parts if p != q))
    require(all(rows[a[0]] >> b[0] & 1 for i, a in enumerate(ans) for b in ans[i + 1:]), 'SCC ordering')
    return ans

def normalize(rows):
    cc = components(rows)
    active = [i for i, c in enumerate(cc) if len(c) > 1]
    require(bool(active), 'nontransitive source')
    lo, hi = active[0], active[-1]
    prefix = tuple(c[0] for c in cc[:lo])
    suffix = tuple(c[0] for c in cc[hi + 1:])
    middle = tuple(v for c in cc[lo:hi + 1] for v in c)
    raw = induced(rows, middle)
    code, order = canonical(raw)
    core = induced(raw, order)
    source_order = prefix + tuple(middle[i] for i in order) + suffix
    return (len(core), code), len(prefix), len(suffix), core, source_order

def join(left, right):
    a, b = len(left), len(right)
    return tuple(r | (((1 << b) - 1) << a) for r in left) + tuple(r << a for r in right)

def transitive(n):
    return tuple(((1 << n) - 1) ^ ((1 << (i + 1)) - 1) for i in range(n))

def golden_target(word):
    """Literal inherited source join, nonspine complement, label reversal."""
    source = ()
    for letter in reversed(word):
        n = 1 if letter == 's' else 3 if letter == 'c' else 4
        mask = int(letter[1:]) if letter.startswith('q') else 0
        rows = [0] * n; k = 0
        for i in range(n):
            for j in range(i + 1, n):
                forward = j == i + 1 or mask >> k & 1
                if j != i + 1:
                    k += 1
                a, b = (i, j) if forward else (j, i)
                rows[a] |= 1 << b
        source = join(source, tuple(rows))
    rows = list(source)
    for i in range(len(rows)):
        for j in range(i + 2, len(rows)):
            rows[i] ^= 1 << j
            rows[j] ^= 1 << i
    return induced(tuple(rows), tuple(reversed(range(len(rows)))))

def birth(rows, cut):
    n = len(rows)
    return tuple(r | (0 if cut >> i & 1 else 1 << n) for i, r in enumerate(rows)) + (cut,)

def attachment_profile(rows):
    """Independent split-at-root compiler, returning all 2^n attachment counts."""
    n = len(rows); full = (1 << n) - 1
    ends = [[0] * n for _ in range(1 << n)]
    starts = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        ends[1 << i][i] = starts[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for i in range(n):
            if not (mask >> i & 1) or mask == 1 << i:
                continue
            prior = mask ^ (1 << i)
            ends[mask][i] = sum(ends[prior][j] for j in range(n) if prior >> j & 1 and rows[j] >> i & 1)
            starts[mask][i] = sum(starts[prior][j] for j in range(n) if prior >> j & 1 and rows[i] >> j & 1)
    covers = [[0] * n for _ in range(n)]
    for mask in range(1, full):
        other = full ^ mask
        for i, a in enumerate(ends[mask]):
            if a:
                for j, b in enumerate(starts[other]):
                    covers[i][j] += a * b
    values = []
    for cut in range(1 << n):
        incoming = [i for i in range(n) if not (cut >> i & 1)]
        outgoing = [j for j in range(n) if cut >> j & 1]
        values.append(sum(ends[full][i] for i in incoming) + sum(starts[full][j] for j in outgoing) + sum(covers[i][j] for i in incoming for j in outgoing))
    return values

def lift(core, record, a, b):
    old_a, old_b = record['padding']
    require(a >= old_a and b >= old_b, 'padding lift order')
    cut = record['cut'] << (a - old_a)
    # Extra prefix vertices beat the new root; the root beats extra suffix.
    cut |= ((1 << (b - old_b)) - 1) << (a + len(core) + old_b)
    source = join(join(transitive(a), core), transitive(b))
    return source, cut

def main():
    data = json.loads(BANK.read_text())
    target_sets = {55: (167, 169), 105: (317, 319)}
    source_values = {3, 5, 11, 35, 55, 105}
    selected_targets = {v for pair in target_sets.values() for v in pair}
    selected_targets |= {107, 109, 111}  # Permanent-gap and golden-bank controls.
    inventories = defaultdict(Counter)
    by_pair = defaultdict(lambda: defaultdict(list))
    all_strong_sources = defaultdict(dict)
    checked_classes = checked_deletions = independent_path_checks = 0
    for ns, entries in data['orders'].items():
        n = int(ns)
        for raw_code, H, self_converse, aut in entries:
            if H not in source_values | selected_targets:
                continue
            rows = decode(n, int(raw_code))
            require(len(components(rows)) == 1, ('strong', n, raw_code))
            require(hdp(rows) == H, ('bank decode/count', n, raw_code))
            require(sum(1 for _ in words(rows)) == H, ('DFS count', n, raw_code))
            independent_path_checks += 1
            if H in source_values:
                inventories[H][n] += 1
                if H in target_sets:
                    key, a, b, core, order = normalize(rows)
                    require((a, b) == (0, 0), 'strong source has no padding')
                    all_strong_sources[H][key] = core
            if H not in selected_targets:
                continue
            checked_classes += 1
            for root in range(n):
                old_labels = tuple(i for i in range(n) if i != root)
                source = induced(rows, old_labels)
                h = hdp(source)
                checked_deletions += 1
                if (h, H) not in {(55, 167), (55, 169), (105, 317), (105, 319), (35, 107), (35, 109), (55, 111)}:
                    continue
                require(sum(1 for _ in words(source)) == h, 'source independent count')
                independent_path_checks += 1
                key, a, b, core, order = normalize(source)
                actual_order = tuple(old_labels[i] for i in order)
                cut = sum(((rows[root] >> v) & 1) << j for j, v in enumerate(actual_order))
                normalized_source = join(join(transitive(a), core), transitive(b))
                require(induced(source, order) == normalized_source, 'source normalization')
                require(hdp(birth(normalized_source, cut)) == H, 'normalized cut transport')
                by_pair[h, H][key].append({'padding': [a, b], 'cut': cut, 'child_bank': [n, raw_code, root]})

    report = {
        'status': 'FINITE-EXACT relative to inherited audited complete strong H<375 inventory; all-order consequence uses coprime-core and padding lemmas',
        'bank': str(BANK.relative_to(ROOT)),
        'bank_sha256': sha256(BANK.read_bytes()).hexdigest(),
        'checked_child_classes': checked_classes,
        'checked_vertex_deletions': checked_deletions,
        'independent_path_checks': independent_path_checks,
        'strong_source_counts_by_order': {str(h): dict(sorted(c.items())) for h, c in sorted(inventories.items())},
        'families': {},
    }
    for h, targets in target_sets.items():
        maps = [by_pair[h, N] for N in targets]
        common = set(maps[0]) & set(maps[1])
        kernels = set(maps[0]) | set(maps[1]) | set(all_strong_sources[h])
        table = []
        for key in sorted(kernels):
            core = decode(*key)
            cc = components(core)
            row = {'key': list(key), 'rows': list(core), 'SCC_counts': [hdp(induced(core, c)) for c in cc], 'branches': {}}
            for N, mp in zip(targets, maps):
                points = sorted({tuple(w['padding']) for w in mp.get(key, [])})
                minimal = [p for p in points if not any(q != p and q[0] <= p[0] and q[1] <= p[1] for q in points)]
                row['branches'][str(N)] = {'all_core_padding_pairs': [list(p) for p in points], 'minimal_padding_pairs': [list(p) for p in minimal]}
                if points:
                    row['branches'][str(N)]['witness'] = min(mp[key], key=lambda w: (sum(w['padding']), w['padding'], w['cut']))
            if key in common:
                w1, w2 = min(((w1, w2) for w1 in maps[0][key] for w2 in maps[1][key]), key=lambda pair: (max(pair[0]['padding'][0], pair[1]['padding'][0]) + max(pair[0]['padding'][1], pair[1]['padding'][1]), pair[0]['padding'], pair[1]['padding'], pair[0]['cut'], pair[1]['cut']))
                a = max(w1['padding'][0], w2['padding'][0])
                b = max(w1['padding'][1], w2['padding'][1])
                source, cut1 = lift(core, w1, a, b)
                source2, cut2 = lift(core, w2, a, b)
                require(source == source2 and hdp(source) == h, 'common actual source')
                require([hdp(birth(source, cut)) for cut in (cut1, cut2)] == list(targets), 'common branch outputs')
                require([sum(1 for _ in words(birth(source, cut))) for cut in (cut1, cut2)] == list(targets), 'independent common output count')
                row['common_witness'] = {'padding': [a, b], 'source_rows': list(source), 'cuts': [cut1, cut2], 'original_witnesses': [w1, w2]}
            table.append(row)
        report['families'][str(h)] = {'targets': list(targets), 'branch_kernel_counts': [len(m) for m in maps], 'common_kernel_count': len(common), 'kernels': table}
    require(not by_pair[35, 107] and not by_pair[35, 109], 'permanent 35 hostile')
    report['controls'] = {'35_to_107': 'absent', '35_to_109': 'absent', '55_to_111_kernel_count': len(by_pair[55, 111]), '55_to_111_kernels': [list(k) for k in sorted(by_pair[55, 111])]}
    # A finite hostile to deleting the named source: equal H does not imply
    # equal kernels. Add a source/sink and verify exact common-source lifting.
    padding_controls = 0
    for h, family in report['families'].items():
        for row in family['kernels']:
            if 'common_witness' not in row:
                continue
            w = row['common_witness']
            core = tuple(row['rows'])
            a, b = w['padding']
            for delta in ((1, 0), (0, 2), (2, 1)):
                for N, cut in zip(family['targets'], w['cuts']):
                    src, lifted = lift(core, {'padding': [a, b], 'cut': cut}, a + delta[0], b + delta[1])
                    require(hdp(src) == int(h) and hdp(birth(src, lifted)) == N, 'exterior padding control')
                    padding_controls += 1
    report['padding_controls'] = padding_controls
    profile_checks = 0
    for h, family in report['families'].items():
        for row in family['kernels']:
            core = tuple(row['rows'])
            if len(components(core)) != 1:
                continue
            values = attachment_profile(core)
            require(values[0] == values[-1] == int(h), 'uniform cuts preserve H')
            for N in family['targets']:
                actual_cuts = [i for i, value in enumerate(values) if value == N]
                row['branches'][str(N)]['all_unpadded_cuts'] = actual_cuts
                recorded_unpadded = [0, 0] in row['branches'][str(N)]['all_core_padding_pairs']
                require(bool(actual_cuts) == recorded_unpadded, 'independent whole source cut profile')
            if h == '55':
                row['cuts_to_111'] = [i for i, value in enumerate(values) if value == 111]
                require(bool(row['cuts_to_111']) == ((len(core), canonical(core)[0]) in by_pair[55, 111]), '111 source kernel control')
            profile_checks += len(values)
    report['independent_full_cut_profiles'] = {'strong_sources': sum(len(v) for v in all_strong_sources.values()), 'cuts': profile_checks}
    golden_file = ROOT / '05-knowledge/results/golden_hamilton_transfer_20261103/verification.json'
    golden = json.loads(golden_file.read_text())
    report['golden_source_file_sha256'] = sha256(golden_file.read_bytes()).hexdigest()
    report['golden_source_controls'] = []
    for item in golden['H55_sources']:
        source = golden_target(item['word'])
        key, a, b, core, order = normalize(source)
        values = attachment_profile(source)
        require(values == item['all_attachment_counts'], 'golden literal source and cut orientation')
        require(hdp(source) == 55, 'golden source value')
        require(not any(N in values for N in (111, 167, 169)), 'golden source hostile')
        report['golden_source_controls'].append({'key': list(key), 'padding': [a, b], 'no_children': [111, 167, 169]})
    print(json.dumps(report, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
