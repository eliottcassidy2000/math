#!/usr/bin/env python3
"""Complete two-point sections; chart-first and record-first constructions."""
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
RAW = json.loads((HERE / 'inputs.json').read_text())
CORE = RAW['cores']
DATA = {int(s): d for s, d in RAW['records'].items()}
R = [s for s in CORE if DATA[s]['cost'] == 2]
E = sorted(set(CORE) - set(R))
N = 2048
FULL = (1 << N) - 1


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def add(a, b):
    return ((((a >> 4) + (b >> 4)) % 128) << 4) | ((a ^ b) & 15)


def minus(a, b):
    return ((((a >> 4) - (b >> 4)) % 128) << 4) | ((a ^ b) & 15)


def charge(selected):
    q = 0
    for s in selected:
        q = add(q, DATA[s]['label'])
    return q


def chart_first(grade, costly_pad):
    charts = {}
    pads = [()] if not costly_pad else [(s,) for s in E]
    for epad in pads:
        rgrade = grade - len(epad)
        for fixed in combinations(R, rgrade - 3):
            table = {}
            for moving in combinations([s for s in R if s not in fixed], 3):
                state = tuple(sorted(epad + fixed + moving))
                q = charge(state)
                table.setdefault(q, state)
            charts[epad + fixed] = table
    return charts


def record_first(grade, costly_pad):
    """Enumerate actual full records before selecting their fixed coordinates."""
    result = defaultdict(set)
    pads = [()] if not costly_pad else [(s,) for s in E]
    for epad in pads:
        rgrade = grade - len(epad)
        for selected in combinations(R, rgrade):
            q = charge(epad + selected)
            for fixed in combinations(selected, rgrade - 3):
                result[epad + fixed].add(q)
    return result


def coverage(charts):
    rows = [0] * N
    for table in charts.values():
        image = sum(1 << q for q in table)
        for q in table:
            rows[q] |= image
    return rows


def witness(charts, q, target):
    for fixed, table in charts.items():
        if q in table and target in table:
            a, b = table[q], table[target]
            changed = set(a) ^ set(b)
            need(charge(a) == q and charge(b) == target, 'BOTH actual charges, not only their difference')
            need(len(changed) <= 6 and sum(DATA[s]['cost'] for s in changed) <= 12,
                 'common costly/cheap padding is truly unchanged')
            return dict(first=list(a), second=list(b), common_padding=list(fixed),
                        changed=sorted(changed), clock_work=sum(DATA[s]['cost'] for s in changed))
    raise RuntimeError('requested pair absent from complete chart family')


def main():
    families = {}
    checks = 0
    for key, grade, costly in [('cheap5', 5, False), ('padded5', 5, True), ('cheap6', 6, False)]:
        charts = chart_first(grade, costly)
        independent = record_first(grade, costly)
        need({fixed: set(table) for fixed, table in charts.items()} == dict(independent),
             'chart-first and independently full-record-first images agree')
        for fixed, table in charts.items():
            for q, state in table.items():
                need(len(state) == grade and len(set(state)) == grade and set(fixed) <= set(state)
                     and charge(state) == q, 'every charge image entry has an actual full record witness')
                checks += 1
        families[key] = charts
    first = coverage(families['cheap5'])
    repaired = coverage({**families['cheap5'], **families['padded5']})
    six = coverage(families['cheap6'])
    missing = sum(N - row.bit_count() for row in first)
    need(missing == 1008, 'exact hostile: cheap-only grade5 does not have universal two-point coverage')
    need(not (first[2] >> 333 & 1), 'specific cheap-only missing ordered pair (2,333)')
    need(all(row == FULL for row in repaired), 'all 2048 squared grade5 ordered charge pairs covered')
    need(all(row == FULL for row in six), 'all 2048 squared grade6 ordered charge pairs covered')
    # An actual maximum-cost helper swap, rather than an arbitrary debt pair.
    before_helpers, after_helpers = [556, 2125], [556, 753]
    routes = []
    for grade, charts in [(5, {**families['cheap5'], **families['padded5']}), (6, families['cheap6'])]:
        for q in range(N):
            target = add(q, minus(charge(before_helpers), charge(after_helpers)))
            w = witness(charts, q, target)
            if w['clock_work'] == 12:
                break
        need(w['clock_work'] == 12, 'literal six-cheap-slot bridge with maximum15-clock helper exchange')
        total_work = w['clock_work'] + sum(DATA[s]['cost'] for s in set(before_helpers) ^ set(after_helpers))
        need(total_work == 27, 'literal27-clock bridge')
        routes.append(dict(grade=grade, first_core=w['first'], second_core=w['second'],
                          first_helpers=before_helpers, second_helpers=after_helpers,
                          actual_virtual_charge=add(q, charge(before_helpers)),
                          work=total_work, support=8))
    out = dict(status='FINITE-EXACT two-point compatible sections',
               checked_section_witnesses=checks,
               chart_counts={k: len(v) for k, v in families.items()},
               charge_image_entries={k: sum(map(len, v.values())) for k, v in families.items()},
               cheap_only_grade5_missing_ordered_pairs=missing,
               grade5_ordered_pairs_covered=N * N, grade6_ordered_pairs_covered=N * N,
               repaired_hostile=witness(families['padded5'], 2, 333), literal_bridges=routes)
    print(json.dumps(out, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
