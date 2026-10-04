"""Exact coverage refinement of the frozen 48 rational-anchor charts.

Finite measures: all admissible m <= 30. Endpoint routes: the first three
admissible m values and lifts 0,1,2 for each chart. All-height separation
uses finite periodic-word and old-bank parent-cylinder certificates.
"""
from collections import Counter
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import argparse
import hashlib
import json
import re
import runpy

ROOT = Path(__file__).resolve().parents[2]
MAX_M = 30
ENDPOINT_HEIGHTS = 3
ENDPOINT_LIFTS = (0, 1, 2)
HUBS = (1, 5, 17, 89)


def need(ok, message):
    if not ok:
        raise ValueError(message)


def intersects(a, b):
    return (a[0]-b[0]) % (1 << min(a[1], b[1])) == 0


def contains(a, b):
    """Whether the whole cylinder b lies in a."""
    return a[1] <= b[1] and intersects(a, b)


def normalize(cylinders):
    """Disjoint antichain; sibling coalescing is unnecessary for exact mass."""
    out = []
    for r, k in sorted(set(cylinders), key=lambda c: (c[1], c[0])):
        need(k >= 0 and 0 <= r < 1 << k, "canonical dyadic cylinder")
        if not any(contains(c, (r, k)) for c in out):
            out.append((r, k))
    return out


def subtract_one(source, removed):
    if not intersects(source, removed):
        return [source]
    if contains(removed, source):
        return []
    # Walk toward the removed descendant, retaining its sibling at each split.
    s, h = removed
    return [(s % (1 << j) + (1-((s >> j) & 1))*(1 << j), j+1)
            for j in range(source[1], h)]


def difference(sources, removed):
    out = normalize(sources)
    for cylinder in normalize(removed):
        out = [child for source in out for child in subtract_one(source, cylinder)]
    return normalize(out)


def mass(cylinders):
    return sum((F(1, 1 << k) for _, k in normalize(cylinders)), F())


def odd_step(n):
    need(n > 0 and n % 2, "positive odd state")
    z = 3*n+1
    a = (z & -z).bit_length()-1
    return z >> a, a


def v2(n):
    need(n != 0, "nonzero valuation argument")
    n = abs(n)
    return (n & -n).bit_length()-1


def dyadic_tests():
    # Independent residue-set oracle over every pair of cylinders of depth <=5.
    cylinders = [(r, k) for k in range(6) for r in range(1 << k)]
    universe = set(range(32))
    count = 0
    for a in cylinders:
        aa = {x for x in universe if x % (1 << a[1]) == a[0]}
        for b in cylinders:
            bb = {x for x in universe if x % (1 << b[1]) == b[0]}
            result = difference([a], [b])
            observed = {x for x in universe
                        if any(x % (1 << k) == r for r, k in result)}
            need(observed == aa-bb, "independent exact difference set")
            need(mass(result) == F(len(observed), 32), "difference mass")
            need(mass([a, b]) == F(len(aa | bb), 32), "exact union mass")
            count += 1
    hostile = difference([(1, 1)], [(1, 3)])
    need(set(hostile) == {(3, 2), (5, 3)} and mass(hostile) == F(3, 8),
         "partial overlap must retain the uncovered siblings")
    return count


def read_inputs():
    path = ROOT / '05-knowledge/results/denominator_arithmetic_filter_20261004.out'
    atlas, _ = json.JSONDecoder().raw_decode(path.read_text(encoding='utf-8-sig'))
    rows = atlas['negative_anchor_gate']['compiler_instances']
    need(len(rows) == 48, "frozen complete rational-only negative chart population")
    rows = sorted(rows, key=lambda r: (len(r['word']), tuple(r['word'])))
    for row in rows:
        row['id'] = 'w' + '_'.join(map(str, row['word']))
    need(len({r['id'] for r in rows}) == 48, "distinct marked words")
    bank_path = ROOT / '05-knowledge/results/reset_20260926_swaplift.out'
    raw = []
    for line in bank_path.read_text().splitlines():
        if re.fullmatch(r'\d+(?: \d+){8}', line):
            values = list(map(int, line.split()))
            raw.append(tuple(values[-2:]))
    need(len(raw) == 171, "frozen old-bank rows")
    bank = normalize(raw)
    need(len(bank) == 65, "frozen disjoint old-bank union")
    need(mass(bank) == F(6985206796614369409, 1 << 65), "old-bank exact density")
    old = runpy.run_path(str(ROOT/'04-computation/experiments/reset_20260926_swaplift.py'))
    regenerated = [old['core_certificate'](q) for q in range(1, 342, 2)]
    need(raw == [(r['residue'], r['K']) for r in regenerated],
         "independent old-bank regeneration agrees row by row")
    return rows, bank, {'atlas_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                        'bank_sha256': hashlib.sha256(bank_path.read_bytes()).hexdigest()}


def pair_separation(rows, compiler):
    proofs = []
    for i, a in enumerate(rows):
        w = tuple(a['word'])
        need(all(w != w[:j]*(len(w)//j) for j in range(1, len(w)) if len(w) % j == 0),
             "primitive marked word")
        for b in rows[i+1:]:
            v = tuple(b['word'])
            period = lcm(len(w), len(v))
            mismatch = next((j for j in range(period)
                             if w[j % len(w)] != v[j % len(v)]), None)
            need(mismatch is not None, "distinct infinite periodic nominal words")
            # Equal first-descent time T is a positive multiple of period.
            # Exact valuations agree before T's final, enlarged valuation.
            # A mismatch j therefore requires T <= j+1, leaving at most T=period.
            m, n = period//len(w), period//len(v)
            exceptional = mismatch == period-1 and m >= a['m'] and n >= b['m']
            if exceptional:
                ca = compiler['cylinder'](w, m)[:2]
                cb = compiler['cylinder'](v, n)[:2]
                need(not intersects(ca, cb), "possible final-letter exception is dyadically disjoint")
            proofs.append(dict(first=a['id'], second=b['id'], comparison_period=period,
                               first_mismatch_zero_based=mismatch,
                               exceptional_time_admissible=exceptional))
    need(len(proofs) == 1128 and not any(p['exceptional_time_admissible'] for p in proofs),
         "all 1128 pairs separate before any admissible exceptional descent")
    return proofs


def baseline_separation(rows, compiler, mixed):
    cycle17 = (1, 1, 1, 2, 1, 1, 4)
    baselines = [('minus5', (), (1, 2), 1), ('minus17', (), cycle17, 1),
                 ('head12_minus17', (1, 2), cycle17, 1),
                 ('head1_minus17', (1,), cycle17, 1),
                 ('rational112', (), (1, 1, 2), 2)]
    proofs = []
    for row in rows:
        w = tuple(row['word'])
        for name, head, cycle, k0 in baselines:
            H, q = len(head), len(cycle)
            bound = H+lcm(len(w), q)
            def value(j):
                return head[j] if j < H else cycle[(j-H) % q]
            mismatch = next((j for j in range(bound)
                             if w[j % len(w)] != value(j)), None)
            if mismatch is None:
                need(w == (1, 1, 2) and name == 'rational112', "only inherited identical family")
                proofs.append(dict(chart=row['id'], baseline=name, relation='identical'))
                continue
            possible = []
            for m in range(row['m'], (mismatch+1)//len(w)+1):
                T = len(w)*m
                if T > H and (T-H) % q == 0 and (T-H)//q >= k0:
                    k = (T-H)//q
                    a = compiler['cylinder'](w, m)[:2]
                    if head:
                        b_row = mixed['compile_family'](head, cycle, k)
                        b = (b_row['residue'], b_row['K'])
                    else:
                        b = compiler['cylinder'](cycle, k)[:2]
                    need(not intersects(a, b), "exceptional baseline overlap excluded exactly")
                    possible.append([m, k])
            proofs.append(dict(chart=row['id'], baseline=name, relation='disjoint',
                               comparison_length=bound, first_mismatch_zero_based=mismatch,
                               possible_early_repetition_pairs=possible))
    need(sum(p['relation'] == 'identical' for p in proofs) == 1,
         "47 genuinely distinct families against all named infinite baselines")
    need(not any(p.get('possible_early_repetition_pairs') for p in proofs),
         "no early clock-compatible baseline exception in this atlas")
    return proofs


def bank_parent(row, bank):
    h, d, A = row['h'], row['d'], sum(row['word'])
    depth = 1
    while True:
        residue = -h*pow(d, -1, 1 << depth) % (1 << depth)
        if all(not intersects((residue, depth), b) for b in bank):
            break
        need(depth <= max(k for _, k in bank), "chart center must escape this finite bank")
        depth += 1
    first_m = max(row['m'], (depth+A-1)//A)
    return dict(residue=residue, exponent=depth, from_m=first_m)


def coverage(rows, bank, compiler):
    ranking, cells = [], []
    for row in rows:
        w = tuple(row['word'])
        P, Q, h, d = compiler['anchor'](w)
        need((h, d) == (row['h'], row['d']), "independently composed anchor")
        m0 = row['m']
        need(Q**m0 > h and (m0 == 1 or Q**(m0-1) <= h), "least admissible repetition")
        parent = bank_parent(row, bank)
        lower = F()
        statuses = Counter()
        excluded = []
        chart_cells = []
        for m in range(m0, MAX_M+1):
            r, K, t = compiler['cylinder'](w, m)
            c = (r, K)
            uncovered = difference([c], bank)
            old_added = mass(uncovered)
            status = 'covered' if not uncovered else 'new' if old_added == F(1, 1 << K) else 'partial'
            statuses[status] += 1
            if status == 'covered':
                excluded.append(m)
            if m >= parent['from_m']:
                need(contains((parent['residue'], parent['exponent']), c) and status == 'new',
                     "all-height parent certificate agrees with finite descendants")
            inherited = w == (1, 1, 2)
            added_cells = [] if inherited else uncovered
            lower += mass(added_cells)
            chart_cells.extend(added_cells)
            cells.append(dict(chart=row['id'], m=m, residue=r, exponent=K, budget=t,
                              old_bank_status=status, added_cells=[list(v) for v in added_cells]))
        tail = F() if w == (1, 1, 2) else F(1, (P-1)*P**MAX_M)
        ranking.append(dict(chart=row['id'], word=list(w), anchor=f'-{h}/{d}',
                            golden_denominator=row['golden_denominator'], minimum_m=m0,
                            missed_bank_parent=parent, finite_old_bank_status=dict(statuses),
                            old_bank_excluded_repetitions=excluded,
                            added_density_lower=str(lower), strict_tail_bound=str(tail),
                            added_density_upper=str(lower+tail),
                            finite_added_cell_count=len(chart_cells),
                            inherited_rational112=(w == (1, 1, 2))))
    ranking.sort(key=lambda r: (-F(r['added_density_lower']), tuple(r['word'])))
    total_cells = [tuple(c) for row in cells for c in row['added_cells']]
    lower = sum((F(r['added_density_lower']) for r in ranking), F())
    need(mass(total_cells) == lower, "independent dyadic union equals all-family disjoint sum")
    need(max(r['missed_bank_parent']['from_m'] for r in ranking) == 3,
         "every chart tail is new to the old bank by repetition three")
    need(not any(r['finite_old_bank_status'].get('partial') for r in ranking),
         "this universe has no partial bank overlap, despite a general subtraction engine")
    tail = sum((F(r['strict_tail_bound']) for r in ranking), F())
    top = [r for r in ranking if tuple(r['word']) in ((1, 1, 1, 3), (1, 1, 2, 2))]
    others = [r for r in ranking if r not in top]
    need(len(top) == 2 and min(F(r['added_density_lower']) for r in top) >
         max(F(r['added_density_upper']) for r in others), "top two dominate every other all-height interval")
    return ranking, cells, dict(lower=str(lower), strict_tail_bound=str(tail), upper=str(lower+tail),
                               finite_union_cell_count=len(normalize(total_cells)))


def endpoint_routes(rows, compiler):
    # A node stores an actual edge and its exact remaining first-hit distance.
    dag = {1: {'next': None, 'valuation': None, 'rank': 0}}
    cap = 10000
    def certify(start):
        n, trail, seen = start, [], set()
        while n not in dag:
            need(len(trail) < cap and n.bit_length() <= cap, "declared endpoint route resource cap")
            need(n not in seen, "no unsealed repeated state")
            seen.add(n)
            out, a = odd_step(n)
            trail.append((n, out, a))
            n = out
        reused_at, new_edges = n, len(trail)
        for n, out, a in reversed(trail):
            dag[n] = dict(next=out, valuation=a, rank=dag[out]['rank']+1)
        return reused_at, new_edges
    for hub in HUBS:
        certify(hub)
    reports = []
    for row in rows:
        w = tuple(row['word'])
        for m in range(row['m'], row['m']+ENDPOINT_HEIGHTS):
            residue, K, t = compiler['cylinder'](w, m)
            for lift in ENDPOINT_LIFTS:
                source = residue+lift*(1 << K)
                endpoint = compiler['verify'](w, m, source)
                need(1 <= endpoint < source and source > 1, "retained source and positive endpoint")
                reused_at, new_edges = certify(endpoint)
                n, first_hit, actual = source, 0, []
                while n != 1:
                    need(first_hit < cap, "independent full-source replay cap")
                    n, a = odd_step(n)
                    actual.append(a)
                    first_hit += 1
                need(first_hit == m*len(w)+dag[endpoint]['rank'], "independent first-hit certificate length")
                nominal = w*m
                need(tuple(actual[:len(nominal)-1]) == nominal[:-1], "literal repeated prefix")
                extra = actual[len(nominal)-1]-nominal[-1]
                need(extra >= t, "literal terminal budget")
                reports.append(dict(chart=row['id'], m=m, lift=lift, source=str(source),
                                    first_descent_step=m*len(w), endpoint=str(endpoint),
                                    extra_terminal_valuation=extra, suffix_pointer=str(endpoint),
                                    suffix_first_hit_steps=dag[endpoint]['rank'],
                                    full_first_hit_steps=first_hit,
                                    preexisting_suffix_join=str(reused_at), new_suffix_edges=new_edges))
    need(len(reports) == 432, "complete finite endpoint universe")
    visits = Counter()
    for report in reports:
        n = int(report['endpoint'])
        while n != 1:
            visits[n] += 1
            edge = dag[n]
            need(edge['rank'] == dag[edge['next']]['rank']+1, "well-founded shared suffix rank")
            need(odd_step(n) == (edge['next'], edge['valuation']), "all retained suffix edges are literal")
            n = edge['next']
    stats = dict(source_count=len(reports), distinct_endpoints=len({r['endpoint'] for r in reports}),
                 unique_suffix_edges=len(dag)-1,
                 unshared_suffix_edge_count=sum(r['suffix_first_hit_steps'] for r in reports),
                 maximum_suffix_first_hit_steps=max(r['suffix_first_hit_steps'] for r in reports),
                 maximum_full_first_hit_steps=max(r['full_first_hit_steps'] for r in reports),
                 maximum_source_bit_length=max(int(r['source']).bit_length() for r in reports),
                 top_shared_nodes=[dict(state=str(n), routes=count, suffix_length=dag[n]['rank'])
                                   for n, count in visits.most_common(12)],
                 odd_step_cap=cap, bit_length_cap=cap)
    exported = {str(n): dict(next=None if edge['next'] is None else str(edge['next']),
                            valuation=edge['valuation'], rank=edge['rank'])
                for n, edge in sorted(dag.items())}
    return reports, exported, stats


def sealed_hub_families(rows, compiler, dag):
    """Symbolic all-height seals; expand only when the exponent is <=4096."""
    seals, literal = [], 0
    for row in rows:
        w, m = tuple(row['word']), row['m']
        P, Q, h, d = compiler['anchor'](w)
        _, _, budget = compiler['cylinder'](w, m)
        power, modulus = len(w)*m, P**m
        period = 2*3**(power-1)
        for hub in HUBS:
            target = -h*pow(d*hub, -1, modulus) % modulus
            t0 = compiler['log2_mod3'](target, power)
            t = t0 + max(0, (budget-t0+period-1)//period)*period
            need(t >= budget and (h+d*hub*pow(2, t, modulus)) % modulus == 0,
                 "complete ternary seal guard and original-source exit budget")
            need(pow(2, period, modulus) == 1 and pow(2, period//2, modulus) != 1,
                 "even part of exact order")
            if power > 1:
                need(pow(2, period//3, modulus) != 1, "ternary part of exact order")
            entry = dict(chart=row['id'], m=m, hub=hub, least_exponent=t,
                         exponent_period=period, minimum_budget=budget,
                         hub_suffix_first_hit_steps=dag[str(hub)]['rank'],
                         source_expression=f'((({h}+{d}*{hub}*2^tau)/{P}^{m})*{Q}^{m}-{h})/{d}',
                         exponents='tau=least_exponent+exponent_period*k, k>=0')
            if t <= 4096:
                b = (h+d*hub*(1 << t))//modulus
                numerator = b*Q**m-h
                need(numerator % d == 0, "literal sealed source integral")
                source = numerator//d
                need(compiler['verify'](w, m, source) == hub, "literal sealed first descent to shared hub")
                entry['expanded_source_bit_length'] = source.bit_length()
                literal += 1
            seals.append(entry)
    need(len(seals) == 192, "all 48 charts times four known hubs")
    return seals, literal


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', type=Path)
    args = parser.parse_args()
    compiler = runpy.run_path(str(ROOT/'04-computation/experiments/rational_anchor_returns_20261004.py'))
    mixed = runpy.run_path(str(ROOT/'04-computation/experiments/collatz_mixed_return_compiler_20261003.py'))
    checks = dyadic_tests()
    rows, bank, provenance = read_inputs()
    pairs = pair_separation(rows, compiler)
    baseline = baseline_separation(rows, compiler, mixed)
    ranking, cells, total = coverage(rows, bank, compiler)
    routes, dag, route_stats = endpoint_routes(rows, compiler)
    seals, literal_seals = sealed_hub_families(rows, compiler, dag)
    for row in ranking:
        subset = [r for r in routes if r['chart'] == row['chart']]
        options = [s for s in seals if s['chart'] == row['chart']]
        best = min(options, key=lambda s: (s['least_exponent'], s['hub']))
        row['finite_suffix_steps_sum'] = sum(r['suffix_first_hit_steps'] for r in subset)
        row['finite_suffix_steps_maximum'] = max(r['suffix_first_hit_steps'] for r in subset)
        row['best_sealed_hub'] = best['hub']
        row['best_seal_exponent'] = best['least_exponent']
    report = dict(status='PROVED separation and completion mechanisms; FINITE-EXACT declared atlas; OPEN global coverage',
                  finite_universe=dict(charts=48, max_coverage_m=MAX_M,
                                       route_heights_per_chart=ENDPOINT_HEIGHTS, route_lifts=list(ENDPOINT_LIFTS)),
                  provenance=provenance, independent_dyadic_pair_tests=checks,
                  old_bank=dict(raw_rows=171, disjoint_cylinders=65, density=str(mass(bank))),
                  all_height_chart_pair_certificates=pairs,
                  all_height_baseline_certificates=baseline,
                  total_added_density=total, ranking=ranking, finite_cylinders=cells,
                  endpoint_route_statistics=route_stats, completed_routes=routes, shared_suffix_dag=dag,
                  symbolic_hub_seals=seals, literal_hub_seals_checked=literal_seals)
    if args.json:
        args.json.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print('PROVED: 48 pairwise disjoint all-height chart families; 47 disjoint from all five named prior families.')
    print('FINITE-EXACT:', checks, 'independent dyadic pair tests; 171 old rows regenerate to65 disjoint cylinders.')
    print('Complete coverage universe: all admissible m<=30; cylinder count:', len(cells))
    print('All chart tails avoid the old bank from m=3 onward; parent depths:',
          dict(sorted(Counter(r['missed_bank_parent']['exponent'] for r in ranking).items())))
    print('Early old-bank exclusions:', [(r['chart'], r['old_bank_excluded_repetitions'])
                                         for r in ranking if r['old_bank_excluded_repetitions']])
    print('Exact added density lower:', total['lower'])
    print('Strict omitted tail upper:', total['strict_tail_bound'])
    print('Illustrative decimal lower:', format(float(F(total['lower'])), '.17g'))
    print('Finite exact ranking: chart, added mass through30, best sealed hub/exponent, 9-source suffix steps sum/max')
    for row in ranking:
        print(row['chart'], row['added_density_lower'], row['best_sealed_hub'], row['best_seal_exponent'],
              row['finite_suffix_steps_sum'], row['finite_suffix_steps_maximum'])
    print('Two top all-height density intervals exceed all remaining upper bounds: w1_1_1_3,w1_1_2_2')
    print('Endpoint routes:', json.dumps(route_stats, sort_keys=True))
    print('Symbolic sealed hub families:', len(seals), '; literal exponent<=4096 replays:', literal_seals)
    print('PASS: all checks remain active under optimized Python; no global coverage assertion.')


if __name__ == '__main__':
    main()
