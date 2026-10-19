#!/usr/bin/env python3
"""Source-authenticated Nyquist cancellation; no entry or optimality oracle.

The shared route DAG ends at the fixed 627-point frontier.  Literal source
rows include zero weights.  Generation discovers first arrivals under a hard
cap; verification reproduces every edge and every retained source weight.
"""
from collections import Counter, defaultdict, deque
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
import argparse
import json
from pathlib import Path

import root_clock_channels_20261018o as inherited

HERE = Path(__file__).resolve().parent
RESULT = HERE.parent.parent / '05-knowledge/results/native_cancellation_105_20261018g.json'
FRONTIER = HERE / 'first_arrival_20261018/V_closed_Y1000.txt'
SCALES = (4096, 5997, 11002, 16384)
CAP = 10000
CHECKS = 0


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def integer(x, minimum=0):
    if type(x) is not int or x < minimum:
        raise ValueError('exact integer outside domain')
    return x


def primitive(s):
    integer(s, 1)
    if not s & 1 or s % 8 == 5:
        raise ValueError('not a canonical primitive source')
    return s


def digest(obj):
    return sha256(json.dumps(obj, separators=(',', ':'), sort_keys=True).encode()).hexdigest()


def add(v, w):
    return v[0] + w[0], v[1] + w[1]


def negate(v):
    return -v[0], -v[1]


def cross(v, w):
    a, b = v
    c, d = w
    return 2*a*c - a*d - b*c + 2*b*d


ALPHABET = tuple((a, b) for a in range(-2, 3) for b in range(-2, 3)
                 if inherited.norm((a, b)) in (1, 3))
TRIPLES = tuple(t for t in combinations_with_replacement(ALPHABET, 3)
                if tuple(map(sum, zip(*t))) == (0, 0))


def charge(s, n):
    primitive(s)
    integer(n, 1)
    if s >= 4*n:
        raise ValueError('source outside its retained scale')
    return inherited.charge(s, n)


def independent_charge(s, n):
    """Direct shell-coordinate sum, independent of the three-real-row decoder."""
    lower = ((1, 2), (-2, -1), (1, -1))
    upper = ((0, 1), (-1, -1), (1, 0))
    result = (0, 0)
    while s < 4*n:
        result = add(result, (lower if s < 2*n else upper)[s % 3])
        s = 4*s + 1
    return result


class Routes:
    def __init__(self, frontier):
        if type(frontier) is not frozenset or len(frontier) != 627:
            raise ValueError('fixed frontier required')
        for s in frontier:
            primitive(s)
            if inherited.edge(s)[0] not in frontier:
                raise ValueError('frontier not closed')
        self.frontier = frontier
        self.memo = {s: (s, 0, 0) for s in frontier}
        self.edges = {}

    def route(self, s):
        primitive(s)
        x, path, seen = s, [], set()
        while x not in self.memo:
            if x in seen or len(path) >= CAP:
                raise RuntimeError(('UNRESOLVED first arrival: no certificate', s, x))
            seen.add(x)
            y, dt = inherited.edge(x)
            primitive(y)
            integer(dt, 1)
            path.append((x, y, dt))
            x = y
        terminal, clock, depth = self.memo[x]
        if len(path) + depth > CAP:
            raise RuntimeError(('UNRESOLVED cap after cached suffix', s))
        for x, y, dt in reversed(path):
            clock += dt
            depth += 1
            self.edges[x] = (y, dt)
            self.memo[x] = (terminal, clock, depth)
        return self.memo[s]


def source_row(s, n, routes):
    a, b = charge(s, n)
    u, t, d = routes.route(s)
    # Source, both compatible registers, terminal, integer clock, depth,
    # unsigned Eisenstein coordinates, and the original source norm.
    return [s, s % 15, s % 35, u, t, d, a, b, inherited.norm((a, b))]


def signed(row):
    return (row[6], row[7]) if row[4] % 2 == 0 else (-row[6], -row[7])


def audit_row(row, n, routes):
    if type(row) is not list or len(row) != 9 or any(type(x) is not int for x in row):
        raise ValueError('exact source row required')
    if row != source_row(row[0], n, routes):
        raise ValueError('source, charge, route or compatible register changed')
    return row


def packet_data(source_ids, rows):
    if type(source_ids) is not list or len(source_ids) not in (2, 3):
        raise ValueError('a pair or triple required')
    if any(type(s) is not int for s in source_ids) or len(set(source_ids)) != len(source_ids):
        raise ValueError('distinct exact source identities required')
    if any(s not in rows for s in source_ids):
        raise ValueError('source absent from this physical row')
    selected = [rows[s] for s in source_ids]
    if len({r[3] for r in selected}) != 1 or any(r[8] == 0 for r in selected):
        raise ValueError('nonzero sources must have the same actual terminal')
    vectors = [signed(r) for r in selected]
    if tuple(map(sum, zip(*vectors))) != (0, 0):
        raise ValueError('packet is not Nyquist zero')
    mass = sum(r[8] for r in selected)
    debt = sum(cross(v, w) for v, w in combinations(vectors, 2))
    if debt != -mass:
        raise ValueError('cross-term debt mismatch')
    return selected[0][3], mass, debt


def greedy_packets(rows):
    buckets = defaultdict(lambda: defaultdict(deque))
    for row in rows:
        if row[8]:
            buckets[row[3]][signed(row)].append(row[0])
    packets = []
    for u in sorted(buckets):
        q = buckets[u]
        for v in ALPHABET:
            w = negate(v)
            if v < w:
                while q[v] and q[w]:
                    packets.append(sorted([q[v].popleft(), q[w].popleft()]))
        for triple in TRIPLES:
            needed = Counter(triple)
            while all(len(q[v]) >= count for v, count in needed.items()):
                packets.append(sorted(q[v].popleft() for v in triple))
    return packets


def terminal_vectors(rows, omitted=frozenset()):
    ans = defaultdict(lambda: (0, 0))
    for row in rows:
        if row[0] not in omitted:
            ans[row[3]] = add(ans[row[3]], signed(row))
    return {u: v for u, v in ans.items() if v != (0, 0)}


def audit_scale(n, rows, packets, routes):
    integer(n, 1)
    expected = list(range(1, 4*n, 2))
    expected = [s for s in expected if s % 8 != 5]
    if type(rows) is not list or [r[0] for r in rows] != expected:
        raise ValueError('complete ordered original source universe required')
    for row in rows:
        audit_row(row, n, routes)
        check((row[6], row[7]) == independent_charge(row[0], n), 'independent shell charges')
    bysource = {r[0]: r for r in rows}
    removed, mass = set(), 0
    histogram = Counter()
    frequency_hostile = None
    for ids in packets:
        u, mu, debt = packet_data(ids, bysource)
        if removed.intersection(ids):
            raise ValueError('a source was spent twice')
        removed.update(ids)
        mass += mu
        histogram[len(ids)] += 1
        check(debt == -mu, 'paid packet cross debt')
        if frequency_hostile is None:
            zero = tuple(sum(bysource[s][6+j] for s in ids) for j in range(2))
            if zero != (0, 0):
                frequency_hostile = {'sources': ids, 'terminal': u,
                    'zero_frequency_vector': list(zero), 'zero_frequency_norm': inherited.norm(zero)}
    before, after = terminal_vectors(rows), terminal_vectors(rows, removed)
    check(before == after, 'every terminal vector retained exactly')
    energy = sum(inherited.norm(v) for v in before.values())
    source_mass = sum(r[8] for r in rows)
    check(energy == sum(inherited.norm(v) for v in after.values()), 'Nyquist energy retained')
    check(sum(r[8] for r in rows if r[0] not in removed) == source_mass-mass, 'source capacity ledger')
    check(frequency_hostile is not None, 'cancellation is not frequency-blind')
    return {'N': n, 'sources': len(rows), 'positive_sources': sum(r[8] != 0 for r in rows),
            'source_norm': source_mass, 'terminal_energy': energy,
            'pairs': histogram[2], 'triples': histogram[3], 'removed_sources': len(removed),
            'removed_source_norm': mass, 'cross_debt': -mass,
            'residual_source_norm': source_mass-mass,
            'nonzero_terminal_count': len(before),
            'max_depth': max(r[5] for r in rows), 'max_clock': max(r[4] for r in rows),
            'source_rows_sha256': digest(rows), 'packets_sha256': digest(packets),
            'terminal_vectors': [[u, *v] for u, v in sorted(before.items())],
            'frequency_hostile': frequency_hostile}


def expected_rejection(fn, why):
    try:
        fn()
    except (ValueError, RuntimeError):
        check(True, why)
    else:
        check(False, why)


def run():
    global CHECKS
    CHECKS = 0
    ranks = {int(line.split()[0]): int(line.split()[1]) for line in FRONTIER.read_text().splitlines()
             if line.strip() and not line.startswith('#')}
    frontier = frozenset(ranks)
    routes = Routes(frontier)
    check(ranks[1] == 0, 'frontier ROOT rank')
    for s in frontier - {1}:
        check(ranks[inherited.edge(s)[0]] < ranks[s], 'inherited frontier has a strict ROOT rank')
    check(len(ALPHABET) == 12, 'the twelve nonzero Eisenstein shell charges')
    reference = json.loads((RESULT.parent/'root_clock_channels_20261018o.json').read_text())
    reference_rows = {r['N']: r for r in reference['rows']}
    scales = []
    for n in SCALES:
        rows = [source_row(s, n, routes) for s in range(1, 4*n, 2) if s % 8 != 5]
        packets = greedy_packets(rows)
        stats = audit_scale(n, rows, packets, routes)
        check(stats['source_norm'] == reference_rows[n]['source_norm'], 'inherited original mass')
        check(stats['terminal_energy'] == reference_rows[n]['frequency_norms']['pi'], 'inherited Nyquist energy')
        # All original sources, not only selected packet members, survive scale transport.
        for row in rows:
            check(charge(row[0], 64*n) == (row[6], row[7]), 'exact scale64 source transport')
        scales.append({'statistics': stats, 'source_rows': rows, 'packets': packets})
        print(n, 'pairs/triples', stats['pairs'], stats['triples'],
              'paid norm', stats['removed_source_norm'], 'energy', stats['terminal_energy'], flush=True)

    # Small complete integer-scale universe, independent direct coordinates.
    four_hostile = None
    for n in range(1, 129):
        for s in range(1, 4*n, 2):
            if s % 8 == 5:
                continue
            v = independent_charge(s, n)
            check(v == independent_charge(s, 64*n), 'all small scale64 controls')
            if four_hostile is None and v != independent_charge(s, 4*n):
                four_hostile = {'N': n, 'source': s, 'A_N': list(v),
                                'A_4N': list(independent_charge(s, 4*n))}
    check(four_hostile is not None, 'scale4 does not transport charges')
    packet_four_hostile = None
    for item in scales:
        n = item['statistics']['N']
        bysource = {r[0]: r for r in item['source_rows']}
        for ids in item['packets']:
            new = (0, 0)
            for s in ids:
                v = charge(s, 4*n)
                new = add(new, v if bysource[s][4] % 2 == 0 else negate(v))
            if new != (0, 0):
                packet_four_hostile = {'N': n, 'sources': ids,
                    'terminal': bysource[ids[0]][3], 'signed_sum_at_4N': list(new),
                    'norm_at_4N': inherited.norm(new)}
                break
        if packet_four_hostile is not None:
            break
    check(packet_four_hostile is not None, 'scale4 can destroy an actual zero packet')
    # Symbolic-sized scales are still read by integer B-chain cutoffs.
    large = 10**40 + 7
    for s in (1, 3, 7, 4*large-3, 4*large-1):
        if s % 8 != 5:
            for k in (1, 2, 5):
                check(charge(s, large) == charge(s, large*64**k), 'iterated large-scale fixed source')

    # Public-boundary mutation hostiles, including a whole-packet duplicate.
    first = scales[0]
    rows, packets = first['source_rows'], first['packets']
    good = rows[0]
    for index, value in ((0, True), (1, float(good[1])), (2, (good[2]+1) % 35),
                         (3, good[3]+2), (4, good[4]+2), (6, good[6]+1), (8, good[8]+1)):
        bad = good[:]
        bad[index] = value
        expected_rejection(lambda b=bad: audit_row(b, SCALES[0], routes), 'modified source row rejected')
    bysource = {r[0]: r for r in rows}
    bad = packets[0][:]
    bad[0] = float(bad[0])
    expected_rejection(lambda: packet_data(bad, bysource), 'float source ID rejected')
    expected_rejection(lambda: packet_data([packets[0][0]]*2, bysource), 'same source twice rejected')
    expected_rejection(lambda: audit_scale(SCALES[0], rows, [packets[0], packets[0]], routes),
                       'double-spend packet rejected')
    expected_rejection(lambda: charge(1, True), 'Boolean scale rejected')
    expected_rejection(lambda: charge(5, 2), 'nonprimitive braid source rejected')
    expected_rejection(lambda: charge(9, 2), 'source outside scale rejected')

    # Every retained route node has one authenticated successor and strictly smaller
    # remaining depth; this independently verifies the memoized first-arrival DAG.
    edge_rows = []
    for s, (v, dt) in sorted(routes.edges.items()):
        check((v, dt) == inherited.edge(s), 'actual edge replay')
        u, t, d = routes.memo[s]
        uu, tt, dd = routes.memo[v]
        check(s not in frontier and (u, t, d) == (uu, tt+dt, dd+1), 'first-arrival suffix and rank')
        edge_rows.append([s, v, dt, u, t, d])
    return {'status': 'PROVED scoped scale transport; FINITE-EXACT native packets; universal entry OPEN',
            'checks': CHECKS, 'scales': scales, 'cap': CAP, 'unresolved_routes': 0,
            'frontier': sorted(frontier), 'frontier_sha256': digest(sorted(frontier)),
            'frontier_ranks': [[s, ranks[s]] for s in sorted(frontier)],
            'producer_sha256': sha256((HERE/'root_clock_channels_20261018o.py').read_bytes().replace(b'\r\n', b'\n')).hexdigest(),
            'row_schema': ['source', 'mod15', 'mod35', 'terminal', 'clock', 'depth', 'a', 'b', 'source_norm'],
            'route_schema': ['source', 'next', 'edge_clock', 'terminal', 'remaining_clock', 'remaining_depth'],
            'route_dag': edge_rows, 'route_dag_sha256': digest(edge_rows),
            'alphabet': ALPHABET, 'zero_triple_types': TRIPLES, 'scale4_hostile': four_hostile,
            'packet_scale4_hostile': packet_four_hostile}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    payload = run()
    text = json.dumps(payload, sort_keys=True, separators=(',', ':')) + '\n'
    if args.verify:
        if RESULT.read_text(encoding='utf-8').replace('\r\n', '\n') != text:
            raise RuntimeError('saved receipt differs')
    else:
        RESULT.write_text(text, encoding='utf-8', newline='\n')
    print('checks', CHECKS, 'LF_SHA256', sha256(text.encode()).hexdigest())


if __name__ == '__main__':
    main()
