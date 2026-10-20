"""Independently reconstruct the blue universe and replay every residual cycle.

No producer modules are imported. The finite universe is all tournaments
on [n] whose standard Hamiltonian path has positional anti-reversal,
for n=4,6,8. The retained filter is strong and exactly two outdegree values.
"""
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def arc(a, u, v):
    return bool(a[u] >> v & 1)


def strong(a):
    n = len(a)
    # Independent two-search criterion, rather than SCC decomposition.
    for backward in (False, True):
        reached = {0}
        stack = [0]
        while stack:
            u = stack.pop()
            for v in range(n):
                if v not in reached and (arc(a, v, u) if backward else arc(a, u, v)):
                    reached.add(v)
                    stack.append(v)
        if len(reached) != n:
            return False
    return True


def universe(n):
    pairs = set(combinations(range(n), 2))
    orbits = []
    while pairs:
        i, j = min(pairs)
        orb = {(i, j), (n-1-j, n-1-i)}
        pairs.difference_update(orb)
        orbits.append(sorted(orb))
    forced = [o for o in orbits if any(j == i+1 for i, j in o)]
    free = [o for o in orbits if o not in forced]
    need(len(free) == (n//2)*(n//2-1), 'wrong reflection-pair dimension')
    for mask in range(1 << len(free)):
        a = [0]*n
        for orb in forced:
            for u, v in orb:
                a[u] |= 1 << v
        for bit, orb in enumerate(free):
            for u, v in orb:
                if mask >> bit & 1:
                    a[u] |= 1 << v
                else:
                    a[v] |= 1 << u
        yield tuple(a)


def main():
    receipt_path = HERE/'moves_residual.json'
    raw = receipt_path.read_bytes()
    data = json.loads(raw)
    need([row['n'] for row in data] == [4, 6, 8], 'wrong receipt orders')
    expected = {4: (4, 3), 6: (64, 21), 8: (4096, 537)}
    rows = []
    for row in data:
        n = row['n']
        generated = set(universe(n))
        retained = {a for a in generated
                    if len({x.bit_count() for x in a}) == 2 and strong(a)}
        matrices = [tuple(c['a']) for c in row['certificates']]
        need(len(matrices) == len(set(matrices)), 'duplicate matrix certificate')
        need(set(matrices) == retained, 'certificate set does not equal full retained universe')
        need((len(generated), len(retained)) == expected[n], 'independent universe mismatch')
        need(row['blue_universe'] == len(generated), 'reported universe differs')
        need(row['strong_two_score_count'] == len(retained), 'reported filter count differs')
        imbalances = {}
        for cert in row['certificates']:
            a, p = tuple(cert['a']), cert['cycle']
            need(sorted(p) == list(range(n)), 'cycle repeats or omits vertices')
            edges = set(zip(p, p[1:]+p[:1]))
            need(all(arc(a, u, v) for u, v in edges), 'missing cycle arc')
            u, v = cert['edge']
            need((u, v) in edges, 'claimed unbalanced edge is outside the cycle')
            imbalance = a[u].bit_count()+a[v].bit_count()-(n-1)
            need(imbalance != 0 and imbalance == cert['imbalance'], 'degree certificate fails')
            imbalances[str(imbalance)] = imbalances.get(str(imbalance), 0)+1
        rows.append(dict(n=n, all_blue_presentations=len(generated),
                         complete_retained_matrices=len(retained),
                         cycle_certificates_replayed=len(matrices),
                         imbalance_histogram=imbalances))
    report = dict(status='PASS: independent complete-universe reconstruction and certificate replay',
                  producer_receipt_sha256=sha256(raw).hexdigest(),
                  audit_source_sha256_lf=sha256(Path(__file__).read_text().encode()).hexdigest(),
                  total_certificates=sum(r['cycle_certificates_replayed'] for r in rows),
                  census=rows)
    output = json.dumps(report, indent=2)+'\n'
    Path(__file__).with_suffix('.json').write_text(output, encoding='utf-8')
    print(output, end='')


if __name__ == '__main__':
    main()
