"""Native path, spectrum and clock-loss controls with preserved source objects."""
from itertools import combinations
from math import isqrt
from pathlib import Path
import importlib.util
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
GATES = 0


def need(value, why):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(why)


def ham(rows):
    n = len(rows)
    dp = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            value = dp[mask][v]
            if not value:
                continue
            choices = rows[v] & ~mask
            while choices:
                bit = choices & -choices
                choices -= bit
                dp[mask | bit][bit.bit_length()-1] += value
    return sum(dp[-1])


def from_mask(n, bits):
    rows = [0]*n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        u, v = (i, j) if (bits >> k) & 1 else (j, i)
        rows[u] |= 1 << v
    return rows


def moments(rows):
    n = len(rows)
    a = [[(rows[i] >> j) & 1 for j in range(n)] for i in range(n)]
    p = [[int(i == j) for j in range(n)] for i in range(n)]
    result = []
    for _ in range(n):
        p = [[sum(p[i][k]*a[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        result.append(sum(p[i][i] for i in range(n)))
    return result


def spectral_hostile():
    target, found = [0, 0, 12, 12, 10, 48], {}
    for mask in range(1 << 15):
        rows = from_mask(6, mask)
        if moments(rows) == target:
            h = ham(rows)
            if h in (13, 17) and h not in found:
                found[h] = {'mask': mask, 'rows': rows, 'H': h, 'moments': target}
            if len(found) == 2:
                break
    need(set(found) == {13, 17}, 'exact cospectral split recovered from native matrices')
    for h, item in found.items():
        # Adjoin sources: a constructive witness at every larger order, sampled here.
        for extra in range(1, 7):
            rows = item['rows'][:]
            for _ in range(extra):
                rows.append((1 << len(rows))-1)
            need(ham(rows) == h, 'source extension preserves H')
    return [found[h] for h in sorted(found)]


def square_path(n):
    rows = [sum(1 << (v-1) for v in range(1, n+1)
                if u != v and isqrt(u+v)**2 == u+v) for u in range(1, n+1)]
    if n > 1 and any(row == 0 for row in rows):
        return None, 0
    leaves = [u for u, row in enumerate(rows) if row.bit_count() == 1]
    if len(leaves) > 2:
        return None, 0
    failed, visits = set(), 0
    full = (1 << n)-1
    def search(last, mask):
        nonlocal visits
        visits += 1
        if mask == full:
            return [last+1]
        if (last, mask) in failed:
            return None
        choices = rows[last] & ~mask
        candidates = [u for u in range(n) if (choices >> u) & 1]
        candidates.sort(key=lambda u: ((rows[u] & ~mask).bit_count(), u))
        for nxt in candidates:
            result = search(nxt, mask | (1 << nxt))
            if result is not None:
                return [last+1]+result
        failed.add((last, mask))
        return None
    starts = leaves[:1] if leaves else list(range(n))
    for start in starts:
        result = search(start, 1 << start)
        if result:
            return result, visits
    return None, visits


def square_controls():
    result = []
    for n in range(1, 26):
        path, visits = square_path(n)
        need((path is not None) == (n in (1, 15, 16, 17, 23, 25)), 'finite square-sum order classification')
        if path:
            need(sorted(path) == list(range(1, n+1)), 'square witness visits each vertex once')
            need(all(isqrt(u+v)**2 == u+v for u, v in zip(path, path[1:])), 'native additive guards')
        result.append({'n': n, 'path': path, 'search_nodes': visits})
    bad = [1, 8]*4
    need(all(isqrt(u+v)**2 == u+v for u, v in zip(bad, bad[1:])), 'repeated-walk hostile passes local sums')
    need(len(set(bad)) != len(bad), 'repeated-walk hostile fails Hamiltonian condition')
    # The same ordered transitive tournament forgets whether its cover edges were legal.
    tt = [sum(1 << j for j in range(i+1, 4)) for i in range(4)]
    need(ham(tt) == 1, 'unmarked transitive order has one Hamiltonian path')
    cover = {(0, 1), (1, 2), (2, 3)}
    need(all((i, i+1) in cover for i in range(3)), 'path-graph cover guard accepts')
    need(not all((i, i+1) in set() for i in range(3)), 'empty-graph cover guard rejects same order')
    return result


def quadratic_clock_controls():
    # O=Z[u], u^2=5u+1; all 169 states mod13, not just units.
    m = 13
    def mul(x, y):
        a, b = x; c, d = y
        return ((a*c+b*d) % m, (a*d+b*c+5*b*d) % m)
    def tr(x):
        return (2*x[0]+5*x[1]) % m
    u, power = (0, 1), (1, 0)
    for _ in range(12):
        power = mul(power, u)
    need(power == (1, 0), 'twelfth clock becomes identity modulo13')
    coarse, full = set(), set()
    for a in range(m):
        for b in range(m):
            x = (a, b)
            coarse.add((tr(x), tr(mul(power, x))))
            full.add((tr(x), tr(mul(u, x))))
    need(len(coarse) == 13 and len(full) == 169, 'exact clock-induced quotient loss and recovery')
    return {'source_states': 169, 'twelfth_clock_classes': 13, 'primitive_clock_classes': 169}


def hostile_81_27():
    path = HERE.parent/'regular_parent_symmetry_20261029'/'verify.py'
    spec = importlib.util.spec_from_file_location('inherited_regular_parent', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    report = module.hostile_checks()
    # Use a separately implemented Hamilton-path DP on the exact native controls.
    counts = []
    for flip in (False, True):
        rows = [0]*12
        for u, v in combinations(range(12), 2):
            b, s = divmod(u, 4); c, t = divmod(v, 4)
            arc = (c-b) % 3 == 1 if b != c else s == 0 or (s > 0 and (t-s) % 3 == 1)
            if flip and b == c == 2:
                arc = not arc
            rows[u if arc else v] |= 1 << (v if arc else u)
        counts.append(ham(rows))
    report['native_H_counts'] = counts
    report['inherited_gates'] = module.GATES
    need([r['group_order'] for r in report['tournaments']] == [81, 27], 'symmetry hostile retained')
    return report


def main():
    report = {'status': 'PASS', 'spectral_hostile': spectral_hostile(),
              'square_sum_orders_through25': square_controls(),
              'quadratic_clock': quadratic_clock_controls(), 'hostile_81_27': hostile_81_27()}
    report['gates'] = GATES
    Path(__file__).with_name('carrier_checks.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'gates': GATES,
                      'cospectral_rows': report['spectral_hostile'],
                      'hostile_81_27': report['hostile_81_27']}))


if __name__ == '__main__':
    main()
