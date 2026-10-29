"""Source-preserving half-orbit recurrence for fixed three-path covers."""
from itertools import permutations
import json
from random import Random

from atom_precision import hostile
from complement_kernel import HERE, need
from defect_transfer import interval_target


def fixed_query(rows, rho, starts):
    n = len(rows)
    need(n % 2 == 0 and sorted(rho) == list(range(n))
         and all(rho[rho[v]] == v and rho[v] != v for v in range(n)),
         'fixed-point-free involution supplied')
    arc = lambda u, v: bool(rows[u] >> v & 1)
    need(all(arc(u, v) != arc(v, u) and arc(u, v) == arc(rho[v], rho[u])
             for u in range(n) for v in range(u+1, n)),
         'literal tournament anti-involution guard')
    need(len(starts) == 3 and len(set(starts) | {rho[v] for v in starts}) == 6,
         'three distinct involution orbits of source ports')
    used = set(starts) | {rho[v] for v in starts}
    free = [(v, rho[v]) for v in range(n) if v < rho[v] and v not in used]
    h, s0, s1, s2 = len(free), *starts
    full = (1 << h)-1

    def half_paths(start, outer):
        dp = [{} for _ in range(1 << h)]
        dp[0][start] = 1
        result = []
        for mask in range(1 << h):
            result.append(sum(weight for v, weight in dp[mask].items()
                              if arc(v, rho[s2] if outer else rho[v])))
            for v, weight in dp[mask].items():
                for i, pair in enumerate(free):
                    if mask >> i & 1:
                        continue
                    for w in pair:
                        if arc(v, w):
                            nxt = dp[mask | 1 << i]
                            nxt[w] = nxt.get(w, 0)+weight
        return result

    outer, middle = half_paths(s0, True), half_paths(s1, False)
    return sum(outer[mask]*middle[full ^ mask] for mask in range(full+1))


def brute_fixed(rows, rho, starts):
    ends = tuple(rho[v] for v in reversed(starts))
    unused = set(range(len(rows)))-set(starts)-set(ends)
    h = len(unused)
    count, total = 0, 0
    for middle in permutations(sorted(unused)):
        for a in range(h+1):
            for b in range(a, h+1):
                paths = ((starts[0],)+middle[:a]+(ends[0],),
                         (starts[1],)+middle[a:b]+(ends[1],),
                         (starts[2],)+middle[b:]+(ends[2],))
                if not all(rows[u] >> v & 1 for path in paths
                           for u, v in zip(path, path[1:])):
                    continue
                total += 1
                reflected = tuple(tuple(rho[v] for v in reversed(path)) for path in reversed(paths))
                count += reflected == paths
    return total, count


def random_anti(n, rng):
    rho = [n-1-v for v in range(n)]
    arcs = {}
    for u in range(n):
        for v in range(u+1, n):
            if (u, v) in arcs:
                continue
            value = rng.randrange(2)
            for a, b, bit in ((u, v, value), (v, u, 1-value),
                              (rho[v], rho[u], value), (rho[u], rho[v], 1-value)):
                need((a, b) not in arcs or arcs[a, b] == bit, 'consistent anti-orbit orientation')
                arcs[a, b] = bit
    return [sum(arcs.get((u, v), 0) << v for v in range(n)) for u in range(n)], rho


def run():
    rng = Random(202610291340)
    controls = []
    for n in (6, 8, 10, 12):
        for _ in range(12):
            rows, rho = random_anti(n, rng)
            orbit_names = rng.sample(range(n//2), 3)
            starts = tuple(v if rng.randrange(2) else rho[v] for v in orbit_names)
            total, brute = brute_fixed(rows, rho, starts)
            result = fixed_query(rows, rho, starts)
            need(result == brute and result % 2 == total % 2,
                 'independent full covers, fixed action and parity')
            controls.append(dict(rows=rows, rho=rho, starts=starts, total=total, fixed=result))
    periodic = []
    old = json.loads((HERE/'terminal_selector.json').read_text())['reflection']
    for n in range(8, 33, 2):
        _, _, defects = interval_target([2]*(n//2))
        rows = [sum(1 << v for v in range(n) if u != v
                    and ((u < v) != ((min(u, v), max(u, v)) in defects)))
                for u in range(n)]
        value = fixed_query(rows, [n-1-v for v in range(n)], (0, 1, 2))
        if n <= 12:
            need(value == next(r['fixed'] for r in old if r['n'] == n),
                 'previous full enumeration specializes')
        periodic.append(dict(n=n, fixed=value))
    transitive = [sum(1 << j for j in range(i+1, 6)) for i in range(6)]
    need(fixed_query(transitive, [5-i for i in range(6)], (0, 1, 2)) == 1,
         'anti-involution alone does not force evenness')
    hostile_guards = []
    for flip in (False, True):
        rows = hostile(flip)
        scores = sorted(row.bit_count() for row in rows)
        converse_scores = sorted(11-d for d in scores)
        need(scores != converse_scores, '81/27 sources are not self-converse')
        rejected = False
        try:
            fixed_query(rows, [11-i for i in range(12)], (0, 1, 2))
        except ValueError as error:
            rejected = str(error) == 'literal tournament anti-involution guard'
        need(rejected, 'do not force a symmetry map onto either hostile')
        hostile_guards.append(dict(flip=flip, scores=scores, converse_scores=converse_scores))
    return dict(random_controls=controls, periodic=periodic,
                transitive_six_fixed=1, hostile_guards=hostile_guards)


if __name__ == '__main__':
    import sys
    receipt = json.loads(json.dumps(run()))
    path = HERE/'reflection_query.json'
    if '--write' in sys.argv:
        path.write_text(json.dumps(receipt, separators=(',', ':'))+'\n')
    else:
        need(receipt == json.loads(path.read_text()), 'unchanged half-orbit recurrence receipt')
    print('PASS half-orbit fixed-cover recurrence; 48 independent controls and 81/27 guards')
