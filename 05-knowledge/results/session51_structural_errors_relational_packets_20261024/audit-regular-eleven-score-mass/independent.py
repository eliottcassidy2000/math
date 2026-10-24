#!/usr/bin/env python3
"""Independent fixed-labelled score counting and regular11 orbit-mass closure.
No full order-eleven or almost-regular-ten atlas is read. The primary recurrence
uses individual neighbour subsets with unit weights, not grouped choices.
"""
import collections
import functools
import hashlib
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
NONRIGID_PIN = '3fd0d159b9e356fb5d86906e657e9ef992061408f2391419ceaf5d248813acaa'

def require(condition, description):
    if not condition:
        raise ValueError(description)

def feasible(d):
    """Necessary conditions only; neither sufficiency nor Landau is assumed."""
    n = len(d)
    if sum(d) != n * (n - 1) // 2 or any(x < 0 or x >= n for x in d):
        return False
    # Any k vertices contain k(k-1)/2 internal edges.
    return all(sum(d[:k]) >= k * (k - 1) // 2 for k in range(1, n))

@functools.cache
def literal_count(d):
    """Input sorted for memoization, answer for that ONE assigned label vector."""
    if not d:
        return 1
    if not feasible(d):
        return 0
    pivot_degree = d[-1]  # distinguish the last (maximum-degree) vertex
    other = d[:-1]
    total = 0
    # Each chosen subset denotes distinct labelled outneighbours of the pivot.
    for out in itertools.combinations(range(len(other)), pivot_degree):
        out = set(out)
        residual = tuple(sorted(r if j in out else r - 1 for j, r in enumerate(other)))
        total += literal_count(residual)
    return total

@functools.cache
def grouped_count(d):
    """Separate minimum-degree pivot and grouped binomial recurrence."""
    if not d:
        return 1
    n = len(d)
    # Deliberately do not share feasibility pruning with the primary method.
    if any(r < 0 or r >= n for r in d) or sum(d) != n * (n - 1) // 2:
        return 0
    a = d[0]
    groups = list(collections.Counter(d[1:]).items())
    answer = 0
    def visit(i, need, residual, weight):
        nonlocal answer
        if i == len(groups):
            if need == 0:
                answer += weight * grouped_count(tuple(sorted(residual)))
            return
        r, m = groups[i]
        for x in range(min(need, m) + 1):
            visit(i + 1, need - x, residual + [r] * x + [r - 1] * (m - x), weight * math.comb(m, x))
    visit(0, a, [], 1)
    return answer

def brute_score_counts(n):
    """Every labelled tournament, direct edge orientation, no isomorphism."""
    edges = list(itertools.combinations(range(n), 2))
    counts = collections.Counter()
    for mask in range(1 << len(edges)):
        d = [0] * n
        for bit, (u, v) in enumerate(edges):
            d[u if mask & (1 << bit) else v] += 1
        counts[tuple(d)] += 1
    return counts

def exhaustive_small():
    receipts = []
    for n in range(1, 7):
        counts = brute_score_counts(n)
        checks = 0
        for d in itertools.combinations_with_replacement(range(n), n):
            # counts[d] means the exact fixed degree assignment, NOT its orbit.
            expected = counts[d]
            require(literal_count(d) == expected, 'literal versus exhaustive fixed score ' + str(d))
            require(grouped_count(d) == expected, 'grouped versus exhaustive fixed score ' + str(d))
            checks += 1
        for assigned, expected in counts.items():
            key = tuple(sorted(assigned))
            require(literal_count(key) == expected and grouped_count(key) == expected, 'every realized labelled assignment')
        require(sum(counts.values()) == 2 ** (n * (n - 1) // 2), 'all labelled orientations')
        receipts.append(dict(order=n,orientations=sum(counts.values()),nondecreasing_vectors_tested=checks,realized_fixed_assignments_tested=len(counts)))
    return receipts

def validate_rows(rows):
    n = len(rows)
    require(n == 11, 'eleven rows')
    require(all(type(x) is int and 0 <= x < (1 << n) and not (x >> i) & 1 for i,x in enumerate(rows)), 'bounded loopless rows')
    require(all(((rows[i] >> j) & 1) + ((rows[j] >> i) & 1) == 1 for i in range(n) for j in range(i)), 'one directed edge per pair')
    require(all(x.bit_count() == 5 for x in rows), 'regular degrees')

def all_isomorphisms(left, right, stop_after_one=False):
    """Independent full ordinary-isomorphism DFS; signatures only prune."""
    n = len(left)
    if len(right) != n:
        return []
    ld = [r.bit_count() for r in left]
    rd = [r.bit_count() for r in right]
    if sorted(ld) != sorted(rd):
        return []
    lc = [[(left[i] & left[j]).bit_count() for j in range(n)] for i in range(n)]
    rc = [[(right[i] & right[j]).bit_count() for j in range(n)] for i in range(n)]
    lf = [sorted((ld[j],lc[i][j],(left[i]>>j)&1) for j in range(n) if j != i) for i in range(n)]
    rf = [sorted((rd[j],rc[i][j],(right[i]>>j)&1) for j in range(n) if j != i) for i in range(n)]
    initial = [sum(1 << j for j in range(n) if ld[i] == rd[j] and lf[i] == rf[j]) for i in range(n)]
    assigned = {}
    answer = []
    def walk(used):
        if len(assigned) == n:
            answer.append(tuple(assigned[i] for i in range(n)))
            return stop_after_one
        best = None
        for i in range(n):
            if i in assigned:
                continue
            bits = initial[i] & ~used
            candidates = []
            while bits:
                b = bits & -bits
                bits -= b
                j = b.bit_length() - 1
                if all(((left[i] >> old) & 1) == ((right[j] >> image) & 1) and lc[i][old] == rc[j][image] for old,image in assigned.items()):
                    candidates.append(j)
            if not candidates:
                return False
            if best is None or len(candidates) < len(best[1]):
                best = (i, candidates)
        i, candidates = best
        for j in candidates:
            assigned[i] = j
            if walk(used | (1 << j)):
                return True
            del assigned[i]
        return False
    walk(0)
    return answer

def orbit_sizes(group):
    remaining = set(range(11))
    out = []
    while remaining:
        v = min(remaining)
        orb = {p[v] for p in group}
        require(orb <= remaining, 'group orbits partition')
        remaining -= orb
        out.append(len(orb))
    return sorted(out)

def verify_manifest(folder, expected=None):
    raw = (folder / 'MANIFEST.json').read_bytes()
    if expected is not None:
        require(hashlib.sha256(raw).hexdigest() == expected, 'manifest pin')
    manifest = json.loads(raw)
    seen = set()
    for record in manifest['files']:
        name = record['path']
        require(isinstance(name,str) and name not in seen, 'unique manifest path')
        seen.add(name)
        p = Path(name)
        require(not p.is_absolute() and '..' not in p.parts, 'contained manifest path')
        data = (folder / p).read_bytes()
        require(len(data) == record['bytes'], 'byte count: ' + name)
        require(hashlib.sha256(data).hexdigest() == record['sha256'], 'file digest: ' + name)
    return hashlib.sha256(raw).hexdigest()

def nonrigid_records(data):
    candidates = data['q3_candidates']
    genuine = [x for x in candidates if len(x['vertex_orbit_sizes']) == 3]
    hidden = [x for x in candidates if len(x['vertex_orbit_sizes']) != 3]
    require(len(hidden) == 1, 'one hidden q3 transitive presentation')
    require(sum(bool(all_isomorphisms(hidden[0]['rows'],r['rows'],True)) for r in data['q1_classes']) == 1, 'hidden presentation counted once')
    return data['q1_classes'] + genuine + data['q7_classes'] + data['q9_classes']

def audit_nonrigid(data):
    records = nonrigid_records(data)
    certified = []
    for record in records:
        rows = tuple(record['rows'])
        validate_rows(rows)
        g = all_isomorphisms(rows, rows)
        sizes = orbit_sizes(g)
        require(len(g) == record['full_automorphism_order'], 'full ordinary automorphism order')
        require(sizes == record['vertex_orbit_sizes'], 'full vertex orbit sizes')
        require(len(g) > 1 and math.factorial(11) % len(g) == 0, 'nonrigid group divides symmetric group')
        require(record['error'] == 11 * (len(sizes) - 1), 'nonrigid error from full orbits')
        require(all(not all_isomorphisms(rows, prior['rows'], True) for prior in certified), 'distinct nonrigid ordinary classes')
        certified.append(dict(rows=rows,aut=len(g),q=len(sizes),labelled_mass=math.factorial(11)//len(g)))
    pointed = sum(x['q'] for x in certified)
    error = sum(11 * (x['q'] - 1) for x in certified)
    require(len(certified) == data['nonrigid_classes'], 'source class count')
    require(pointed == data['nonrigid_pointed_states'], 'source pointed mass')
    require(error == data['nonrigid_error'], 'source error')
    hist = collections.Counter((x['aut'],x['q']) for x in certified)
    return dict(classes=len(certified),pointed_classes=pointed,error=error,unpointing_deficit=11*len(certified)-pointed,
                labelled_mass=sum(x['labelled_mass'] for x in certified),
                group_and_orbit_spectrum=[dict(aut_order=a,vertex_orbits=q,classes=k,labelled_mass=k*math.factorial(11)//a) for (a,q),k in sorted(hist.items())],
                full_groups_recomputed=len(certified),distinct_pairs=len(certified)*(len(certified)-1)//2)

def derive_closure(total, nonrigid):
    f = math.factorial(11)
    remainder = total - nonrigid['labelled_mass']
    require(remainder >= 0, 'nonnegative rigid mass')
    require(remainder % f == 0, 'rigid mass divisible by 11 factorial')
    rigid = remainder // f
    classes = rigid + nonrigid['classes']
    parents = 11 * rigid + nonrigid['pointed_classes']
    error = 110 * rigid + nonrigid['error']
    require(11 * classes - parents == nonrigid['unpointing_deficit'], 'pointing deficit identity')
    require(error == 11 * (parents - classes), 'pointed error identity')
    return dict(labelled_regular11=total,eleven_factorial=f,nonrigid_labelled_mass=nonrigid['labelled_mass'],
                rigid_labelled_mass=remainder,rigid_regular11_classes=rigid,regular11_classes=classes,
                almost_regular10_classes=parents,regular11_error=error)

def main():
    source_hash = verify_manifest(HERE/'nonrigid', NONRIGID_PIN)
    literal_count.cache_clear()
    grouped_count.cache_clear()
    regular = []
    for n in (1,3,5,7,9,11):
        d = ((n-1)//2,) * n
        x = literal_count(d)
        y = grouped_count(d)
        require(x == y, 'independent recurrences at regular order ' + str(n))
        regular.append(dict(order=n,labelled_count=x))
    fixed_parent = literal_count((4,)*5 + (5,)*5)
    require(fixed_parent == grouped_count((4,)*5 + (5,)*5), 'fixed parent count')
    total = regular[-1]['labelled_count']
    require(total == math.comb(10,5)*fixed_parent, 'labelled extension bijection')
    small = exhaustive_small()
    nonrigid = audit_nonrigid(json.loads((HERE/'nonrigid'/'normal.json').read_text()))
    closure = derive_closure(total, nonrigid)
    output = dict(status='PASS',method='independent literal labelled-neighbour subsets plus grouped/binomial crosscheck',
                  source_manifest_sha256=source_hash,regular_labelled_counts=regular,exhaustive_fixed_score_checks=small,
                  fixed_assigned_parent_count=fixed_parent,score_assignment_multiplier=math.comb(10,5),
                  nonrigid=nonrigid,closure=closure,
                  boundary='Score-count and full automorphism-mass sidecar. No bare-metagraph sufficiency and no unrestricted order-eleven class count claim.')
    print(json.dumps(output,sort_keys=True,indent=2))

if __name__ == '__main__':
    main()
