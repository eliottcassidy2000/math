"""Independent exact referee: affine-coset repair compiler and native barriers.

No producer import. Geometry is counted from complete integer grid-line masks;
the producer used triangle determinants and anchored direction counts. Quotient
edges here exchange physical values. Reciprocal events are rebuilt literally.
The analytic all-prime conclusions are proved in the companion report.
"""
from collections import Counter, defaultdict, deque
from itertools import combinations, permutations
from math import comb, factorial, gcd, isqrt
from pathlib import Path
import argparse
import hashlib
import json
import sys

sys.stdout.reconfigure(newline="\n")
HERE = Path(__file__).resolve().parent
STEM = "continuing5_20260906_wildcard_swap_compiler"
parser = argparse.ArgumentParser()
parser.add_argument("--certificate", type=Path)
args = parser.parse_args()
cert_path = args.certificate or HERE / (STEM + "_certificate.json")
if not cert_path.exists():
    cert_path = HERE.parent / "05-knowledge" / "results" / cert_path.name
source_path = HERE / (STEM + ".py")
GATES = 0


def require(condition, label):
    global GATES
    if not condition:
        raise RuntimeError(label)
    GATES += 1


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


require(digest(source_path) == "d15877188ab7714390cf9334838151c9b20be39ca553d333530a14ce9cdccad1", "frozen producer source")
require(digest(cert_path) == "cb41dab50c3505c48dac732f67013337bb936bddf1f0d762a08e48e3a337a982", "frozen certificate")
require(b"\r" not in cert_path.read_bytes(), "literal LF certificate")
data = json.loads(cert_path.read_text(encoding="utf-8"))


def det(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1]) - (c[0]-a[0])*(b[1]-a[1])


def key_line(a, b):
    aa, bb = b[1]-a[1], a[0]-b[0]
    cc = aa*a[0] + bb*a[1]
    div = gcd(aa, bb)
    aa, bb, cc = aa//div, bb//div, cc//div
    if aa < 0 or (aa == 0 and bb < 0):
        aa, bb, cc = -aa, -bb, -cc
    return aa, bb, cc


class Grid:
    def __init__(self, p):
        self.p = p
        points = [(i, j) for i in range(p) for j in range(p)]
        masks = defaultdict(int)
        for i, j in combinations(range(p*p), 2):
            masks[key_line(points[i], points[j])] |= (1 << i) | (1 << j)
        self.lines = [mask for mask in masks.values() if mask.bit_count() >= 3]
        self.choose = [comb(i, 3) if i >= 3 else 0 for i in range(p*p+1)]

    def mask(self, points):
        out = 0
        for x, y in points:
            out |= 1 << (self.p*x+y)
        return out

    def count(self, points):
        mask = self.mask(points)
        return sum(self.choose[(mask & line).bit_count()] for line in self.lines)

    def unsafe(self, points):
        mask = self.mask(points)
        return any((mask & line).bit_count() >= 3 for line in self.lines)


def cycle_board(p, col, row=None):
    row = tuple(range(p)) if row is None else row
    return [(row[i], col[j]) for i in range(p) for j in (i, (i+1) % p)]


def normalize(pi):
    p = len(pi)
    scale = pow((pi[1]-pi[0]) % p, -1, p)
    return tuple((x-pi[0])*scale % p for x in pi)


def af(pi, a, b):
    return tuple((a*x+b) % len(pi) for x in pi)


def swap_values(pi, a, b):
    return tuple(b if x == a else a if x == b else x for x in pi)


def word_length(pi):
    unseen = set(range(len(pi)))
    cycles = 0
    while unseen:
        j = next(iter(unseen))
        while j in unseen:
            unseen.remove(j)
            j = pi[j]
        cycles += 1
    return len(pi)-cycles


def distance_formula(sig, pi):
    p = len(sig)
    inv = {v:i for i, v in enumerate(sig)}
    return min(word_length(tuple(inv[(a*v+b) % p] for v in pi))
               for a in range(1, p) for b in range(p))


def bfs(adj, start, scores=None):
    levels = {start:0}
    todo = deque([start])
    while todo:
        u = todo.popleft()
        for v in adj[u]:
            # Reverse traversal of the nonincreasing directed graph.
            if v not in levels and (scores is None or scores[v] >= scores[u]):
                levels[v] = levels[u]+1
                todo.append(v)
    return levels


def cdict(counter):
    return {str(k):v for k, v in counter.items()}


stored = {}
geometry_cases = 0
for p in (3, 5, 7):
    grid = Grid(p)
    groups = defaultdict(list)
    counts = {}
    assisted = {}
    rho = (5, 1, 2, 3, 4, 0, 6) if p == 7 else None
    for pi in permutations(range(p)):
        r = normalize(pi)
        groups[r].append(pi)
        counts[pi] = grid.count(cycle_board(p, pi))
        geometry_cases += 1
        if p == 7:
            assisted[pi] = grid.count(cycle_board(p, pi, rho))
            geometry_cases += 1
        # Complete finite check of the modular slope lemma, every anchor.
        for i in range(p):
            slopes = Counter((pi[j]-pi[i])*pow((j-i) % p, -1, p) % p
                             for j in range(p) if i != j)
            require(max(slopes.values()) >= 2, "odd-prime modular matching at every point")
    require(len(groups) == factorial(p-2), "full coset universe")
    require(all(len(group) == p*(p-1) for group in groups.values()), "coset cardinality")
    reps = sorted(groups)
    adj = {r:{normalize(swap_values(r, a, b)) for a, b in combinations(range(p), 2)}-{r} for r in reps}
    start = tuple(range(p))
    levels = bfs(adj, start)
    qdata = data["quotients"][str(p)]
    require(cdict(Counter(levels.values())) == qdata["radius_profile"], "independent quotient radius")
    scores = {r:min(counts[pi] for pi in groups[r]) for r in reps}
    require(len(qdata["records"]) == len(reps), "complete quotient records")
    seen = set()
    for rec in qdata["records"]:
        r = tuple(rec["representative"])
        require(r in groups and r not in seen, "unique quotient record")
        seen.add(r)
        require(distance_formula(start, r) == levels[r] == rec["distance"], "cycle distance versus physical-value BFS")
        require(rec["minimum"] == scores[r], "native orbit minimum")
        require(rec["histogram"] == cdict(Counter(counts[pi] for pi in groups[r])), "complete native orbit histogram")
        masks = []
        for a in range(1, p):
            mask = 0
            for b in range(p):
                if counts[af(r, a, b)]:
                    mask |= 1 << b
            masks.append(mask)
        require(masks == rec["forbidden_masks"], "every saved affine forbidden mask")
    require(qdata["full_histogram"] == cdict(Counter(counts.values())), "full permutation histogram")
    require(qdata["minimum_profile"] == cdict(Counter(scores.values())), "orbit minimum profile")
    weak = [r for r in reps if scores[r] > 0 and all(scores[v] >= scores[r] for v in adj[r])]
    strict = [r for r in reps if scores[r] > 0 and all(scores[v] > scores[r] for v in adj[r])]
    require(qdata["weak_bad_minima"] == len(weak), "all weak minima")
    require(qdata["strict_bad_minima"] == len(strict), "all strict minima")
    # All roots of the quotient graph, not just the designated identity.
    for i, r in enumerate(reps):
        shifted = bfs(adj, r)
        require(Counter(shifted.values()) == Counter(levels.values()), "vertex-transitive radius control")
        target = reps[(i*31+7) % len(reps)]
        require(distance_formula(r, target) == shifted[target], "nonidentity coset orientation")
    require(max(levels.values()) <= max(0, p-4), "general compiler bound in finite universes")
    stored[p] = grid, groups, adj, counts, scores, assisted
    print(f"p={p}: {factorial(p)} column permutations; {len(reps)} cosets; radius {dict(sorted(Counter(levels.values()).items()))}.")

grid, groups, adj, counts, scores, assisted = stored[7]
require(Counter(counts.values())[0] == 0 and min(counts.values()) == 1, "p7 full fixed-row failure")
require(Counter(stored[5][3].values())[0] == 4, "p5 unrestricted positive control")
loc = data["strict_local_minimum"]
r = tuple(loc["representative"])
require(scores[r] == 1 and cdict(Counter(scores[v] for v in adj[r])) == loc["neighbors"], "strict local minimum neighborhood")
require(counts[af(r, *loc["attaining_affine"])] == 1, "strict local minimum actual realization")
repaired = data["row_assisted_compiler"]
rho = tuple(repaired["rho"])
scores2 = {r:min(assisted[pi] for pi in groups[r]) for r in groups}
good = [r for r in groups if scores2[r] == 0]
require(len(good) == 1, "unique row-assisted successful orbit")
monotone = bfs(adj, good[0], scores2)
require(Counter(monotone.values()) == Counter({0:1, 1:21, 2:91, 3:7}), "directed monotone distance profile")
require(cdict(Counter(monotone.values())) == repaired["monotone_distance_profile"], "saved monotone profile")
require(cdict(Counter(scores2.values())) == repaired["minimum_profile"], "row-assisted score profile")
stalls = {r for r in groups if scores2[r] and all(scores2[v] >= scores2[r] for v in adj[r])}
require(len(stalls) == 9 and stalls == {tuple(r) for r in repaired["strict_descent_stalls"]}, "all nine strict-descent stalls")
records = {tuple(rec["representative"]):rec for rec in repaired["records"]}
require(set(records) == set(groups) and len(repaired["records"]) == 120, "complete 120-record compiler")
bare_increases = 0
for r, rec in records.items():
    best = min((assisted[af(r, a, b)], a, b) for a in range(1, 7) for b in range(7))
    require(best == (rec["minimum"], *rec["best_affine"]), "literal row-assisted best affine pair")
    require(rec["monotone_distance"] == monotone[r], "saved directed distance")
    if monotone[r]:
        nxt = tuple(rec["successor"])
        current = af(r, *rec["best_affine"])
        i, j = rec["source_column_swap"]
        require(0 <= i < j < 7, "valid source swap")
        # Execute the physical exchange of current values, not source labels.
        intermediate = swap_values(current, current[i], current[j])
        final = af(intermediate, *rec["following_affine"])
        desired = af(nxt, *records[nxt]["best_affine"])
        require(final == desired and normalize(intermediate) == nxt, "concrete physical successor word")
        require(monotone[nxt] == monotone[r]-1 and assisted[final] <= assisted[current], "completed move monotonicity and progress")
        if assisted[intermediate] > assisted[current]:
            bare_increases += 1
        if r in stalls:
            require(scores2[nxt] == scores2[r], "required first plateau")
    # Follow every entire word independently to its literal successful endpoint.
    at = r
    for _ in range(monotone[r]):
        at = tuple(records[at]["successor"])
    require(at == good[0] and assisted[af(at, *records[at]["best_affine"]) ] == 0, "full terminating successful path")
require(bare_increases > 0, "bare intermediate is not universally monotone")
print(f"All 120 row-assisted compilers pass; {bare_increases} saved bare swaps increase X before their affine completion.")

repair = data["repair"]
pi = tuple(range(7))
for i, j in repair["column_swaps"]:
    pi = swap_values(pi, pi[i], pi[j])
pi = af(pi, *repair["affine"])
require(pi == tuple(repair["pi"]), "declared successful column word")
require(cycle_board(7, pi, rho) == [tuple(pt) for pt in repair["points"]], "declared exact point list")
require(assisted[pi] == 0 and counts[pi] == repair["before_row_swap_triples"] == 6, "load-bearing row addresses")
graph = defaultdict(set)
for x, y in cycle_board(7, pi, rho):
    graph[x].add(7+y)
    graph[7+y].add(x)
require(len(graph) == 14 and all(len(v) == 2 for v in graph.values()), "degree-two incidence")
require(len(bfs(graph, 0)) == 14, "one C14 incidence component")
words = [tuple(range(7))] + [swap_values(tuple(range(7)), i, j) for i, j in combinations(range(7), 2)]
small_bank = 0
for row in words:
    for col in words:
        for a in range(1, 7):
            for b in range(7):
                require(grid.unsafe(cycle_board(7, af(col, a, b), row)), "literal one-per-shore bank fails")
                small_bank += 1
require(small_bank == 20328, "explicit restricted operation universe")
print("All 20328 one-row/one-column/affine-column words fail; the recorded larger word succeeds.")

# Independent tangent search by growing diamonds, rather than a square search.
layer_cases = 0
safe_subsets = 0
for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
    for slope in range(1, p):
        chosen = None
        for radius in range(1, isqrt(2*p)+1):
            for u in range(-radius, radius+1):
                for v in {radius-abs(u), abs(u)-radius}:
                    if (v-slope*u) % p == 0:
                        chosen = (u, v)
                        break
                if chosen is not None:
                    break
            if chosen is not None:
                break
        require(chosen is not None, "diamond tangent radius")
        u, v = chosen
        r = abs(u)+abs(v)
        require(gcd(u, v) == 1 and r < p, "primitive shortest tangent")
        for b in range(p):
            pts = [(x, (slope*x+b) % p) for x in range(p)]
            layers = Counter(v*x-u*y for x, y in pts)
            cap = sum(min(2, n) for n in layers.values())
            require(len(layers) <= r and cap <= 2*r, "native layer capacity")
            require(p-cap == sum(max(0, n-2) for n in layers.values()), "typed final-line deletion deficit")
            parities = Counter((x % 2, y % 2) for x, y in pts)
            require(sum(comb(n, 2) for n in parities.values()) > 0, "zero-swap midpoint obstruction")
            if p <= 7:
                for mask in range(1 << p):
                    subset = [pts[i] for i in range(p) if (mask >> i) & 1]
                    if not any(det(*tri) == 0 for tri in combinations(subset, 3)):
                        require(len(subset) <= cap, "every small actual-safe retained subset")
                        safe_subsets += 1
            layer_cases += 1
require(layer_cases == data["layer_controls"], "complete inherited line-control universe")
print(f"{layer_cases} modular lines, {safe_subsets} actual-safe subsets, all native capacity controls pass.")

reciprocal_cases = 0
modular_false_positives = 0
for p in (3, 5, 7, 11, 13):
    inverses = [0] + [next(y for y in range(1, p) if x*y % p == 1) for x in range(1, p)]
    pts = list(enumerate(inverses))
    native = [tri for tri in combinations(pts, 3) if det(*tri) == 0]
    require(native == [((0, 0), (1, 1), (p-1, p-1))], "one native reciprocal triple")
    modular = [tri for tri in combinations(pts, 3) if det(*tri) % p == 0]
    require(len(modular) == (p-1)//2 and all((0, 0) in tri for tri in modular), "reciprocal modular classification")
    require(p == 3 or not any(all((y-a*x-b) % p == 0 for x, y in pts) for a in range(1, p) for b in range(p)), "no complete line in reciprocal graph")
    wins = []
    pdata = data["reciprocal"][str(p)]
    require([rec["a"] for rec in pdata["records"]] == list(range(1, p)), "complete zero-swap bank")
    for rec in pdata["records"]:
        a = rec["a"]
        b = inverses.index(a)
        output = list(enumerate(swap_values(tuple(inverses), 0, a)))
        P, Q = (0, a), (b, 0)
        literal = {tuple(tri) for tri in combinations(output, 3) if det(*tri) == 0}
        counts_by_type = [0, 0, 0]
        for tri in combinations(output, 3):
            actual = det(*tri) == 0
            if P in tri and Q in tri:
                x, y = next(pt for pt in tri if pt not in (P, Q))
                decoded = a*x+b*y == a*b
                typ = 2
            elif P in tri or Q in tri:
                isP = P in tri
                moved = P if isP else Q
                (x, y), (xx, yy) = [pt for pt in tri if pt != moved]
                if not isP:
                    x, y, xx, yy = y, x, yy, xx
                intercept = a if isP else b
                carry = (y+yy-intercept)//p
                decoded = ((y+yy == intercept and x*y == xx*yy) or
                           (y+yy == intercept+p and x*(p-y) == xx*(p-yy)))
                if det(*tri) % p == 0 and not actual:
                    modular_false_positives += 1
                require(not decoded or carry in (0, 1), "exact reciprocal carry range")
                typ = 0 if isP else 1
            else:
                decoded, typ = False, None
            require(decoded == actual, "literal individual reciprocal event iff")
            if actual:
                counts_by_type[typ] += 1
        require(counts_by_type == rec["counts"], "three disjoint event multiplicities")
        require(literal == {tuple(tuple(pt) for pt in tri) for tri in rec["triples"]}, "entire saved reciprocal triple set")
        if not literal:
            wins.append(a)
        reciprocal_cases += 1
    require(wins == pdata["successful_zero_swaps"] == {3:[1], 5:[1,4], 7:[1], 11:[], 13:[1]}[p], "exact successful zero-swap values")
    print(f"Reciprocal p={p}: successful zero-output swaps {wins}.")
require(modular_false_positives > 0, "native integer factors cannot be replaced by congruences")

# F4 multiplication table, independent of the producer's polynomial loop.
table = ((0,0,0,0), (0,1,2,3), (0,2,3,1), (0,3,1,2))
frobenius = (0,1,3,2)
require(all(table[x][x] == frobenius[x] for x in range(4)), "F4 Frobenius table")
for i, j, k in combinations(range(4), 3):
    require(table[i ^ j][frobenius[i] ^ frobenius[k]] != table[i ^ k][frobenius[i] ^ frobenius[j]], "odd-characteristic firewall")

print(f"Independent grid-line geometry cases: {geometry_cases}; reciprocal zero swaps: {reciprocal_cases}.")
print("All-prime proofs audited separately; conditional existence and completed-move monotonicity retained.")
print(f"PASS: {GATES} always-active exact gates.")
