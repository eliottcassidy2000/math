#!/usr/bin/env python3
"""Exact bounded controls for the third-direction repair theorem. No packages."""
from collections import Counter, defaultdict, deque
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, permutations, product
import hashlib
import json
import math
from pathlib import Path
import sys

sys.stdout.reconfigure(newline="\n")
STEM = Path(__file__).stem
GATES = Counter()


def check(name, condition):
    GATES[name] += 1
    if not condition:
        raise RuntimeError(name)


def boards(n):
    pairs = tuple(combinations(range(n), 2))
    degree = [0] * n
    rows = []

    def visit(r):
        if r == n:
            if all(d == 2 for d in degree):
                yield tuple(rows)
            return
        for pair in pairs:
            if any(degree[c] == 2 for c in pair):
                continue
            for c in pair:
                degree[c] += 1
            if all(2 - d <= n - r - 1 for d in degree):
                rows.append(pair)
                yield from visit(r + 1)
                rows.pop()
            for c in pair:
                degree[c] -= 1
    yield from visit(0)


def cells(rows):
    return tuple((r, c) for r, pair in enumerate(rows) for c in pair)


def families(points):
    answer = []
    for slope in (1, -1, 2):
        line = defaultdict(int)
        for i, (r, c) in enumerate(points):
            line[c - slope * r] |= 1 << i
        answer.append(dict(line))
    return answer


def exact_deletions(masks):
    triples = set()
    for mask in masks:
        bits = [1 << i for i in range(mask.bit_length()) if mask >> i & 1]
        triples.update(sum(t) for t in combinations(bits, 3))

    @lru_cache(None)
    def hit(remaining):
        if not remaining:
            return 0
        first = remaining[0]
        best = None
        while first:
            bit = first & -first
            first -= bit
            candidate = bit | hit(tuple(t for t in remaining if not t & bit))
            if best is None or candidate.bit_count() < best.bit_count():
                best = candidate
        return best
    deletion = hit(tuple(sorted(triples)))
    return deletion.bit_count(), deletion


def boolean_dual(fams):
    active = [m for f in fams for m in f.values() if m.bit_count() > 2]
    best = 0
    for selected in range(1 << len(active)):
        union = 0
        number = 0
        for j, mask in enumerate(active):
            if selected >> j & 1:
                union |= mask
                number += 1
        best = max(best, union.bit_count() - 2 * number)
    return best


def compiled_dual(fams):
    active = [m for f in fams[:2] for m in f.values() if m.bit_count() > 2]
    best = 0
    for selected in range(1 << len(active)):
        union = 0
        number = 0
        for j, mask in enumerate(active):
            if selected >> j & 1:
                union |= mask
                number += 1
        credit = sum(max(0, (m & ~union).bit_count() - 2) for m in fams[2].values())
        best = max(best, union.bit_count() - 2 * number + credit)
    return best


def flow_deletions(points):
    # Independent augmenting paths, with each original cell represented by a
    # unit-capacity edge; line capacities alone would be an invalid relaxation.
    ds = sorted({c-r for r, c in points})
    ss = sorted({c+r for r, c in points})
    dn = {v: i+1 for i, v in enumerate(ds)}
    sn = {v: i+1+len(ds) for i, v in enumerate(ss)}
    sink = 1+len(ds)+len(ss)
    cap = [[0]*(sink+1) for _ in range(sink+1)]
    for d in ds:
        cap[0][dn[d]] = 2
    for s in ss:
        cap[sn[s]][sink] = 2
    for r, c in points:
        cap[dn[c-r]][sn[c+r]] = 1
    total = 0
    while True:
        prev = {0: None}
        queue = deque([0])
        while queue and sink not in prev:
            u = queue.popleft()
            for v, capacity in enumerate(cap[u]):
                if capacity and v not in prev:
                    prev[v] = u
                    queue.append(v)
        if sink not in prev:
            return len(points)-total
        v = sink
        while v:
            u = prev[v]
            cap[u][v] -= 1
            cap[v][u] += 1
            v = u
        total += 1


def inspect(rows, certify=True):
    points = cells(rows)
    fams = families(points)
    old = [m for f in fams[:2] for m in f.values()]
    tau2, delete2 = exact_deletions(old)
    tau3, delete3 = exact_deletions(old + list(fams[2].values()))
    active = 0
    for mask in old:
        if mask.bit_count() > 2:
            active |= mask
    safe = ((1 << len(points))-1) & ~active
    bonus = sum(max(0, (m & safe).bit_count()-2) for m in fams[2].values())
    isolated = sum(m.bit_count() == 3 and all(
        fams[k][points[i][1] - slope*points[i][0]].bit_count() == 1
        for i in range(len(points)) if m >> i & 1
        for k, slope in enumerate((1, -1))) for m in fams[2].values())
    beta2 = boolean_dual(fams[:2])
    beta3 = boolean_dual(fams)
    if certify:
        check("unit_flow_vs_exact", flow_deletions(points) == tau2)
        check("two_direction_duality", beta2 == tau2)
        check("three_direction_compiler", compiled_dual(fams) == beta3)
        check("independent_bonus", tau3 >= beta3 >= tau2+bonus >= tau2+isolated)
    return dict(rows=rows, tau2=tau2, tau3=tau3, J3=bonus, I3=isolated,
                beta3=beta3, delete2=delete2, delete3=delete3,
                slope2_excess=sum(max(0,m.bit_count()-2) for m in fams[2].values()),
                overfull=[[[key, m] for key, m in f.items() if m.bit_count()>2] for f in fams])


counts = {}
for n in range(2, 6):
    stats = Counter()
    for rows in boards(n):
        data = inspect(rows)
        stats["boards"] += 1
        stats["strict_gain"] += data["tau3"] > data["tau2"]
        stats["isolated_positive"] += data["I3"] > 0
    counts[n] = dict(stats)
    check("complete_census_cardinality", stats["boards"] == {2:1,3:6,4:90,5:2040}[n])
    check("no_isolated_before_six", stats["isolated_positive"] == 0)
check("complete_n5_strict_gain", counts[5]["strict_gain"] == 142)

named_rows = {
    "naive_addition_refuted": ((0,1),(0,1),(2,3),(2,4),(3,4)),
    "dual_beats_safe_bonus": ((0,1),(0,3),(2,4),(3,4),(1,2)),
    "safe_bonus_two": ((0,1),(2,3),(0,4),(1,2),(3,4)),
    "triangle_integrality_gap": ((0,3),(1,2),(0,4),(3,4),(1,2)),
    "minimal_isolated": ((0,1),(1,3),(0,5),(2,3),(2,4),(4,5)),
}
witnesses = {name: inspect(rows) for name, rows in named_rows.items()}
check("naive_addition_hostile", witnesses["naive_addition_refuted"]["tau3"] == 4 <
      witnesses["naive_addition_refuted"]["tau2"]+witnesses["naive_addition_refuted"]["slope2_excess"])
check("stronger_dual_witness", tuple(witnesses["dual_beats_safe_bonus"][k] for k in ("tau2","tau3","J3","beta3")) == (1,2,0,2))
check("safe_two_witness", tuple(witnesses["safe_bonus_two"][k] for k in ("tau2","tau3","J3")) == (0,2,2))
check("isolated_witness", tuple(witnesses["minimal_isolated"][k] for k in ("tau2","tau3","J3","I3")) == (3,4,1,1))
triangle = witnesses["triangle_integrality_gap"]
check("triangle_integer_values", (triangle["beta3"],triangle["tau3"]) == (1,2))
constraints = [mask for family in triangle["overfull"] for _, mask in family]
check("triangle_three_triples", len(constraints) == 3 and all(m.bit_count()==3 for m in constraints))
fractional = [Q(1,2) if i in (0,5,6) else Q(0) for i in range(10)]
check("triangle_fractional_primal", sum(fractional)==Q(3,2) and all(sum(fractional[i] for i in range(10) if m>>i&1)>=1 for m in constraints))
check("triangle_fractional_dual", all(sum(Q(1,2) for m in constraints if m>>i&1)<=1 for i in range(10)))
matrix = [[int(m>>i&1) for i in (0,5,6)] for m in constraints]
det = sum(matrix[0][j] * (matrix[1][(j+1)%3]*matrix[2][(j+2)%3]-matrix[1][(j+2)%3]*matrix[2][(j+1)%3]) for j in range(3))
check("nonunimodular_minor", abs(det)==2)

# All single row/column transpositions of the five explicit boards.
for rows in named_rows.values():
    points = set(cells(rows))
    baseline = inspect(rows, False)["tau3"]
    for i,j in combinations(range(len(rows)),2):
        for axis in (0,1):
            transformed = set()
            for point in points:
                p = list(point)
                if p[axis] == i: p[axis] = j
                elif p[axis] == j: p[axis] = i
                transformed.add(tuple(p))
            fams = families(tuple(sorted(transformed)))
            value,_ = exact_deletions([m for f in fams for m in f.values()])
            check("transposition_lipschitz", abs(value-baseline)<=len(points-transformed)<=4)


def physical_skeleton(n, kind, order):
    native = [((r,r^1) if kind=="squares" else (r,(r+1)%n)) for r in range(n)]
    return tuple(tuple(sorted(native[old])) for old in order)


def eligible(target_rows, rowcols, colrows):
    targets = set(target_rows)
    for r in target_rows:
        for column in rowcols[r]:
            v = next(x for x in colrows[column] if x != r)
            bad = {v, 2*r-v}
            if (2*r+v)%3 == 0:
                bad.add((2*r+v)//3)
            if (targets-{r}) & bad:
                return False
    return True


def seven_hit(r,c,q,target_rows):
    return c-2*r==q or any(c-r==q+t or c+r==q+3*t for t in target_rows)


minimal_rows = named_rows["minimal_isolated"]
minimal_columns = {c:tuple(r for r in range(6) if c in minimal_rows[r]) for c in range(6)}
check("realized_event_need_not_be_eligible", not eligible((0,1,2),minimal_rows,minimal_columns))
event_stats = Counter()
for n in (6,8,10,12):
    orders = (tuple(range(n)), tuple(range(0,n,2))+tuple(range(1,n,2)))
    for kind, order in product(("cycle","squares"),orders):
        rowcols = physical_skeleton(n,kind,order)
        colrows = {c:tuple(r for r in range(n) if c in rowcols[r]) for c in range(n)}
        event_stats["skeleton_row_orders"] += 1
        checked_choices = 0
        for q in range(-2*(n-1), n):
            line_rows = tuple(r for r in range(n) if 0<=q+2*r<n)
            bad_count = 0
            for targets in combinations(line_rows,3):
                if not eligible(targets,rowcols,colrows):
                    bad_count += 1
                    continue
                event_stats["eligible_triples"] += 1
                # Every eligible triple is checked for all eight forced choices.
                for chosen in product(*(rowcols[r] for r in targets)):
                    labels = tuple(q+2*r for r in targets)
                    check("forced_distinct", len(set(chosen))==3)
                    forced = dict(zip(chosen,labels))
                    check("companions_avoid", all(not seven_hit(v,forced[c],q,targets)
                        for r,c in zip(targets,chosen) for v in colrows[c] if v!=r))
                    event_stats["forced_choices"] += 1
                    # Bounded full forbidden-matrix bank per skeleton/order.
                    if checked_choices >= 16:
                        continue
                    checked_choices += 1
                    source_left = [c for c in range(n) if c not in forced]
                    labels_left = [c for c in range(n) if c not in labels]
                    forbidden = {(c,y) for c in source_left for y in labels_left
                                 if any(seven_hit(r,y,q,targets) for r in colrows[c])}
                    total_length = sum(0<=q+2*r<n for r in range(n))
                    total_length += sum(sum(c-r==q+t for r in range(n) for c in range(n))
                                      +sum(c+r==q+3*t for r in range(n) for c in range(n)) for t in targets)
                    check("forbidden_row_degree", all(sum(c==a for a,b in forbidden)<=14 for c in source_left))
                    check("forbidden_column_degree", all(sum(y==b for a,b in forbidden)<=14 for y in labels_left))
                    check("forbidden_size", len(forbidden)<=2*total_length<=13*n+1)
                    event_stats["full_forbidden_matrices"] += 1
                    if n<=8:
                        for image in permutations(labels_left):
                            assignment = dict(forced)
                            assignment.update(zip(source_left,image))
                            avoid = all((c,assignment[c]) not in forbidden for c in source_left)
                            board = tuple((r,assignment[c]) for r in range(n) for c in rowcols[r])
                            target_points = {(r,q+2*r) for r in targets}
                            literal = {p for p in board if seven_hit(*p,q,targets)} == target_points
                            check("conditional_event_equivalence", avoid==literal)
                            event_stats["conditional_permutations"] += 1
            check("eligibility_uniform_count", bad_count <= 6*len(line_rows)*max(0,len(line_rows)-2))


def pmul(a,b):
    out = [Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j] += x*y
    return out


def power(p,k):
    out = [Q(1)]
    for _ in range(k): out = pmul(out,p)
    return out


def integral(p,a,b):
    return sum(x*(b**(i+1)-a**(i+1))/Q(i+1) for i,x in enumerate(p))


def abs_integral(a,b,lo,hi):
    cuts = [lo,hi]
    root = -b/a
    if lo<root<hi: cuts.append(root)
    cuts.sort()
    return sum(abs(a*(v*v-u*u)/2+b*(v-u)) for u,v in zip(cuts,cuts[1:]))


pieces = [(-2,-1,[Q(1),Q(1,2)],[Q(5,3),Q(7,6),Q(1,6)]),
          (-1,0,[Q(1,2)],[Q(2,3),Q(-1,3),Q(-1,3)]),
          (0,1,[Q(1,2),Q(-1,2)],[Q(2,3),Q(-5,6),Q(1,6)])]
L3 = sum(integral(power(L,3),Q(a),Q(b)) for a,b,L,J in pieces)
L4 = sum(integral(power(L,4),Q(a),Q(b)) for a,b,L,J in pieces)
L2J = sum(integral(pmul(power(L,2),J),Q(a),Q(b)) for a,b,L,J in pieces)
check("line_length_integrals", (L3,L4,L2J)==(Q(3,16),Q(7,80),Q(187,720)))
average_cost = 2*(L4+3*L2J)/L3
check("jensen_cost", average_cost == Q(416,45))
for numerator in range(-120,61):
    q = Q(numerator,60)
    lo,hi = max(Q(0),-q/2),min(Q(1),(1-q)/2)
    actual = 2*(hi-lo)-abs_integral(Q(1),q,lo,hi)-abs_integral(Q(3),q-1,lo,hi)
    J = next(J for a,b,L,J in pieces if a<=q<=b)
    predicted = sum(c*q**i for i,c in enumerate(J))
    check("independent_absolute_integral", actual==predicted)
for n in range(1,81):
    count_lines = sum(math.comb(sum(0<=q+2*r<n for r in range(n)),3)
                      for q in range(-2*(n-1),n))
    count_spans = sum((n-h)*(h-1)*(n-2*h) for h in range(2,(n-1)//2+1))
    check("independent_triple_count", count_lines==count_spans)

certificate = dict(status="FINITE-EXACT bounded controls; analytic theorem proved in paired report",
    universe="All 2137 simple two-regular boards for n=2..5; five named controls (one at n=6); 16 skeleton/row-order event banks",
    census=counts, witnesses=witnesses, event_bank=dict(event_stats), gates=dict(GATES),
    gate_total=sum(GATES.values()), triangle_lp="3/2", triangle_integer="2",
    integral_values={"L3":str(L3),"L4":str(L4),"L2J":str(L2J),"cost":str(average_cost)},
    delta3_formula="exp(-416/45)/4", delta3_decimal=math.exp(-416/45)/4,
    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
outdir = Path(__file__).parent
if outdir.name == "04-computation":
    outdir = outdir.parent/"05-knowledge"/"results"
output = outdir/(STEM+"_certificate.json")
output.write_bytes((json.dumps(certificate,indent=2,sort_keys=True)+"\n").encode("utf-8"))
print("PASS continuing11 third-direction repair controls")
print("Complete census n=2..5:", json.dumps(counts,sort_keys=True))
print("Event bank:", json.dumps(dict(event_stats),sort_keys=True))
print("Triangle: Boolean bound 1; fractional optimum 3/2; exact repair 2")
print("Uniform extra mean coefficient: exp(-416/45)/4 = %.17g" % certificate["delta3_decimal"])
print("Exact integrals:", json.dumps(certificate["integral_values"],sort_keys=True))
print("Gates:", sum(GATES.values()))
print("Certificate:", output.name)
