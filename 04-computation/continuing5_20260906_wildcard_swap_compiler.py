"""Exact operation controls; standard library only, no repository imports."""
from collections import Counter, deque
from itertools import combinations, permutations
from math import factorial, gcd, isqrt
from pathlib import Path
import hashlib
import json
import sys

sys.stdout.reconfigure(newline="\n")
GATES = 0
BASE = Path(__file__).resolve().parent
DEST = BASE.parent / "05-knowledge" / "results" if BASE.name == "04-computation" else BASE


def gate(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def compose(a, b):
    return tuple(a[x] for x in b)


def inverse(a):
    inv = [0]*len(a)
    for i, j in enumerate(a):
        inv[j] = i
    return tuple(inv)


def affine(p, a, b):
    return tuple((a*i+b) % p for i in range(p))


def norm(pi):
    p = len(pi)
    a = pow((pi[1]-pi[0]) % p, -1, p)
    return tuple((y-pi[0])*a % p for y in pi)


def swap(pi, i, j):
    result = list(pi)
    result[i], result[j] = result[j], result[i]
    return tuple(result)


def cycles(pi):
    visited = set()
    result = 0
    for i in range(len(pi)):
        if i in visited:
            continue
        result += 1
        while i not in visited:
            visited.add(i)
            i = pi[i]
    return result


def det(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1])-(c[0]-a[0])*(b[1]-a[1])


def literal(points):
    return sum(det(*tri) == 0 for tri in combinations(points, 3))


def safe(points):
    return not any(det(*tri) == 0 for tri in combinations(points, 3))


def direction_count(points):
    total = 0
    for x, y in points:
        cs = Counter()
        for xx, yy in points:
            dx, dy = xx-x, yy-y
            if dx == dy == 0:
                continue
            g = gcd(dx, dy)
            dx, dy = dx//g, dy//g
            if dx < 0 or (dx == 0 and dy < 0):
                dx, dy = -dx, -dy
            cs[dx, dy] += 1
        total += sum(k*(k-1)//2 for k in cs.values())
    gate(total % 3 == 0, "anchored triple multiplicity three")
    return total//3


def board(p, pi, rho=None):
    if rho is None:
        rho = tuple(range(p))
    return [(rho[i], pi[j]) for i in range(p) for j in (i, (i+1) % p)]


def cycle_profile(points, p):
    adj = {i: [] for i in range(2*p)}
    for r, c in points:
        adj[r].append(p+c)
        adj[p+c].append(r)
    gate(all(len(v) == 2 for v in adj.values()), "saturated degrees retained")
    todo = set(adj)
    profile = []
    while todo:
        start = min(todo)
        seen = {start}
        stack = [start]
        for u in stack:
            for v in adj[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        todo -= seen
        profile.append(len(seen))
    return tuple(sorted(profile))


def interval_counts(points, p):
    """Native signed-lift decoder, independent of literal determinant outputs."""
    out = [[0]*p for _ in range(p)]
    for tri in combinations(points, 3):
        (x0,y0),(x1,y1),(x2,y2) = sorted(tri)
        if len({x0,x1,x2}) < 3 or len({y0,y1,y2}) < 3:
            continue
        if det(*tri) % p:
            continue
        u = (y1-y0)*pow(x1-x0, -1, p) % p
        g = gcd(x1-x0, x2-x0)
        span = (x2-x0)//g
        K = (p-1)//span
        inv = pow(u*g % p, -1, p)
        for k in range(-K, K+1):
            if not k:
                continue
            a = k*inv % p
            for initial in range(max(0, -k*span), min(p-1, p-1-k*span)+1):
                b = (initial-a*y0) % p
                out[a][b] += 1
    return out


def distance_by_cycles(r, affines):
    p = len(r)
    return min(p-cycles(compose(h, r)) for h in affines)


data = {"quotients": {}, "repair": {}, "reciprocal": {}}
for p in (3, 5, 7):
    ident = tuple(range(p))
    affines = [affine(p, a, b) for a in range(1, p) for b in range(p)]
    reps = [(0, 1)+r for r in permutations(range(2, p))]
    gate(len(reps) == factorial(p-2), "free affine left coset count")
    adj = {r: {norm(swap(r, i, j)) for i, j in combinations(range(p), 2)}-{r} for r in reps}
    dist = {ident: 0}
    queue = deque([ident])
    while queue:
        u = queue.popleft()
        for v in adj[u]:
            if v not in dist:
                dist[v] = dist[u]+1
                queue.append(v)
    gate(len(dist) == len(reps), "connected quotient graph")
    expected = {3: {0: 1}, 5: {0: 1, 1: 5}, 7: {0: 1, 1: 21, 2: 98}}[p]
    gate(dict(Counter(dist.values())) == expected, "exact radius profile")
    for r in reps:
        gate(dist[r] == distance_by_cycles(r, affines), "BFS versus minimum cycle formula")
        if p >= 5:
            gate(dist[r] <= p-4, "classical modular triple compiler bound")
        if p >= 7:
            gate(len(adj[r]) == p*(p-1)//2, "all transposition neighbors distinct")
        # Every point has a repeated modular slope; product is 1, never -1.
        for i in range(p):
            slopes = [(r[j]-r[i])*pow((j-i) % p, -1, p) % p for j in range(p) if i != j]
            product = 1
            for slope in slopes:
                product = product*slope % p
            gate(product == 1 and len(set(slopes)) < p-1, "modular slope product / repeated direction")
    score = {}
    records = []
    full_hist = Counter()
    for r in reps:
        initial_board = board(p, r)
        prediction = interval_counts(initial_board, p)
        hist = Counter()
        for a in range(1, p):
            for b in range(p):
                pi = compose(affine(p, a, b), r)
                points = board(p, pi)
                count = literal(points)
                gate(count == prediction[a][b], "exact affine event decoder")
                gate(count == direction_count(points), "independent primitive-direction triple count")
                gate(cycle_profile(points, p) == (2*p,), "incidence cycle unchanged")
                hist[count] += 1
        score[r] = min(hist)
        full_hist.update(hist)
        records.append({"representative": r, "distance": dist[r], "minimum": min(hist),
                        "histogram": dict(sorted(hist.items())),
                        "forbidden_masks": [sum(1 << b for b in range(p) if prediction[a][b]) for a in range(1, p)]})
    gate(sum(full_hist.values()) == factorial(p), "complete labelled column universe")
    if p == 5:
        gate(full_hist[0] == 4 and all(len(a) == 5 for a in adj.values()), "p5 K6 and four successes")
    if p == 7:
        gate(full_hist[0] == 0 and min(full_hist) == 1, "all p7 fixed-row columns fail")
        gate(all(all(mask == (1 << p)-1 for mask in row["forbidden_masks"]) for row in records),
             "all 120 exact affine interval unions cover")
        r = (0,1,2,5,3,6,4)
        gate(score[r] == 1, "strict bad local minimum")
        neighbor_profile = dict(sorted(Counter(score[v] for v in adj[r]).items()))
        gate(neighbor_profile == {2: 8, 3: 7, 4: 3, 5: 3}, "every adjacent orbit worsens")
        data["strict_local_minimum"] = {"representative": r, "neighbors": neighbor_profile,
                                         "attaining_affine": [3,5]}
        gate(literal(board(7, compose(affine(7,3,5), r))) == 1, "actual representative attains local minimum")
    weak_minima = [r for r in reps if score[r] > 0 and all(score[w] >= score[r] for w in adj[r])]
    strict_minima = [r for r in weak_minima if all(score[w] > score[r] for w in adj[r])]
    if p == 7:
        gate(len(weak_minima) == 22 and len(strict_minima) == 2, "complete quotient local minima")
    data["quotients"][str(p)] = {"radius_profile": expected, "records": records,
                                "full_histogram": dict(sorted(full_hist.items())),
                                "minimum_profile": dict(sorted(Counter(score.values()).items())),
                                "weak_bad_minima": len(weak_minima), "strict_bad_minima": len(strict_minima)}
    print("p="+str(p)+": quotient distance profile "+str(expected)+"; full successes "+str(full_hist[0]))

# Exact p7 positive control after one row move, with a certified two-column word.
p = 7
rho = swap(tuple(range(p)), 0, 5)
word = swap(swap(tuple(range(p)), 0, 1), 1, 2)
pi = compose(affine(p, 6, 3), word)
gate(pi == (2,1,3,0,6,5,4), "physical two-column word and final affine")
points = board(p, pi, rho)
gate(safe(points) and direction_count(points) == 0, "row-assisted actual no-three repair")
gate(cycle_profile(points, p) == (14,), "repair retains full C14 skeleton")
gate(literal(board(p, pi)) > 0, "row-address change is material")
data["repair"] = {"row_swap": [0,5], "column_swaps": [[0,1],[1,2]], "affine": [6,3],
                   "pi": pi, "points": points, "before_row_swap_triples": literal(board(p, pi))}

# Complete constructive row class: permit equal-score moves and retain an actual lift.
family_score = {}
family_best = {}
for r in reps:
    prediction = interval_counts(board(p,r,rho),p)
    best = min((prediction[a][b],a,b) for a in range(1,p) for b in range(p))
    family_score[r],aa,bb = best
    family_best[r] = (aa,bb)
    for a in range(1,p):
        for b in range(p):
            pts = board(p,compose(affine(p,a,b),r),rho)
            count = literal(pts)
            gate(count == prediction[a][b], "row-assisted full affine decoder")
            gate(direction_count(pts) == count, "row-assisted independent native directions")
good = [r for r in reps if family_score[r] == 0]
gate(len(good) == 1, "row-assisted unique successful affine orbit")
monotone_dist = {r:0 for r in good}
successor = {}
queue = deque(good)
while queue:
    u = queue.popleft()
    for v in sorted(adj[u]):
        if v not in monotone_dist and family_score[v] >= family_score[u]:
            monotone_dist[v] = monotone_dist[u]+1
            successor[v] = u
            queue.append(v)
gate(dict(Counter(monotone_dist.values())) == {0:1,1:21,2:91,3:7}, "complete monotone compiler distance profile")
weak = [r for r in reps if family_score[r] > 0 and all(family_score[w] >= family_score[r] for w in adj[r])]
gate(len(weak) == 9, "strict descent stalls at nine bad states")
family_records = []
for r in reps:
    record = {"representative":r, "minimum":family_score[r], "best_affine":family_best[r],
              "monotone_distance":monotone_dist[r]}
    if r in successor:
        nxt = successor[r]
        i,j = next((i,j) for i,j in combinations(range(p),2) if norm(swap(r,i,j)) == nxt)
        current = compose(affine(p,*family_best[r]),r)
        desired = compose(affine(p,*family_best[nxt]),nxt)
        moved = swap(current,i,j)
        a = (desired[1]-desired[0])*pow((moved[1]-moved[0]) % p,-1,p) % p
        b = (desired[0]-a*moved[0]) % p
        gate(compose(affine(p,a,b),moved) == desired, "successor has an actual one-swap affine lift")
        gate(literal(board(p,desired,rho)) <= literal(board(p,current,rho)), "actual successor integer triples do not increase")
        gate(monotone_dist[nxt] == monotone_dist[r]-1, "successor terminates at a true safe board")
        record.update({"successor":nxt,"source_column_swap":[i,j],"following_affine":[a,b]})
    family_records.append(record)
data["row_assisted_compiler"] = {"rho":rho,"records":family_records,
                                  "minimum_profile":dict(sorted(Counter(family_score.values()).items())),
                                  "monotone_distance_profile":dict(sorted(Counter(monotone_dist.values()).items())),
                                  "strict_descent_stalls":weak}
print("p7 rho=(0 5): all 120 orbits have a nonincreasing actual-triple compiler with <=3 swaps; nine need a plateau move.")
row_words = [tuple(range(p))] + [swap(tuple(range(p)), i, j) for i,j in combinations(range(p),2)]
column_words = row_words
control_count = 0
for rr in row_words:
    for cc in column_words:
        for a in range(1,p):
            for b in range(p):
                gate(not safe(board(p, compose(affine(p,a,b), cc), rr)), "at most one row and one column swap never sufficient")
                control_count += 1
gate(control_count == 20328, "complete small repair-budget universe")
print("p7 row swap (0,5), column swaps (0,1),(1,2), affine (6,3): exact success; smaller one-per-shore box empty.")

# Lattice-layer sidecar underlying the all-prime support barrier.
layer_controls = 0
for p in (5,7,11,13,17,19,23,29,31):
    h = isqrt(2*p)
    for slope in range(1,p):
        candidates = [(abs(u)+abs(v),u,v) for u in range(-h,h+1) for v in range(-h,h+1)
                      if (u or v) and (v-slope*u) % p == 0]
        r,u,v = min(candidates)
        gate(r <= h and gcd(u,v) == 1, "short primitive tangent direction")
        for intercept in range(p):
            line = [(x,(slope*x+intercept) % p) for x in range(p)]
            layers = Counter(v*x-u*y for x,y in line)
            capacity = sum(min(2,k) for k in layers.values())
            gate(len(layers) <= r and capacity <= 2*r, "exact occupied-layer capacity")
            gate(sum(max(0,k-2) for k in layers.values()) == p-capacity, "exact directional deletion deficit")
            if p <= 7:
                # All subsets, including hostile ones: layer feasibility is necessary only.
                for mask in range(1 << p):
                    retained = [line[j] for j in range(p) if mask >> j & 1]
                    if safe(retained):
                        gate(len(retained) <= capacity, "every actual safe subset respects capacity")
            layer_controls += 1
    print("p="+str(p)+": all-prime total-swap lower bound "+str(max(1,(p-2*h+1)//2))+" for any board containing a full modular line.")
data["layer_controls"] = layer_controls

# A completed reciprocal graph evades the full-line obstruction but needs exact carries.
def reciprocal_points(p, a=None):
    inv = [0]+[pow(i,-1,p) for i in range(1,p)]
    if a is None:
        return list(enumerate(inv))
    b = pow(a,-1,p)
    inv[0], inv[b] = inv[b], inv[0]
    return list(enumerate(inv))


def reciprocal_packet(p, a):
    b = pow(a,-1,p)
    retained = [(x,pow(x,-1,p)) for x in range(1,p) if x != b]
    through_P, through_Q, through_both = [], [], []
    for (x,y),(xx,yy) in combinations(retained,2):
        if (y+yy == a and x*y == xx*yy) or (y+yy == a+p and x*(p-y) == xx*(p-yy)):
            through_P.append(((0,a),(x,y),(xx,yy)))
        if (x+xx == b and x*y == xx*yy) or (x+xx == b+p and y*(p-x) == yy*(p-xx)):
            through_Q.append(((b,0),(x,y),(xx,yy)))
    for x,y in retained:
        if a*x+b*y == a*b:
            through_both.append(((0,a),(b,0),(x,y)))
    return through_P, through_Q, through_both


for p in (3,5,7,11,13):
    origin = reciprocal_points(p)
    exact_bad = [tri for tri in combinations(origin,3) if det(*tri) == 0]
    gate(exact_bad == [((0,0),(1,1),(p-1,p-1))], "one native reciprocal triple")
    modular = [tri for tri in combinations(origin,3) if det(*tri) % p == 0]
    gate(len(modular) == (p-1)//2 and all((0,0) in tri for tri in modular), "classical reciprocal modular triple family")
    records = []
    wins = []
    for a in range(1,p):
        output = reciprocal_points(p,a)
        actual = {tuple(sorted(tri)) for tri in combinations(output,3) if det(*tri) == 0}
        packets = reciprocal_packet(p,a)
        decoded = {tuple(sorted(tri)) for part in packets for tri in part}
        gate(actual == decoded, "factor-layer decoder iff for every individual triple")
        gate(sum(map(len,packets)) == len(actual), "three disjoint event types")
        gate(direction_count(output) == len(actual), "reciprocal independent native directions")
        if not actual:
            wins.append(a)
        records.append({"a": a, "counts": list(map(len,packets)), "triples": sorted(actual)})
    gate(wins == {3:[1],5:[1,4],7:[1],11:[],13:[1]}[p], "complete bounded reciprocal swap bank")
    data["reciprocal"][str(p)] = {"successful_zero_swaps": wins, "records": records}
    print("reciprocal p="+str(p)+": successful swaps of output0 with a: "+str(wins))

# Characteristic two hostile to the classical slope-product step.
def mul4(a,b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        a <<= 1
        if a & 4:
            a ^= 7
        b >>= 1
    return out
f4 = tuple(mul4(a,a) for a in range(4))
gate(sorted(f4) == list(range(4)), "F4 Frobenius is a permutation")
for i,j,k in combinations(range(4),3):
    gate(mul4(j^i,f4[k]^f4[i]) != mul4(k^i,f4[j]^f4[i]), "F4 graph has no modular collinear triple")

cert_path = DEST / (Path(__file__).stem+"_certificate.json")
cert_path.write_text(json.dumps(data,sort_keys=True,indent=2)+"\n",encoding="utf-8",newline="\n")
print("Certificate SHA256: "+hashlib.sha256(cert_path.read_bytes()).hexdigest())
print("No unrestricted all-prime success or algorithmic running-time claim.")
print("PASS: "+str(GATES)+" always-active exact gates.")
