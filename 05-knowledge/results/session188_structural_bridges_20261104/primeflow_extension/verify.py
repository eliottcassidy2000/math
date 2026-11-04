"""Finite-lattice encoding and coupled arithmetic certificate checks.

Standard library only; this script imports no sibling implementation.
"""
import hashlib
import itertools
import json
import random
from math import isqrt, lcm, prod


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def val(n, r):
    need(n > 0, "positive valuation argument required")
    k = 0
    while n % r == 0:
        n //= r
        k += 1
    return k


def order(p, r):
    a = 1
    for k in range(1, r):
        a = a*p % r
        if a == 1:
            return k
    raise RuntimeError("bad order")


def direct(p, n, targets):
    value = (p**n-1)//(p-1)
    return tuple(val(value, r) for r in targets)


class Source:
    def __init__(self, p, targets, euler=False):
        need(prime(p) and p % 2, "odd source prime")
        need(not euler or p % 4 == 1, "Euler role")
        need(all(prime(r) and r % 2 and r != p for r in targets), "target primes")
        self.p, self.targets, self.euler = p, tuple(targets), euler
        self.d = tuple(order(p, r) for r in targets)
        self.c = tuple(0 if d == 1 else val(p**d-1, r) for d, r in zip(self.d, targets))

    def compile(self, lower):
        length = 2 if self.euler else 1
        for r, d, c, k in zip(self.targets, self.d, self.c, lower):
            if not k:
                continue
            if (self.euler and d % 4 == 0) or (not self.euler and d % 2 == 0):
                return None
            rank = r**k if d == 1 else d*r**max(0, k-c)
            length = lcm(length, rank)
        profile = tuple(val(length, r) if d == 1 else c+val(length, r) if length % d == 0 else 0
                        for r, d, c in zip(self.targets, self.d, self.c))
        return length, profile

    def domain(self, lower, upper):
        result = []
        for row in itertools.product(*(range(a, b+1) for a, b in zip(lower, upper))):
            compiled = self.compile(row)
            if compiled is None or compiled[1] != row:
                continue
            witness = compiled[0]
            if witness == 1:
                witness = max((*self.targets, *self.d, 2))+1
                while not prime(witness):
                    witness += 1
            need(direct(self.p, witness, self.targets) == row, "domain witness")
            result.append(row)
        return result


def propagate(sources, initial_lower, initial_upper, totals):
    lower, upper = [list(r) for r in initial_lower], [list(r) for r in initial_upper]
    trace = []
    rows, columns = len(sources), len(totals)
    while True:
        previous = ([r[:] for r in lower], [r[:] for r in upper])
        for i, source in enumerate(sources):
            result = source.compile(lower[i])
            if result is None:
                trace.append({"rule": "disallowed_row", "source": source.p})
                return False, lower, upper, trace
            length, closed = result
            if any(a > b for a, b in zip(closed, upper[i])):
                trace.append({"rule": "row_overflow", "source": source.p,
                              "lower": lower[i][:], "upper": upper[i][:],
                              "length": length, "closure": list(closed)})
                return False, lower, upper, trace
            if tuple(lower[i]) != closed:
                trace.append({"rule": "row_closure", "source": source.p,
                              "old_lower": lower[i][:], "new_lower": list(closed),
                              "length": length})
                lower[i] = list(closed)
            # Exact row coordinate maxima (Step 5), using monotone thresholds.
            for j in range(columns):
                lo, hi = lower[i][j], upper[i][j]
                while lo < hi:
                    midpoint = (lo+hi+1)//2
                    demand = lower[i][:]
                    demand[j] = midpoint
                    completion = source.compile(demand)
                    feasible = completion is not None and all(a <= b for a, b in zip(completion[1], upper[i]))
                    if feasible:
                        lo = midpoint
                    else:
                        hi = midpoint-1
                if lo != upper[i][j]:
                    trace.append({"rule": "row_upper_hull", "source": source.p,
                                  "target": source.targets[j], "old": upper[i][j], "new": lo})
                    upper[i][j] = lo
        for j, total in enumerate(totals):
            low_sum, high_sum = sum(r[j] for r in lower), sum(r[j] for r in upper)
            if not low_sum <= total <= high_sum:
                trace.append({"rule": "column_overflow", "column": j,
                              "total": total, "lower_sum": low_sum, "upper_sum": high_sum})
                return False, lower, upper, trace
            for i in range(rows):
                lo = max(lower[i][j], total-sum(upper[k][j] for k in range(rows) if k != i))
                hi = min(upper[i][j], total-sum(lower[k][j] for k in range(rows) if k != i))
                need(lo <= hi, "column update unexpectedly crossed after range check")
                if (lo, hi) != (lower[i][j], upper[i][j]):
                    trace.append({"rule": "column_residual", "column": j,
                                  "source": sources[i].p, "total": total,
                                  "old": [lower[i][j], upper[i][j]], "new": [lo, hi]})
                    lower[i][j], upper[i][j] = lo, hi
        if previous == (lower, upper):
            return True, lower, upper, trace
        old_width = sum(b-a for row_a, row_b in zip(*previous) for a, b in zip(row_a, row_b))
        new_width = sum(b-a for row_a, row_b in zip(lower, upper) for a, b in zip(row_a, row_b))
        need(new_width < old_width, "propagation did not progress")


def lattice_checks():
    primes = [3, 5, 7, 11, 13, 17]
    counts, pair_checks, interval_checks = {}, 0, 0
    digest = hashlib.sha256()
    for n in range(1, 7):
        count = 0
        pairs = list(itertools.combinations(range(n), 2))
        for mask in range(1 << len(pairs)):
            le = [[i == j for j in range(n)] for i in range(n)]
            for bit, (i, j) in enumerate(pairs):
                if mask >> bit & 1:
                    le[i][j] = True
            if any(le[i][j] and le[j][k] and not le[i][k]
                   for i in range(n) for j in range(n) for k in range(n)):
                continue
            meet, join = {}, {}
            for x in range(n):
                for y in range(n):
                    lows = [z for z in range(n) if le[z][x] and le[z][y]]
                    highs = [z for z in range(n) if le[x][z] and le[y][z]]
                    best_lows = [z for z in lows if all(le[t][z] for t in lows)]
                    best_highs = [z for z in highs if all(le[z][t] for t in highs)]
                    if len(best_lows) == len(best_highs) == 1:
                        meet[x, y], join[x, y] = best_lows[0], best_highs[0]
            if len(meet) != n*n or len(join) != n*n:
                continue
            # Natural ordering forces the lattice bottom and top to be 0,n-1.
            need(all(le[0][x] and le[x][n-1] for x in range(n)), "bounded lattice")
            labels = [prod(primes[t] for t in range(n-1) if not le[x][t]) for x in range(n)]
            constants = [1+x % 3 for x in range(1, n)]
            image = [tuple(c if le[x][y] else 0 for x, c in zip(range(1, n), constants)) for y in range(n)]
            need(len(set(image)) == n, "encoding injection")
            for x in range(n):
                for y in range(n):
                    need((labels[y] % labels[x] == 0) == le[x][y], "order divisibility")
                    need(lcm(labels[x], labels[y]) == labels[join[x,y]], "join LCM")
                    need(tuple(min(a,b) for a,b in zip(image[x],image[y])) == image[meet[x,y]], "ambient meet")
                    active = [z for z in range(1,n) if image[x][z-1] or image[y][z-1]]
                    length = lcm(*(labels[z] for z in active)) if active else 1
                    closed = tuple(c if length % labels[z] == 0 else 0 for z,c in zip(range(1,n),constants))
                    need(closed == image[join[x,y]], "ambient join")
                    pair_checks += 1
            fixed = []
            for bits in range(1 << (n-1)):
                row = tuple(c if bits >> (z-1) & 1 else 0 for z,c in zip(range(1,n),constants))
                length = lcm(*(labels[z] for z in range(1,n) if row[z-1]))
                closed = tuple(c if length % labels[z] == 0 else 0 for z,c in zip(range(1,n),constants))
                if closed == row:
                    fixed.append(row)
                interval_checks += 1
            need(set(fixed) == set(image), "interval contains extra profile")
            digest.update(json.dumps([n,mask,labels],separators=(',',':')).encode())
            count += 1
        counts[str(n)] = count
    return {"naturally_ordered_lattices_by_size": counts,
            "encoding_pair_checks": pair_checks, "interval_binary_checks": interval_checks,
            "encoding_sha256": digest.hexdigest(),
            "scope": "Symbolic lattice encodings, not arithmetic witnesses for every enumerated lattice."}


def main():
    universal = lattice_checks()
    penta = Source(18481, (67,331,50821,4621))
    need(penta.p % 1155 == 1, "source label avoidance")
    need(penta.d == (33,165,105,1155), "N5 exact orders")
    need(penta.c == (1,1,1,1), "N5 base valuations")
    lengths = [3,33,165,105,1155]
    penta_profiles = [direct(penta.p, n, penta.targets) for n in lengths]
    expected = [(0,0,0,0),(1,0,0,0),(1,1,0,0),(0,0,1,0),(1,1,1,1)]
    need(penta_profiles == expected, "N5 profiles")
    need(set(penta.domain((0,0,0,0),(1,1,1,1))) == set(expected), "N5 exact interval")

    sources = [Source(3,(11,23)), Source(5,(11,23),True)]
    lower, upper, totals = [[0,0],[0,0]], [[2,1],[1,1]], [2,2]
    domains = [s.domain(a,b) for s,a,b in zip(sources,lower,upper)]
    need(all(domains), "each row feasible")
    for j in range(2):
        need(any(sum(row[j] for row in choice) == totals[j]
                 for choice in itertools.product(*domains)), "each column separately feasible")
    choices = list(itertools.product(*domains))
    need(not any(all(sum(row[j] for row in choice) == totals[j] for j in range(2))
                 for choice in choices), "coupled problem unexpectedly feasible")
    coupled = propagate(sources,lower,upper,totals)
    need(not coupled[0], "propagation missed coupled contradiction")
    need(direct(3,5,(11,23)) == (2,0) and direct(5,2,(11,23)) == (0,0), "column 11 witnesses")
    need(direct(3,11,(11,23)) == (0,1) and direct(5,22,(11,23)) == (0,1), "column 23 witnesses")

    parity_sources = [Source(p,(13,)) for p in (191,653,991,1543)]
    need(all(s.d == (3,) and s.c == (2,) for s in parity_sources), "parity arithmetic")
    parity_domains = [s.domain((0,),(2,)) for s in parity_sources]
    need(all(set(d) == {(0,),(2,)} for d in parity_domains), "parity domains")
    parity = propagate(parity_sources,[[0]]*4,[[2]]*4,[5])
    need(parity[0] and parity[1] == [[0]]*4 and parity[2] == [[2]]*4, "parity fixed point")
    need(not any(sum(row[0] for row in choice) == 5
                 for choice in itertools.product(*parity_domains)), "odd sum from even rows")
    need(all(direct(s.p,3,(13,)) == (2,) and direct(s.p,5,(13,)) == (0,)
             for s in parity_sources), "parity direct witnesses")

    # Soundness of all propagation rules against exact finite row allocations.
    rng = random.Random(188)
    pool = [Source(3,(11,23)),Source(5,(11,23),True),Source(13,(11,23))]
    random_counts = {"cases": 0, "feasible": 0, "rejected": 0,
                     "infeasible_nonrejecting": 0, "preserved_witnesses": 0}
    for case in range(400):
        ss = pool[:2+(case % 2)]
        uu = [[rng.randrange(3) for _ in range(2)] for _ in ss]
        ll = [[rng.randrange(b+1) if case % 3 else 0 for b in row] for row in uu]
        dd = [s.domain(a,b) for s,a,b in zip(ss,ll,uu)]
        target = [rng.randrange(1+sum(row[j] for row in uu)) for j in range(2)]
        if case % 4 == 0 and all(dd):
            chosen = [rng.choice(d) for d in dd]
            target = [sum(r[j] for r in chosen) for j in range(2)]
        solutions = [choice for choice in itertools.product(*dd)
                     if all(sum(r[j] for r in choice) == target[j] for j in range(2))]
        passed, final_l, final_u, _ = propagate(ss,ll,uu,target)
        need(passed or not solutions, "false rejection")
        for sol in solutions:
            need(all(final_l[i][j] <= sol[i][j] <= final_u[i][j]
                     for i in range(len(ss)) for j in range(2)), "removed feasible allocation")
        random_counts["cases"] += 1
        random_counts["feasible"] += bool(solutions)
        random_counts["rejected"] += not passed
        random_counts["infeasible_nonrejecting"] += passed and not solutions
        random_counts["preserved_witnesses"] += len(solutions)

    return {"status":"PASS", "lattice_universality_checks":universal,
            "pentagon":{"source":penta.p,"targets":penta.targets,"orders":penta.d,
                        "constants":penta.c,"physical_lengths":lengths,"profiles":penta_profiles},
            "coupled_certificate":{"domains":domains,"column_totals":totals,"trace":coupled[3]},
            "parity_certificate":{"sources":[s.p for s in parity_sources],"target":13,
                                  "domains":parity_domains,"total":5,"modulus":2,
                                  "propagation_fixed_point":True,"trace":parity[3]},
            "propagation_soundness_checks":random_counts,
            "scope":"Local geometry and conditional partial-row allocation, not a complete OPN support or global nonexistence theorem."}


if __name__ == '__main__':
    print(json.dumps(main(),sort_keys=True,indent=2))
