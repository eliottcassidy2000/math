"""Exact finite controls for geometric_carrier.md; not a scheme-proof engine."""
import itertools
import json
from collections import Counter
from pathlib import Path


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def unit(n, i):
    return tuple(int(j == i) for j in range(n))


def norm(v, prime):
    first = next(x for x in v if x)
    inv = pow(first, -1, prime)
    return tuple(x * inv % prime for x in v)


def projective(n, prime):
    # Every line has a unique representative with first nonzero entry one.
    return tuple(v for v in itertools.product(range(prime), repeat=n)
                 if any(v) and next(x for x in v if x) == 1)


def C(x, prime):
    return tuple(norm(v[::-1], prime) for v in x)


def Q(x):
    return x[0], x[2], x[1]


def main():
    ranks = {}
    for k in range(2, 13):
        divisors = list(itertools.product(range(k + 1), range(2), range(2)))
        proper = [v for v in divisors if v not in [(0, 0, 0), (k, 1, 1)]]
        points = {(unit(k + 1, a), unit(2, b), unit(2, c)): (a, b, c)
                  for a, b, c in proper}
        require(len(points) == 4 * k + 2, "coordinate injection")
        fixed = [0, 0, 0]
        for x, (a, b, c) in points.items():
            require(points[C(x, 2)] == (k - a, 1 - b, 1 - c), "C")
            require(points[Q(x)] == (a, c, b), "Q")
            for i, image in enumerate((C(x, 2), Q(x), C(Q(x), 2))):
                fixed[i] += image == x
        require(fixed == [0, 2 * k, 2 if k % 2 == 0 else 0], "fixed counts")
        if k == 2:
            rank = Counter(sum(v) for v in proper)
            sign = lambda z: (z > 0) - (z < 0)
            f = Counter(sign(2*a-k*(b+c)) for a,b,c in proper)
            g = Counter(sign(b-c) for a,b,c in proper)
            require([rank[i] for i in (1, 2, 3)] == [3, 4, 3], "rank")
            require([g[i] for i in (-1, 0, 1)] == [3, 4, 3], "g rows")
            require([f[i] for i in (-1, 0, 1)] == [4, 2, 4], "f columns")
            ranks = {"rank": dict(rank), "sign_g": dict(g), "sign_f": dict(f)}

    fields = []
    for k, prime in [(k, p) for k in range(2, 6) for p in (2, 3)] + [(2, 5)]:
        factors = [projective(n, prime) for n in (k + 1, 2, 2)]
        all_points = set(itertools.product(*factors))
        extremes = {(unit(k + 1, 0), unit(2, 0), unit(2, 0)),
                    (unit(k + 1, k), unit(2, 1), unit(2, 1))}
        opened = all_points - extremes
        polynomial = sum(prime ** i for i in range(k + 1)) * (1 + prime) ** 2 - 2
        require(len(opened) == polynomial, "finite-field open count")
        for x in opened:
            cx = C(x, prime)
            require(cx in opened and Q(x) in opened, "invariant open")
            require(C(cx, prime) == x and Q(Q(x)) == x, "involutions")
            require(C(Q(x), prime) == Q(cx), "commutation")
        # Count the standard product cells independently by last nonzero index.
        dimensions = Counter(sum(max(i for i,z in enumerate(v) if z) for v in x)
                             for x in all_points)
        coefficients = Counter(a+b+c for a,b,c in
                               itertools.product(range(k+1), range(2), range(2)))
        require(all(dimensions[d] == coefficients[d] * prime ** d for d in coefficients),
                "standard cell dimensions")
        fields.append({"k": k, "prime": prime, "open_points": len(opened)})

    supports = [v for v in itertools.product((0, 1), repeat=3) if any(v)]
    require(len(supports) == 7, "seven literal Hom points")
    reversal_fixed = [v for v in supports if v == v[::-1]]
    require(len(reversal_fixed) == 3, "coordinate reversal fixed supports")
    torus_ranks = Counter()
    for a,b,c in itertools.product(supports, projective(2,2), projective(2,2)):
        if (a,b,c) in {(unit(3,0),unit(2,0),unit(2,0)),
                       (unit(3,2),unit(2,1),unit(2,1))}:
            continue
        torus_ranks[sum(map(sum, (a,b,c)))-3] += 1
    require([torus_ranks[i] for i in range(5)] == [10,24,19,7,1], "torus strata")
    require(sum(torus_ranks.values()) == 61, "literal open Hom count")
    require(sum(torus_ranks[i]*2**i for i in torus_ranks) == 206, "q=3 torus count")
    require(3*2+4*2**2+3*2**3 == 46 and 61-46 == 2**4-1, "cell deletion hostile")

    result = {"status": "FINITE-EXACT", "checks": CHECKS,
              "universe": "rank-zero action k=2..12; normalized prime-field projective vectors as listed",
              "finite_field_controls": fields, "ten_point_partitions": ranks,
              "U2_torus_rank_counts": dict(torus_ranks),
              "seven_coordinate_reversal_fixed": reversal_fixed,
              "limits": "Coordinate/action/count checks; geometric existence and all-k assertions use written proofs."}
    output = json.dumps(result, indent=2, sort_keys=True) + "\n"
    Path(__file__).with_name("geometry_receipt.json").write_text(output, encoding="utf-8")
    print(output, end="")


if __name__ == "__main__":
    main()
