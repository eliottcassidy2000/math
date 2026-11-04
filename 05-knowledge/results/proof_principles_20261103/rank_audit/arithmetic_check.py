"""Independent arithmetic audit: literal geometric sums and exact periods.

No arithmetic producer is imported. Values and budgets are independently
computed from integer geometric sums, including every finite-box vector.
"""
from argparse import ArgumentParser
from hashlib import sha256
from itertools import product
from math import lcm, prod
from pathlib import Path
import json

SOURCE = Path(__file__)
SOURCE_HASH = sha256(SOURCE.read_bytes()).hexdigest()
CHECKS = 0


def require(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def val(x, p):
    answer = 0
    while x % p == 0:
        answer += 1
        x //= p
    return answer


def geometric(p, length):
    return (p ** length - 1) // (p - 1)


def order(p, r):
    return next(d for d in range(1, r) if pow(p, d, r) == 1)


def row_realization(p, targets, desired, role):
    # Assemble divisibility requirements, then check the candidate by a
    # direct geometric sum rather than the producer's LTE evaluation code.
    n = 2 if role else 1
    for r, t in zip(targets, desired):
        if t:
            d = order(p, r)
            base = 0 if d == 1 else val(p ** d - 1, r)
            if t < base:
                return None
            n = lcm(n, (r ** t if d == 1 else d * r ** (t - base)))
    if (role and n % 4 != 2) or (not role and n % 2 != 1):
        return None
    if n == 1:
        # Only targets 5 and 11 are used below, all orders <=10.
        n = 13
    actual = tuple(val(geometric(p, n), r) for r in targets)
    return n if actual == desired else None


def exact_periods():
    targets = (5, 11)
    result = []
    for p in (3, 7, 13):
        period = lcm(4, *(order(p, r) * r ** 3 for r in targets))
        seen = [set(), set()]
        powers = [1, 1]
        sums = [0, 0]
        for n in range(1, period + 1):
            for j, r in enumerate(targets):
                modulus = r ** 3
                sums[j] = (sums[j] + powers[j]) % modulus
                powers[j] = powers[j] * p % modulus
            if all(sums):
                exact = tuple(val(x, r) for x, r in zip(sums, targets))
                if n % 2:
                    seen[0].add(exact)
                elif n % 4 == 2:
                    seen[1].add(exact)
        require(sums == [0, 0] and powers == [1, 1], "full modular period")
        for role in (0, 1):
            predicted = {t for t in product(range(3), repeat=2)
                         if row_realization(p, targets, t, role) is not None}
            require(seen[role] == predicted, ("exact row iff", p, role))
        result.append({"source": p, "targets": targets, "period": period,
                       "feasible": [sorted(x) for x in seen]})
    return result


def complete_support(q):
    support = (3, 5, 7, 11, 13)
    h = len(support) - 1
    baselines = []
    permitted = set()
    for r in support:
        c = 0
        for p in support:
            if p == r:
                continue
            d = order(p, r)
            if (d % 2 == 1) or (p == q and d % 4 == 2):
                permitted.add((p, r))
                if d > 1:
                    c += val(p ** d - 1, r)
        baselines.append(c)
    C = max(baselines)
    B = max(2, 2 * h, C + 1)
    while 3 ** (B - C) <= (B + 1) ** h:
        B += 1
    require((C, B) == ((3, 13) if q == 5 else (4, 14)), "box endpoint")
    require(3 ** (B - C) > (B + 1) ** h, "strict endpoint certificate")
    # Integer certificate for the Bernoulli monotonicity tail at its worst
    # boundary; the proof establishes persistence for all larger integers.
    require(2 * h < B + 2, "Bernoulli ratio exceeds 3/2")
    for p in support:
        reached = {p}
        while True:
            extended = reached | {v for u, v in permitted if u in reached}
            if extended == reached:
                break
            reached = extended
        require(reached == set(support), "graph admits whole support")
    domains = [range(1, B, 4) if p == q else range(2, B, 2) for p in support]
    direct_rows = {(p, a): tuple(val(geometric(p, a + 1), r) for r in support)
                   for p, exponents in zip(support, domains) for a in exponents}
    count = balanced = perfect = 0
    first_failure = {str(p): 0 for p in support}
    for exponents in product(*domains):
        count += 1
        budget = tuple(sum(direct_rows[p, a][j] for p, a in zip(support, exponents))
                       for j in range(len(support)))
        mismatch = next((j for j in range(len(support)) if budget[j] != exponents[j]), None)
        if mismatch is None:
            balanced += 1
        else:
            first_failure[str(support[mismatch])] += 1
        # Direct final equality is independently checked for all candidates,
        # not merely those passing the budget filter.
        N = prod(p ** a for p, a in zip(support, exponents))
        total_sigma = prod(geometric(p, a + 1) for p, a in zip(support, exponents))
        perfect += total_sigma == 2 * N
        require(val(total_sigma, 2) == 1, "Euler parity supplies exactly one factor 2")
        if q == 5:
            require(budget[4] <= 1 < exponents[4], "q=5 hand contradiction")
        else:
            require(budget[1] <= 1 < exponents[1], "q=13 hand contradiction")
    require(count == (3888 if q == 5 else 5184), "box cardinality")
    require(balanced == perfect == 0, "complete support exclusion")
    return {"euler": q, "support": support, "C": C, "B_exclusive": B,
            "baselines": baselines, "vectors": count, "balanced": balanced,
            "perfect": perfect, "first_failure": first_failure}


def run():
    periods = exact_periods()
    boxes = [complete_support(q) for q in (5, 13)]
    sharp = geometric(3, 55)
    require((val(sharp, 11), val(sharp, 23)) == (3, 1), "sharp joint transfer")
    require(val(geometric(3, 5), 11) == 2 and val(geometric(3, 5), 23) == 0,
            "first antecedent alone")
    require(val(geometric(3, 11), 11) == 0 and val(geometric(3, 11), 23) == 1,
            "second antecedent alone")
    require(geometric(3, 3) == 13 and geometric(3, 5) == 121,
            "increasing exponent can remove an edge")
    require(val(geometric(23, 3), 11) == 0 and val(geometric(23, 11), 11) == 1,
            "order-one denominator exception")
    require(sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_HASH, "source stable")
    return {"status": "FINITE-EXACT independent arithmetic audit passed",
            "source_sha256": SOURCE_HASH, "checks": CHECKS,
            "independent_periods": periods, "direct_sigma_boxes": boxes,
            "sharp_source3_length55_valuations": [3, 1],
            "scope": "Complete modular periods for two targets with exact valuations0..2, plus every9072 vector in the two proved fixed-support boxes. Universal claims audited separately."}


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    text = json.dumps(run(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8", newline="\n")
    else:
        print(text, end="")
