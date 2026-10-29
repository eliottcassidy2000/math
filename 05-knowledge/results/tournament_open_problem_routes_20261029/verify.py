"""Independent finite guards for the proofs; no finite test proves an asymptotic."""
from collections import Counter
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent

def require(test, message):
    if not test:
        raise RuntimeError(message)

def rows(n, code):
    a = [0] * n
    for bit, (i, j) in enumerate(combinations(range(n), 2)):
        if code >> bit & 1:
            a[i] |= 1 << j
        else:
            a[j] |= 1 << i
    return a

def transitive(a, mask):
    # Independent source-peeling test, not a triangle-pattern test.
    while mask:
        for v in range(len(a)):
            if mask >> v & 1 and (a[v] & mask).bit_count() == mask.bit_count() - 1:
                mask ^= 1 << v
                break
        else:
            return False
    return True

def triangles(a):
    return [sum(1 << i for i in s) for s in combinations(range(len(a)), 3)
            if all((a[i] & sum(1 << j for j in s)).bit_count() == 1 for i in s)]

def vc_at_most_one(a):
    return all(len({x & ((1 << i) | (1 << j)) for x in a}) < 4
               for i, j in combinations(range(len(a)), 2))

def classified(a):
    full = (1 << len(a)) - 1
    if transitive(a, full):
        return True
    for c in triangles(a):
        rest = full ^ c
        if transitive(a, rest) and all(a[i] & rest == rest for i in range(len(a)) if c >> i & 1):
            return True
    return False

def ham(a):
    n = len(a)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not dp[mask][v]:
                continue
            for w in range(n):
                if a[v] >> w & 1 and not mask >> w & 1:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[-1])

def largest_transitive(a):
    return max(s.bit_count() for s in range(1 << len(a)) if transitive(a, s))

def substitute(q, blocks):
    sizes = [len(b) for b in blocks]
    offsets = [sum(sizes[:i]) for i in range(len(sizes))]
    a = [0] * sum(sizes)
    for i, block in enumerate(blocks):
        for v, r in enumerate(block):
            a[offsets[i] + v] = r << offsets[i]
            for j in range(len(q)):
                if q[i] >> j & 1:
                    a[offsets[i] + v] |= ((1 << sizes[j]) - 1) << offsets[j]
    return a

def main():
    census = []
    total = 0
    for n in range(1, 7):
        low = 0
        for code in range(1 << (n * (n - 1) // 2)):
            a = rows(n, code)
            small = vc_at_most_one(a)
            require(small == classified(a), ("VC classification", n, code))
            low += small
            cs = len(triangles(a))
            h = ham(a)
            tr = largest_transitive(a)
            require(h >= 1 + 2 * cs and 2 * tr >= 2 * n - h + 1, ("H bound", n, code))
            total += 1
        expected = factorial(n) if n < 3 else 4 * factorial(n) // 3
        require(low == expected, ("labelled count", n, low, expected))
        census.append({"n": n, "tournaments": 1 << (n * (n - 1) // 2), "VC_le_1": low})

    # Weighted quotient law checked against literal induced-subset search.
    blocks = [[0], rows(2, 1), [2, 4, 1]]
    queries = 0
    for code in range(8):
        q = rows(3, code)
        for choices in product(range(3), repeat=3):
            bs = [blocks[i] for i in choices]
            actual = largest_transitive(substitute(q, bs))
            ws = [largest_transitive(b) for b in bs]
            predicted = max(sum(ws[i] for i in range(3) if s >> i & 1)
                            for s in range(8) if transitive(q, s))
            require(actual == predicted, ("substitution", code, choices))
            queries += 1

    # The definition in HYP-9375, independently built from Euler's criterion.
    primes = [n for n in range(3, 1000, 4) if all(n % d for d in range(2, int(n**0.5) + 1))]
    prime_rows = []
    for p in primes:
        prime_rows.append(sum(1 << j for j, q in enumerate(primes) if p != q and pow(p, (q-1)//2, q) == 1))
    require(all(((prime_rows[i] >> j) ^ (prime_rows[j] >> i)) & 1
                for i in range(len(primes)) for j in range(i)), "reciprocity")
    # Literal exhaustive cut maximization at the first 12 primes: for a fixed A,
    # the maximizing B is the set of positive column sums (or negative ones).
    n = 12
    s = [[0 if i == j else (1 if prime_rows[i] >> j & 1 else -1) for j in range(n)] for i in range(n)]
    cut = 0
    for mask in range(1 << n):
        sums = [sum(s[i][j] for i in range(n) if mask >> i & 1) for j in range(n)]
        cut = max(cut, sum(max(x, 0) for x in sums), -sum(min(x, 0) for x in sums))
    # Unlabelled four-class counts: degree multisets distinguish all four classes.
    four = Counter()
    first = primes[:24]
    for vs in combinations(range(len(first)), 4):
        a = [sum(1 << j for j, v in enumerate(vs) if prime_rows[u] >> v & 1) for u in vs]
        signature = tuple(sorted(x.bit_count() for x in a))
        four[str(signature)] += 1
    require(sum(four.values()) == 10626 and len(four) == 4, "four-class census")
    # Non-tournament hostile for allowing primes 1 mod 4.
    require(pow(5, 14, 29) == 1 and pow(29, 2, 5) == 1, "mixed-prime hostile")
    result = {"status": "FINITE-EXACT checks; analytic proof uses cited large sieve, not these data",
              "all_labelled_tournaments_tested": total, "vc_census": census,
              "weighted_substitution_queries": queries,
              "prime_prefix": first, "four_class_counts": dict(sorted(four.items())),
              "first12_max_absolute_signed_cut": cut,
              "large_sieve_constant_estimated": False,
              "independent_referee": False}
    (HERE / "verification.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
