"""Independent literal-query audit. Imports no response implementation."""
import itertools
import json
from functools import lru_cache


def check(ok, detail):
    if not ok:
        raise RuntimeError(detail)


def graph(n, code):
    rows = [0] * n
    for bit, (i, j) in enumerate(itertools.combinations(range(n), 2)):
        u, v = (i, j) if code & (1 << bit) else (j, i)
        rows[u] |= 1 << v
    return tuple(rows)


def hamilton(rows):
    full = (1 << len(rows)) - 1

    @lru_cache(None)
    def walk(last, unused):
        if not unused:
            return 1
        candidates = rows[last] & unused
        total = 0
        while candidates:
            bit = candidates & -candidates
            candidates -= bit
            total += walk(bit.bit_length() - 1, unused ^ bit)
        return total

    return sum(walk(v, full ^ (1 << v)) for v in range(len(rows)))


def response(rows, a, c):
    n = len(rows)
    out = list(rows) + [1 << (n + 1), 0]
    for v in range(n):
        if a & (1 << v):
            out[v] |= 1 << n
        else:
            out[n] |= 1 << v
        if c & (1 << v):
            out[v] |= 1 << (n + 1)
        else:
            out[n + 1] |= 1 << v
    return hamilton(tuple(out))


def coefficient(rows, A, B, cache):
    n = len(rows)
    positions = [i for i in range(n) if A >> i & 1]
    positions += [n+i for i in range(n) if B >> i & 1]
    total = 0
    for bits in range(1 << len(positions)):
        assignment = sum(1 << p for k, p in enumerate(positions) if bits >> k & 1)
        a, c = assignment & ((1 << n) - 1), assignment >> n
        if (a, c) not in cache:
            cache[a, c] = response(rows, a, c)
        total += (-1) ** (len(positions)-bits.bit_count()) * cache[a, c]
    return total


def cover(rows, names, end, start):
    # Literal disjoint two-path covers; no endpoint dynamic programming.
    total = 0
    for word in itertools.permutations(names):
        for cut in range(1, len(names)):
            if word[cut-1] != end or word[cut] != start:
                continue
            if all(k == cut or rows[word[k-1]] >> word[k] & 1
                   for k in range(1, len(names))):
                total += 1
    return total


def three_cover_sum(rows, A, B):
    total = 0
    n = len(rows)
    for word in itertools.permutations(range(n)):
        for first, second in itertools.combinations(range(1, n), 2):
            if second == first + 1:
                continue
            boundaries = [{word[first-1], word[first]},
                          {word[second-1], word[second]}]
            if boundaries != [set(A), set(B)] and boundaries != [set(B), set(A)]:
                continue
            if all(k in (first, second) or rows[word[k-1]] >> word[k] & 1
                   for k in range(1, n)):
                total += 1
    return total


def split(transitive_size, root_mask):
    m = transitive_size
    # Root=0, exterior=1..m, second twin=m+1.
    rows = [0] * (m + 2)
    for i in range(1, m + 1):
        for j in range(i + 1, m + 1):
            rows[i] |= 1 << j
        if root_mask >> (i-1) & 1:
            rows[0] |= 1 << i
            rows[m+1] |= 1 << i
        else:
            rows[i] |= (1 << 0) | (1 << (m+1))
    rows[0] |= 1 << (m+1)
    flipped = list(rows)
    flipped[0] ^= 1 << (m+1)
    flipped[m+1] |= 1
    return tuple(rows), tuple(flipped)


def is_strong(rows):
    full = (1 << len(rows))-1
    for v in range(len(rows)):
        seen, frontier = 1 << v, [v]
        while frontier:
            u = frontier.pop()
            for w in range(len(rows)):
                if rows[u] >> w & 1 and not seen >> w & 1:
                    seen |= 1 << w
                    frontier.append(w)
        if seen != full:
            return False
    return True


def main():
    counts = {"tournaments": 0, "literal_queries": 0, "quartics": 0,
              "shared_name": 0, "disjoint_pairs": 0, "double_shared_zero": 0,
              "transitive_twin_cases": 0}
    fixtures = [graph(n, code) for n in range(1, 5)
                for code in range(1 << (n*(n-1)//2))]
    fixtures += [graph(5, (179*k+23) % 1024) for k in range(16)]
    fixtures += [graph(6, (8191*k+7) % 32768) for k in range(4)]
    hostile, hostile_flip = split(3, 1)
    fixtures += [hostile, hostile_flip]
    for rows in fixtures:
        n = len(rows)
        cache = {}
        for A in itertools.combinations(range(n), 2):
            for B in itertools.combinations(range(n), 2):
                coefficient_got = coefficient(rows, sum(1 << i for i in A),
                                              sum(1 << i for i in B), cache)
                shared = set(A) & set(B)
                if len(shared) == 2:
                    want = 0
                    counts["double_shared_zero"] += 1
                elif len(shared) == 1:
                    i = next(iter(shared))
                    j = next(iter(set(A)-shared))
                    k = next(iter(set(B)-shared))
                    names = [v for v in range(n) if v != i]
                    want = cover(rows, names, j, k) + cover(rows, names, k, j)
                    counts["shared_name"] += 1
                else:
                    want = three_cover_sum(rows, A, B)
                    counts["disjoint_pairs"] += 1
                check(coefficient_got == want, ("quartic", rows, A, B, coefficient_got, want))
                if n <= 4:
                    universal = 0 if len(shared) == 2 else 2 if shared else 4
                    check(coefficient_got == universal, "small-order universal quartic")
                counts["quartics"] += 1
        counts["literal_queries"] += len(cache)
        counts["tournaments"] += 1

    twin_records = []
    for m, bits in ([(m, bits) for m in range(2, 6) for bits in range(1 << m)]
                    + [(m, 1) for m in range(6, 9)]):
        rows, flipped = split(m, bits)
        A, B = (1 << (m+1)) | 2, 1 | (1 << m)
        ca, cb = {}, {}
        got = coefficient(rows, A, B, ca) - coefficient(flipped, A, B, cb)
        exterior = list(range(1, m+1))
        derivative = cover(rows, exterior, m, 1) - cover(rows, exterior, 1, m)
        want = 0 if m == 2 else 2 ** (m-2)
        check(got == derivative == want, ("twin", m, bits, got, derivative, want))
        if bits == 1:
            check(is_strong(rows) and is_strong(flipped), "strong family")
        counts["transitive_twin_cases"] += 1
        counts["literal_queries"] += len(ca)+len(cb)
        twin_records.append({"m": m, "root_mask": bits, "coefficient": got})
    check(hostile == (18,12,25,17,2), "hostile adjacency")
    check(hostile_flip == (2,12,25,17,3), "hostile flip adjacency")
    return {"status": "PASS", "counts": counts, "twin_records": twin_records,
            "claim_review": "Repeated-label formulas, universal order-four ceiling, sharp strong order-five witness, and all-order transitive twin proof pass.",
            "limitations": "No full-response collision, named reconstruction, or hidden prime-core exclusion claimed."}


if __name__ == "__main__":
    print(json.dumps(main(), sort_keys=True, indent=2))
