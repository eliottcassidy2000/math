"""Exact controls for GROWTH.md; no external packages, no assert statements."""
from functools import lru_cache
from itertools import permutations
from pathlib import Path
import hashlib
import json


def require(test, message):
    if not test:
        raise RuntimeError(message)


def rows_from_bits(n, bits):
    rows = [0] * n
    pos = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> pos & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
            pos += 1
    return rows


def all_counts(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    counts = [0] * (1 << n)
    counts[0] = 1
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        counts[mask] = sum(dp[mask])
        for v, value in enumerate(dp[mask]):
            if value:
                available = rows[v] & ~mask
                while available:
                    bit = available & -available
                    available -= bit
                    dp[mask | bit][bit.bit_length() - 1] += value
    return counts


def cyclic_mask(rows, subset):
    vertices = [v for v in range(len(rows)) if subset >> v & 1]
    result = 0
    for i in vertices:
        for j in vertices:
            if not rows[i] >> j & 1:
                continue
            for k in vertices:
                if rows[j] >> k & 1 and rows[k] >> i & 1:
                    result |= (1 << i) | (1 << j) | (1 << k)
    return result


@lru_cache(None)
def core_budget(h):
    best = 0
    pow5, b = 1, 0
    while pow5 <= h:
        pow3, a = pow5, 0
        while pow3 <= h:
            best = max(best, 3 * a + 4 * b)
            pow3 *= 3
            a += 1
        pow5 *= 5
        b += 1
    return best


def threshold_word(h):
    backwards = []
    while h > 1:
        if h % 3 == 0:
            backwards.append("I")
            h //= 3
        else:
            backwards.append("U")
            h -= 2
    return backwards[::-1]


def threshold_rows(word):
    rows = []
    for op in word:
        n = len(rows)
        if op == "U":
            for i in range(n):
                rows[i] |= 1 << n
            rows.append((1 << n) - 1)
        else:
            rows.append(0)
    return rows


def independent_weight(rows):
    @lru_cache(None)
    def calc(mask):
        if not mask:
            return 1
        v = mask.bit_length() - 1
        rest = mask ^ (1 << v)
        return calc(rest) + 2 * calc(rest & ~rows[v])
    return calc((1 << len(rows)) - 1)


def main():
    census = []
    comparisons = 0
    independent = 0
    for n in range(7):
        total = 1 << (n * (n - 1) // 2)
        max_positive = 0
        for bits in range(total):
            rows = rows_from_bits(n, bits)
            hs = all_counts(rows)
            prefixes = [(1 << j) - 1 for j in range(n + 1)]
            cores = [cyclic_mask(rows, mask) for mask in prefixes]
            positives = [int(hs[prefixes[j]] > hs[prefixes[j - 1]])
                         for j in range(1, n + 1)]
            max_positive = max(max_positive, sum(positives))
            for i in range(n + 1):
                for j in range(i, n + 1):
                    positive = sum(positives[i:j])
                    require(cores[i] & ~cores[j] == 0, "cyclic vertex lost")
                    require(cores[i].bit_count() + positive <= cores[j].bit_count(),
                            "G1 persistent cyclic-vertex inequality")
                    require(cores[j].bit_count() <= core_budget(hs[prefixes[j]]),
                            "inherited exact core budget")
                    comparisons += 1
            if n <= 5:
                brute = sum(all(rows[p[i]] >> p[i + 1] & 1 for i in range(n - 1))
                            for p in permutations(range(n)))
                require(brute == hs[-1], "independent full path count")
                independent += 1
        census.append({"order": n, "labelled_tournaments": total,
                       "largest_positive_prefix_depth": max_positive})
    targets = list(range(1, 10002, 2)) + [10**30 + 1, 3**80, 2**127 - 1]
    examples = []
    for h in targets:
        word = threshold_word(h)
        require(independent_weight(threshold_rows(word)) == h, "threshold compiler")
        power, floorlog = 1, 0
        while power * 3 <= h:
            floorlog += 1
            power *= 3
        require(len(word) <= 3 * floorlog + 2, "threshold word bound")
        if h in (5, 7, 17, 21) or h > 10001:
            examples.append({"h": h, "vertices": len(word), "word": "".join(word)})
    report = {"status": "VERIFIED finite controls; universal proof in GROWTH.md",
              "universe": "all labelled tournaments of order 0 through 6; no filters",
              "tournaments": sum(c["labelled_tournaments"] for c in census),
              "all_start_prefix_comparisons": comparisons,
              "permutation_cross_checks_through_order5": independent,
              "census": census, "threshold_targets_checked": len(targets),
              "threshold_examples": examples,
              "source_sha256_lf": hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n', b'\n')).hexdigest()}
    destination = Path(__file__).with_name("growth_receipt.json")
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
