"""Full-word precedence, simple-cycle supports and flip parity, by enumeration."""

from itertools import combinations
from pathlib import Path

BANK = Path(__file__).resolve().parents[3] / "session212" / "cycle-levels"


def simple_paths(rows):
    def extend(word, used):
        yield word
        choices = rows[word[-1]] & ~used
        while choices:
            bit = choices & -choices
            choices ^= bit
            yield from extend(word + (bit.bit_length() - 1,), used | bit)

    for v in range(len(rows)):
        yield from extend((v,), 1 << v)


def inspect(rows):
    n = len(rows)
    pairs = list(combinations(range(n), 2))
    pair_bits = {pair: 1 << i for i, pair in enumerate(pairs)}
    direct = sum(pair_bits[a, b] for a, b in pairs if rows[a] >> b & 1)
    full = (1 << n) - 1
    order = [0] * (full + 1)
    cycles = [0] * (full + 1)
    detours = 0
    for path in simple_paths(rows):
        size = len(path)
        if size < 2:
            continue
        support = sum(1 << v for v in path)
        order[support] ^= sum(pair_bits[a, b] for a, b in combinations(path, 2) if a < b)
        u, v = path[0], path[-1]
        if size >= 3 and u == min(path) and rows[v] >> u & 1:
            cycles[support] ^= 1
        if size >= 3 and (size - 1) % 4 in (1, 2):
            detours ^= pair_bits[tuple(sorted((u, v)))]

    for a, b in pairs:
        bit = pair_bits[a, b]
        containing = [s for s in range(full + 1) if s >> a & s >> b & 1]
        assert bool((direct ^ order[full]) & bit) == (sum(cycles[s] for s in containing) % 2)
        if n >= 3:
            assert (sum(bool(order[s] & bit) for s in containing) % 2) == cycles[full]

    # A full word contributes to the 'smaller before larger' mask. Its reverse
    # pair indicators yield its inversions, since the total word count is odd.
    all_pairs = (1 << len(pairs)) - 1
    odd_fiber = (order[full] ^ all_pairs).bit_count() % 2
    return odd_fiber, detours ^ all_pairs


def main():
    for n, expected in enumerate((0, 1, 1, 2, 4, 12, 56, 456)):
        if n < 2:
            continue
        pairs = list(combinations(range(n), 2))
        seen = 0
        flip_tests = 0
        cache = {}
        for line in (BANK / f"classes{n}.rows").read_text().splitlines():
            rows = tuple(map(int, line.split()))
            odd, derivatives = inspect(rows)
            seen += 1
            if n <= 6:
                for i, (a, b) in enumerate(pairs):
                    neighbor = list(rows)
                    neighbor[a] ^= 1 << b
                    neighbor[b] ^= 1 << a
                    neighbor = tuple(neighbor)
                    if neighbor not in cache:
                        cache[neighbor] = inspect(neighbor)[0]
                    assert odd ^ cache[neighbor] == derivatives >> i & 1
                    flip_tests += 1
        assert seen == expected
        print(f"order {n}: {seen} classes, {seen * len(pairs)} marked pairs, {flip_tests} direct flips")


if __name__ == "__main__":
    main()
