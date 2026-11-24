"""Explicit integer witnesses with transitive odd-precedence tournament."""
from math import gcd

# Label each seed by its precedence order. These labels need not form a path.
SEEDS = {
    5: (6, 8, 10, 1),
    17: (30, 60, 56, 48, 32, 1),
    37: (30, 44, 48, 52, 34, 1),
    77: (222, 252, 248, 240, 96, 129, 160, 16),
    129: (254, 508, 504, 496, 480, 448, 384, 256, 1),
    751: (150, 504, 506, 257, 40, 73, 153, 312, 113),
}



# Successive orders supply the three residues still missing at nine.
LATE_RESIDUES = {
    1563: (782, 508, 504, 432, 353, 193, 393, 17, 160, 510),
    1243: (838, 508, 504, 369, 161, 449, 400, 265, 16, 510),
    5363: (1806, 1532, 504, 1456, 353, 193, 393, 17, 160, 1534, 500),
}


def count_and_order(rows):
    n = len(rows)
    predecessors = [0] * n
    full = (1 << n) - 1
    assert all(rows[v] >> n == 0 and not rows[v] & (1 << v) for v in range(n))
    assert all(((rows[a] >> b) ^ (rows[b] >> a)) & 1 for a in range(n) for b in range(a))

    def extend(used, choices):
        if used == full:
            return 1
        paths = 0
        while choices:
            bit = choices & -choices
            choices ^= bit
            vertex = bit.bit_length() - 1
            available = full ^ (used | bit)
            count = extend(used | bit, rows[vertex] & available)
            if count % 2:
                predecessors[vertex] ^= used
            paths += count
        return paths

    height = extend(0, full)
    order = sorted(range(n), key=lambda v: predecessors[v].bit_count())
    assert [predecessors[v].bit_count() for v in order] == list(range(n))
    return height, order


def product_residues(modulus):
    # Exact breadth-first inverse; every new residue remembers an actual height factor.
    parent = {1: None}
    queue = [1]
    for source in queue:
        for factor in SEEDS:
            dest = source * factor % modulus
            if dest not in parent:
                parent[dest] = (source, factor)
                queue.append(dest)
    return parent


def recover(parent, target):
    word = []
    while parent[target] is not None:
        target, factor = parent[target]
        word.append(factor)
    return word[::-1]


def main():
    for h, rows in SEEDS.items():
        assert count_and_order(rows) == (h, list(range(len(rows))))
    for h, rows in LATE_RESIDUES.items():
        assert count_and_order(rows)[0] == h
    assert {h % 24 for h in LATE_RESIDUES} == {3, 11, 19}
    modulus = 2**6 * 3**2 * 11**2
    parent = product_residues(modulus)
    assert len(parent) == modulus // 2
    euclid = 2**60 * (2**61 - 1)
    assert all(gcd(euclid + 1, h) == 1 for h in SEEDS)
    print("height bank", sorted(SEEDS), "odd residues", len(parent), "modulus", modulus)
    for target in (euclid - 1, euclid + 1):
        word = recover(parent, target % modulus)
        value = 1
        for factor in word:
            value *= factor
        assert value % modulus == target % modulus
        print("target residue", target % modulus, "component heights", word, "exact height", value)


if __name__ == "__main__":
    main()
