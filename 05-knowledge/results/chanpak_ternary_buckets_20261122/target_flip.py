"""Look up a prescribed even midpoint and recount both tournaments of its one-arc flip."""
import argparse
from pathlib import Path

TABLE = Path(__file__).with_name("target_flips.rows")
MIDPOINTS = range(264, 1503, 2)


def load_table():
    certificates = {}
    for line_number, line in enumerate(TABLE.read_text().splitlines(), 1):
        target, u, v, *rows = map(int, line.split())
        n = len(rows)
        if target in certificates or not 0 <= u < n or not 0 <= v < n:
            raise ValueError(f"duplicate midpoint or invalid arc on line {line_number}")
        if any(row < 0 or row >= 1<<n or row>>i & 1 for i, row in enumerate(rows)):
            raise ValueError(f"invalid row on line {line_number}")
        for i in range(n):
            for j in range(i):
                if (rows[i]>>j & 1) == (rows[j]>>i & 1):
                    raise ValueError(f"invalid pair on line {line_number}")
        if not rows[u]>>v & 1:
            raise ValueError(f"missing directed arc on line {line_number}")
        certificates[target] = u, v, rows
    if set(certificates) != set(MIDPOINTS):
        raise ValueError("the certificate midpoints must be exactly 264,266,...,1502")
    return certificates


def height(rows):
    # Extend actual paths one outgoing arc at a time. No profile or flip identities.
    full = (1<<len(rows))-1
    ways = [[0]*len(rows) for _ in range(full+1)]
    for i in range(len(rows)):
        ways[1<<i][i] = 1
    for used, counts in enumerate(ways):
        for i, count in enumerate(counts):
            if not count:
                continue
            available = rows[i] & ~used
            while available:
                bit = available & -available
                available ^= bit
                ways[used|bit][bit.bit_length()-1] += count
    return sum(ways[full])


def flip_and_check(target, certificate):
    u, v, lower = certificate
    upper = lower.copy()
    upper[u] ^= 1<<v
    upper[v] ^= 1<<u
    endpoints = height(lower), height(upper)
    if endpoints != (target-1, target+1):
        raise ValueError(f"midpoint {target} has wrong endpoint heights {endpoints}")
    return upper


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.add_mutually_exclusive_group(required=True)
    args.add_argument("midpoint", type=int, nargs="?")
    args.add_argument("--all", action="store_true", help="recount all 620 certificates")
    options = parser.parse_args()
    certificates = load_table()
    if options.all:
        for target, certificate in certificates.items():
            flip_and_check(target, certificate)
        print(f"Recounted {len(certificates)} one-arc certificates, midpoints 264 through 1502")
        return
    target = options.midpoint
    if target not in certificates:
        parser.error("midpoint must be even and between 264 and 1502")
    u, v, lower = certificates[target]
    upper = flip_and_check(target, certificates[target])
    print(f"Midpoint {target}, order {len(lower)}, reverse {u} -> {v}")
    print(f"H={target-1} rows", *lower)
    print(f"H={target+1} rows", *upper)


if __name__ == "__main__":
    main()
