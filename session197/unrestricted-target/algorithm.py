"""Exact unrestricted target decision, using Busch's published strong-order floor.

This is a finite exhaustive algorithm, not a proof of the proposed {7,21} spectrum.
The binary-input running time is 2**O(b*b), not polynomial in b.
"""
from itertools import combinations


def strong_floor(n):
    if n < 3:
        raise ValueError('nontrivial strong order must be at least three')
    if n % 3 == 0:
        return 3 * 5 ** ((n - 3) // 3)
    if n % 3 == 1:
        return 5 ** ((n - 1) // 3)
    return 9 * 5 ** ((n - 5) // 3)


def reduced_floor(n):
    if n == 0:
        return 1  # bookkeeping for the empty component list, not an HP convention
    if n < 3:
        return None
    if n == 5:
        return 9
    if n % 3 == 0:
        return 3 ** (n // 3)
    if n % 3 == 1:
        return 5 * 3 ** ((n - 4) // 3)
    return 25 * 3 ** ((n - 8) // 3)


def reduced_bound(h):
    """B=max{n:3**n<=h**3}; callers handle h=1 separately."""
    if h < 1:
        raise ValueError('positive target required')
    n, power = 0, 1
    while 3 * power <= h ** 3:
        n += 1
        power *= 3
    return n


def strong_bound(h):
    n = 3
    while strong_floor(n) <= h:
        n += 1
    return n - 1 if n > 3 else 0


def tournaments(n):
    pairs = list(combinations(range(n), 2))
    for mask in range(1 << len(pairs)):
        rows = [0] * n
        for bit, (u, v) in enumerate(pairs):
            if mask >> bit & 1:
                u, v = v, u
            rows[u] |= 1 << v
        yield tuple(rows)


def is_strong(rows):
    if not rows:
        return False
    n = len(rows)
    full = (1 << n) - 1
    for start in range(n):
        seen = 1 << start
        while True:
            after = seen
            for v in range(n):
                if seen >> v & 1:
                    after |= rows[v]
            if after == seen:
                break
            seen = after
        if seen != full:
            return False
    return True


def hp(rows, cap=None):
    """Forward endpoint DP; cap is the overflow value, usually target+1."""
    n = len(rows)
    if not n:
        return 0
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(full + 1)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, full + 1):
        for v, value in enumerate(dp[mask]):
            if not value:
                continue
            choices = rows[v] & (full ^ mask)
            while choices:
                bit = choices & -choices
                choices -= bit
                target = dp[mask | bit]
                j = bit.bit_length() - 1
                target[j] += value
                if cap is not None:
                    target[j] = min(cap, target[j])
    value = sum(dp[full])
    return value if cap is None else min(cap, value)


def strong_bank(limit):
    """All strong path counts <=limit; one physical representative per count."""
    bank = {}
    for n in range(3, strong_bound(limit) + 1):
        for rows in tournaments(n):
            if is_strong(rows):
                h = hp(rows, limit + 1)
                if h <= limit and h not in bank:
                    bank[h] = rows
    return bank


def ordered_join(blocks):
    n = sum(map(len, blocks))
    rows, offset = [], 0
    for block in blocks:
        later = ((1 << n) - 1) ^ ((1 << (offset + len(block))) - 1)
        rows.extend((r << offset) | later for r in block)
        offset += len(block)
    return tuple(rows)


def from_bank(h, bank):
    if h == 1:
        return (0,)
    if h < 1 or h % 2 == 0:
        return None
    # Every intermediate value divides the target; all nontrivial factors >=3.
    reached = {1: []}
    pending = [1]
    candidates = sorted(c for c in bank if c > 1 and h % c == 0)
    while pending:
        value = pending.pop()
        for c in candidates:
            child = value * c
            if child <= h and h % child == 0 and child not in reached:
                reached[child] = reached[value] + [c]
                if child == h:
                    return ordered_join([bank[x] for x in reached[child]])
                pending.append(child)
    return None


def decide(h):
    if not isinstance(h, int):
        raise TypeError('integer target required')
    if h == 1:
        return (0,)
    if h < 1 or h % 2 == 0:
        return None
    return from_bank(h, strong_bank(h))


if __name__ == '__main__':
    import argparse, json
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('target', type=int)
    args = p.parse_args()
    rows = decide(args.target)
    print(json.dumps({'target': args.target, 'attainable': rows is not None,
                      'rows': rows}, sort_keys=True))
