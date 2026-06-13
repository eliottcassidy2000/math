"""Fast class-SIZE-histogram counter for signed-LRC collisions (monad-explorer-S708).
Same Phi^2 hashing as s707g but records the multiset of class sizes, so we can read off
the |G_eps| distribution that drives the count law. Validated vs exact (s708) for C<=39.
"""
import numpy as np, math
from collections import Counter


def is_prime(C):
    return C > 1 and all(C % d for d in range(2, int(C**0.5) + 1))


def factor(C):
    f, m, d = [], C, 2
    while d * d <= m:
        while m % d == 0:
            f.append(d); m //= d
        d += 1
    if m > 1:
        f.append(m)
    return "x".join(map(str, f))


def size_hist(n, batch=1 << 18, ndec=7):
    C = 2 * n - 1
    n1 = n - 1
    H = np.arange(1, n1 + 1)
    S = np.sin(2 * math.pi * np.outer(H, H) / C)
    nfree = n1 - 1
    total = 1 << nfree
    counts = Counter()                      # signature bytes -> #cuts
    for start in range(0, total, batch):
        end = min(start + batch, total)
        idx = np.arange(start, end, dtype=np.uint64)
        eps = np.ones((n1, end - start), dtype=np.int8)
        for b in range(nfree):
            bit = ((idx >> np.uint64(b)) & np.uint64(1)).astype(np.int8)
            eps[b + 1, :] = np.where(bit == 1, 1, -1)
        Phi = S @ eps
        Phi2 = np.round(Phi * Phi, ndec)
        for row in Phi2.T:
            counts[row.tobytes()] += 1
    sizehist = Counter(counts.values())     # class-size -> #classes
    orbit = len(counts)
    defic = total - orbit
    return C, total, orbit, defic, dict(sorted(sizehist.items()))


if __name__ == "__main__":
    print(f"{'n':>3} {'C':>4} {'fact':>9} {'2^(n-2)':>9} {'orbit':>9} {'defic':>7}  size-hist(size:#classes)")
    for n in range(5, 26):
        C = 2 * n - 1
        if is_prime(C) or (n - 2) > 23:
            continue
        Cc, tot, orb, df, sh = size_hist(n)
        print(f"{n:>3} {C:>4} {factor(C):>9} {tot:>9} {orb:>9} {df:>7}  {sh}", flush=True)
