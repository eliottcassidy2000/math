"""Exact isomorphism-census boundary for the two-reservoir triangle cut."""
import sys
from pathlib import Path
from math import isqrt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from fermat_profiles import canon, extend, strong, profile, ends


def glue(U, W):
    m, n = len(U), len(W)
    return (tuple(row | (((1 << n) - 1) << m) for row in U)
            + tuple((row << m) | (1 << (m + n)) for row in W)
            + ((1 << m) - 1,))


def value(U, W):
    h, q = U
    k, r = W
    return (q - h) * (r - k) + 2 * h * k


def prime(n):
    return n > 1 and all(n % d for d in range(2, isqrt(n) + 1))


def main():
    classes = {(0,)}
    P = {}
    least = {3: (3, 9), 4: (5, 27), 5: (9, 75), 6: (15, 201),
             7: (25, 523), 8: (45, 1371)}
    for n in range(2, 9):
        classes = {canon(extend(A, c)) for A in classes for c in range(1 << (n - 1))}
        if n < 3:
            continue
        P[n] = {}
        for A in sorted(classes):
            if strong(A):
                P[n].setdefault(profile(A), A)
        H, K = least[n]
        if (min(h for h, q in P[n]) != H or min(q - h for h, q in P[n]) != K
                or (H, H + K) not in P[n]):
            raise RuntimeError(("joint floor", n))
    if len(classes) != 6880:
        raise RuntimeError("wrong class count")

    bounds = {(4, 7): 14371, (5, 6): 15345, (4, 8): 37467, (5, 7): 39675,
              (6, 6): 40851, (6, 7): 105873, (7, 7): 274779}
    for (m, n), N in bounds.items():
        h, K = least[m]
        k, L = least[n]
        U, W = P[m][(h, h + K)], P[n][(k, k + L)]
        if value((h, h + K), (k, k + L)) != N or sum(ends(glue(U, W))[-1]) != N:
            raise RuntimeError(("literal lower witness", m, n, N))
        print(m, n, N, flush=True)

    for m, n, expected in [(4, 7, {29167, 29477, 30347, 31657, 32837, 32911}),
                           (5, 6, set())]:
        h, a = least[m]
        k, b = least[n]
        q, r = h + a, k + b
        too_big = [value((h, 4*q), (k, r)), value((h, q), (k, 4*r)),
                   value((h, 2*q), (k, 2*r)), value((3*h, 12*q), (k, r)),
                   value((h, q), (3*k, 12*r))]
        if min(too_big) <= 34000:
            raise RuntimeError(("unaccounted pads or extra components", too_big))
        found = {}
        for U in P[m]:
            for W in P[n]:
                for x, y in [(0, 0), (0, 1), (1, 0)]:
                    N = value((U[0], U[1] << x), (W[0], W[1] << y))
                    if 29000 <= N <= 34000 and (prime(N) or (N % 7 == 0 and prime(N // 7))):
                        found[N] = U, W, x, y
        if set(found) != expected:
            raise RuntimeError(("nearby prime and 7p outputs", found))
        for N, (U, W, x, y) in sorted(found.items()):
            A = extend(P[m][U], 0) if x else P[m][U]
            B = extend(P[n][W], 0) if y else P[n][W]
            if sum(ends(glue(A, B))[-1]) != N:
                raise RuntimeError(("nearby literal witness", N))
            print("atom", m, n, N, U, W, x, y, A, B, flush=True)



if __name__ == "__main__":
    main()
