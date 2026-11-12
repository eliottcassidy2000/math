#!/usr/bin/env python3
"""Lane thinorbit (S93), finishing pass.  The invariant plane h = O of the port-chain block library.

structure_scan.py (4) found invertible blocks sharing the invariant plane with normal (1, 0, -1), i.e. blocks N with
(h' - O') = lam (h - O).  The start state (1, 0, 1) lies on the plane h = O, so the blocks preserving it generate a
genuinely 2 x 2 system on the coordinates (h, E) (O = h):  (h, E) -> A (h, E),
    A = [[N11 + N13, N12], [N21 + N23, N22]].
This script (exact integer arithmetic) lists, for a library file with .reps:
 (1) every block with (1, 0, -1) N = lam (1, 0, -1), its lam, A, det A, and for |det A| = 1 the continued-fraction
     word (a_1, ..., a_r) with A = M(a_1) ... M(a_r), M(a) = [[a, 1], [1, 0]] (so that h is a continuant);
 (2) every block mapping all of Z^3 into the plane (lam = 0);
 (3) the orbit of (h, E) = (1, 0) under the plane-preserving blocks: number of distinct states and of odd values
     with h <= 2^k (growth exponent and coverage), for all such blocks and for the unimodular ones only.
Usage: finish_plane_scan.py MATFILE [KMAX_BITS]
"""
import sys


def cf_word(A):
    """A = [[p, q], [r, s]] nonnegative, det +-1: word (a_1..a_r), a_i >= 0, with A = M(a_1) ... M(a_r); None if stuck."""
    p, q, r, s = A
    word = []
    guard = 0
    while (p, q, r, s) != (1, 0, 0, 1):
        guard += 1
        if guard > 200:
            return None
        # A = M(a) A'  with A' = [[r, s], [p - a r, q - a s]] nonnegative
        if r == 0 and s == 0:
            return None
        cands = []
        if r > 0:
            cands.append(p // r)
        if s > 0:
            cands.append(q // s)
        a = min(cands)
        # keep A' unimodular nonnegative and not degenerate: if the remainder row would be (0, 0) reduce a
        while a >= 0 and (p - a * r, q - a * s) == (0, 0):
            a -= 1
        if a < 0:
            return None
        word.append(a)
        p, q, r, s = r, s, p - a * r, q - a * s
        if min(p, q, r, s) < 0:
            return None
    return word


def main():
    path = sys.argv[1]
    kbits = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    mats, reps = [], []
    for line in open(path + ".reps"):
        a, b = line.split("|")
        mats.append(tuple(map(int, a.split())))
        reps.append(b.strip())
    print(f"{len(mats)} block matrices from {path}")
    keep = []
    into = []
    for i, N in enumerate(mats):
        f = (N[0] - N[6], N[1] - N[7], N[2] - N[8])            # (1,0,-1) N
        if f[1] != 0 or f[0] != -f[2]:
            continue
        lam = f[0]
        det = (N[0] * (N[4] * N[8] - N[5] * N[7]) - N[1] * (N[3] * N[8] - N[5] * N[6]) + N[2] * (N[3] * N[7] - N[4] * N[6]))
        A = (N[0] + N[2], N[1], N[3] + N[5], N[4])
        dA = A[0] * A[3] - A[1] * A[2]
        if lam == 0:
            into.append(i)
        keep.append((i, lam, A, dA, det))
    print(f"(1) blocks with (1,0,-1) N = lam (1,0,-1): {len(keep)}  (lam = 0, image inside the plane: {len(into)})")
    seenA = {}
    for (i, lam, A, dA, det) in keep:
        seenA.setdefault(A, []).append((i, lam, det))
    print(f"    distinct induced 2 x 2 matrices A on (h, E): {len(seenA)}")
    uni = []
    for A in sorted(seenA, key=lambda A: (abs(A[0] * A[3] - A[1] * A[2]) != 1, A)):
        dA = A[0] * A[3] - A[1] * A[2]
        tag = ""
        if abs(dA) == 1:
            uni.append(A)
            tag = f"  continuant word {cf_word(A)}"
        ex = seenA[A][0]
        print(f"    A = {list(A)} det {dA:5d} trace {A[0] + A[3]:4d}{tag}   lam/det(N) of the blocks: "
              f"{sorted(set((l, d) for (_, l, d) in seenA[A]))}  e.g. block [{reps[ex[0]]}] (x{len(seenA[A])})")
    print(f"    unimodular induced matrices: {len(uni)}")
    # (3) orbits
    LIM = 1 << kbits
    for name, gens in (("all plane-preserving blocks", sorted(seenA)), ("unimodular induced matrices only", uni)):
        if not gens:
            print(f"(3) {name}: none")
            continue
        seen = {(1, 0)}
        stack = [(1, 0)]
        while stack:
            h, e = stack.pop()
            for (p, q, r, s) in gens:
                h2 = p * h + q * e
                if h2 > LIM:
                    continue
                st = (h2, r * h + s * e)
                if st not in seen:
                    seen.add(st)
                    stack.append(st)
                    if len(seen) > 30_000_000:
                        stack = []
                        break
        cnt = [0] * (kbits + 2)
        vals = set()
        for (h, e) in seen:
            k = max(0, (h - 1).bit_length())
            cnt[k] += 1
            vals.add(h)
        cum = 0
        row = []
        for k in range(kbits + 1):
            cum += cnt[k]
            row.append(cum)
        print(f"(3) {name} ({len(gens)} generators): distinct states (h, E) with h <= 2^k, k = 0..{kbits}:")
        print("    " + " ".join(map(str, row)))
        cov = []
        for k in range(1, kbits):
            lo, hi = 1 << k, (2 << k) - 1
            c = sum(1 for v in vals if lo < v <= hi and v % 2 == 1)
            cov.append(f"{k}:{c}/{1 << (k - 1)}")
        print("    odd values reached in (2^k, 2^(k+1)): " + " ".join(cov))
        print("    even values reached:", sum(1 for v in vals if v % 2 == 0))


if __name__ == "__main__":
    main()
