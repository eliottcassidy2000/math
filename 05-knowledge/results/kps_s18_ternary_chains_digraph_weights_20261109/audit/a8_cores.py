#!/usr/bin/env python3
"""AUDIT: THM-5071 (5). Strong loopless semicomplete digraphs of order <= 5 with H in {7, 21}: isomorphism classes
and their weighted packing decompositions (THM-5071 (2)), recomputed independently.

H by counting permutations. Cycles of D by brute force over vertex subsets and cyclic orders; kinds: 'B' (uses a
reciprocal pair, weight 1), 'To' (tournament pairs only, odd, weight 2), 'Te' (tournament pairs only, even, weight 0).
For semicomplete D the complement has no cycle through a missing pair, so these are all admissible cycles.
Each class prints b, the cycle-kind counts, and H as a sum over packings grouped by the multiset of kinds.
Compares with the packet's out/core_weights.out (matched by canonical form).
Usage: a8_cores.py PACKET_CORE_WEIGHTS_OUT
"""
import itertools
import sys


def H(n, A):
    return sum(1 for p in itertools.permutations(range(n)) if all(A[p[i]][p[i + 1]] for i in range(n - 1)))


def strong(n, A):
    def reach(s, fwd):
        seen, st = {s}, [s]
        while st:
            v = st.pop()
            for t in range(n):
                if t not in seen and (A[v][t] if fwd else A[t][v]):
                    seen.add(t)
                    st.append(t)
        return len(seen) == n
    return reach(0, True) and reach(0, False)


def canon(n, A):
    return min(tuple(A[p[i]][p[j]] for i in range(n) for j in range(n)) for p in itertools.permutations(range(n)))


def cycles(n, A):
    out = []
    for k in range(2, n + 1):
        for sub in itertools.combinations(range(n), k):
            for perm in itertools.permutations(sub[1:]):
                c = (sub[0],) + perm
                if all(A[c[i]][c[(i + 1) % k]] for i in range(k)):
                    usesB = any(A[c[(i + 1) % k]][c[i]] for i in range(k))
                    kind = "B" if usesB else ("To" if k % 2 else "Te")
                    out.append((sum(1 << v for v in c), kind))
    return out


def decomposition(n, A):
    cyc = [(m, kd) for (m, kd) in cycles(n, A)]
    w = {"B": 1, "To": 2, "Te": 0}
    terms = {}

    def rec(i, used, prod, kinds):
        if i == len(cyc):
            if prod:
                key = "+".join(sorted(kinds)) or "empty"
                terms[key] = terms.get(key, 0) + prod
            return
        rec(i + 1, used, prod, kinds)
        m, kd = cyc[i]
        if not m & used and w[kd]:
            rec(i + 1, used | m, prod * w[kd], kinds + [kd])

    rec(0, 0, 1, [])
    counts = {}
    for _, kd in cyc:
        counts[kd] = counts.get(kd, 0) + 1
    return terms, counts


def main():
    theirs = {}
    for line in open(sys.argv[1]):
        parts = line.split()
        n = int(parts[0][2:])
        rows = parts[3][5:].split("/")
        A = [[int(c) for c in r] for r in rows]
        theirs[canon(n, A)] = line.strip()
    found = 0
    agree = 0
    for n in range(2, 6):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        classes = {}
        for code in range(3 ** len(pairs)):
            A = [[0] * n for _ in range(n)]
            c = code
            for (i, j) in pairs:
                s = c % 3 + 1
                c //= 3
                A[i][j], A[j][i] = s & 1, s >> 1
            h = H(n, A)
            if h in (7, 21) and strong(n, A):
                key = canon(n, A)
                classes.setdefault(key, (h, [r[:] for r in A]))
        tally = {}
        for key, (h, A) in sorted(classes.items(), key=lambda t: t[1][0]):
            terms, counts = decomposition(n, A)
            assert sum(terms.values()) == h
            b = sum(1 for i in range(n) for j in range(i + 1, n) if A[i][j] and A[j][i])
            tally[(h, b)] = tally.get((h, b), 0) + 1
            found += 1
            mine = f"n={n} H={h} b={b} kinds={dict(sorted(counts.items()))} H = " + " + ".join(
                f"{v}[{k}]" for k, v in sorted(terms.items(), key=lambda t: (t[0].count("+"), t[0])))
            print(mine)
            if key in theirs:
                t = theirs[key]
                # compare numbers: their cycle counts and term values
                ok = (f"b={b}" in t and f"'B': {counts.get('B', 0)}" in t if counts.get('B') else True)
                ok = ok and all(f"{v}[" in t for v in terms.values())
                agree += ok
                if not ok:
                    print("   DISAGREES WITH PACKET:", t)
            else:
                print("   NOT IN PACKET LIST")
        print(f"order {n}: classes by (H, b): {dict(sorted(tally.items()))}")
    print(f"total classes {found}; packet lines {len(theirs)}; matched with consistent numbers {agree}")


if __name__ == "__main__":
    main()
