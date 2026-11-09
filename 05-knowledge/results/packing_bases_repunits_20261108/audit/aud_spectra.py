"""aud_spectra.py - independent audit analysis for THM-4940 / HYP-9530 (audit of opus-2026-11-08-S86).

Written from scratch for the audit (does not import or copy the session's base_spectra.py / carry_free_layer.py).

Inputs (written by aud_fewcycle.c, the audit's own enumerator):
  FC   : every strong class with alpha_1 <= 80 (mode 0), lines "n a1 .. a7 H_ocf H_dp sc code"
  IC   : every strong intercyclic class with alpha_1 <= 1023 (mode 1), same format
Both runs ended with an EMPTY level, so by Moon's deletion lemma (+ monotonicity / heredity) they are complete at all
orders without using Moon's cycle-count floor.

Method for the spectra: H_b is multiplicative over order joins, and the strong components of a tournament with
H_b = 1 + b m (m <= 80) all have alpha_1 <= 80.  For every base b the attained alpha-sums m <= 80 are obtained by a
dynamic programme over the distinct strong values v = H_b(C) (an unbounded multiplicative knapsack, values capped at
1 + 80 b).  Every product is realised by an ordinal sum, so the set is exact.
usage: python aud_spectra.py FC IC > aud_spectra.out
"""
import sys
from collections import Counter, defaultdict

MB = 80


def load(path):
    rows = []
    with open(path) as f:
        for line in f:
            t = line.split()
            n = int(t[0])
            a = tuple(int(x) for x in t[1:8])
            H = int(t[8])
            sc = int(t[10])
            assert H == 1 + sum(c << (k + 1) for k, c in enumerate(a))
            rows.append((n, a, sc, t[11]))
    return rows


def hb(a, b):
    return 1 + sum(c * b ** (k + 1) for k, c in enumerate(a))


def attained(strong_vecs, b, mb=MB):
    """alpha-sums m <= mb attained in base b: 1 + b m is a product of strong values H_b(C) (multiplicative
    unbounded knapsack over m; every quotient of numbers = 1 mod b by a number = 1 mod b is again = 1 mod b)"""
    top = 1 + b * mb
    vals = sorted({hb(a, b) for a in strong_vecs if hb(a, b) <= top})
    ok = {0}
    for m in range(1, mb + 1):
        x = 1 + b * m
        for v in vals:
            if v > x:
                break
            if x % v == 0:
                y = x // v
                if (y - 1) % b == 0 and (y - 1) // b in ok:
                    ok.add(m)
                    break
    return ok


def one_apex_values(maxlen, cap):
    """odd-cycle counts of all one-apex tournaments (transitive chain v_1 -> ... -> v_m plus apex x; bit w_i = 1
    means x -> v_i).  Every odd cycle uses x: x -> v_i -> (increasing chain path) -> v_j -> x with w_i = 1, w_j = 0
    and an even number of chain vertices.  E[d] (number of such chain paths from i to j = i + d) is computed by a
    path DP, not taken from the closed form c(d)."""
    E = [0] * (maxlen + 1)
    for d in range(1, maxlen + 1):
        even = [0] * (d + 1)
        odd = [0] * (d + 1)
        odd[0] = 1                      # the path consisting of v_i alone has one chain vertex
        for j in range(1, d + 1):
            even[j] = sum(odd[k] for k in range(j))
            odd[j] = sum(even[k] for k in range(j))
        E[d] = even[d]
    vals = {0}
    for m in range(1, maxlen + 1):
        for mask in range(1 << m):
            ones = [i for i in range(m) if (mask >> i) & 1]
            zeros = [j for j in range(m) if not (mask >> j) & 1]
            tot = 0
            for i in ones:
                for j in zeros:
                    if j > i:
                        tot += E[j - i]
            if tot < cap:
                vals.add(tot)
    return vals, E


def main():
    FC = load(sys.argv[1])
    IC = load(sys.argv[2])
    print(f"strong classes with alpha_1 <= {MB} (own enumerator): {len(FC)}; per order:",
          sorted(Counter(n for n, *_ in FC).items()))
    strong_vecs = sorted({a for n, a, sc, c in FC})
    print(f"distinct strong packing vectors: {len(strong_vecs)}")
    mins = defaultdict(lambda: 10 ** 9)
    for n, a, sc, c in FC:
        mins[n] = min(mins[n], a[0])
    moon = {n: sum(n - l + 1 for l in range(3, n + 1, 2)) for n in range(3, 20)}
    print("min alpha_1 per order vs Moon floor m(n):", [(n, mins[n], moon[n]) for n in sorted(mins)],
          "; m(19) =", moon[19])

    # ---- A below 81 from FC, below 1024 from IC, vs Session 16's formula (one-apex DP + rays)
    A81 = {0} | {a[0] for n, a, sc, c in FC if a[1] == 0}
    A1024 = {0} | {a[0] for n, a, sc, c in IC}
    assert all(a[1] == 0 for n, a, sc, c in IC)
    assert A81 == {x for x in A1024 if x <= MB}, "the two inventories disagree on A below 81"
    oa, E = one_apex_values(14, 1024)
    print("chain-path counts E[d] (DP):", E[1:], " closed form c(d) = 1, 2^(d-2):",
          E[1:] == [1] + [2 ** (d - 2) for d in range(2, 15)])
    rays = {3 * 2 ** s for s in range(1, 12) if 3 * 2 ** s < 1024} | {9 * 2 ** s - 2 for s in range(0, 12) if 9 * 2 ** s - 2 < 1024}
    S16 = oa | rays
    print(f"A below 1024 (own intercyclic enumerator): {len(A1024)} values; Session 16 formula: {len(S16)}; equal: {A1024 == S16}")
    holes81 = [x for x in range(MB + 1) if x not in A81]
    holes1024 = [x for x in range(1024) if x not in A1024]
    print(f"holes of A in [0, 80] ({len(holes81)}): {holes81}")
    print(f"holes of A below 1024: {len(holes1024)}; first 45: {holes1024[:45]}")

    # ---- singletons of the carry-free layer
    mult = Counter(a[0] for n, a, sc, c in IC)
    scs = defaultdict(list)
    for n, a, sc, c in IC:
        scs[a[0]].append((n, sc))
    singles = sorted(x for x in mult if mult[x] == 1)
    print(f"carry-free singletons below 1024: {len(singles)}: {singles}")
    print("  orders:", [scs[x][0][0] for x in singles])
    print("  all self-converse:", all(scs[x][0][1] == 1 for x in singles))
    print("  multiplicity histogram:", sorted(Counter(mult.values()).items()))

    # ---- gap tables
    att = {}
    bases = list(range(2, 91)) + [100, 127, 128, 500, 1000, 10 ** 6]
    print()
    print(f"gap tables: alpha-sums m <= {MB} NOT attained in base b")
    for b in bases:
        att[b] = attained(strong_vecs, b)
        g = [m for m in range(MB + 1) if m not in att[b]]
        if b <= 20 or b in (80, 81, 85, 90) or b > 90:
            print(f"  b = {b}: gaps {g}" + ("   (= holes of A)" if g == holes81 else ""))
    print("  bases b in [81, 90] U {100, 127, 128, 500, 1000, 10^6} with gap set == holes of A below 81:",
          all([m for m in range(MB + 1) if m not in att[b]] == holes81 for b in bases if b >= 81))
    print("  bases b in [2, 80] whose gap set equals the holes of A below 81:",
          [b for b in range(2, 81) if [m for m in range(MB + 1) if m not in att[b]] == holes81])

    # ---- rescue sets
    print()
    print("rescue sets R(a) = {b <= 90 : 1 + a b attained} for the holes a <= 80 of A:")
    always23 = True
    for a in holes81:
        R = [b for b in range(2, 91) if a in att[b]]
        if a not in (3, 10) and not (2 in R and 3 in R):
            always23 = False
        if a in (3, 10, 14, 17):
            print(f"  a = {a}: {R}")
    print("  every hole a <= 80 other than 3, 10 is rescued in both bases 2 and 3:", always23)

    # ---- decompositions of alpha-sums 3 and 10 (packing vectors of ALL tournaments with alpha_1 <= 10)
    small = {(1,)}
    polys = [(1,) + tuple(a) for a in strong_vecs if a[0] <= 10]
    frontier = set(small)
    while frontier:
        new = set()
        for p in frontier:
            for q in polys:
                r = [0] * (len(p) + len(q) - 1)
                for i, x in enumerate(p):
                    for j, y in enumerate(q):
                        r[i + j] += x * y
                while r[-1] == 0:
                    r.pop()
                r = tuple(r)
                if r[1] <= 10 and r not in small:
                    new.add(r)
        small |= new
        frontier = new
    vecs = {tuple(list(p[1:]) + [0] * (7 - len(p[1:]))) for p in small}
    print()
    print(f"packing vectors of all tournaments with alpha_1 <= 10: {len(vecs)}")
    for msum in (3, 10):
        print(f"  alpha-sum {msum}:")
        for b in range(2, msum + 2):
            decs = []

            def rec(k, rem, cur):
                if k == 7:
                    if rem == 0:
                        decs.append(tuple(cur))
                    return
                w = b ** k
                for c in range(rem // w + 1):
                    rec(k + 1, rem - c * w, cur + [c])
            rec(0, msum, [])
            chain = [d for d in decs if all(d[0] >= k + 1 for k in range(1, 7) if d[k])]
            out = []
            for d in chain:
                out.append(f"{d[:3]}{' REALISED' if d in vecs else ''}")
            print(f"    b = {b}: {'; '.join(out)}")

    # ---- repunits 111_b: criterion and the attained/not list for b + 1 <= 80 from the tables
    print()
    ok = all(((b + 1) in att[b]) == ((b + 1) in A81) for b in range(2, 80))
    print("repunit criterion '111_b attained iff b + 1 in A' agrees with the exact tables for 2 <= b <= 79:", ok)
    fail = [b for b in range(2, 1023) if (b + 1) not in A1024]
    print("failing bases b (b + 1 a hole of A), b <= 60:", [b for b in fail if b <= 60])

    def isprime(x):
        if x < 2:
            return False
        i = 2
        while i * i <= x:
            if x % i == 0:
                return False
            i += 1
        return True
    rp = [b for b in range(2, 1023) if isprime(b * b + b + 1)]
    print("bases b <= 100 with 111_b prime:", [b for b in rp if b <= 100])
    print("  of these, failing (b + 1 not in A), b <= 100:", [(b, b * b + b + 1) for b in rp if b <= 100 and b in fail])
    print("  failing repunit primes with b <= 1022:", len([b for b in rp if b in fail]))

    # ---- witnesses
    print()
    w61 = [(n, c) for n, a, sc, c in FC if a[:3] == (6, 1, 0)]
    print("strong classes with alpha = (6, 1):", [(n, sc) for n, a, sc, c in FC if a[:3] == (6, 1, 0)], "H_4 =", hb((6, 1), 4), "H_2 =", hb((6, 1), 2))
    print("C3 => C3: alpha = (2, 1); H_8 =", hb((2, 1), 8), " H_4 =", hb((2, 1), 4), " I(Omega,1) =", hb((2, 1), 1))


if __name__ == "__main__":
    main()
