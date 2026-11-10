#!/usr/bin/env python3
"""THM-5125: the quadratic-residue (Redei) tournament model.  Exact checks.

Orientation used here:  p -> q  iff  p is a quadratic residue modulo q   (p, q distinct primes = 3 mod 4).
Quadratic reciprocity for two primes = 3 mod 4 makes this a tournament.

(1) Universality (constructive check of the Dirichlet/CRT argument): every tournament on at most 5 vertices is the
    residue tournament of a set of primes = 3 mod 4 (search), and likewise with primes = 3 mod 8.
(2) Heegner domination: for q in {3, 7, 11, 19, 43, 67, 163} and every prime p = 3 mod 4 with p < (q+1)/4, q -> p.
(3) The Heegner tournament on {3, 7, 11, 19, 43, 67, 163}: arcs, scores, H of every sub-tournament, prefix values,
    Monsky nullity for the six primes = 3 mod 8.
(4) Least realisers: a(N) = least squarefree n with all prime factors = 3 mod 4 whose residue tournament has
    exactly N Hamiltonian paths ("N residue chains"), by exhaustive search inside a stated box.
Run:  python3 redei_model.py > redei_model.out
"""
import itertools
import sys
from math import prod


def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def qr_tournament(primes):
    n = len(primes)
    beats = [0] * n
    for i, p in enumerate(primes):
        for j, q in enumerate(primes):
            if i != j and legendre(p, q) == 1:
                beats[i] |= 1 << j
    for i in range(n):
        for j in range(i + 1, n):
            if ((beats[i] >> j) & 1) + ((beats[j] >> i) & 1) != 1:
                raise RuntimeError("reciprocity must give a tournament")
    return beats


def hs_table(n, beats):
    size = 1 << n
    cnt = [[0] * n for _ in range(size)]
    for v in range(n):
        cnt[1 << v][v] = 1
    Hs = [0] * size
    Hs[0] = 1
    for mask in range(1, size):
        row = cnt[mask]
        tot = 0
        for v in range(n):
            c = row[v]
            if c:
                tot += c
                nxt = beats[v] & ~mask
                while nxt:
                    w = nxt & -nxt
                    cnt[mask | w][w.bit_length() - 1] += c
                    nxt ^= w
        Hs[mask] = tot
    return Hs


def H(n, beats):
    return hs_table(n, beats)[(1 << n) - 1]


def canon(n, beats):
    best = None
    for perm in itertools.permutations(range(n)):
        code = 0
        k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (beats[perm[i]] >> perm[j]) & 1:
                    code |= 1 << k
                k += 1
        if best is None or code < best:
            best = code
    return best


def rank_f2(rows, n):
    rows = list(rows)
    r = 0
    for col in range(n):
        piv = None
        for i in range(r, len(rows)):
            if rows[i] >> col & 1:
                piv = i
                break
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> col & 1:
                rows[i] ^= rows[r]
        r += 1
    return r


def mul_f2(A, B, n):
    C = []
    for i in range(n):
        row = 0
        a = A[i]
        k = 0
        while a:
            if a & 1:
                row ^= B[k]
            a >>= 1
            k += 1
        C.append(row)
    return C


def monsky_nullity(primes):
    """nullity over F_2 of M^2 + M + I, M = adjacency + diag(row parities), with i -> j iff p_j is a quadratic
    residue modulo p_i (Legendre symbol +1).  The cited note 05-knowledge/results/bsd_selmer_growth_20261013a.md
    section 3 (primes = 3 mod 8) writes its condition with an additive bracket symbol; the nullity is the same for
    either orientation (THM-5125 (4c); checked by the audit on every tournament of orders 3-6)."""
    n = len(primes)
    A = [0] * n
    for i, p in enumerate(primes):
        for j, q in enumerate(primes):
            if i != j and legendre(q, p) == 1:
                A[i] |= 1 << j
    M = [A[i] ^ ((bin(A[i]).count("1") & 1) << i) for i in range(n)]
    M2 = mul_f2(M, M, n)
    return n - rank_f2([M2[i] ^ M[i] ^ (1 << i) for i in range(n)], n)


def main():
    ok = True
    # ---------------- (1) universality, constructive ----------------
    for modulus, res in ((4, 3), (8, 3)):
        pool = [p for p in range(3, 6000) if p % modulus == res and is_prime(p)]
        for n in range(1, 6):
            classes = set()
            pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
            for bits in range(1 << len(pairs)):
                beats = [0] * n
                for k, (i, j) in enumerate(pairs):
                    if bits >> k & 1:
                        beats[i] |= 1 << j
                    else:
                        beats[j] |= 1 << i
                classes.add(canon(n, beats))
            found = set()

            def extend(chosen):
                if len(found) == len(classes):
                    return
                if len(chosen) == n:
                    found.add(canon(n, qr_tournament(chosen)))
                    return
                start = pool.index(chosen[-1]) + 1 if chosen else 0
                seen_patterns = set()
                for p in pool[start:start + 500]:
                    pat = tuple(legendre(p, q) for q in chosen)
                    if pat in seen_patterns:
                        continue          # one representative prime per residue pattern suffices
                    seen_patterns.add(pat)
                    extend(chosen + [p])
                    if len(seen_patterns) == 1 << len(chosen):
                        break

            extend([])
            good = (found == classes)
            ok &= good
            print(f"(1) primes = {res} mod {modulus}: order {n}: {len(found)} of {len(classes)} isomorphism classes realised",
                  "OK" if good else "FAIL")

    # ---------------- (2) Heegner domination ----------------
    heeg = [3, 7, 11, 19, 43, 67, 163]
    for q in heeg:
        small = [p for p in range(3, q) if p % 4 == 3 and is_prime(p) and 4 * p < q + 1]
        good = all(legendre(q, p) == 1 for p in small)
        ok &= good
        tag = " (prime)" if is_prime((q + 1) // 4) else ""
        print(f"(2) q = {q}: (q+1)/4 = {(q + 1) // 4}{tag}; primes p = 3 mod 4 below it: {small}; q is a residue modulo each:",
              "OK" if good else "FAIL")

    # ---------------- (3) the Heegner tournament ----------------
    n = len(heeg)
    beats = qr_tournament(heeg)
    print("(3) Heegner tournament, p -> q iff p is a quadratic residue mod q:")
    for i, p in enumerate(heeg):
        print(f"      {p:3d} -> {[heeg[j] for j in range(n) if beats[i] >> j & 1]}")
    Hs = hs_table(n, beats)
    print("    scores:", [bin(b).count("1") for b in beats], " H =", Hs[-1])
    print("    prefix values H({3}), H({3,7}), ...:", [Hs[(1 << k) - 1] for k in range(1, n + 1)])
    by = {}
    for m in range(1, 1 << n):
        by.setdefault(bin(m).count("1"), {}).setdefault(Hs[m], []).append(m)
    for k in sorted(by):
        print(f"    sub-tournaments of size {k}: " + ", ".join(f"H={h} x{len(v)}" for h, v in sorted(by[k].items())))
    six = [p for p in heeg if p % 8 == 3]
    m6 = sum(1 << i for i, p in enumerate(heeg) if p % 8 == 3)
    print(f"    the six primes = 3 mod 8 {six}: product {prod(six)}, H = {Hs[m6]}, Monsky nullity = {monsky_nullity(six)}"
          " (pure 2-Selmer dimension of y^2 = x^3 - n^2 x by the cited formula)")
    print("    controls: 627 = 3*11*19 ->", H(3, qr_tournament([3, 11, 19])), "paths, Monsky nullity", monsky_nullity([3, 11, 19]),
          "; 1947 = 3*11*59 ->", H(3, qr_tournament([3, 11, 59])), "paths, Monsky nullity", monsky_nullity([3, 11, 59]))
    ok &= (Hs[-1] == 71 and Hs[m6] == 29 and monsky_nullity([3, 11, 19]) == 2 and monsky_nullity([3, 11, 59]) == 0)

    # ---------------- (4) least realisers ----------------
    BOUND = 3 * 10 ** 9
    pool = [p for p in range(3, 600) if p % 4 == 3 and is_prime(p)]
    best = {}

    def search(start, chosen, product):
        if chosen:
            h = H(len(chosen), qr_tournament(chosen))
            if h not in best or product < best[h][0]:
                best[h] = (product, list(chosen))
        if len(chosen) >= 7:
            return
        for idx in range(start, len(pool)):
            p = pool[idx]
            if product * p > BOUND:
                break
            search(idx + 1, chosen + [p], product * p)

    search(0, [], 1)
    # least order at which each count occurs (exhaustive, orders <= 6)
    least_order = {}
    for n in range(1, 7):
        prs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        for bits in range(1 << len(prs)):
            bt = [0] * n
            for k, (i, j) in enumerate(prs):
                if bits >> k & 1:
                    bt[i] |= 1 << j
                else:
                    bt[j] |= 1 << i
            least_order.setdefault(H(n, bt), n)
    smallest = [p for p in range(3, 200) if p % 4 == 3 and is_prime(p)]
    print(f"(4) least n < {BOUND} (squarefree, prime factors = 3 mod 4 and < 600, at most 7 of them) with exactly N residue chains.")
    print("    'global' = certified least over ALL such n: the count needs k = (number of prime factors found) vertices")
    print("    (exhaustive order spectra, orders <= 6) and n < 600 * (product of the k-1 smallest primes = 3 mod 4),")
    print("    so every smaller candidate lies inside the search box; 'box' = least inside the box only.")
    for N in sorted(best):
        if N <= 45:
            n_found, fac = best[N]
            k = len(fac)
            need = least_order.get(N, 7)
            cert = (need == k) and (n_found < 600 * prod(smallest[:k - 1]))
            print(f"      N = {N:3d}: n = {n_found} = {'*'.join(map(str, fac))}   [{'global' if cert else 'box'}]")
    missing = [N for N in range(1, 46, 2) if N not in best]
    print("    odd N <= 45 without a realiser in this search box:", missing, "(7 and 21 are never values, THM-1370)")
    ok &= (7 in missing and 21 in missing)
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
