#!/usr/bin/env python3
"""Independent audit of THM-5125 statements 1-3 (residue tournament, Heegner domination, Heegner tournament,
least realisers).  Own code: Jacobi symbols by the binary reciprocity algorithm (the packet uses Euler's criterion),
Hamiltonian paths by depth-first enumeration and, for the Heegner tournament, by brute-force permutations.

Orientation as in the packet:  p -> q  iff  p is a quadratic residue mod q  (p, q primes = 3 mod 4).

(1) extension property, spot check: for the anchors (3, 7, 11, 19) every one of the 16 arc patterns is realised by a
    prime r = 3 mod 8 (least witness listed); the CRT class used in the proof is checked to be coprime to the modulus.
(2) Heegner domination: q in {3, 7, 11, 19, 43, 67, 163}; every prime p < (q+1)/4 is inert ((-q/p) = -1, p odd;
    p = 2: q = 3 mod 8); for p = 3 mod 4 this is q -> p.  Control: the statement fails for every non-Heegner prime
    q = 3 mod 4, 23 <= q < 2000, that has a prime p = 3 mod 4 below (q+1)/4 ... (listed count).
(3) Heegner tournament: arcs, scores, triangles, H = 71 (permutations), prefixes, six primes = 3 mod 8, value set of
    sub-tournaments, Monsky nullity in both orientations and from the full 2t x 2t Monsky matrix of the cited note.
(4) least realisers a(N), N <= 45: (a) fully exhaustive over EVERY squarefree n <= 400,000 with prime factors
    = 3 mod 4 (no box, no least-order input); (b) for each N, every n < a(N) with at least k(N) prime factors
    (k(N) = least order realising N, from an own exhaustive spectrum of orders <= 6) is enumerated with no bound on
    the size of the primes and shown not to have N residue chains.
"""
import itertools
import sys
import time
from math import gcd, prod


def jacobi(a, n):
    assert n > 0 and n % 2 == 1
    a %= n
    result = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def primes_upto(N):
    sieve = bytearray([1]) * (N + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i:: i] = bytearray(len(sieve[i * i:: i]))
    return [i for i in range(N + 1) if sieve[i]]


def arcs(primes):
    """out[i] = bitmask of j with p_i -> p_j (p_i a QR mod p_j)."""
    k = len(primes)
    out = [0] * k
    for i in range(k):
        for j in range(k):
            if i != j and jacobi(primes[i], primes[j]) == 1:
                out[i] |= 1 << j
    for i in range(k):
        for j in range(i + 1, k):
            assert ((out[i] >> j) & 1) + ((out[j] >> i) & 1) == 1, "not a tournament"
    return out


def hp_dfs(out, mask=None):
    """Hamiltonian paths of the sub-tournament on mask, depth-first."""
    k = len(out)
    if mask is None:
        mask = (1 << k) - 1
    if mask == 0:
        return 1
    total = 0

    def rec(v, used):
        nonlocal total
        if used == mask:
            total += 1
            return
        c = out[v] & mask & ~used
        while c:
            low = c & -c
            c ^= low
            rec(low.bit_length() - 1, used | low)

    m = mask
    while m:
        low = m & -m
        m ^= low
        rec(low.bit_length() - 1, low)
    return total


def hp_perm(out):
    k = len(out)
    return sum(1 for p in itertools.permutations(range(k)) if all((out[p[i]] >> p[i + 1]) & 1 for i in range(k - 1)))


def rank_gf2(rows):
    rows = [r for r in rows]
    rank = 0
    while rows:
        piv = rows.pop()
        if piv == 0:
            continue
        rank += 1
        low = piv & -piv
        rows = [r ^ piv if r & low else r for r in rows]
    return rank


def mat_mul_gf2(A, B):
    out = []
    for a in A:
        r = 0
        k = 0
        while a:
            if a & 1:
                r ^= B[k]
            a >>= 1
            k += 1
        out.append(r)
    return out


def phi3_nullity(M):
    n = len(M)
    M2 = mat_mul_gf2(M, M)
    return n - rank_gf2([M2[i] ^ M[i] ^ (1 << i) for i in range(n)])


def monsky_variants(primes):
    """(nullity with A_ij = [(p_j/p_i) = +1], nullity with A_ij = [(p_j/p_i) = -1] (cited note), nullity of the
    full 2t x 2t Monsky matrix [[A + D_2, D_2], [D_2, A + D_-2]] with the cited note's A)."""
    t = len(primes)
    res = []
    for sign in (1, -1):
        A = [0] * t
        for i in range(t):
            for j in range(t):
                if i != j and jacobi(primes[j], primes[i]) == sign:
                    A[i] |= 1 << j
            if bin(A[i]).count("1") & 1:
                A[i] |= 1 << i
        res.append(phi3_nullity(A))
        if sign == -1:
            D2 = [1 if jacobi(2, p) == -1 else 0 for p in primes]
            Dm2 = [1 if jacobi(-2 % p, p) == -1 else 0 for p in primes]
            big = []
            for i in range(t):
                big.append((A[i] ^ (D2[i] << i)) | ((D2[i] << i) << t))
            for i in range(t):
                big.append((D2[i] << i) | ((A[i] ^ (Dm2[i] << i)) << t))
            res.append(2 * t - rank_gf2(big))
    return tuple(res)


def main():
    t0 = time.time()
    ok = True

    # sanity of the Jacobi routine against Euler's criterion on a few hundred pairs
    ps = [p for p in primes_upto(400) if p > 2]
    bad = sum(1 for p in ps for a in range(1, 60) if a % p and jacobi(a, p) != (1 if pow(a, (p - 1) // 2, p) == 1 else -1))
    print("(0) Jacobi routine vs Euler criterion:", "OK" if bad == 0 else "FAIL")
    ok &= bad == 0

    # ---- (1) extension property, spot check
    anchors = [3, 7, 11, 19]
    P3mod8 = [p for p in primes_upto(200000) if p % 8 == 3]
    found = {}
    for r in P3mod8:
        if r in anchors:
            continue
        pat = tuple(jacobi(r, u) for u in anchors)     # +1: r -> u, -1: u -> r
        found.setdefault(pat, r)
        if len(found) == 16:
            break
    good = len(found) == 16
    print(f"(1) anchors {anchors}: all 16 arc patterns realised by primes r = 3 mod 8:", "OK" if good else "FAIL",
          "; least witnesses", sorted(found.values()))
    ok &= good
    # the CRT class: r = c_w mod w (c_w a unit), r = 3 mod 8 -> coprime to 8 prod w
    M = 8 * prod(anchors)
    for pat in found:
        r = found[pat]
        assert gcd(r % M, M) == 1 and r % 8 == 3 and r % 4 == 3
        for u, s in zip(anchors, pat):
            assert jacobi(u, r) == -s          # reciprocity for two primes = 3 mod 4

    # ---- (2) Heegner domination
    heeg = [3, 7, 11, 19, 43, 67, 163]
    allp = primes_upto(2100)
    for q in heeg:
        bound = (q + 1) // 4
        small = [p for p in allp if p < bound]
        inert = all((q % 8 == 3) if p == 2 else jacobi(-q % p, p) == -1 for p in small)
        dom = [p for p in small if p % 4 == 3]
        good = inert and all(jacobi(q, p) == 1 for p in dom)
        print(f"(2) q = {q}: (q+1)/4 = {bound}; primes below it {small}: all inert: {inert}; q -> p for p = 3 mod 4 in "
              f"{dom}:", "OK" if good else "FAIL")
        ok &= good
    viol = 0
    tested = 0
    for q in allp:
        if q % 4 == 3 and q not in heeg and q >= 23:
            dom = [p for p in allp if p % 4 == 3 and 4 * p < q + 1]
            if dom:
                tested += 1
                if not all(jacobi(q, p) == 1 for p in dom):
                    viol += 1
    print(f"    control: among the {tested} non-Heegner primes q = 3 mod 4, 23 <= q < 2100, domination fails for {viol}")

    # ---- (3) Heegner tournament
    out = arcs(heeg)
    k = len(heeg)
    print("(3) arcs:", "; ".join(f"{heeg[i]} -> {[heeg[j] for j in range(k) if out[i] >> j & 1]}" for i in range(k)))
    claimed = {3: [11], 7: [3, 19], 11: [7, 19, 43], 19: [3, 67], 43: [3, 7, 19, 163], 67: [3, 7, 11, 43],
               163: [3, 7, 11, 19, 67]}
    good = all([heeg[j] for j in range(k) if out[i] >> j & 1] == claimed[heeg[i]] for i in range(k))
    scores = [bin(o).count("1") for o in out]
    tri = sum(1 for a, b, c in itertools.combinations(range(k), 3)
              if hp_dfs(out, (1 << a) | (1 << b) | (1 << c)) == 3)
    H = hp_perm(out)
    prefixes = [hp_dfs(out, (1 << j) - 1) for j in range(1, k + 1)]
    six = [i for i in range(k) if heeg[i] % 8 == 3]
    m6 = sum(1 << i for i in six)
    h6 = hp_dfs(out, m6)
    values = sorted({hp_dfs(out, m) for m in range(1, 1 << k)})
    good &= (scores == [1, 2, 3, 2, 4, 4, 5] and tri == 8 and H == 71 and hp_dfs(out) == 71
             and prefixes == [1, 1, 3, 5, 9, 33, 71] and h6 == 29
             and values == [1, 3, 5, 9, 13, 15, 17, 23, 29, 33, 71])
    print(f"    scores {scores}, cyclic triangles {tri}, H = {H} (permutations), prefixes {prefixes}, six primes = 3 mod 8 "
          f"(product {prod(heeg[i] for i in six)}): H = {h6}; sub-tournament values {values}", "OK" if good else "FAIL")
    ok &= good
    # {3,7,11} cyclic, {3,7,11,19} strong with H = 5
    forced = [(19, 3), (43, 3), (43, 7), (67, 3), (67, 7), (67, 11), (163, 3), (163, 7), (163, 11), (163, 19)]
    good = all(jacobi(a, b) == 1 for a, b in forced) and len(forced) == 10
    n_forced = sum(1 for q in heeg for p in heeg if p % 4 == 3 and 4 * p < q + 1)
    print(f"    arcs forced by domination inside the Heegner set: {n_forced} of 21:", "OK" if good and n_forced == 10 else "FAIL")
    ok &= good and n_forced == 10
    six_primes = [heeg[i] for i in six]
    mv = monsky_variants(six_primes)
    c1 = monsky_variants([3, 11, 19])
    c2 = monsky_variants([3, 11, 59])
    good = mv == (4, 4, 4) and c1 == (2, 2, 2) and c2 == (0, 0, 0)
    print(f"    Monsky nullity of {six_primes}: packet convention {mv[0]}, cited-note convention {mv[1]}, full 12 x 12 "
          f"Monsky matrix {mv[2]}; controls 627: {c1}, 1947: {c2}", "OK" if good else "FAIL")
    ok &= good

    # ---- (4) least realisers
    # (a) fully exhaustive below X1
    X1 = 400000
    P3 = [p for p in primes_upto(X1) if p % 4 == 3]
    best = {}
    count = 0

    def dfs(start, chosen, product, outm):
        nonlocal count
        if chosen:
            count += 1
            h = hp_dfs(outm) if len(chosen) >= 3 else 1
            if h not in best or product < best[h]:
                best[h] = product
        for idx in range(start, len(P3)):
            p = P3[idx]
            if product * p > X1:
                break
            new = 0
            newout = list(outm)
            kk = len(chosen)
            for i, q in enumerate(chosen):
                if jacobi(q, p) == 1:
                    newout[i] |= 1 << kk
                else:
                    new |= 1 << i
            newout.append(new)
            dfs(idx + 1, chosen + [p], product * p, newout)

    dfs(0, [], 1, [])
    print(f"(4a) every squarefree n <= {X1} with all prime factors = 3 mod 4 ({count} integers, no box): least n by "
          f"number of residue chains: {dict(sorted(best.items()))}", f"[{time.time()-t0:.0f}s]")
    good = (best.get(1) == 3 and best.get(3) == 231 and best.get(5) == 4389 and best.get(9) == 188727
            and best.get(11) == 100947 and best.get(13) == 206283 and best.get(15) == 311619)
    print("     a(1), a(3), a(5), a(9), a(11), a(13), a(15) = 3, 231, 4389, 188727, 100947, 206283, 311619:",
          "OK" if good else "FAIL")
    ok &= good

    # (b) own spectra of orders <= 6 -> least order k(N)
    least_order = {}
    for n in range(1, 7):
        prs = [(i, j) for j in range(1, n) for i in range(j)]
        for code in range(1 << len(prs)):
            o = [0] * n
            for t, (i, j) in enumerate(prs):
                if code >> t & 1:
                    o[i] |= 1 << j
                else:
                    o[j] |= 1 << i
            least_order.setdefault(hp_dfs(o), n)
    spec6 = sorted(h for h in least_order)
    print(f"(4b) spectrum of orders <= 6: {spec6}; least orders:",
          {h: least_order[h] for h in spec6}, f"[{time.time()-t0:.0f}s]")
    table = {1: 3, 3: 231, 5: 4389, 9: 188727, 11: 100947, 13: 206283, 15: 311619, 17: 48894657, 19: 8378601,
             23: 4340721, 25: 9115953, 27: 26233053, 29: 11035101, 31: 17121489, 33: 7741041, 35: 1901942427,
             37: 4744509, 39: 1209233487, 41: 7167237, 43: 9660189, 45: 3129357}
    XMAX = max(table.values())
    # primes needed: a k-factor n <= XMAX has largest prime <= XMAX / (product of the k-1 smallest primes of P3)
    small = [3, 7, 11, 19, 23, 31, 43, 47]
    allok = True
    for N in sorted(table):
        target = table[N]
        kmin = least_order.get(N, 7)
        pmax = target // prod(small[:kmin - 1]) + 1
        Pl = [p for p in primes_upto(max(pmax, 50)) if p % 4 == 3]
        smaller_hits = []
        seen_target = False
        nodes = 0

        def dfs2(start, chosen, product, outm):
            nonlocal seen_target, nodes
            kk = len(chosen)
            if kk >= kmin:
                nodes += 1
                h = hp_dfs(outm)
                if h == N:
                    if product < target:
                        smaller_hits.append(product)
                    elif product == target:
                        seen_target = True
            for idx in range(start, len(Pl)):
                p = Pl[idx]
                # need at least kmin factors in total: the cheapest completion uses the next primes after p
                need = max(kmin - kk - 1, 0)
                if idx + need >= len(Pl) and need > 0:
                    break
                lower = product * p
                for t in range(1, need + 1):
                    lower *= Pl[idx + t]
                if lower > target:
                    break
                newout = list(outm)
                new = 0
                for i, q in enumerate(chosen):
                    if jacobi(q, p) == 1:
                        newout[i] |= 1 << kk
                    else:
                        new |= 1 << i
                newout.append(new)
                dfs2(idx + 1, chosen + [p], product * p, newout)

        dfs2(0, [], 1, [])
        good = seen_target and not smaller_hits
        allok &= good
        fac = []
        m = target
        for p in Pl:
            if m % p == 0:
                fac.append(p)
                m //= p
        print(f"     N = {N:2d}: claimed a(N) = {target} = {'*'.join(map(str, fac))}{'' if m == 1 else ' * ?'}; least order "
              f"{'>= 7' if N not in least_order else kmin}; candidates n <= a(N) with >= {kmin} prime factors examined "
              f"(primes up to {Pl[-1]}): {nodes}; realised at a(N): {seen_target}; smaller realisers: {smaller_hits}",
              "OK" if good else "FAIL", f"[{time.time()-t0:.0f}s]", flush=True)
    ok &= allok
    print("AUDIT REDEI:", "ALL CHECKS PASSED" if ok else "SOME CHECK FAILED", f"[{time.time()-t0:.0f}s]")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
