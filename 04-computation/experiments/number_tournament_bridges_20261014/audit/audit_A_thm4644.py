#!/usr/bin/env python3
r"""Independent audit of THM-4644 (referee script; does not import the author's code).

A1  explicit orbit decomposition of x -> g x on Z_(p^2) \ {0} for every odd prime p < 4000:
    orbit count N, negation pairing (no orbit closed under x -> -x), E = N/2, vs the stated formula and vs Wieferich.
A2  brute force: count ALL circulant tournaments S on Z_(n) (one element from each pair {x, -x}) with g S = S,
    n = 9, 25, 49 (base 2) and the Wieferich-regime analogue n = 49 with g = 18 (18^3 = 1 mod 49).
A3  orders: ord_3511(2), ord_1093(2), ord_1093(3), (3^7-1)/2.
A4  rigidity: multiplier stabilizers of all circulant tournaments for p <= 31; brute-force automorphism counts
    (own backtracking) for p <= 17.
"""
import sys, time
import numpy as np
import numba as nb


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def mult_order(a, m):
    """multiplicative order of a mod m by repeated multiplication (no library)"""
    a %= m
    x = a
    k = 1
    while x != 1:
        x = x * a % m
        k += 1
        if k > m:
            raise ValueError('not a unit')
    return k


@nb.njit(cache=True)
def orbit_census(n, g):
    """orbits of x -> g x mod n on 1..n-1; returns (number of orbits, number of orbits closed under negation)"""
    seen = np.zeros(n, dtype=np.int32)      # orbit id + 1
    norb = 0
    for x in range(1, n):
        if seen[x] != 0:
            continue
        norb += 1
        y = x
        while seen[y] == 0:
            seen[y] = norb
            y = (y * g) % n
    closed = 0
    # an orbit O is closed under negation iff -x in O for its representative
    done = np.zeros(norb + 1, dtype=np.uint8)
    for x in range(1, n):
        o = seen[x]
        if done[o] == 0:
            done[o] = 1
            if seen[n - x] == o:
                closed += 1
    return norb, closed


def part_A1(out):
    out.append('== A1. exponent E for every odd prime p < 4000 (explicit orbit decomposition on Z_(p^2) \\ {0})')
    jumps = []
    bad = []
    rows = {}
    for p in range(3, 4000):
        if not is_prime(p):
            continue
        d = mult_order(2, p)
        dp = d
        while dp % 2 == 0:
            dp //= 2
        g = pow(2, d // dp, p * p)
        og = mult_order(g, p * p)
        N, closed = orbit_census(p * p, g)
        if closed != 0 or N % 2 != 0:
            bad.append((p, N, closed))
            continue
        E = N // 2
        wief = pow(2, p - 1, p * p) == 1
        formula = (p + 1) * (p - 1) // (2 * dp) if wief else (p - 1) // dp
        generic = (p - 1) // dp
        if E != formula:
            bad.append((p, 'formula', E, formula))
        if E != generic:
            jumps.append((p, E, generic, wief))
        if (og % 2 == 0) or (og != dp * (1 if wief else p)):
            bad.append((p, 'order', og, dp))
        if p in (3, 5, 7, 31, 127, 1093, 3511):
            rows[p] = (d, dp, og, E, generic, wief)
    for p, (d, dp, og, E, generic, wief) in sorted(rows.items()):
        out.append(f'   p = {p:5d}: ord_p(2) = {d}, d\' = {dp}, ord_(p^2)(g) = {og}, E = {E} (non-Wieferich value {generic}), Wieferich {wief}')
    out.append(f'   primes p < 4000 where E differs from (p-1)/d\': {jumps}')
    out.append(f'   anomalies (orbit closed under negation / odd N / formula mismatch / even order): {bad}')
    return jumps, bad


def enumerate_invariant(n, g):
    """brute force over all 2^((n-1)/2) connection sets; count those with g S = S (bitmask representation)"""
    pairs = [(x, n - x) for x in range(1, (n - 1) // 2 + 1)]
    m = len(pairs)
    # precompute the image of each element under multiplication by g
    img = [(x * g) % n for x in range(n)]
    return _enum_inv(np.array([a for a, b in pairs], dtype=np.int64), np.array([b for a, b in pairs], dtype=np.int64),
                     np.array(img, dtype=np.int64), n)


@nb.njit(cache=True)
def _enum_inv(pa, pb, img, n):
    m = pa.shape[0]
    cnt = 0
    for bits in range(1 << m):
        # S as membership array via bitmask over residues (n <= 64 here)
        smask = np.int64(0)
        for i in range(m):
            if (bits >> i) & 1:
                smask |= np.int64(1) << pb[i]
            else:
                smask |= np.int64(1) << pa[i]
        ok = True
        for x in range(1, n):
            if (smask >> x) & 1:
                if not ((smask >> img[x]) & 1):
                    ok = False
                    break
        if ok:
            cnt += 1
    return cnt


def part_A2(out):
    out.append('\n== A2. brute force: g-invariant circulant tournaments on Z_n (all connection sets enumerated)')
    res = {}
    for (p, g, label) in ((3, 4, 'base 2: g = 2^2 = 4'), (5, 16, 'base 2: g = 2^4 = 16'), (7, 2, 'base 2: g = 2'),
                          (7, 18, 'Wieferich-regime analogue: g = 18, 18^3 = 1 mod 49'),
                          (5, 7 ** 4 % 25, 'base 7 (7^4 = 1 mod 25): g = 7^4 mod 25 = 1')):
        n = p * p
        t0 = time.time()
        cnt = enumerate_invariant(n, g)
        og = mult_order(g, n) if g != 1 else 1
        dprime = mult_order(g % p, p) if g % p != 1 else 1
        wief_like = (og == dprime)
        N = (p * (p - 1)) // og + (p - 1) // dprime
        pred = 2 ** (N // 2)
        out.append(f'   Z_{n:2d}, {label}: ord_(p^2)(g) = {og}, ord_p(g) = {dprime}; invariant tournaments (brute force) = {cnt},'
                   f' orbit formula 2^(N/2) = {pred};  match {cnt == pred};  {time.time() - t0:.1f} s')
        res[(p, g)] = (cnt, pred)
    return res


@nb.njit(cache=True)
def count_aut(adj, n):
    """number of automorphisms of a tournament given by boolean adjacency (adj[i, j] = i -> j), by backtracking"""
    perm = -np.ones(n, dtype=np.int64)
    used = np.zeros(n, dtype=np.uint8)
    # out-degrees must be preserved
    outdeg = np.zeros(n, dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if adj[i, j]:
                outdeg[i] += 1
    count = 0
    # iterative DFS over positions 0..n-1
    cand = np.zeros(n, dtype=np.int64)    # next candidate to try at each depth
    depth = 0
    cand[0] = 0
    while depth >= 0:
        if depth == n:
            count += 1
            depth -= 1
            used[perm[depth]] = 0
            perm[depth] = -1
            continue
        found = False
        c = cand[depth]
        while c < n:
            if used[c] == 0 and outdeg[c] == outdeg[depth]:
                ok = True
                for k in range(depth):
                    if adj[k, depth] != adj[perm[k], c]:
                        ok = False
                        break
                if ok:
                    found = True
                    break
            c += 1
        if found:
            perm[depth] = c
            used[c] = 1
            cand[depth] = c + 1
            depth += 1
            if depth < n:
                cand[depth] = 0
        else:
            cand[depth] = 0
            depth -= 1
            if depth >= 0:
                used[perm[depth]] = 0
                perm[depth] = -1
    return count


def circulant_sets(p):
    pairs = [(x, p - x) for x in range(1, (p - 1) // 2 + 1)]
    for bits in range(1 << len(pairs)):
        yield frozenset(pairs[i][(bits >> i) & 1] for i in range(len(pairs)))


def part_A4(out):
    out.append('\n== A4. rigidity of circulant tournaments of prime order')
    for p in (3, 5, 7, 11, 13, 17):
        sizes = {}
        for S in circulant_sets(p):
            adj = np.zeros((p, p), dtype=np.bool_)
            for i in range(p):
                for j in range(p):
                    if i != j and (j - i) % p in S:
                        adj[i, j] = True
            a = count_aut(adj, p)
            sizes[a] = sizes.get(a, 0) + 1
        out.append(f'   p = {p:2d}: brute-force |Aut| over all {2 ** ((p - 1) // 2)} circulant tournaments: {dict(sorted(sizes.items()))}')
    rigid = []
    for p in range(3, 32):
        if not is_prime(p):
            continue
        allrigid = True
        for S in circulant_sets(p):
            M = [a for a in range(2, p) if frozenset(a * x % p for x in S) == S]
            if M:
                allrigid = False
                break
        if allrigid:
            rigid.append(p)
    out.append(f'   primes p <= 31 at which every circulant tournament has trivial multiplier stabilizer: {rigid}')
    # free action count at Fermat primes
    out.append('   rotational classes at Fermat primes 2^((p-1)/2)/(p-1): ' +
               str({p: (2 ** ((p - 1) // 2)) // (p - 1) if p <= 17 else f'2^{(p - 1) // 2 - (p - 1).bit_length() + 1}' for p in (3, 5, 17, 257, 65537)}))


def part_A3(out):
    out.append('\n== A3. orders')
    out.append(f'   ord_3511(2) = {mult_order(2, 3511)}; (3511-1)/2 = {3510 // 2}; 3511 mod 4 = {3511 % 4}, mod 8 = {3511 % 8}')
    out.append(f'   ord_1093(2) = {mult_order(2, 1093)}; ord_1093(3) = {mult_order(3, 1093)}; (3^7-1)/2 = {(3 ** 7 - 1) // 2}')
    out.append(f'   2^1092 mod 1093^2 = {pow(2, 1092, 1093 ** 2)}; 2^3510 mod 3511^2 = {pow(2, 3510, 3511 ** 2)}')
    # <2> orbit counts on units (THM-4631's c_1, c_2) vs <g> orbit counts (THM-4644's unit orbits)
    for p in (1093, 3511):
        d = mult_order(2, p)
        dp = d
        while dp % 2 == 0:
            dp //= 2
        L2 = mult_order(2, p * p)
        g = pow(2, d // dp, p * p)
        og = mult_order(g, p * p)
        out.append(f'   p = {p}: <2>-orbits on units: level 1 {(p - 1) // d}, level 2 {p * (p - 1) // L2};'
                   f'  <g>-orbits on units: level 1 {(p - 1) // dp}, level 2 {p * (p - 1) // og}')


if __name__ == '__main__':
    out = []
    t0 = time.time()
    part_A1(out)
    part_A2(out)
    part_A3(out)
    part_A4(out)
    out.append(f'\n(total {time.time() - t0:.0f} s)')
    print('\n'.join(out))
