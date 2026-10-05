#!/usr/bin/env python3
"""
The stitched sheets at the centre of the multiplication table, the A/B families, the atom operators, and the
5-tournament tiles (opus, 2026-10-05).  Exact checks of the owner's observations.

1. Centre of the odd multiplication table mod m = 2k+1: the block [[k^2, k(k+1)], [k(k+1), (k+1)^2]] is
   [[d, -d], [-d, d]] with d = 4^-1 mod m, and 3d - 1 = -d, 3(-d) + 1 = d (mod m): the two sheets swap the two
   central values; this is the rational 2-cycle {1/4, -1/4} of x -> 3x-1 / x -> 3x+1, i.e. the sign conjugacy
   U_-(x) = -U_+(-x).  Even modulus 2k with k odd: k^2 = k, (k+-1)^2 = k+1, (k-1)(k+1) = k-1 (no stitching).
2. Fundamental domains: the triangle 0 <= x <= y <= k is a fundamental domain of {0, +-1, ..., +-k}^2 under the
   group of order 8 (swap, two sign flips), with T_(k+1) tiles (15 for k = 4: mod 9 fully, mod 10 without its
   centre cross); the values transform by xy -> +-xy.
3. Triples F_n = (n, ceil(n/2), 2n+1): B_k = (2k-1, k, 4k-1), A_k = (2k, k, 4k+1).  Exact: 4 y = -+1 mod z
   (y = 4^-1 mod z on B, y = -4^-1 on A), 3y - 1 = -y mod z on B and 3y + 1 = -y mod z on A; x(B_k) + x(A_k) = z(B_k);
   the operator p(I,J) = (x(I) + x(J), z(I) x(J)) sends (A_k, B_k), (B_k, A_k) to the x-pair of the atom
   k(4k-1) = y z; the operator q(I,J) = (2(I+J)-1, 4J) on consecutive integers is (4I+1, 4I+4), whose odd member is
   the Collatz sibling of I (U(4I+1) = U(I)).  The mod-4 classes of z are the sheet-pair cases of the sibling
   automaton: v_2(3z+1) = 1 exactly on B, v_2(3z-1) = 1 exactly on A.
4. Primes among z(B_k) = 4k-1 and z(A_k) = 4k+1 to 10^6 (Chebyshev bias), and the independence of the 2-adic
   valuation profile from primality (geometric law in both).
5. Tournaments on 5 vertices: 1024 labelled, iso classes, strong classes; the fixed-Hamiltonian-path tiling has
   6 off-path arcs in diagonals 3, 2, 1 (64 tilings), the selfie version 10 arcs + 5 loops = 15 tiles.
"""
import sys, math, itertools, time
from collections import Counter

def main():
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    t0 = time.time()
    # ---- 1
    ok = 0; tot = 0
    for m in range(3, 402, 2):
        k = (m - 1) // 2; d = pow(4, -1, m)
        tot += 1
        if (k * k - d) % m == 0 and (k * (k + 1) + d) % m == 0 and (3 * d - 1 + d) % m == 0 and (3 * (-d) + 1 - d) % m == 0:
            ok += 1
    P(f"1. odd m = 2k+1 <= 401: centre block = [[d,-d],[-d,d]] with d = 4^-1, and 3d-1 = -d, 3(-d)+1 = d (mod m): {ok}/{tot}")
    for m in (9, 11, 13, 15, 21):
        k = (m - 1) // 2; d = pow(4, -1, m)
        P(f"   m={m}: k^2={k*k % m}, k(k+1)={k*(k+1) % m}, (k+1)^2={(k+1)**2 % m}; 4^-1={d}, -4^-1={(-d) % m}; 3*{d}-1={(3*d-1) % m}, 3*{(-d)%m}+1={(3*(-d)+1) % m}")
    oke = 0; tote = 0
    for m in range(6, 402, 4):          # m = 2k, k odd
        k = m // 2; tote += 1
        if (k * k - k) % m == 0 and ((k - 1) ** 2 - (k + 1)) % m == 0 and ((k + 1) ** 2 - (k + 1)) % m == 0 and ((k - 1) * (k + 1) - (k - 1)) % m == 0:
            oke += 1
    P(f"   even m = 2k, k odd, m <= 398: k^2 = k, (k+-1)^2 = k+1, (k-1)(k+1) = k-1 (mod m): {oke}/{tote}; mod 10: 5,6,6,4 (the central 5 with 4 and 6 diagonalized); 3*4+1={13%10}, 3*6-1={17%10}: no stitching")
    # ---- 2
    def orbits(cells, m):
        seen = set(); reps = []
        for (x, y) in cells:
            if (x, y) in seen: continue
            orb = set()
            for sx in (1, -1):
                for sy in (1, -1):
                    for sw in (False, True):
                        a, b = (sx * x) % m, (sy * y) % m
                        if sw: a, b = b, a
                        orb.add((a, b))
            seen |= orb; reps.append((x, y))
        return reps
    for m in (9, 10, 11):
        cells = [(x, y) for x in range(m) for y in range(m)]
        reps = orbits(cells, m)
        k = (m - 1) // 2
        block = [(x, y) for x in range(m) for y in range(m) if (m % 2 == 1) or (x != m // 2 and y != m // 2)]
        reps_block = orbits(block, m)
        P(f"2. mod {m}: orbits of the full table under swap x sign flips: {len(reps)}; block without the centre cross: {len(reps_block)} = T_{k+1 if m%2 else m//2} = {(k+1)*(k+2)//2 if m%2 else (m//2)*(m//2+1)//2}; representatives 0<=x<=y<=k: {sum(1 for (x,y) in reps_block if x <= y <= (m-1)//2)}")
    # ---- 3
    ok3 = 0; tot3 = 0
    for n in range(1, 2001):
        x, y, z = n, (n + 1) // 2, 2 * n + 1
        tot3 += 1
        if n % 2 == 1:   # B
            if (4 * y - 1) % z == 0 and (3 * y - 1 + y) % z == 0: ok3 += 1
        else:            # A
            if (4 * y + 1) % z == 0 and (3 * y + 1 + y) % z == 0: ok3 += 1
    P(f"3. triples (n, ceil(n/2), 2n+1), n <= 2000: 4y = +1 mod z and 3y-1 = -y on B (n odd), 4y = -1 and 3y+1 = -y on A (n even): {ok3}/{tot3}")
    okp = 0
    for k in range(1, 1001):
        B = (2 * k - 1, k, 4 * k - 1); A = (2 * k, k, 4 * k + 1)
        pBA = (B[0] + A[0], B[2] * A[0]); pAB = (A[0] + B[0], A[2] * B[0])
        N = k * (4 * k - 1)
        if pBA[0] == B[2] and pAB[0] == B[2] and {pBA[1], pAB[1]} == {2 * N - 1, 2 * N} and N == B[1] * B[2]:
            okp += 1
    P(f"   p(B_k,A_k) = (z(B_k), 2N), p(A_k,B_k) = (z(B_k), 2N-1) with N = k(4k-1) = y z: {okp}/1000 (the x-pair of the atom y*z)")
    def U(n):
        n = 3 * n + 1
        while n % 2 == 0: n //= 2
        return n
    okq = 0
    for I in range(1, 10001):
        O, E = 2 * (I + I + 1) - 1, 4 * (I + 1)
        if O == 4 * I + 1 and E == 4 * I + 4 and (I % 2 == 0 or U(O) == U(I)): okq += 1
    P(f"   q(I, I+1) = (4I+1, 4I+4); for odd I the odd output is the sibling, U(4I+1) = U(I): {okq}/10000")
    okv = 0
    for k in range(1, 10001):
        zB, zA = 4 * k - 1, 4 * k + 1
        if (3 * zB + 1) % 4 == 2 and (3 * zB - 1) % 4 == 0 and (3 * zA - 1) % 4 == 2 and (3 * zA + 1) % 4 == 0: okv += 1
    P(f"   sheet split: v2(3z+1) = 1 and v2(3z-1) >= 2 on B; v2(3z-1) = 1 and v2(3z+1) >= 2 on A: {okv}/10000 (the P/M cases of the sibling automaton)")
    # ---- 4
    N4 = 1_000_000
    sieve = bytearray([1]) * (N4 + 1); sieve[0] = sieve[1] = 0
    for i in range(2, int(N4 ** 0.5) + 1):
        if sieve[i]: sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    pB = sum(1 for z in range(3, N4 + 1, 4) if sieve[z]); pA = sum(1 for z in range(5, N4 + 1, 4) if sieve[z])
    P(f"4. primes <= 10^6: in B (4k-1) {pB}, in A (4k+1) {pA} (Chebyshev bias toward B: {pB - pA})")
    def vprof(zs, sheet):
        c = Counter()
        for z in zs:
            v = 3 * z + sheet
            a = 0
            while v % 2 == 0: v //= 2; a += 1
            c[min(a, 6)] += 1
        n = sum(c.values()); return [round(c[a] / n, 4) for a in range(1, 7)]
    primesA = [z for z in range(5, N4 + 1, 4) if sieve[z]]; allA = list(range(5, N4 + 1, 4))
    primesB = [z for z in range(3, N4 + 1, 4) if sieve[z]]; allB = list(range(3, N4 + 1, 4))
    P(f"   valuation profile of 3z+1 (a=1..5, >=6) on A: primes {vprof(primesA, 1)} vs all {vprof(allA, 1)}; on B of 3z-1: primes {vprof(primesB, -1)} vs all {vprof(allB, -1)} (geometric(1/2) beyond the forced first bit: primality does not see the sheets)")
    # ---- 5
    n = 5; pairs = list(itertools.combinations(range(n), 2))
    perms = list(itertools.permutations(range(n)))
    def canon(T):
        best = None
        for p in perms:
            key = tuple(sorted(((p[i], p[j]) if T[(i, j)] else (p[j], p[i])) for (i, j) in pairs))
            if best is None or key < best: best = key
        return best
    classes = {}; strong = 0
    for bits in range(1 << len(pairs)):
        T = {pairs[t]: (bits >> t) & 1 for t in range(len(pairs))}
        key = canon(T)
        if key not in classes:
            # strong?  reachability
            adj = {i: set() for i in range(n)}
            for (i, j) in pairs:
                if T[(i, j)]: adj[i].add(j)
                else: adj[j].add(i)
            def reach(s):
                seen = {s}; st = [s]
                while st:
                    u = st.pop()
                    for w in adj[u]:
                        if w not in seen: seen.add(w); st.append(w)
                return seen
            st_ = all(len(reach(s)) == n for s in range(n))
            classes[key] = st_
    offpath = [(i, j) for (i, j) in pairs if j - i >= 2]
    diag = Counter(j - i for (i, j) in offpath)
    P(f"5. tournaments on 5 vertices: {1 << len(pairs)} labelled, {len(classes)} isomorphism classes, {sum(classes.values())} strong; fixed Hamiltonian path 1->2->3->4->5 leaves {len(offpath)} off-path arcs in diagonals (by length 2,3,4): {[diag[2], diag[3], diag[4]]}, so 2^{len(offpath)} = {1 << len(offpath)} tilings; a selfie 5-tournament has 10 arcs + 5 loops = 15 tiles = T_5, and the frame around the 6 off-path tiles is the 4 path arcs plus the 5 loops = 9")
    P(f"  [{time.time()-t0:.0f}s]")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
