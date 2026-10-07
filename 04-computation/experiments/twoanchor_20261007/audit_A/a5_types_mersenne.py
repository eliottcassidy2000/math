#!/usr/bin/env python3
"""Audit A, item 5: THM-4601 (v) run-type reduction and (vi) Mersenne shells, independent code.

(v)  end-of-run D-chain states (Terras time 2J after x) for random sources (own generator), D = 1..8:
     do they depend only on (type, D)?  (types: j = 3, 4, 5, 6, and j >= 7 split by r)
(vi) (a) q -> w'(q) = 3^(J-3) (3^(2^m q) - 1)/2^(m+2) mod 2^L is well defined on q mod 2^L and a bijection of odd residues
         (m = 4..10, L = 1..14);  (b) for actual Mersenne exponents K = 1 + 2^m q the integer Y_end - 1 equals 2^r w'(q);
     (c) the two examples 2^1889 - 1 ~> 2^1886 - 1 and 2^5249 - 1 ~> 2^5246 - 1 by direct U-orbits, with the stated words;
     (d) shell fractions for K < 12000 at depth 200 recomputed (own chain code).
"""
import random, sys
from collections import defaultdict, Counter
from a2_barrier import make_source, T

NCHK = 0
def check(c, m):
    global NCHK
    if not c: raise AssertionError(m)
    NCHK += 1

def v2(n):
    return (n & -n).bit_length() - 1

def uword(z, L):
    out = []
    for _ in range(L):
        y = 3*z + 1; a = v2(y); out.append(a); z = y >> a
    return out, z

def end_state(x, D, J):
    y = (x + 1)//3**D - 1
    u, v, k = x, y, D
    for _ in range(2*J):
        k += (u & 1) - (v & 1); u, v = T(u), T(v)
    return (k, u - 3**k*v) if k >= 0 else (k, None)

def part_v():
    rnd = random.Random(555)
    table = defaultdict(set)
    for j in (3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 25, 26):
        J = (j - 1)//2
        for _ in range(60):
            K = rnd.randint(9, 60)
            n, t, x = make_source(rnd, K, j, 100)
            for D in range(1, 9):
                table[(j, D)].add(end_state(x, D, J))
    for j in (3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 25, 26):
        row = []
        for D in range(1, 9):
            st = table[(j, D)]
            check(len(st) == 1, f"state determined by (j, D) for j={j}, D={D}")
            row.append(next(iter(st)))
        print(f"   j={j:2d} (J={(j-1)//2}, r={j - 2*((j-1)//2)}): " + "  ".join(f"D{D}:{s}" for D, s in zip(range(1, 9), row)))
    # within the type 'j >= 7, r', does the state depend on J?
    for r in (1, 2):
        js = [j for j in (7, 8, 9, 10, 11, 12, 13, 14, 25, 26) if j - 2*((j-1)//2) == r]
        for D in range(1, 9):
            sts = {next(iter(table[(j, D)])) for j in js}
            print(f"   type j>=7, r={r}, D={D}: {len(sts)} distinct end states over j in {js}")
    # r-independence: same J, different r -> same state
    for J in (1, 2, 3, 4, 5, 6, 12):
        for D in range(1, 9):
            check(table[(2*J + 1, D)] == table[(2*J + 2, D)], "end state depends on J only, not on r")
    print("   end-of-run state depends on (J, D) only (identical for r = 1, 2 at equal J)")

def wprime(m, q, L):
    J = 1 + m//2
    mod = 1 << (L + m + 2)
    a = (pow(3, (1 << m)*q, mod) - 1) % mod
    check(a % (1 << (m + 2)) == 0 and (a >> (m + 2)) & 1 == 1, "v2(3^(2^m q) - 1) = m+2")
    u = (a >> (m + 2)) % (1 << L)
    return (pow(3, J - 3, 1 << L) * u) % (1 << L) if J >= 3 else None

def part_vi_ab():
    for m in range(4, 11):
        for L in range(1, 15):
            imgs = {}
            for q in range(1, 1 << (L + 1), 2):            # two periods
                w = wprime(m, q, L)
                check(w & 1, "odd image")
                key = q % (1 << L)
                if key in imgs:
                    check(imgs[key] == w, "depends only on q mod 2^L")
                imgs[key] = w
            check(sorted(imgs.values()) == list(range(1, 1 << L, 2)), f"bijection on odd residues m={m} L={L}")
    print("   (vi a) q -> w'(q) mod 2^L well defined and bijective on odd residues for m = 4..10, L = 1..14")
    # (b) actual integers
    n = 0
    for K in range(17, 2600, 2):
        m = v2(K - 1)
        if m < 4: continue
        J = 1 + m//2; r = 1 + (m % 2)
        x = 2*3**(K - 1) - 1
        check(v2(x - 1) == 3 + m, "j = 3 + m")
        y3 = (x + 1)//27 - 1
        u, v, k = x, y3, 3
        for _ in range(2*J):
            k += (u & 1) - (v & 1); u, v = T(u), T(v)
        check(k == 3 and u - 1 == 27*(v - 1), "universal state on the Mersenne line")
        wtrue = 3**(J - 3) * ((3**(K - 1) - 1) >> (m + 2))
        check(v - 1 == (wtrue << r), "Y_end - 1 = 2^r w'")
        n += 1
    print(f"   (vi b) {n} Mersenne exponents K < 2600 with m >= 4: universal state and Y_end = 1 + 2^r w'(q) exact")

def part_vi_c():
    for K, m, J, c in ((1889, 5, 3, 3), (5249, 7, 4, 1)):
        check(v2(K - 1) == m and 1 + m//2 == J and 1 + (m % 2) == 2, "shell data")
        n = (1 << K) - 1; h = (1 << (K - 3)) - 1
        src_exp = [1]*(K - 1) + [2]*J + [10] + [c]
        chd_exp = [1]*(K - 4) + [4, 1, 1] + [2]*(J - 3) + [3, 1, 1, 3] + [c + 2]
        ws, zs = uword(n, len(src_exp))
        wc, zc = uword(h, len(chd_exp))
        check(ws == src_exp, f"source word of 2^{K}-1")
        check(wc == chd_exp, f"child word of 2^{K-3}-1")
        check(zs == zc, "same odd endpoint")
        # Terras-time bookkeeping: equal number of odd steps, Terras times differ by D = 3
        check(len(ws) == len(wc) and sum(ws) == sum(wc) + 3, "equal odd steps; Terras times differ by 3")
        print(f"   (vi c) 2^{K}-1 ~> 2^{K-3}-1: words 1^{K-1} 2^{J} (10) {c}  and  1^{K-4} (4,1,1) 2^{J-3} (3,1,1,3) {c+2}; "
              f"common odd endpoint has {zs.bit_length()} bits; odd steps {len(ws)}")

def part_vi_d(KMAX=12000, S=200):
    shell = defaultdict(lambda: [0, 0])
    for K in range(17, KMAX, 2):
        m = v2(K - 1)
        if m < 4: continue
        J = 1 + m//2
        x = 2*3**(K - 1) - 1; y = (x + 1)//27 - 1
        u, v, k = x, y, 3
        for _ in range(2*J):
            k += (u & 1) - (v & 1); u, v = T(u), T(v)
        got = False
        for s in range(S + 1):
            if u == v and k == 0:
                got = True; break
            k += (u & 1) - (v & 1); u, v = T(u), T(v)
        shell[m][0] += 1; shell[m][1] += got
    print("   (vi d) shells K < %d, depth %d: " % (KMAX, S) +
          ", ".join(f"m={m}: {a}/{c}={a/c:.3f}" for m, (c, a) in sorted(shell.items())))
    p1 = [sum(shell[m][i] for m in shell if m % 2 == 0) for i in (0, 1)]
    p2 = [sum(shell[m][i] for m in shell if m % 2 == 1) for i in (0, 1)]
    print(f"          pooled r=1 (m even) {p1[1]}/{p1[0]} = {p1[1]/p1[0]:.3f}; r=2 (m odd) {p2[1]}/{p2[0]} = {p2[1]/p2[0]:.3f}")

if __name__ == "__main__":
    print("(v) end-of-run states")
    part_v()
    print("(vi) Mersenne shells")
    part_vi_ab()
    part_vi_c()
    if len(sys.argv) > 1 and sys.argv[1] == 'shells':
        part_vi_d()
    print(f"ALL (v)/(vi) CHECKS PASSED ({NCHK} assertions)")
