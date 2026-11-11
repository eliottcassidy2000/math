# summary_tables.py -- auditor proof-5134-AE.  Tables from my own censuses census_22/24/26/27/28/30.json.
#  T1: exact remainder of THM-5134 (B) for the twelve moduli at X = 2^22 .. 2^30 (2^30 is beyond the session range)
#  T2: P4-type counts: in how many (modulus, sheet) pairs does the remainder decrease from one size to the next
#  T3: for every group g0 > 1: is the measured class factor closer to the first-order value Phi(g0) (measured types)
#      or to the limit L_m(g0)?   (the session registered this comparison, P3, for m = 9, 27, 25 only)
import json, math
from fractions import Fraction as Fr

HERE = "C:/tmp/wt-s92/05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/proof-5134-AE/"
MODS = [9, 27, 81, 25, 125, 49, 45, 75, 63, 225, 15, 21]


def primes_of(n):
    ps = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            ps.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        ps.append(n)
    return ps


def c_d(m, d):
    v = Fr(1, d * d)
    for p in primes_of(m // d):
        v *= (1 - Fr(1, p * p))
    return v


def stats(L):
    res = json.load(open(HERE + "census_%d.json" % L))
    out = {}
    for m in MODS:
        dat = res["data"][str(m)]
        divs = dat["divisors"]
        for sh in range(2):
            off = dat["off"][sh]
            dia = dat["diag"][sh]
            CM = [2 * sum(off[t][r] for t in range(len(divs))) + sum(dia[t][r] for t in range(len(divs))) for r in range(m)]
            T = sum(CM)
            D = {d: 2 * sum(off[t]) for t, d in enumerate(divs)}
            rem = Fr(0)
            grp = {}
            for r in range(m):
                g0 = math.gcd(m, r) if r else m
                pred = sum(Fr(d, m) * D[d] for d in divs if g0 % d == 0)
                rem = max(rem, abs(CM[r] - pred))
                grp.setdefault(g0, []).append(m * CM[r] / T)
            fac = {}
            for g0, v in grp.items():
                fac[g0] = (sum(v) / len(v), float(sum(d * Fr(D[d], T) for d in divs if g0 % d == 0)), float(sum(d * c_d(m, d) for d in divs if g0 % d == 0)))
            out[(m, sh)] = (float(m * rem / T), fac, T)
    return out


def main():
    sizes = [22, 24, 26, 28, 30]
    st = {L: stats(L) for L in sizes}
    print("T1. exact remainder of (B), sheet 0 / sheet 1   (own census; 2^30 is new)")
    print("   m   " + "   ".join("      2^%d      " % L for L in sizes))
    for m in MODS:
        print(" %4d  " % m + "   ".join("%.5f / %.5f" % (st[L][(m, 0)][0], st[L][(m, 1)][0]) for L in sizes))
    print("   #T  " + "   ".join("%d / %d" % (st[L][(9, 0)][2], st[L][(9, 1)][2]) for L in sizes))
    print()
    print("T2. number of (modulus, sheet) pairs in which the remainder decreases")
    for a, b in zip(sizes[:-1], sizes[1:]):
        six = [9, 27, 81, 25, 125, 49]
        n6 = sum(1 for m in six for sh in range(2) if st[b][(m, sh)][0] < st[a][(m, sh)][0])
        n12 = sum(1 for m in MODS for sh in range(2) if st[b][(m, sh)][0] < st[a][(m, sh)][0])
        grew = ["%d/s%d" % (m, sh) for m in MODS for sh in range(2) if not st[b][(m, sh)][0] < st[a][(m, sh)][0]]
        print("   2^%d -> 2^%d: %2d of 12 (registered six prime powers 9, 27, 81, 25, 125, 49);  %2d of 24 (all twelve moduli);  grew at: %s" % (a, b, n6, n12, ", ".join(grew)))
    print()
    print("T3. groups g0 > 1: measured factor closer to the first-order value Phi (with the measured types) than to the limit L_m ?")
    for L in sizes:
        row = []
        for m in MODS:
            c = t = 0
            for sh in range(2):
                for g0, (mean, phi, lim) in st[L][(m, sh)][1].items():
                    if g0 > 1:
                        t += 1
                        c += abs(mean - phi) < abs(mean - lim)
            row.append("%d: %d/%d" % (m, c, t))
        print("   X = 2^%d:  " % L + "   ".join(row))
    print("   (at the multiples of 7 the limit is the closer one in most groups: the remainder of (B) there is still larger than the first-order correction)")
    print()
    for L in (28, 30):
        print("   detail X = 2^%d, m = 49: " % L + ";  ".join("sheet %d g0=%d: measured %.4f, Phi %.4f, limit %.4f" % (sh, g0, v[0], v[1], v[2]) for sh in range(2) for g0, v in sorted(st[L][(49, sh)][1].items()) if g0 > 1))


if __name__ == "__main__":
    main()
