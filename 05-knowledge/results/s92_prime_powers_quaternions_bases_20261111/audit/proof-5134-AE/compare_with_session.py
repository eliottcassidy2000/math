# compare_with_session.py -- auditor proof-5134-AE.
# Compares every number of the session outputs cum_22/24/26/28.out (remainder table, noise, class factors:
# measured mean [min, max] | first-order | limit) with the values recomputed from my own census (census_L.json).
import sys, json, math, re
from fractions import Fraction as Fr

SESS = "C:/tmp/wt-s92/04-computation/experiments/s92_prime_powers_20261111/"
HERE = "C:/tmp/wt-s92/05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/proof-5134-AE/"


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


def mine(L):
    res = json.load(open(HERE + "census_%d.json" % L))
    out = {}
    for ms, dat in res["data"].items():
        m = int(ms)
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
                fac[g0] = (sum(v) / len(v), min(v), max(v), float(sum(d * Fr(D[d], T) for d in divs if g0 % d == 0)),
                           float(sum(d * c_d(m, d) for d in divs if g0 % d == 0)))
            out[(m, sh)] = (T, float(m * rem / T), math.sqrt(m / T), fac)
    return out


def main():
    allok = True
    for L in (22, 24, 26, 28):
        my = mine(L)
        txt = open(SESS + "cum_%d.out" % L).read().splitlines()
        nrem = nfac = 0
        maxdiff = 0.0
        hdr = re.search(r"sheet 0: (\d+), sheet 1: (\d+)", txt[0])
        T0, T1 = int(hdr.group(1)), int(hdr.group(2))
        ok = (T0 == my[(9, 0)][0] and T1 == my[(9, 1)][0])
        for ln in txt:
            mm = re.match(r"^\s+(\d+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*$", ln)
            if mm:
                m = int(mm.group(1))
                for sh in range(2):
                    d = abs(float(mm.group(2 + sh)) - my[(m, sh)][1])
                    maxdiff = max(maxdiff, d)
                    if d > 6e-6:
                        ok = False
                        print("  MISMATCH remainder m=%d sheet %d: session %s mine %.5f" % (m, sh, mm.group(2 + sh), my[(m, sh)][1]))
                    nrem += 1
                # noise: the session prints one noise column; compare with both sheets
                nz = float(mm.group(4))
                if min(abs(nz - my[(m, 0)][2]), abs(nz - my[(m, 1)][2])) > 6e-6:
                    print("  note: noise column m=%d: session %.5f, mine sheet0 %.5f sheet1 %.5f" % (m, nz, my[(m, 0)][2], my[(m, 1)][2]))
            mm = re.match(r"^\s+m =\s+(\d+) sheet (\d): (.*)$", ln)
            if mm:
                m, sh = int(mm.group(1)), int(mm.group(2))
                for part in mm.group(3).split(";"):
                    pm = re.match(r"\s*g0=(\d+): ([0-9.]+) \[([0-9.]+), ([0-9.]+)\] \| ([0-9.]+) \| ([0-9.]+)", part)
                    g0 = int(pm.group(1))
                    vals = [float(pm.group(i)) for i in range(2, 7)]
                    minev = my[(m, sh)][3][g0]
                    for a, b in zip(vals, minev):
                        maxdiff = max(maxdiff, abs(a - b))
                        if abs(a - b) > 6e-5:
                            ok = False
                            print("  MISMATCH factor m=%d sheet %d g0=%d: session %s mine %s" % (m, sh, g0, vals, ["%.4f" % t for t in minev]))
                    nfac += 5
        print("X = 2^%d: triple counts equal: %s; %d remainders and %d class-factor entries compared, largest difference %.1e  -> %s"
              % (L, (T0, T1) == (my[(9, 0)][0], my[(9, 1)][0]), nrem, nfac, maxdiff, "ALL AGREE" if ok else "DISAGREEMENT"))
        allok = allok and ok
    print("OVERALL:", "session outputs cum_22/24/26/28.out reproduced" if allok else "DISAGREEMENTS FOUND")


if __name__ == "__main__":
    main()
