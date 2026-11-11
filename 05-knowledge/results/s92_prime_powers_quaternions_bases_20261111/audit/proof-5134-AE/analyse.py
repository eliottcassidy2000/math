# analyse.py -- auditor proof-5134-AE.  Reads the exact counts of census_f1.py and evaluates THM-5134 (B), (C), (D).
# usage: python3 analyse.py CENSUS.json [BRUTE.json]
import sys, json, math
from fractions import Fraction as Fr


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


def L_m(m, e):
    g = math.gcd(m, e % m) if e % m else m
    return sum(d * c_d(m, d) for d in range(1, m + 1) if g % d == 0)


def main():
    res = json.load(open(sys.argv[1]))
    L = res["log2X"]
    X = 1 << L
    print("independent census (auditor proof-5134-AE), X = 2^%d" % L)
    brute = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else None
    for ms, dat in res["data"].items():
        m = int(ms)
        divs = dat["divisors"]
        assert abs(sum(c_d(m, d) for d in divs) - 1) == 0
        print()
        print("=" * 100)
        print("modulus %d   divisors %s   c_d = %s" % (m, divs, ["%.5f" % float(c_d(m, d)) for d in divs]))
        for sh in range(2):
            off = dat["off"][sh]
            dia = dat["diag"][sh]
            csig = -1 if sh == 0 else 0
            # exact structure: type d  =>  d | M   (M = N - c_sigma)
            viol = 0
            for t, d in enumerate(divs):
                for r in range(m):
                    if (off[t][r] or dia[t][r]) and r % d != 0:
                        viol += off[t][r] + dia[t][r]
            CM = [2 * sum(off[t][r] for t in range(len(divs))) + sum(dia[t][r] for t in range(len(divs))) for r in range(m)]
            T = sum(CM)
            ndiag = sum(sum(row) for row in dia)
            D = {d: 2 * sum(off[t]) for t, d in enumerate(divs)}          # triples h1 != h2 by type of the smaller word (theorem definition)
            Dd = {d: D[d] + sum(dia[t]) for t, d in enumerate(divs)}      # diagonal counted by the type of its state
            assert sum(D.values()) + ndiag == T
            if brute is not None:
                b = brute["data"][ms][sh]
                assert brute["tot"][sh] == T, (brute["tot"][sh], T)
                bo = {(a[0], a[1]): a[2] for a in b[0]}
                bd = {(a[0], a[1]): a[2] for a in b[1]}
                for t, d in enumerate(divs):
                    for r in range(m):
                        assert bo.get((d, r), 0) == 2 * off[t][r], (m, sh, d, r)
                        assert bd.get((d, r), 0) == dia[t][r], (m, sh, d, r)
                print("  sheet %d: fast census == brute-force ordered enumeration (every (type, class) cell)" % sh)
            rem = Fr(0)
            remd = Fr(0)
            worst = None
            for r in range(m):
                g0 = math.gcd(m, r) if r else m
                pred = sum(Fr(d, m) * D[d] for d in divs if g0 % d == 0)
                predd = sum(Fr(d, m) * Dd[d] for d in divs if g0 % d == 0)
                if abs(CM[r] - pred) > rem:
                    rem = abs(CM[r] - pred)
                    worst = r
                remd = max(remd, abs(CM[r] - predd))
            print("  sheet %d: #T = %d  diagonal = %d  type-violations (d does not divide N - c_sigma) = %d" % (sh, T, ndiag, viol))
            print("     remainder of (B)  m max_c |#{N=c} - sum (d/m) D^(d)| / #T :  %.5f   (D without diagonal, as defined)   %.5f (diagonal counted in D)   worst class M = %d (N = %d mod %d)"
                  % (float(m * rem / T), float(m * remd / T), worst, (worst + csig) % m, m))
            print("     sampling noise sqrt(m/#T) = %.5f" % math.sqrt(m / T))
            print("     type shares F_d = D^(d)/#T :  " + "  ".join("d=%d: %.5f (c_d %.5f, ratio %.4f)" % (d, D[d] / T, float(c_d(m, d)), (D[d] / T) / float(c_d(m, d))) for d in divs))
            print("     (F_d - c_d) #T/(2X)        :  " + "  ".join("d=%d: %+.5f" % (d, (D[d] - float(c_d(m, d)) * T) / (2 * X)) for d in divs))
            # group factors
            grp = {}
            for r in range(m):
                g0 = math.gcd(m, r) if r else m
                grp.setdefault(g0, []).append(CM[r])
            line = []
            for g0 in sorted(grp):
                fac = m * sum(grp[g0]) / (len(grp[g0]) * T)
                phi = sum(d * D[d] / T for d in divs if g0 % d == 0)
                lim = float(L_m(m, g0))
                mx = max(abs(m * v / T - phi) for v in grp[g0])
                line.append("g0=%d: factor %.4f  Phi(g0) %.4f  limit L_m %.4f  max|class - Phi| %.4f" % (g0, fac, phi, lim, mx))
            print("     class groups by g0 = gcd(m, N - c_sigma):")
            for t in line:
                print("        " + t)


if __name__ == "__main__":
    main()
