# window_check.py -- auditor proof-5134-AE.
# The pre-registered window of odd N in [2^27, 2^28), recomputed from my own cumulative censuses
# (census_28.json minus census_27.json; N is odd, so the window is exactly the difference).
# Evaluates P1, P2, P3 against the registered file prereg_window_28_27.json and prints, next to the registered
# a-priori prediction, the first-order value computed from the MEASURED type shares of the window.
import json, math
from fractions import Fraction as Fr

HERE = "C:/tmp/wt-s92/05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/proof-5134-AE/"
SESS = "C:/tmp/wt-s92/04-computation/experiments/s92_prime_powers_20261111/"


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


def main():
    hi = json.load(open(HERE + "census_28.json"))
    lo = json.load(open(HERE + "census_27.json"))
    reg = json.load(open(SESS + "prereg_window_28_27.json"))
    closer = []
    worst_model = 0.0
    okP1 = okP2 = True
    for m in reg["moduli"]:
        ms = str(m)
        divs = hi["data"][ms]["divisors"]
        devs = []
        flats = []
        for sh in range(2):
            off = [[hi["data"][ms]["off"][sh][t][r] - lo["data"][ms]["off"][sh][t][r] for r in range(m)] for t in range(len(divs))]
            dia = [[hi["data"][ms]["diag"][sh][t][r] - lo["data"][ms]["diag"][sh][t][r] for r in range(m)] for t in range(len(divs))]
            CM = [2 * sum(off[t][r] for t in range(len(divs))) + sum(dia[t][r] for t in range(len(divs))) for r in range(m)]
            T = sum(CM)
            D = {d: 2 * sum(off[t]) for t, d in enumerate(divs)}
            dev = 0.0
            flat = 0.0
            cells = []
            for g0 in divs:
                cl = [r for r in range(m) if (math.gcd(m, r) if r else m) == g0]
                fac = [m * CM[r] / T for r in cl]
                mean = sum(fac) / len(fac)
                pr = reg["phi"]["%d,%d,%d" % (sh, m, g0)]
                phi = sum(d * D[d] / T for d in divs if g0 % d == 0)
                lim = float(sum(d * c_d(m, d) for d in divs if g0 % d == 0))
                dev = max(dev, abs(mean - pr))
                flat = max(flat, max(fac) - mean, mean - min(fac))
                worst_model = max(worst_model, abs(pr - phi))
                if g0 > 1 and m in (9, 27, 25):
                    closer.append(abs(mean - pr) < abs(mean - lim))
                cells.append("g0=%d: %.4f | reg %.4f | meas-types %.4f | lim %.4f" % (g0, mean, pr, phi, lim))
            print("m = %3d sheet %d (#T = %d): " % (m, sh, T) + ";  ".join(cells))
            devs.append(dev)
            flats.append(flat)
            okP1 = okP1 and dev <= reg["tol"][ms]
            okP2 = okP2 and flat <= reg["flat"][ms]
        print("    m = %3d: largest |group mean - registered prediction| %.4f / %.4f (tol %.3f)   flatness %.4f / %.4f (bound %.3f)"
              % (m, devs[0], devs[1], reg["tol"][ms], flats[0], flats[1], reg["flat"][ms]))
    print("P1:", "CONFIRMED" if okP1 else "FAILED", "  P2:", "CONFIRMED" if okP2 else "FAILED", "  P3: %d of %d" % (sum(closer), len(closer)))
    print("largest |registered a-priori prediction - first-order value from the measured window types| over all groups: %.5f" % worst_model)


if __name__ == "__main__":
    main()
