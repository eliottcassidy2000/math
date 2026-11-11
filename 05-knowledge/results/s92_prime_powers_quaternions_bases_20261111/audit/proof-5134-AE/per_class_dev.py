# per_class_dev.py -- auditor proof-5134-AE.
# Signed per-class deviations m (#{M = c'} - sum_{d | gcd(m, c')} (d/m) D^(d)) / #T of the first-order law (B),
# by size, to see whether the remainder is a decaying oscillation (as the proof says) or a persistent offset.
# Also the RMS over the classes.   M = N - c_sigma, so c' = 0 is the distinguished class.
import json, math
from fractions import Fraction as Fr

HERE = "C:/tmp/wt-s92/05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/proof-5134-AE/"


def main():
    sizes = [22, 24, 26, 27, 28, 30]
    data = {L: json.load(open(HERE + "census_%d.json" % L)) for L in sizes}
    for m in (9, 25, 15, 49):
        for sh in range(2):
            print("modulus %d sheet %d: m (count - first-order)/#T by class c' = M mod m" % (m, sh))
            for L in sizes:
                dat = data[L]["data"][str(m)]
                divs = dat["divisors"]
                off = dat["off"][sh]
                dia = dat["diag"][sh]
                CM = [2 * sum(off[t][r] for t in range(len(divs))) + sum(dia[t][r] for t in range(len(divs))) for r in range(m)]
                T = sum(CM)
                D = {d: 2 * sum(off[t]) for t, d in enumerate(divs)}
                dev = []
                for r in range(m):
                    g0 = math.gcd(m, r) if r else m
                    pred = sum(Fr(d, m) * D[d] for d in divs if g0 % d == 0)
                    dev.append(float(m * (CM[r] - pred) / T))
                rms = math.sqrt(sum(v * v for v in dev) / m)
                if m <= 15:
                    print("   2^%d: " % L + " ".join("%+.4f" % v for v in dev) + "   rms %.4f  max %.4f" % (rms, max(abs(v) for v in dev)))
                else:
                    # group summary for a prime square p^2 (25, 49): classes by c' mod p
                    p = math.isqrt(m)
                    byp = [sum(dev[r] for r in range(m) if r % p == a) / p for a in range(p)]
                    print("   2^%d: mean over the classes c' = a mod %d, a = 0..%d: " % (L, p, p - 1) + " ".join("%+.4f" % v for v in byp) + "   rms %.4f  max %.4f" % (rms, max(abs(v) for v in dev)))


if __name__ == "__main__":
    main()
