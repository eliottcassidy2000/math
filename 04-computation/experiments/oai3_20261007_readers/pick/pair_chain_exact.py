"""Exact Haar probabilities for the pair chain by dynamic programming over states.
Lag-1 Mersenne switch: (k,e) = (0,1), v = p - 1 = 16 m (4 forced zero bits), then fair bits.
Prints exact P(merged by T-step n) (as a dyadic fraction and float) and the number of
distinct live states.  usage: python3 pair_chain_exact.py NMAX
"""
import sys
from fractions import Fraction
from pair_chain_verify import step

def exact(NMAX=60, k0=0, num0=1, forced_zero=4, report=None):
    dist = {(k0, num0): 1}          # weights are integers; total weight 2^(n - forced)
    merged = 0                      # merged weight (scaled to current denominator)
    denom_exp = 0
    rows = []
    for n in range(NMAX):
        new = {}
        if n < forced_zero:
            for st, w in dist.items():
                s2 = step(st[0], st[1], 0)[:2]
                new[s2] = new.get(s2, 0) + w
        else:
            merged *= 2; denom_exp += 1
            for st, w in dist.items():
                for b in (0, 1):
                    s2 = step(st[0], st[1], b)[:2]
                    new[s2] = new.get(s2, 0) + w
        m = new.pop((0, 0), 0)
        merged += m
        dist = new
        P = Fraction(merged, 2 ** denom_exp)
        rows.append((n + 1, P, len(dist)))
        if report and (n + 1) in report:
            print(f"n={n+1:4d}  P(merged by n) = {float(P):.6f}  live states = {len(dist)}", flush=True)
    return rows

if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    rep = set(range(5, NMAX + 1, 5)) | {27, 28, 29, 30, 31}
    rows = exact(NMAX, report=rep)
    for n, P, ns in rows:
        if n in (27, 28, 29, 30, 31) or n == NMAX:
            print(f"exact n={n}: P = {P.numerator}/{P.denominator}" if P.denominator < 2**40 else f"exact n={n}: P ~ {float(P):.10f}")
