#!/usr/bin/env python3
"""
Parallel worker for the LRC sigma_2 scan. Identical math to lrc_sigma2_scan.run,
but the DFS is restricted to k-subsets whose smallest elements form the given
PREFIX (a strictly increasing tuple of the smallest chosen speeds).
Usage: python lrc_sigma2_worker.py k out_file p1 [p2 ...]
  where p1<p2<... is the fixed prefix (the smallest elements of the subset).
Writes a line: "prefix=<tuple> sigma2=<frac or None> wit=<tuple or None> surv=<n>"
"""
import sys, os
from fractions import Fraction
from math import gcd
from functools import reduce

def frac_norm_num_den(num, den):
    f = num % den
    if 2 * f <= den:
        return Fraction(f, den)
    return Fraction(den - f, den)

def minnorm_at(S, t):
    num = t.numerator; den = t.denominator
    m = None
    for v in S:
        d = frac_norm_num_den(v * num, den)
        if m is None or d < m:
            m = d
            if m == 0:
                return Fraction(0)
    return m

def exact_M(S):
    ds = set()
    Sl = list(S); n = len(Sl)
    for i in range(n):
        for j in range(i + 1, n):
            ds.add(Sl[i] + Sl[j])
            d = abs(Sl[i] - Sl[j])
            if d: ds.add(d)
    ds = {d for d in ds if d > 0}
    best = Fraction(0)
    for d in ds:
        for m in range(1, d):
            val = minnorm_at(S, Fraction(m, d))
            if val > best:
                best = val
    return best

def setgcd(S):
    return reduce(gcd, S)

def main():
    k = int(sys.argv[1]); out = sys.argv[2]
    prefix = [int(x) for x in sys.argv[3:]]
    BOX = 3 * k
    floor = Fraction(1, k + 1)
    mediant = Fraction(2, 2 * k + 1)
    QMAX = 3 * k

    probes = []
    for q in range(2, QMAX + 1):
        for m in range(1, q):
            if gcd(m, q) == 1:
                probes.append(Fraction(m, q))
    P = len(probes)

    kill = {}
    for v in range(1, BOX + 1):
        s = set()
        for idx, t in enumerate(probes):
            if frac_norm_num_den(v * t.numerator, t.denominator) < mediant:
                s.add(idx)
        kill[v] = s
    full = frozenset(range(P))

    union_from = [frozenset()] * (BOX + 2)
    acc = set()
    for v in range(BOX, 0, -1):
        acc |= kill[v]
        union_from[v] = frozenset(acc)

    best_sigma = [None]; best_witness = [None]; survivors = [0]

    def dfs(start, chosen, alive):
        rem = k - len(chosen)
        if rem == 0:
            if not alive:
                S = tuple(chosen)
                if setgcd(S) != 1:
                    return
                survivors[0] += 1
                M = exact_M(S)
                if M == floor:
                    return
                if M < mediant:
                    if best_sigma[0] is None or M < best_sigma[0]:
                        best_sigma[0] = M
                        best_witness[0] = S
            return
        if BOX - start + 1 < rem:
            return
        if not alive <= union_from[start]:
            return
        for v in range(start, BOX + 1):
            if BOX - v + 1 < rem:
                break
            nalive = alive - kill[v]
            if nalive and not (nalive <= union_from[v + 1]):
                continue
            chosen.append(v)
            dfs(v + 1, chosen, nalive)
            chosen.pop()

    # restrict: smallest elements form the given prefix (strictly increasing)
    ok = (len(prefix) <= k and all(prefix[i] < prefix[i+1] for i in range(len(prefix)-1))
          and (not prefix or prefix[-1] + (k - len(prefix)) <= BOX))
    if ok:
        alive = full
        for v in prefix:
            alive = alive - kill[v]
        start = (prefix[-1] + 1) if prefix else 1
        # only descend if remaining coverage still possible
        rem = k - len(prefix)
        if rem == 0:
            # full set already fixed
            dfs(start, list(prefix), alive)
        elif (BOX - start + 1) >= rem and alive <= union_from[start]:
            dfs(start, list(prefix), alive)

    sig = best_sigma[0]; wit = best_witness[0]
    with open(out, "w") as f:
        f.write(f"prefix={tuple(prefix)} sigma2={sig} wit={wit} surv={survivors[0]}\n")
        f.flush()

if __name__ == "__main__":
    main()
