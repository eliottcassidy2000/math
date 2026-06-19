#!/usr/bin/env python3
"""
Exact LRC max-min second-spectrum scan for k=12,13,14.
Self-contained implementation from definitions (kind-pasteur).

M(S) = max_{t in (0,1)} min_{v in S} ||v t||,  ||x|| = dist(x, nearest int).
sigma_2(k) = min M(S) over gcd-1 sets S of k distinct positive ints, S NOT tight
(M != 1/(k+1)).  We hunt for the deepest-dipping (smallest M) such set.

Strategy (DFS with monotone probe pruning):
  floor = 1/(k+1), mediant = 2/(2k+1).
  A set S "survives" (is a candidate for M < mediant) iff for EVERY probe
  time t = m/q (q <= QMAX), min_v ||v t|| < mediant -- i.e. every probe is
  "killed" by some v with ||v t|| < mediant.
  As we add elements to a partial set, the set of killed probes only grows
  (monotone).  So:
    * If a probe is killed by NO available element it can never be killed ->
      but we cannot know "available" globally; instead we DFS and a set is a
      survivor only when ALL probes killed by the chosen elements.
  Pruning used: track alive (un-killed) probes for the partial set.  A full
  k-set survives iff alive == empty.  We also prune by a simple counting
  argument: each remaining slot can kill the alive probes only via the
  elements we add; we just require alive to become empty by depth k.
  (Conservative: we only descend; the killing-set per candidate element is
  precomputed so killing is O(1) set ops.)

For survivors: exact M via crossing candidates t=m/d, d in {v_i +- v_j}.
"""
import sys
from fractions import Fraction
from itertools import combinations
from math import gcd
from functools import reduce

def frac_norm_num_den(num, den):
    # ||num/den|| as Fraction, den>0
    f = num % den            # in [0,den)
    # fractional part = f/den in [0,1)
    if 2 * f <= den:
        return Fraction(f, den)
    else:
        return Fraction(den - f, den)

def minnorm_at(S, t):
    num = t.numerator
    den = t.denominator
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
    Sl = list(S)
    n = len(Sl)
    for i in range(n):
        for j in range(i + 1, n):
            ds.add(Sl[i] + Sl[j])
            ds.add(Sl[i] - Sl[j])  # positive since sorted? not necessarily; abs
            d = abs(Sl[i] - Sl[j])
            if d: ds.add(d)
    ds = {d for d in ds if d > 0}
    best = Fraction(0)
    for d in ds:
        num = None
        for m in range(1, d):
            t = Fraction(m, d)
            val = minnorm_at(S, t)
            if val > best:
                best = val
    return best

def setgcd(S):
    return reduce(gcd, S)

def run(k):
    BOX = 3 * k
    floor = Fraction(1, k + 1)
    mediant = Fraction(2, 2 * k + 1)
    QMAX = 3 * k

    # probes indexed
    probes = []
    for q in range(2, QMAX + 1):
        for m in range(1, q):
            if gcd(m, q) == 1:
                probes.append(Fraction(m, q))
    P = len(probes)

    # For each speed v in 1..BOX, which probes does it kill? (||v t|| < mediant)
    # kill[v] = frozenset of probe indices
    kill = {}
    for v in range(1, BOX + 1):
        s = set()
        for idx, t in enumerate(probes):
            if frac_norm_num_den(v * t.numerator, t.denominator) < mediant:
                s.add(idx)
        kill[v] = s
    full = frozenset(range(P))

    # suffix-union of kills: union_from[s] = union of kill[v] for v in [s, BOX]
    union_from = [frozenset()] * (BOX + 2)
    acc = set()
    for v in range(BOX, 0, -1):
        acc |= kill[v]
        union_from[v] = frozenset(acc)

    best_sigma = [None]
    best_witness = [None]
    survivors = [0]

    def dfs(start, chosen, alive):
        rem = k - len(chosen)
        if rem == 0:
            if not alive:  # all probes killed -> survivor candidate
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
        # need at least rem more speeds available
        if BOX - start + 1 < rem:
            return
        # coverage prune: speeds in [start,BOX] must be able to kill every
        # remaining alive probe; else dead end.
        if not alive <= union_from[start]:
            return
        for v in range(start, BOX + 1):
            if BOX - v + 1 < rem:
                break
            nalive = alive - kill[v]
            # quick prune: remaining speeds [v+1,BOX] must cover nalive
            if nalive and not (nalive <= union_from[v + 1]):
                # even with everything left we can't kill nalive -> skip this v
                # (still try? no: this v chosen, can't recover) -> prune branch
                continue
            chosen.append(v)
            dfs(v + 1, chosen, nalive)
            chosen.pop()

    dfs(1, [], full)
    return best_sigma[0], best_witness[0], survivors[0], floor, mediant

if __name__ == "__main__":
    for k in [12, 13, 14]:
        sig, wit, surv, fl, med = run(k)
        print(f"k={k}: floor=1/{k+1}={fl}={float(fl):.6f}  mediant=2/{2*k+1}={med}={float(med):.6f}")
        print(f"   survivors(all probes killed)={surv}")
        if sig is None:
            print(f"   sigma_2: NONE strictly below mediant -> sigma_2 = mediant = {med} (family bound)")
        else:
            gval = sig - fl
            print(f"   sigma_2(k) = {sig} = {float(sig):.6f}")
            print(f"   witness S = {wit}")
            print(f"   g(k)=sigma_2-floor = {gval} = {float(gval):.6f}")
        sys.stdout.flush()
