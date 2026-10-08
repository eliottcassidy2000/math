#!/usr/bin/env python3
"""Audit E, task E (exact part): shortest merge times from (0,1), exact counts, and exact P(absorbed by L).

Exhaustive breadth-first enumeration of all coin words of length <= L (with multiplicities) in integer coordinates
(k, E), f = E / p^|k|. Pruning is RIGOROUS (unlike shortest_merge.py's |f| > 10^4 cut):
  every step has lambda >= 1/2 and |eta| <= 1/2 (Task A, sharp), so |f'| + 1 >= (|f| + 1)/2; hence a state at time t
  can be absorbed by time L only if |f_t| <= 2^(L-t) - 1 and |k_t| <= L - t.
Thus every absorption at a time <= L is counted exactly; P(absorbed by L) = sum_t count_t 2^-t exactly.
Also checks the reset word for p = 2^a - 1.
Usage: python3 e_shortest_exact.py p L [statecap]
"""
import sys
from fractions import Fraction as Fr

def step(p, k, E, b):
    s = E & 1
    if s == 0:
        if b == 0: return k, E >> 1
        if k >= 0: return k, (p * E + 1 - p ** k) >> 1
        return k, (p * E + p ** (-k) - 1) >> 1
    if b == 0:
        if k >= 0: return k + 1, (p * E + 1) >> 1
        return k + 1, (E + p ** (-k - 1)) >> 1
    if k >= 1: return k - 1, (E - p ** (k - 1)) >> 1
    return k - 1, (p * E - 1) >> 1

def enumerate_absorption(p, L, statecap=6_000_000, start=(0, 1)):
    frontier = {start: 1}
    absorbed = {}
    maxfront = 0
    for t in range(1, L + 1):
        new = {}
        rem = L - t
        lim = (1 << rem) - 1          # |f_t| <= 2^(L-t) - 1  <=>  |E| <= lim * p^|k|
        for (k, E), mult in frontier.items():
            for b in (0, 1):
                k2, E2 = step(p, k, E, b)
                if k2 == 0 and E2 == 0:
                    absorbed[t] = absorbed.get(t, 0) + mult
                    continue
                if abs(k2) > rem: continue
                if abs(E2) > lim * p ** abs(k2): continue
                key = (k2, E2)
                new[key] = new.get(key, 0) + mult
        frontier = new
        maxfront = max(maxfront, len(frontier))
        if len(frontier) > statecap:
            return absorbed, maxfront, t
    return absorbed, maxfront, None

def reset_word_check():
    out = []
    for a in range(2, 8):
        p = 2 ** a - 1
        k, E = 0, 1
        word = [0] + [0] * (a - 1) + [1]
        hit_at = None
        for i, b in enumerate(word, 1):
            k, E = step(p, k, E, b)
            if (k, E) == (0, 0):
                hit_at = i; break
        out.append(f"p={p} (a={a}): word {''.join(map(str, word))} absorbed at time {hit_at} (claimed a+1={a+1})")
    return out

if __name__ == '__main__':
    if sys.argv[1] == 'reset':
        print("\n".join(reset_word_check())); sys.exit()
    p, L = int(sys.argv[1]), int(sys.argv[2])
    cap = int(sys.argv[3]) if len(sys.argv) > 3 else 6_000_000
    absorbed, maxfront, aborted = enumerate_absorption(p, L, cap)
    if aborted:
        print(f"p={p}: ABORTED at t={aborted} (frontier > {cap}); absorptions so far {dict(sorted(absorbed.items()))}", flush=True)
        sys.exit()
    if absorbed:
        t0 = min(absorbed)
        P = sum(Fr(c, 2 ** t) for t, c in absorbed.items())
        print(f"p={p:2d} L={L}: first absorption at t={t0} with {absorbed[t0]} words -> exact bound {absorbed[t0]}*2^-{t0};"
              f" P(absorbed by {L}) = {float(P):.6e} (exact: {P}); absorptions by time: {dict(sorted(absorbed.items()))};"
              f" max frontier {maxfront}", flush=True)
    else:
        print(f"p={p:2d} L={L}: NO absorption at any time <= {L} (rigorous); max frontier {maxfront}", flush=True)
