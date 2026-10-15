#!/usr/bin/env python3
"""The safe Mahler carry language and its slack automaton (exact rationals).

A finite carry word w is safe iff every suffix sum sum_j w_(n+j) (2/3)^(j+1) is < 1 (the finite language of THM-2228/
THM-4072's strict safe shift).  Slack form: put m = min over suffixes of (1 - S_suffix)(3/2)^(length), with the empty
suffix counted (m <= 1).  Appending a letter e maps every suffix slack r to (3/2) r - e and adds a fresh empty suffix
(slack 1), so the whole language is recognised by the one-dimensional automaton

    m  ->  min((3/2) m - e, 1),   valid iff (3/2) m - e > 0,   m_0 = 1.

This script: (1) checks the automaton against direct suffix sums for all words up to length 16; (2) counts safe words
N(l) and the maximal number of ones k_max(l) for l <= LMAX by exact DP over automaton states; (3) the greedy expansion
d of 1 in base 3/2 and its gap sequence; (4) the minimal gap successor rules; (5) the slope bound for Collatz affine maps
along safe words and the invariant interval of the corresponding real IFS.
"""
from fractions import Fraction as F
from itertools import product
import json

LMAX = 60


def step(m, e):
    t = F(3, 2) * m - e
    return None if t <= 0 else min(t, F(1))


def safe_direct(w):
    L = len(w)
    for n in range(L):
        if sum(w[n + j] * F(2, 3) ** (j + 1) for j in range(L - n)) >= 1:
            return False
    return True


def safe_auto(w):
    m = F(1)
    for e in w:
        m = step(m, e)
        if m is None:
            return False
    return True


out = {}
# (1) automaton = direct definition
nchk = 0
for L in range(0, 17):
    for w in product((0, 1), repeat=L):
        assert safe_direct(w) == safe_auto(w), w
        nchk += 1
print(f"(1) slack automaton = suffix-sum definition on all {nchk} words of length <= 16")

# (2) exact DP: state -> (count, max ones)
states = {F(1): (1, 0)}
N, KMAX, NSTATES = [1], [0], [1]
for L in range(1, LMAX + 1):
    new = {}
    for m, (c, k) in states.items():
        for e in (0, 1):
            m2 = step(m, e)
            if m2 is None:
                continue
            c0, k0 = new.get(m2, (0, -1))
            new[m2] = (c0 + c, max(k0, k + e))
    states = new
    N.append(sum(c for c, _ in states.values()))
    KMAX.append(max(k for _, k in states.values()))
    NSTATES.append(len(states))
print("(2) N(l) safe words, l = 0..24:", N[:25])
print("    N(l)/(3/2)^l at l = 10, 20, 30, 40, 50, 60:",
      [round(N[l] / 1.5 ** l, 4) for l in (10, 20, 30, 40, 50, 60)])
print("    k_max(l), l = 0..30:", KMAX[:31])
viol = [l for l in range(LMAX + 1) if KMAX[l] > (l + 2) // 3 + (1 if l % 3 == 0 and l > 0 else 0)]
print("    k_max(l) - ceil(l/3):", [KMAX[l] - (-(-l // 3)) for l in range(LMAX + 1)])
print("    automaton states reached at l = 10, 20, 40, 60:", [NSTATES[l] for l in (10, 20, 40, 60)])
out['N'] = N
out['kmax'] = KMAX

# (3) greedy expansion of 1 in base 3/2: 1 = sum d_j (2/3)^j
r = F(1)
d = []
for j in range(1, 121):
    if F(2, 3) ** j <= r:
        d.append(1)
        r -= F(2, 3) ** j
    else:
        d.append(0)
ones = [j + 1 for j, x in enumerate(d) if x]
gaps = [b - a for a, b in zip(ones, ones[1:])]
print("(3) greedy expansion of 1 in base 3/2, first 60 digits:", ''.join(map(str, d[:60])))
print("    positions of ones:", ones[:16], " gaps:", gaps[:15])
out['d_gaps'] = gaps

# (4) gap successor rules: for each gap g (2..12) after a one placed at the slack-maximal moment, the minimal next gap
#     (over all reachable states), computed exactly by exploring states that occur right after a one.
def run_gaps(gs):
    """is the word 1 0^(g1-1) 1 0^(g2-1) ... 1 safe (as a factor, i.e. from the best possible entry state m = 1)?"""
    w = []
    for g in gs:
        w += [1] + [0] * (g - 1)
    w += [1]
    return safe_auto(w)

print("(4) gap rules (from the slack-maximal entry state):")
print("    gap 1 allowed?", run_gaps([1]))
for g in range(2, 6):
    mins = [h for h in range(1, 30) if run_gaps([g, h])]
    print(f"    after gap {g}: minimal next gap {mins[0] if mins else None}")
mins = [h for h in range(1, 40) if run_gaps([2, 6, h])]
print(f"    after gaps (2, 6): minimal next gap {mins[0]};  after (2, 7): {[h for h in range(1, 40) if run_gaps([2, 7, h])][0]}")
print(f"    after (3, 2): minimal next gap {[h for h in range(1, 40) if run_gaps([3, 2, h])][0]};  "
      f"(3, 3, 2): {[h for h in range(1, 40) if run_gaps([3, 3, 2, h])][0]}")

# (5) slopes and the real IFS of Collatz affine maps along safe words
# along a word w of length l with k ones: T_w(y) = (3^k y + C)/2^l; the maximal slope is 3^kmax(l)/2^l
sl = [3 ** KMAX[l] / 2 ** l for l in range(LMAX + 1)]
print("(5) maximal slope 3^kmax(l)/2^l at l = 1..12:", [round(s, 4) for s in sl[1:13]],
      " ratio to (3^(1/3)/2)^l at l=60:", round(sl[60] / (3 ** (1 / 3) / 2) ** 60, 4))
# sup of T_w(0) = C/2^l over safe words: DP over (state, value) keeping the max value per state
best = {F(1): F(0)}
sup_hist = []
for L in range(1, 41):
    new = {}
    for m, val in best.items():
        for e in (0, 1):
            m2 = step(m, e)
            if m2 is None:
                continue
            v2 = (3 * val + 1) / 2 if e else val / 2
            if m2 not in new or v2 > new[m2]:
                new[m2] = v2
    best = new
    sup_hist.append(max(best.values()))
print("    sup over safe words of T_w(0) (orbit of 0 along w), l = 1..40: max", float(max(sup_hist)),
      " last", float(sup_hist[-1]))
out['sup_T0'] = float(max(sup_hist))
json.dump(out, open('safe_language.json', 'w'))
