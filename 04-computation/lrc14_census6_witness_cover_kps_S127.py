# kind-pasteur-2026-07-09-S127
# The 6-witness greedy cover of the census-branch domain (replaces opus-S200's native_decide census).
#
# Grand-assembly branch order: non-covering(sieve) | NOT GapFamily(spread13) | dominant | repeated
#   | detuned | window-22 / residual.
# window-22 fires only for covering + GapFamily (ratio > 13) families. In [1,22]:
#   min >= 2  =>  max <= 22 <= 26 = 13*2 <= 13*min  =>  ratio <= 13  =>  spread13-dispatched.
# So the census-branch domain = {S subset [1,22], |S|=13, 1 in S, max > 13, covering} = 14002 families.
# CLAIM: 6 explicit rational witnesses cover it -- every family is lonely at one of them
#   (min_i ||v_i * p/Q|| >= 1/14, i.e. Q <= 14*min(r,Q-r) for all v, r = v*p mod Q -- EXACT arithmetic,
#    no native_decide).
from itertools import combinations
from math import gcd

def covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def lonely_at(S, p, Q):
    return all(14 * min((v * p) % Q, Q - (v * p) % Q) >= Q for v in S)

fams = [list(S) for S in combinations(range(1, 23), 13)
        if 1 in S and max(S) > 13 and covering(list(S))]
print(f"census-branch domain (window-22 & GapFamily, [1,22]): {len(fams)} families")

# THE 6 WITNESSES
WIT = [(12, 25), (9, 26), (7, 27), (11, 28), (4, 23), (11, 26)]
print(f"the 6 witnesses: {['%d/%d' % w for w in WIT]}")

# verify the cover is COMPLETE and report per-witness coverage
covered = [False] * len(fams)
for (p, Q) in WIT:
    fresh = 0
    for i, S in enumerate(fams):
        if not covered[i] and lonely_at(S, p, Q):
            covered[i] = True; fresh += 1
    print(f"  tau = {p}/{Q}: +{fresh} newly covered")
nleft = covered.count(False)
print(f"COVER {'COMPLETE (0 uncovered)' if nleft == 0 else 'INCOMPLETE: %d left' % nleft}")

# sanity: confirm 'lonely_at' matches the Mreach>=1/14 loneliness definition on a sample
def mreach_exact(S, p, Q):  # min_i ||v_i p/Q|| as a fraction m/Q
    return min(min((v * p) % Q, Q - (v * p) % Q) for v in S)
print("\nsanity: each family's realized clearance at its covering witness >= Q/14 (i.e. >= 1/14):")
import random; random.seed(0)
bad = 0
for i in random.sample(range(len(fams)), 12):
    S = fams[i]
    w = next((p, Q) for (p, Q) in WIT if lonely_at(S, p, Q))
    p, Q = w; cl = mreach_exact(S, p, Q)
    ok = 14 * cl >= Q
    if not ok: bad += 1
    print(f"  {S}  @ {p}/{Q}: clearance {cl}/{Q} >= 1/14? {ok}")
print(f"violations: {bad}")

# confirm the spread13 boundary claim
allcov = [list(S) for S in combinations(range(1, 23), 13) if covering(list(S))]
disp = sum(1 for S in allcov if min(S) >= 2)
print(f"\n[1,22] covering: {len(allcov)} total; {disp} have min>=2 (spread13-dispatched, ratio<=13); "
      f"{len(allcov) - disp} reach the census (= {len(fams)}). Non-covering are sieve-dispatched.")
