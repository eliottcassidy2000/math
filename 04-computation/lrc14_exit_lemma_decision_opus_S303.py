#!/usr/bin/env python3
"""THM-783: the r=8 visitor-cluster package (opus-2026-07-14-S303).

HISTORICAL/CORRECTION NOTE (codex-S10, THM-784/MISTAKE-147): this script
referees necessary algebraic identities and records a bounded adversarial
census. It is not a decisive experiment and does not prove a universal exit
constant. A fixed slow rainbow plus an arbitrarily fast owner gives unbounded
raw wall runs. THM-783 now states the corrected scopes.

THE ALGEBRA (canonized in THM-783):
  phi-recurrence: starting from a covered simple wall i, if the next global
  event coordinate is also simple, its wall is covered iff
  phi_{i+1} = phi_i + w_{o_i}^{-1} (mod 7). Same-owner steps satisfy it
  identically. The checks below test necessity on observed blocking runs; they
  do not independently test the iff or the needed anchor/simplicity hypotheses.
  PERIOD-SUM LAW: telescoping over one period of any owner b walling twice:
  sum over a != b of n_a * w_a^{-1} = 0 (mod 7), n_a = #a-walls in the period
  (a Beatty count taking two adjacent values).
  FLIP-BREAK: if exactly one owner's count flips between consecutive periods,
  the law fails and the run ends.
  THE LOOPHOLE CANDIDATE: synchronized flips (e.g. owner w=2 makes every odd
  owner's count alternate in lockstep; mod-7 cancellation reduces to
  sum over odd owners of +-w^{-1} = 0 mod 7, ONE tunable condition).

THE CENSUS: 300 sampled tuples (heights to 10^4) plus a 600-step anneal, with
targeted synchronized/near-ratio/mod-7-tuned seeds, in windows of length 1/5.
The observed maximum 6 is telemetry for those banks only. The script does not
test the withdrawn de-phase claim or the corrected conditional extent theorem.
THM-779/MISTAKE-147's coherent seven-owner stalk plus one fast owner has
unbounded raw runs and lies outside these generators.

PARTS: 1 referee vs S302 checker; 2 recurrence necessity on covered runs; 3 the
bounded census + packet experiments; 4 period-sum samples; 5 telemetry summary.
"""
import sys, os, random
from fractions import Fraction as F
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

random.seed(303)

def inv7(w): return pow(w % 7, 5, 7)

def walls_window(W8, lo, hi):
    """all walls (x=(2m+1)/(2w) in (lo,hi)) as sorted exact tuples (2m+1, 2w) -> compare by cross-mult."""
    ev = []
    for w in W8:
        mlo = int(2 * w * lo); mhi = int(2 * w * hi) + 2
        for tm in range(mlo, mhi):          # tm = 2m+1 odd
            if tm % 2 == 0: continue
            # x = tm/(2w) in (lo,hi)?
            if F(tm, 2 * w) <= lo or F(tm, 2 * w) >= hi: continue
            ev.append((F(tm, 2 * w), w, (tm - 1) // 2))
    ev.sort()
    return ev

def token_at(v, x):
    """token of owner v at exact x (None if v is at a wall there)."""
    q = v * x  # Fraction
    n, d = q.numerator, q.denominator
    if d == 2 and n % 2 == 1: return None       # half-integer: wall
    # nearest integer (no ties now)
    r = (2 * n + d) // (2 * d)
    return (-inv7(v) * r) % 7

def wall_ok(W8, x, w):
    """absolute wall condition: the 7 non-walling tokens are exactly F_7."""
    toks = []
    for v in W8:
        if v == w: continue
        t = token_at(v, x)
        if t is None: return False              # simultaneous wall
        toks.append(t)
    return len(set(toks)) == 7

def max_anchored_run(W8, lo, hi):
    """max streak of consecutive walls all satisfying the absolute wall condition."""
    ev = walls_window(W8, lo, hi)
    best = cur = 0
    for x, w, m in ev:
        if wall_ok(W8, x, w):
            cur += 1; best = max(best, cur)
        else:
            cur = 0
    return best, len(ev)

def phi_recurrence_streak(W8, lo, hi):
    """max streak where consecutive-wall phi-recurrence holds (algebraic form)."""
    ev = walls_window(W8, lo, hi)
    best = cur = 0
    for i in range(1, len(ev)):
        x0, w0, m0 = ev[i-1]; x1, w1, m1 = ev[i]
        if x0 == x1: cur = 0; continue
        if (inv7(w1) * m1 - inv7(w0) * (m0 + 1)) % 7 == 0:
            cur += 1; best = max(best, cur)
        else:
            cur = 0
    return best

fails = 0
def check(name, cond):
    global fails
    print(("  [PASS] " if cond else "  [FAIL] ") + name)
    if not cond: fails += 1

print("=" * 88)
print("PART 1+2 -- referee: absolute wall condition vs phi-recurrence consistency")
print("=" * 88)
# On a satisfied streak, recurrence must hold; conversely recurrence streak + anchored
# start = wall-condition streak. Empirical identity check of streak lengths:
agree = True
for trial in range(40):
    W8 = random.sample([w for w in range(1, 300) if w % 7], 8)
    a, _ = max_anchored_run(W8, F(1, 100), F(1, 100) + F(1, 5))
    b = phi_recurrence_streak(W8, F(1, 100), F(1, 100) + F(1, 5))
    # anchored run of length L involves L walls = L-1 recurrence steps + anchor;
    # recurrence streaks can exceed anchored runs (unanchored algebra) but never be shorter than a-1
    if b < a - 1: agree = False
check("phi-recurrence streak >= anchored run - 1 on 40 random instances (necessity)", agree)

print()
print("=" * 88)
print("PART 3 -- THE DECISION CENSUS (anchored runs; heights to 10^4; targeted packets)")
print("=" * 88)
def census(name, tuples, span=F(1, 5)):
    runs = []
    peak = (0, None)
    for W8 in tuples:
        r, nw = max_anchored_run(W8, F(1, 1000), F(1, 1000) + span)
        runs.append(r)
        if r > peak[0]: peak = (r, sorted(W8))
    runs.sort()
    print(f"  {name:<44} n={len(runs):<4} median={runs[len(runs)//2]:<3} "
          f"max={peak[0]:<3} at {peak[1]}")
    return peak

# (a) generic
gen = [random.sample([w for w in range(1, 10000) if w % 7], 8) for _ in range(80)]
p1 = census("generic random (w <= 10^4)", gen)
# (b) w=2 synchronized packets: {2} + 7 odds with sum inv7 == 0 mod 7
sync = []
odds = [w for w in range(3, 4000, 2) if w % 7]
while len(sync) < 80:
    S = random.sample(odds, 7)
    if sum(inv7(w) for w in S) % 7 == 0:
        sync.append([2] + S)
p2 = census("synchronized {2}+odds, sum w^-1 = 0 (7)", sync)
# (c) near-ratio packets: w_a ~ n_a * base + tuned residues
near = []
while len(near) < 80:
    base = random.randint(8, 60)
    if base % 7 == 0: continue
    ns = random.sample(range(2, 40), 7)
    W8 = [base] + [n * base + random.choice([-1, 1]) for n in ns]
    W8 = [w for w in W8 if w % 7 and w > 0]
    if len(set(W8)) == 8: near.append(W8)
p3 = census("near-ratio packets n*base +- 1", near)
# (d) all same residue mod 7
same = []
while len(same) < 60:
    r0 = random.randint(1, 6)
    S = random.sample([w for w in range(1, 6000) if w % 7 == r0], 8)
    same.append(S)
p4 = census("all w = r (mod 7)", same)
# (e) anneal from the best found overall
best = max([p1, p2, p3, p4], key=lambda t: t[0])
W8 = list(best[1]); cur, _ = max_anchored_run(W8, F(1, 1000), F(1, 1000) + F(1, 5))
pool = [w for w in range(1, 10000) if w % 7]
for step in range(600):
    W2 = list(W8); W2[random.randrange(8)] = random.choice(pool)
    if len(set(W2)) < 8: continue
    r2, _ = max_anchored_run(W2, F(1, 1000), F(1, 1000) + F(1, 5))
    if r2 >= cur: W8, cur = W2, r2
print(f"  ANNEALED PEAK: max anchored run = {cur} at {sorted(W8)}")

print()
print("=" * 88)
print("PART 4 -- period-sum law verification on random blocking-run fragments")
print("=" * 88)
law_ok = True; tested = 0
for trial in range(6000):
    W8 = random.sample([w for w in range(1, 400) if w % 7], 8)
    ev = walls_window(W8, F(1, 50), F(1, 50) + F(1, 3))
    # find owner-b consecutive wall pairs whose whole in-between streak satisfies wall_ok
    for i in range(len(ev)):
        x0, w0, m0 = ev[i]
        if not wall_ok(W8, x0, w0): continue
        for j in range(i + 1, len(ev)):
            xj, wj, mj = ev[j]
            if not wall_ok(W8, xj, wj): break
            if wj == w0:
                if mj == m0 + 1:
                    inner = ev[i + 1:j]
                    ssum = sum(inv7(w) for _, w, _ in inner) % 7
                    tested += 1
                    if ssum != 0: law_ok = False
                break
        if tested >= 40: break
    if tested >= 40: break
check(f"period-sum law holds on all {tested} in-run periods found", law_ok)
check(f"period-sum sample size adequate ({tested} periods)", tested >= 20)

print()
print("=" * 88)
print("PART 5 -- the single-visitor break + cluster balance (the unconditional teeth)")
print("=" * 88)
sv = 0; fper = 0; pb_ok = True; pairs = 0
for trial in range(4000):
    W8 = random.sample([w for w in range(1, 500) if w % 7], 8)
    f = max(W8)
    ev = walls_window(W8, F(1, 50), F(1, 50) + F(1, 3))
    i = 0; n = len(ev)
    while i < n:
        if not wall_ok(W8, ev[i][0], ev[i][1]): i += 1; continue
        j = i
        while j + 1 < n and wall_ok(W8, ev[j+1][0], ev[j+1][1]): j += 1
        fidx = [t for t in range(i, j + 1) if ev[t][1] == f]
        for a, b in zip(fidx, fidx[1:]):
            if ev[b][2] != ev[a][2] + 1: continue
            vis = [ev[t][1] for t in range(a + 1, b)]
            fper += 1
            if len(vis) == 1: sv += 1
            if len(vis) == 2:
                pairs += 1
                if (vis[0] + vis[1]) % 7 != 0: pb_ok = False
        i = j + 1
check(f"single-visitor break: 0 single-visitor in-run f-periods (of {fper} examined)", sv == 0)
check(f"pair balance w+w' = 0 mod 7 on all {pairs} two-visitor periods", pb_ok)

print()
print("=" * 88)
print(f"SUMMARY: {'ALL CHECKS PASSED' if fails == 0 else f'{fails} FAILURES'}")
print("=" * 88)
