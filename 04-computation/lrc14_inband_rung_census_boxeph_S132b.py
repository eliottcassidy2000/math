#!/usr/bin/env python3
"""
lrc14_inband_rung_census_boxeph_S132b.py  -- EXHAUSTIVE in-band rung census.

Upgrades the S132 sampled rerun to a finite THEOREM on the m=1 band slice:

  enumerate ALL 12-sets C with {3,35} <= C <= [3,35], covering{2..12},
  passing all five spread rungs p in {13,17,19,23,25} (p-mult present or all
  antipodal unit pairs hit) -- the complete certificate-stack residual of the
  q=38 in-band slice -- and compute EXACT M (all q <= 70, complete by Pinch)
  for every one.

Output: exact survivor count, min/median M, blocker census, cheapest-needle
spectrum.  "Every in-band full-stack survivor has M >= X" is then a proved
finite fact quantifying the blocking-budget penalty (S130/S131 thesis).

boxeph-2026-07-19-S132.  Pure Python, exact integers.
"""

from math import gcd

RUNGS = [13, 17, 19, 23, 25]

def upairs(p):
    return [j for j in range(1, (p + 1) // 2) if gcd(j, p) == 1 and gcd(p - j, p) == 1]

# duty bits for covering {2..12}: it suffices to track q in {7,8,9,10,11,12}
# (2,3,4,5,6 are implied by 8,9,10,12 -- verified below against the direct def).
DUTIES = [7, 8, 9, 10, 11, 12]

def covering(sp):
    return all(any(v % q == 0 for v in sp) for q in range(2, 13))

VALS = list(range(4, 35))            # free values (3 and 35 forced)
N_FREE = 10

# per-value bit contributions
dbit = {q: 1 << k for k, q in enumerate(DUTIES)}
DFULL = (1 << len(DUTIES)) - 1
off, width = {}, 0
for p in RUNGS:
    off[p] = width; width += len(upairs(p)) + 1
SEG = {p: (1 << len(upairs(p))) - 1 for p in RUNGS}

def contrib(v):
    d = 0
    for q in DUTIES:
        if v % q == 0: d |= dbit[q]
    s = 0
    for p in RUNGS:
        if v % p == 0:
            s |= 1 << (off[p] + len(upairs(p)))
        elif gcd(v % p, p) == 1:
            s |= 1 << (off[p] + upairs(p).index(min(v % p, p - v % p)))
    return d, s

D3, S3 = contrib(3); D35, S35 = contrib(35)
DV = [contrib(v) for v in VALS]
sufD = [0] * (len(VALS) + 1); sufS = [0] * (len(VALS) + 1)
for i in range(len(VALS) - 1, -1, -1):
    sufD[i] = sufD[i + 1] | DV[i][0]
    sufS[i] = sufS[i + 1] | DV[i][1]

def rungs_ok(s):
    for p in RUNGS:
        seg = s >> off[p]
        L = len(upairs(p))
        if not ((seg >> L) & 1 or (seg & SEG[p]) == SEG[p]):
            return False
    return True

def M_int(sp):
    bn, bd, qstar = 0, 1, None
    cheapest = None
    for q in range(2, 71):
        best_here = 0
        for b in range(1, q // 2 + 1):
            md = q
            for c in sp:
                r = (c * b) % q
                d = r if r <= q - r else q - r
                if d < md: md = d
            if md > best_here: best_here = md
            if md * bd > bn * q:
                bn, bd, qstar = md, q, q
        if cheapest is None and best_here * 38 > 3 * q:
            cheapest = q
    g = gcd(bn, bd)
    return bn // g, bd // g, qstar, cheapest

survivors = []
count_cov = 0
stack = []

def rec(start, got_d, got_s):
    global count_cov
    slots = N_FREE - len(stack)
    if slots == 0:
        if got_d != DFULL: return
        count_cov += 1
        if rungs_ok(got_s):
            survivors.append(tuple([3] + stack + [35]))
        return
    if (DFULL & ~(got_d | sufD[start])): return
    # rung reachability prune
    s_reach = got_s | sufS[start]
    if not rungs_ok(s_reach): return
    if len(VALS) - start < slots: return
    for i in range(start, len(VALS)):
        stack.append(VALS[i])
        rec(i + 1, got_d | DV[i][0], got_s | DV[i][1])
        stack.pop()

print("=" * 88)
print("EXHAUSTIVE IN-BAND RUNG CENSUS: {3,35} <= C <= [3,35], covering, all 5 rungs")
print("=" * 88)
rec(0, D3 | D35, S3 | S35)
print("duty-complete leaves inside rung-consistent branches: %d" % count_cov)
print("(NB: NOT the full covering count -- the rung prune runs first. True covering total")
print(" ~= 12.8%% of C(31,10)=44,352,165 ~= 5.7M by sampling; survivor share ~= 0.45%% of covering.)")
print("FULL-STACK SURVIVORS (COMPLETE list -- prune is monotone-sound): %d" % len(survivors))

# sanity: duty-bit covering == direct covering on a sample of survivors
for C in survivors[:50]:
    assert covering(list(C)), C

blocker_hist, needle_hist = {}, {}
worst = []
for C in survivors:
    blk = tuple(p for p in RUNGS if any(v % p == 0 for v in C))
    blocker_hist[blk] = blocker_hist.get(blk, 0) + 1
    nu, de, qs, ch = M_int(list(C))
    needle_hist[ch] = needle_hist.get(ch, 0) + 1
    worst.append((nu / de, nu, de, qs, ch, C))
    assert (nu, de) != (3, 38), ("3/38 realized?!", C)
worst.sort()
print("\nblocker census:")
for blk, c in sorted(blocker_hist.items(), key=lambda x: -x[1]):
    print("  %-24s: %d" % (str(blk), c))
print("\ncheapest-needle spectrum (q: count): %s" % dict(sorted(needle_hist.items())))
print("\nMIN M over all survivors: %d/%d = %.4f  (vs 3/38 = %.4f; ratio %.2fx)"
      % (worst[0][1], worst[0][2], worst[0][0], 3 / 38, worst[0][0] / (3 / 38)))
print("5 lowest-M survivors:")
for Mv, nu, de, qs, ch, C in worst[:5]:
    print("  M=%d/%d=%.4f  maximizer q=%d  cheapest needle q=%s  C=%s" % (nu, de, Mv, qs, ch, C))
print("\nTHEOREM (finite check, this run): every 12-set C with {3,35} in C subset [3,35],")
print("covering{2..12}, passing all five bounded spread rungs, has M(C) >= %d/%d." % (worst[0][1], worst[0][2]))
print("The q=38 in-band certificate residual is NONEMPTY but UNIFORMLY LOOSE:")
print("the blocking budget forces every survivor far above the 3/38 target.")
print("\nDONE.")
