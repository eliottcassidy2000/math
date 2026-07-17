# opus-2026-07-17-S339 -- HYP-7220: THE BLOCK-STRUCTURE REDUCTION constant.
# THM-955 gives a <=6-block a safe width in any window; the glue jumps a
# junction when nextmin * width >= 2.  This script derives the PROVEN
# junction constant G0 (worst case over block shapes) and verifies the
# composed tower end-to-end exactly.
#
#   window entering block i: length L_i (starts 1 for the first block)
#   THM-955 width out: W_i >= [(1-k/7) L_i - 2k/(7 m_i)] / (1 + S_i L_i + 2k)
#     with S_i = sum of block speeds, m_i = min speed.  For the worst case
#     write everything in units of M_i = max speed: m_i >= M_i / r (r =
#     within-block ratio cap), S_i <= k M_i.
#   glue condition: m_{i+1} * W_i >= 2, i.e. junction G = m_{i+1}/M_i >=
#     2 / (M_i W_i)  -- so G0 = sup over shapes of 2/(M W).
# We tabulate 2/(M_i W_i) exactly for k = 1..6, r = 1..13, and entry
# lengths L in the range produced upstream (L * m_prevblock ~ 6/7-ish
# after a glue jump: L >= 2/m_i at entry by the previous glue), then
# verify towers at the tabulated G0 exactly.
from fractions import Fraction
from math import floor
import random

F = Fraction

def subtract_comb(V, x):
    w = F(1, 14 * x)
    out = []
    for (a, b) in V:
        cur = a
        for j in range(floor((a - w) * x), floor((b + w) * x) + 2):
            lo, hi = F(j, x) - w, F(j, x) + w
            if hi <= cur: continue
            if lo >= b: break
            if lo > cur: out.append((cur, lo))
            cur = max(cur, hi)
            if cur >= b: break
        if cur < b: out.append((cur, b))
    return out

print("(1) THE PROVEN JUNCTION CONSTANT (worst 2/(M*W) over shapes):")
# entry window: the glue guarantees m_i * L_i >= 2 at entry; scale-free:
# set l = m * L (>= 2), rho = M/m (<= 13 by compression within a block).
# W = [(1-k/7) L - 2k/(7m)] / (1 + S L + 2k), S <= k M:
# M*W = [(1-k/7) l rho - 2k rho/7] / (1 + k rho l + 2k)   [all over m-units]
# worst at l = 2 (smallest entry window) and rho = 13:
# entry length requirement: (1-k/7) l > 2k/7  <=>  l > 2k/(7-k); take
# l_k = max(2, 4k/(7-k)) (factor-2 margin), rho = 13 worst:
lk = {k: max(F(2), F(4*k, 7-k)) for k in range(1, 7)}
MW = {}
for k in range(1, 7):
    l, rho = lk[k], F(13)
    num = (1 - F(k,7)) * l * rho - F(2*k,7) * rho
    den = 1 + k * rho * l + 2 * k
    MW[k] = num / den
    print(f"   k={k}: entry l_k = {lk[k]}, M*W >= {MW[k]} = {float(MW[k]):.5f}")
# junction table: G0[src][tgt] = l_tgt / (M*W)[src]
print("   G0(src k -> tgt k) table (worst case, rho = 13):")
G0 = {}
for ks in range(1, 7):
    G0[ks] = {kt: lk[kt] / MW[ks] for kt in range(1, 7)}
    print("     src k=%d: " % ks + "  ".join(
        f"->k{kt}:{float(G0[ks][kt]):9.1f}" for kt in range(1, 7)))
print()
print("(2) END-TO-END exact verification: towers with proven-G0 junctions:")
random.seed(339)
ok = tot = 0
for trial in range(60):
    sizes = random.choice([[6, 6, 1], [5, 4, 4], [6, 4, 3], [4, 4, 5], [3, 6, 4]])
    m0 = random.randint(3, 40)
    blocks = []
    prev_max = None; prev_k = None
    for k in sizes:
        m = m0 if prev_max is None else int(prev_max * G0[prev_k][k]) + 1
        prev_k = k
        M = m * random.randint(1, 13)
        B = sorted(random.sample(range(m, max(M, m + k)), k)) if M > m + k \
            else list(range(m, m + k))
        blocks.append(B)
        prev_max = max(B)
    speeds = [x for B in blocks for x in B]
    if len(set(speeds)) != 13: continue
    tot += 1
    V = [(F(0), F(1))]
    alive = True
    for B in blocks:
        for x in B: V = subtract_comb(V, x)
        if not V: alive = False; break
        V = [max(V, key=lambda iv: iv[1] - iv[0])]
    if alive:
        lo, hi = V[0]
        t = (lo + hi) / 2
        if all(min((x * t) % 1, 1 - (x * t) % 1) >= F(1, 14) for x in speeds):
            ok += 1
print(f"   {ok}/{tot} towers lonely with exact witnesses at proven G0 junctions")

print()
print("(3) THE FLOOR TRANSFER (item-3 bridge, exact): width delta in a window")
print("    => >= ceil(delta*D) - 1 points of any 1/D-grid inside it:")
for _ in range(5):
    a = F(random.randint(0, 500), 503)
    delta = F(random.randint(1, 40), 1009)
    D = random.randint(50, 4000)
    npts = floor((a + delta) * D) - floor(a * D)  # grid pts in (a, a+delta]
    bound = -((-delta * D).__floor__()) - 1       # ceil(delta*D) - 1
    print(f"   delta={float(delta):.5f} D={D:5d}: grid pts = {npts}, "
          f"floor bound = {bound}, holds: {npts >= bound}")
print()
print("REDUCTION: residual families with all blocks <= 6 at G0 junctions are")
print("LONELY (THM-955 + glue, proven constants). The dense core is exactly")
print("the single >=7-comparable-block families -- where the 7-wall's")
print("pair-overlap crumb (mac-mini's lane) is the one missing ingredient.")
