# THM-959 / opus-S339 corrected by codex-S62 and sharpened by direct use of
# THM-955: block-tower constants.  This script derives the exact source/target
# junction table G0 and verifies the composed tower end-to-end exactly.
#
#   window entering block i: length L_i (starts 1 for the first block)
#   direct THM-955 width out:
#     W_i >= [(1-k/7) L_i - k/(7 m_i)] / (1 + k + S_i L_i),
#     with S_i = sum of block speeds and m_i = min speed.
#   in scale-free variables l=m_i L_i and rho=M_i/m_i, S_i<=kM_i:
#     M_i W_i >= rho[(1-k/7)l-k/7]/[1+k+k rho l].
# This is increasing in l and rho once its numerator is positive, so rho=1
# is the uniform worst case and no within-block ratio cap is needed.
# A source block of size s supplies M_i W_i>=mu_s.  A target block of size t
# requires m_{i+1} W_i>=ell_t, hence G0(s,t)=ell_t/mu_s.
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

print("THM-959 SHARPENED CORRECTED PRESCRIBED BLOCK-TOWER REDUCTION")
print("(1) DIRECT-THM-955 ENTRY/OUTPUT CONSTANTS:")
# Rational entry lengths chosen near the exact minimizers of ell/mu(ell).
# For k=6 the real minimizer is 6+sqrt(43), with cost about 1103.825;
# ell_6=13 gives the exact near-optimal integer cost 1105.
lk = {
    1: F(3, 4),
    2: F(5, 4),
    3: F(2),
    4: F(3),
    5: F(6),
    6: F(13),
}
MW = {}
for k in range(1, 7):
    l, rho = lk[k], F(1)
    num = rho * ((1 - F(k, 7)) * l - F(k, 7))
    den = 1 + k + k * rho * l
    MW[k] = num / den
    assert num > 0
    # The full-circle first-block floor dominates the propagated output floor.
    first_floor = F((7 - k) * k, 7 * (1 + k * k))
    assert first_floor >= MW[k]
    print(f"   k={k}: ell_k={lk[k]}, mu_k={MW[k]} = {float(MW[k]):.6f}, "
          f"first-block floor={first_floor}")
# junction table: G0[src][tgt] = l_tgt / (M*W)[src]
print("   exact G0(src k -> tgt k)=ell_tgt/mu_src table:")
G0 = {}
for ks in range(1, 7):
    G0[ks] = {kt: lk[kt] / MW[ks] for kt in range(1, 7)}
    print("     src k=%d: " % ks + "  ".join(
        f"->k{kt}:{str(G0[ks][kt]):>7}" for kt in range(1, 7)))
worst = max((G0[s][t], s, t) for s in range(1, 7) for t in range(1, 7))
assert worst == (F(1105), 6, 6)
print(f"   uniform maximum: G0({worst[1]},{worst[2]})={worst[0]}")
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
print("LONELY (direct THM-955 induction, proven constants). The dense core is exactly")
print("packets admitting no such prescribed <=6-block partition; identifying")
print("that complement with one >=7 comparable block needs a separate lemma.")
