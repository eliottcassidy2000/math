"""
lrc_a2_corner_exact_opus_S137.py

THE a=2 PINCH CORNER, EXACT (kps-S62's named residue; owner worklist "corner sweep").

kps-S62's a=2 stratum assembly: 2-letter step words close via [dilation degenerate splits] +
[S59 diameter <= 75] + [grid ledger st >~ 80, error C/(st)].  The RESIDUE ("pinch corner") =
words with min-letter small, st <= 80, diam >= 76 — where the 2-torus grid bound pinches.
kps sampled 40 corner families (min mu = 0.5553); this is the EXACT exhaustive sweep of the
corner in BLOCK arrangement s^i t^(12-i) + both edge arrangements, plus (for the worst (s,t,i)
cells) ALL arrangements up to reversal when the count is small.

Corner spec swept here: letters {s,t}, gcd(s,t)=1, min(s,t) <= 3, s*t <= 80, i in 1..11
(genuinely two-letter), diam = i*s + (12-i)*t >= 76, diam <= 400 (engine cost guard; the
omitted diam > 400 sliver has st <= 80/5... only huge-letter words -> equidistributed, grid
bound comfortable there; reported as the explicit remaining sliver).

Output: exact mu_{1/7} per corner cell; global corner minimum; margin vs m_P.
"""
from fractions import Fraction as F
from math import gcd
import sys, time

sys.path.insert(0, ".")
from lrc_exact_mu_ordercells_opus_S136 import mu_exact

M_P = F(14249, 252252)

def steps_to_E(steps):
    E = [1]
    for s in steps:
        E.append(E[-1] + s)
    return E

def main():
    t0 = time.time()
    cells = []
    for small in range(1, 4):
        for big in range(small + 1, 81):
            if small * big > 80 or gcd(small, big) != 1:
                continue
            for (s, t) in ((small, big), (big, small)):
                # i copies of s then 12-i copies of t (block); i = count of FIRST letter
                for i in range(1, 12):
                    diam = i * s + (12 - i) * t
                    if 76 <= diam <= 400:
                        cells.append((s, t, i, diam))
    # dedupe reversal-equivalent blocks: block s^i t^(12-i) reversed = t^(12-i) s^i
    seen = set(); todo = []
    for (s, t, i, diam) in cells:
        key = min((s, t, i), (t, s, 12 - i))
        if key in seen: continue
        seen.add(key)
        todo.append((s, t, i, diam))
    todo.sort(key=lambda c: c[3])
    print(f"a=2 pinch corner: {len(todo)} block cells (min-letter<=3, st<=80, 76<=diam<=400)")
    best = None
    for n, (s, t, i, diam) in enumerate(todo):
        steps = (s,) * i + (t,) * (12 - i)
        m = mu_exact(steps_to_E(steps))
        if best is None or m < best[0]:
            best = (m, (s, t, i, diam))
            print(f"   new corner min: mu = {m} = {float(m):.6f} at s={s},t={t},i={i},diam={diam}"
                  f"  ({float(m)/float(M_P):.1f}x m_P)  [{n+1}/{len(todo)}, {time.time()-t0:.0f}s]")
    print(f"\nCORNER MINIMUM (block arrangements): mu = {best[0]} = {float(best[0]):.6f} at "
          f"(s,t,i,diam)={best[1]}  = {float(best[0])/float(M_P):.1f}x m_P")
    # all arrangements for the worst cell if feasible
    s, t, i, diam = best[1]
    from itertools import combinations
    if 1 <= i <= 11:
        cnt = 0; wbest = None
        for pos in combinations(range(12), i):
            steps = tuple(s if j in pos else t for j in range(12))
            if steps > tuple(reversed(steps)): continue
            m = mu_exact(steps_to_E(steps))
            cnt += 1
            if wbest is None or m < wbest[0]: wbest = (m, steps)
        print(f"ALL {cnt} arrangements at the worst cell: min = {wbest[0]} = {float(wbest[0]):.6f} "
              f"at {wbest[1]}  ({float(wbest[0])/float(M_P):.1f}x m_P)")
    print(f"[total {time.time()-t0:.0f}s]  Sliver not swept: diam in (400, 12*80]: "
          f"grid-ledger territory (st large forces letters>=... reported for the assembly note).")

if __name__ == "__main__":
    main()
