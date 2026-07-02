#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S98 -- HYP-3854 / THM-598: the resolved interval anti-covering
lemma.  Verifies the two pillars and the necessity of the hypotheses:

(a) FORCED PAIR INDEPENDENCE: for a Diophantine-resolved pair (no |m w1 + m' w2|
    small for small m, m'), the window overlap is ~ (2r)^2 L for EVERY phase
    choice -- adversarial minimization cannot beat independence.
(b) NEAR-RESONANT ESCAPE (necessity): a pair like (1000, 2001) inside a short
    window behaves as the frozen pattern (1,2): phases CAN zero the overlap.
(c) THE j=7 FLOOR: Diophantine-resolved 7-clusters at 2r = 1/7: adversarial
    coverage search vs the Bonferroni-2 bound  u/L >= 1 - j/7 + (j-1)/49
    = 6(8-j)/49 -> kappa_7 = 6/49 ~ 0.1224.
(d) THE FROZEN-PATTERN COVER BLOCK: near-equal combs (w, w+1, ..., w+6) over a
    window with L*max-diff << 1 tile like equal speeds: coverage -> ~1
    (why the unresolved column must renormalize instead).
"""
import random, itertools

R = 1 / 14  # arc half-width in phase units; density 2r = 1/7

def covered_fraction(ws, phases, L, grid=30000):
    cov = 0
    for i in range(grid):
        x = (i + 0.5) / grid * L
        for w, ph in zip(ws, phases):
            y = w * x - ph
            y -= round(y)
            if abs(y) < R:
                cov += 1
                break
    return cov / grid

def overlap_fraction(w1, w2, ph1, ph2, L, grid=30000):
    ov = 0
    for i in range(grid):
        x = (i + 0.5) / grid * L
        y1 = w1 * x - ph1; y1 -= round(y1)
        y2 = w2 * x - ph2; y2 -= round(y2)
        if abs(y1) < R and abs(y2) < R:
            ov += 1
    return ov / grid

print("=" * 78)
print("(a) FORCED INDEPENDENCE: Diophantine-resolved pair, window L=0.01")
print("=" * 78)
random.seed(98)
indep = (2 * R) ** 2
for (w1, w2) in [(1009, 1523), (997, 1601), (1013, 2311)]:
    # check: min |m w1 + m' w2| over 1<=|m|,|m'|<=8 is large
    minres = min(abs(m * w1 + mp * w2)
                 for m in range(-8, 9) for mp in range(-8, 9) if (m, mp) != (0, 0))
    worst = 1.0
    for trial in range(2000):
        ph = [random.random(), random.random()]
        f = overlap_fraction(w1, w2, ph[0], ph[1], 0.01)
        worst = min(worst, f)
    print(f"  ({w1},{w2}): min height-8 resonance |m w1+m' w2| = {minres}; "
          f"independence (2r)^2 = {indep:.5f}; adversarial MIN overlap/L = {worst:.5f} "
          f"({'FORCED (>= 0.7x indep)' if worst > 0.7 * indep else 'NOT forced'})")

print()
print("=" * 78)
print("(b) NEAR-RESONANT ESCAPE: (1000, 2001) ~ frozen pattern (1,2), L=0.01")
print("=" * 78)
w1, w2 = 1000, 2001
best_min = 1.0
for trial in range(4000):
    ph = [random.random(), random.random()]
    f = overlap_fraction(w1, w2, ph[0], ph[1], 0.01)
    best_min = min(best_min, f)
print(f"  |1*w2 - 2*w1| = {abs(w2 - 2*w1)} (height-2 resonance, drift over L: "
      f"{abs(w2-2*w1)*0.01:.3f} << 1)")
print(f"  adversarial MIN overlap/L over 4000 draws = {best_min:.5f} "
      f"vs independence {indep:.5f}  ({'ESCAPES (overlap ~ 0)' if best_min < 0.2*indep else 'no escape'})")

print()
print("=" * 78)
print("(c) THE j=7 FLOOR: resolved 7-cluster, adversarial coverage search")
print("=" * 78)
ws7 = [1009, 1523, 1997, 2311, 2749, 3181, 3617]  # pairwise non-resonant primes
minres7 = min(abs(m * a + mp * b)
              for a, b in itertools.combinations(ws7, 2)
              for m in range(-8, 9) for mp in range(-8, 9) if (m, mp) != (0, 0))
print(f"  speeds {ws7}; min pairwise height-8 resonance = {minres7}")
for L in (0.01, 0.005):
    best = 0.0
    best_ph = None
    for trial in range(2500):
        ph = [random.random() for _ in ws7]
        f = covered_fraction(ws7, ph, L, grid=12000)
        if f > best:
            best, best_ph = f, ph
    # coordinate-ascent refinement from the best random start
    step = 0.02
    for sweep in range(60):
        improved = False
        for i in range(len(ws7)):
            for d in (+step, -step):
                ph2 = list(best_ph); ph2[i] = (ph2[i] + d) % 1
                f = covered_fraction(ws7, ph2, L, grid=12000)
                if f > best:
                    best, best_ph = f, ph2; improved = True
        if not improved:
            step /= 2
            if step < 1e-4:
                break
    bound = 1 - 7 / 7 + (7 - 1) / 49  # 1 - 2rj + (2r)^2 (j-1) at j=7
    print(f"  L={L}: adversarial MAX coverage (random + ascent) = {best:.4f}; "
          f"uncovered >= {1-best:.4f} vs Bonferroni-2 floor 6/49 = {bound:.4f}  "
          f"({'BOUND HOLDS' if 1-best >= bound - 0.02 else 'VIOLATION?'})")

print()
print("=" * 78)
print("(d) FROZEN-PATTERN COVER BLOCK: near-equal combs (w..w+6), L small")
print("=" * 78)
w = 3000
ws_eq = [w + i for i in range(7)]
L = 0.0005   # max pairwise drift 6 * L = 0.003 << 2r/(w) arcs align
# equal-speed tiling phases: k/7 offsets
ph_tile = [(k / 7) * 1.0 for k in range(7)]
f_tile = covered_fraction(ws_eq, ph_tile, L, grid=20000)
best_eq = f_tile
for trial in range(1500):
    ph = [random.random() for _ in ws_eq]
    best_eq = max(best_eq, covered_fraction(ws_eq, ph, L, grid=8000))
print(f"  speeds {ws_eq}, L={L} (pairwise drift <= {6*L:.4f}):")
print(f"  tiling phases k/7: coverage = {f_tile:.4f}; best overall = {best_eq:.4f} "
      f"({'TILES (~1) -- must renormalize' if best_eq > 0.97 else 'does not tile'})")
