#!/usr/bin/env python3
r"""
lrc14_triple_atom_law_monad_S11.py   (monad-explorer-2026-07-07-S11, HYP-5157, part 8)

THE TRIPLE ATOM LAW + THE LAYER-CAKE ROUTE TO Sigma_3 AP-EXTREMALITY.

LAW (proof sketch, 5 lines): I(0,p,q; theta) := int_x (theta - span{0, px, qx})_+ dx.
  Substitute y = gx (g = gcd(p,q)) -- measure-preserving mod 1 => WLOG gcd(p',q') = 1.
  Contributions need all three of {0, p'y, q'y} in a closed theta-arc.  Near y = m/q' + delta,
  p'y sits at distance ||p'm/q'|| >= 1/q' from 0 unless p'm == 0 (mod q') <=> m == 0.
  For theta <= 1/q' (q' <= 7 at theta = 1/7), only m = 0 contributes: there the span is
  q'|delta|, giving 2 * int_0^{theta/q'} (theta - q' d) dd = theta^2 / q'.
  HENCE  I = theta^2 * gcd(p,q)/q  whenever q/gcd(p,q) <= 1/theta = 7.
  For q' > 7 the r != 0 windows open and I ~ theta^3 with arithmetic wiggles (exact per case).

COROLLARY (AP closed form):  Sigma_3(AP_k) for k-1 <= 7 is
  theta^2 * sum_{d=2}^{k-1} (k-d) * P(d)/d,   P(d) = Pillai = sum_{j<d} gcd(j,d)  (A018804).
  k=8: theta^2 * [6*1/2 + 5*2/3 + 4*4/4 + 3*4/5 + 2*9/6 + 1*6/7] = 1742/5145.  (P(d): 1,2,4,4,9,6)

LAYER-CAKE DOMINANCE (the route to Sigma_3 <= Sigma_3(AP) for ALL E):
  Sigma_3(E) = theta^2 * sum_{q'=2}^{7} N_{q'}(E)/q' + tail(E),
  N_{q'}(E) = #{triples with reduced max-diff q'}, tail = sum over q' > 7 of exact I.
  Since 1/q' is decreasing, Abel summation turns "M_m(E) := sum_{q'<=m} N_{q'}(E) <= M_m(AP)
  for every m <= 7" into "the theta^2-part of Sigma_3(E) <= that of AP".
  M_2 = #3-APs: the classical "AP maximizes 3-AP count".  TEST the full dominance here.

TESTS:
  1. law verification: exhaustive (p,q) with q <= 30: I == theta^2 g/q iff q/g <= 7
  2. AP closed form at k=8 (exact) + Sigma_3(AP_13) exact (with tail)
  3. M_m dominance over battery + random/structured shapes (k=8, all m = 2..7)
  4. Sigma_4 AP-extremality probe (jump descent)
  5. second-moment universality probe: the two-window S2-truncation as a shape-independent
     upper bound for E[V^2] -- computed exactly for several shapes (are they equal?)
"""
import sys, random
from fractions import Fraction as F
from itertools import combinations
from math import gcd

exec(open('/home/bigo/math/04-computation/lrc14_cubic_moment_gate_monad_S11.py')
     .read().split("if __name__")[0])
_src = open('/home/bigo/math/04-computation/lrc14_window_correlation_calculus_monad_S11.py').read()
_body = _src[:_src.rfind('if __name__')]
_body = _body[_body.index('THETA = F(1, 7)'):]
exec(_body)

THETA = F(1, 7)

def pillai(d):
    return sum(gcd(j, d) for j in range(1, d))

def reduced_pattern(a, b, c):
    p, q = b - a, c - a
    g = gcd(p, q)
    return p // g, q // g, g

if __name__ == "__main__":
    print("=" * 100)
    print("PART 8a -- LAW VERIFICATION: I(0,p,q) == theta^2 * g/q  iff  q/g <= 7   (q <= 30 exhaustive)")
    print("=" * 100)
    ok = viol = 0
    exceptions = []
    for q in range(2, 31):
        for p in range(1, q):
            g = gcd(p, q)
            I = exact_I_subset([0, p, q])
            law = THETA ** 2 * g / q
            if q // g <= 7 and q % g == 0 and (q // g) * g == q:
                pass
            qred = q // g
            if qred <= 7:
                if I == law:
                    ok += 1
                else:
                    viol += 1
                    exceptions.append((p, q, I, law))
            else:
                # record how far from theta^3 * g-scaling
                pass
    print(f"  q' <= 7 cases: {ok} exact matches, {viol} violations")
    for p, q, I, law in exceptions[:8]:
        print(f"    VIOLATION (p,q)=({p},{q}): I = {I}, law = {law}")
    # q' > 7 sample: distribution of I / theta^3
    samp = []
    for q in range(8, 31):
        for p in range(1, q):
            if gcd(p, q) == 1:
                samp.append(float(exact_I_subset([0, p, q]) / THETA ** 3))
    print(f"  q' > 7 (coprime, q <= 30): I/theta^3 in [{min(samp):.4f}, {max(samp):.4f}], "
          f"mean {sum(samp)/len(samp):.4f}  (generic ~= 1)")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART 8b -- AP CLOSED FORMS")
    print("=" * 100)
    s3_ap8 = sigma_m(list(range(8)), 3)
    closed = THETA ** 2 * sum(F(8 - d) * pillai(d) / d for d in range(2, 8))
    print(f"  Sigma_3(AP_8) engine = {s3_ap8}, Pillai closed form = {closed}, "
          f"match = {s3_ap8 == closed}")
    s3_ap13 = sigma_m(list(range(13)), 3)
    closed13_head = THETA ** 2 * sum(F(13 - d) * pillai(d) / d for d in range(2, 8))
    print(f"  Sigma_3(AP_13) engine = {s3_ap13} = {float(s3_ap13):.6f}")
    print(f"    d<=7 Pillai head = {closed13_head} = {float(closed13_head):.6f} "
          f"(reduced-q'<=7 long triples + q'>7 tail make the difference)")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART 8c -- LAYER-CAKE DOMINANCE M_m(E) <= M_m(AP_8), m = 2..7  (k=8)")
    print("=" * 100)
    def M_counts(E):
        cnt = {m: 0 for m in range(2, 8)}
        for a, b, c in combinations(sorted(E), 3):
            _, qr, _ = reduced_pattern(a, b, c)
            if qr <= 7:
                for m in range(max(2, qr), 8):
                    cnt[m] += 1
        return cnt
    ap_cnt = M_counts(range(8))
    print(f"  AP_8 cumulative M_m: {[ap_cnt[m] for m in range(2,8)]}  (m=2..7)")
    rng = random.Random(8118)
    def norm8(E):
        E = sorted(set(E)); E = [e - E[0] for e in E]
        g = 0
        for e in E[1:]:
            g = gcd(g, e)
        return [e // g if g > 1 else e for e in E]
    shapes = [
        [0,1,2,3,4,5,6,8], [0,2,4,6,8,10,11,12], [0,1,2,3,40,41,42,43],
        [0,1,2,4,8,16,32,64], [0,1,3,7,12,20,30,44], [0,2,3,7,9,13,16,17],
        [0,1,2,3,4,5,6,40], [0,3,6,9,12,15,18,22], [0,1,2,3,4,5,12,13],
        [0,1,2,4,5,6,8,9], [0,1,3,4,6,7,9,10], [0,2,4,5,8,10,12,14],
        [0,1,2,3,4,6,8,10], [0,2,4,8,10,14,16,20],
    ]
    for _ in range(4000):
        E = norm8([rng.randint(0, 40) for _ in range(8)])
        if len(E) == 8:
            shapes.append(E)
    viol_dom = []
    for E in shapes:
        c = M_counts(E)
        for m in range(2, 8):
            if c[m] > ap_cnt[m]:
                viol_dom.append((E, m, c[m], ap_cnt[m]))
                break
    print(f"  tested {len(shapes)} shapes: dominance violations = {len(viol_dom)}")
    for E, m, cm, am in viol_dom[:6]:
        print(f"    VIOLATION at m={m}: M_m = {cm} > AP {am} for {E}")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART 8d -- Sigma_4 AP-EXTREMALITY PROBE (jump descent, 80 steps)")
    print("=" * 100)
    s4_ap = sigma_m(list(range(8)), 4)
    print(f"  Sigma_4(AP_8) = {s4_ap} = {float(s4_ap):.6f}")
    best4 = (s4_ap, list(range(8)))
    cache4 = {tuple(range(8)): s4_ap}
    cur, curv = list(range(8)), s4_ap
    for step in range(80):
        E2 = list(cur)
        r = rng.random()
        if r < 0.5:
            E2[rng.randrange(8)] = rng.randint(0, 30)
        else:
            i, j = rng.randrange(8), rng.randrange(8)
            E2[i] = E2[j] + rng.choice([1, -1, 2, -2])
        E2 = norm8(E2)
        if len(E2) != 8:
            continue
        t = tuple(E2)
        if t not in cache4:
            cache4[t] = sigma_m(E2, 4)
        v = cache4[t]
        if v > best4[0]:
            best4 = (v, E2)
            print(f"  step {step}: Sigma_4 = {float(v):.6f} > AP at {E2}  <-- AP NOT the max!")
        if v > curv or rng.random() < 0.3:
            cur, curv = E2, v
    print(f"  max Sigma_4 over {len(cache4)} shapes: {float(best4[0]):.6f} at {best4[1]}"
          f"  (AP {'IS' if best4[1]==list(range(8)) else 'is NOT'} the maximizer)")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART 8e -- SECOND-MOMENT WEIGHT-2 UNIVERSALITY PROBE")
    print("  E[V^2] exact per shape vs the conjectured universal S2-truncation ceiling.")
    print("=" * 100)
    for name, E in [("AP_8", range(8)), ("Sidon", [0,1,3,7,12,20,30,44]),
                    ("geometric", [0,1,2,4,8,16,32,64]), ("two-block", [0,1,2,3,40,41,42,43])]:
        aV, MV, m3V, vmaxV = excess_moments(list(E), [THETA])
        print(f"  {name:12s} E[V] = {float(aV[0]):.6f}  E[V^2] = {float(MV[0][0]):.6f}  "
              f"E[V]^2/E[V^2] = {float(aV[0]**2/MV[0][0]):.4f}")
    print("  (if the S2 ceiling is universal, sup E[V^2] should be attained at max-Sigma_3-ish")
    print("   shapes; the exact ceiling derivation is the next analytic step)")
    sys.stdout.flush()
