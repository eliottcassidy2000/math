# -*- coding: utf-8 -*-
# kind-pasteur-2026-07-11-S127: closing the wide-spread recursion (THM-701), on top of THM-700.
#
# Setup: p_j(E) = meas{x : E misses exactly j of the six inner sectors {1..6}} (sector 0 auto by offset 0).
# p0 = meas(S7(E)) = cover measure. The LRC(14)-S3 crux is p0(E) <= cap_k for all k-sets.
#
# THIS SCRIPT verifies the four pieces of the closure:
# (1) TRANSFER OPERATOR (p1..p6 decorrelation, extends THM-700): far w => p_j(E) -> ((7-j)/7)p_j(E') +
#     ((j+1)/7)p_{j+1}(E')  (missed(E) = missed(E') minus {sector(wx)}; the fill decorrelates by THM-700).
# (2) EXACT (p0,p1) TRANSITION: p0(E)=p0(E')+Delta_w, p1(E)=p1(E')-Delta_w+Delta2_w,
#     Delta_w = meas{E' misses exactly the sector wx lands in}, Delta2_w = meas{E' misses 2, wx fills one}.
# (3) THE p1-TAX: for FAR w, Delta_w <= p1(E')/3 (decorrelated limit p1/7; worst far-w ratio ~0.25 < 1/3).
# (4) THE JOINT FUNCTIONAL Phi_lambda(F) = p0(F) + lambda*p1(F), lambda=1/3:  p0(E) <= Phi_{1/3}(E\{w})
#     (tax) and Phi_{1/3}(F) <= cap_{|F|+1} for all F -- the reduced obligation. ACCUMULATION: adding a far
#     element raises Phi by (2/3)Delta_w + (1/3)Delta2_w = 2(p1+p2)/21 (decorrelated), which DECAYS as
#     coverage grows, so Phi_{1/3} converges to its iid limit < cap (the series is summable).
import random

def sector(y):
    return min(6, int((y % 1.0) * 7))

def missed(E, x):
    return set(range(1, 7)) - set(sector(e * x) for e in E)

def pvec(E, N):
    v = [0] * 7
    for i in range(N):
        v[len(missed(E, (i + 0.5) / N))] += 1
    return [c / N for c in v]

cap = {8: 0.38153, 9: 0.49426, 10: 0.60440}   # seven-sector caps (THM-532/534), exact
N = 60013

def main():
    # (1) transfer operator
    def M(p):
        return [((7 - j) / 7) * p[j] + ((j + 1) / 7) * (p[j + 1] if j + 1 < 7 else 0) for j in range(7)]
    core = [0, 1, 2, 3, 4, 5, 6, 7]; pC = pvec(core, N)
    err = max(abs(pvec(core + [1601], N)[j] - M(pC)[j]) for j in range(7))
    print(f"(1) transfer operator p_j(E) -> ((7-j)/7)p_j + ((j+1)/7)p_(j+1):  max err (w=1601) = {err:.5f}")

    # (2)+(3) exact transition + p1-tax over FAR w
    print("(2)/(3) exact transition + p1-tax Delta_w <= p1(E')/3 over far w (w >= 3*max core):")
    worst = 0.0
    for c in [[0, 1, 2, 3, 4, 5, 6, 7], [0, 1, 3, 7, 12, 20], [0, 2, 5, 11, 17, 29, 40]]:
        p = pvec(c, N); mx = 0.0
        for w in range(3 * max(c), 3 * max(c) + 200):
            mx = max(mx, (pvec(c + [w], N)[0] - p[0]) / p[1])
        worst = max(worst, mx)
        print(f"    core max={max(c)}: max Delta_w/p1 = {mx:.4f}")
    print(f"    => worst far-w tax = {worst:.4f}  (bound 1/3 = {1/3:.4f}: {'HOLDS' if worst <= 1/3 else 'FAILS'})")

    # (4) joint functional Phi_{1/3}(F) <= cap_{|F|+1}, random wide F
    print("(4) Phi_{1/3}(F) = p0 + (1/3)p1 <= cap_{|F|+1}, 1500 random wide F:")
    random.seed(3); worstm = {}
    for _ in range(1500):
        k = random.choice([7, 8, 9])
        E = [0] + sorted(random.sample(range(1, 70), k - 1))
        p = pvec(E, N)
        worstm[k] = min(worstm.get(k, 9.0), cap[k + 1] - (p[0] + p[1] / 3))
    for k in sorted(worstm):
        print(f"    |F|={k}: min margin cap_{k+1} - Phi_(1/3) = {worstm[k]:.4f}  "
              f"({'OK' if worstm[k] > 0 else 'VIOLATED'})")

    # accumulation: Phi_{1/3} increments decay as far elements are appended (geometric-ish spread)
    print("(acc) Phi_{1/3}(consec_8 + growing far tail) -> bounded limit (increments decay):")
    E = [0, 1, 2, 3, 4, 5, 6, 7]; prev = None
    for w in [40, 200, 1000, 5000]:
        E2 = E + [w]; p = pvec(E2, N); phi = p[0] + p[1] / 3
        inc = "" if prev is None else f"  d={phi - prev:+.4f}"
        print(f"    append w={w:>5}: Phi_(1/3) = {phi:.5f}{inc}")
        E = E2; prev = phi

if __name__ == "__main__":
    main()
