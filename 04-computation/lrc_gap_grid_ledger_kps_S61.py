#!/usr/bin/env python3
r"""
lrc_gap_grid_ledger_kps_S61.py   (kind-pasteur-2026-07-07-S61, HYP-4857 part 2)

THE GAP-GRID SUPERSET LEDGER for the spread residual (mac-mini-S41 handoff (b),
extending the kps-S59 subset lemma to rank-2 supersets).

FRAME.  If E ⊆ a + {i*d1 + j*d2 : 0<=i<n1, 0<=j<n2} =: a + G(n1,n2;d1,d2) (a rank-2
GAP cover), the subset lemma gives mu_{1/7}(E) >= mu_{1/7}(grid) pointwise-monotone.
The grid config {frac((i d1 + j d2)x)} = {i*u + j*v mod 1} with (u,v) = (d1 x, d2 x)
sweeping the (d1,d2)-geodesic of T^2.  For coprime d1,d2 the geodesic equidistributes:
    mu_{1/7}(grid(d1,d2)) --> mu2(n1,n2) := meas_{T^2}{(u,v) : maxgap{iu+jv} > 1/7},
a STEP-INDEPENDENT 2-torus constant; the error is a geodesic-discrepancy term
(Fourier modes (m1,m2) with m1 d1 + m2 d2 = 0 -- the same mechanism as the S59
pair bound |joint - product| <= 1/(3 (d1/g)(d2/g))).

So the spread residual's STRUCTURED half reduces to a FINITE table mu2(n1,n2) plus
a step-error bound: any 13-set covered by a grid whose (n1,n2) row clears m_P (with
error margin) inherits the floor -- regardless of how spread its diameter is.

This script:
  PART 1: step-independence, measured -- mu(grid) across many coprime (d1,d2) at
          fixed (n1,n2): spread of values vs 1/(d1*d2) scaling.
  PART 2: the mu2(n1,n2) table (numeric 2-torus MC + fine deterministic grid) for
          n1*n2 <= 80, n2 <= 6 -- vs m_P and vs the k=13 sharp bar 477/1078.
  PART 3: coverability -- which known/structured spread families admit small GAP
          covers (record family = (11,2)-grid subset; interlacings; two-block).
  PART 4: honest error anatomy -- worst |mu(grid) - mu2| vs the pair-bound scale.
"""
import math, random

TH = 1 / 7
MP = 14249 / 252252

def maxgap_pts(pts):
    ph = sorted(pts)
    g = max(ph[i + 1] - ph[i] for i in range(len(ph) - 1)) if len(ph) > 1 else 1.0
    return max(g, ph[0] + 1.0 - ph[-1])

def grid_config(u, v, n1, n2):
    return [ (i * u + j * v) % 1.0 for i in range(n1) for j in range(n2) ]

def mu_grid_1d(d1, d2, n1, n2, res=30000):
    cnt = 0
    for r in range(res):
        x = (r + 0.5) / res
        if maxgap_pts(grid_config((d1 * x) % 1.0, (d2 * x) % 1.0, n1, n2)) > TH:
            cnt += 1
    return cnt / res

def mu2_torus(n1, n2, res=700):
    """deterministic 2-torus grid average (res^2 cells, midpoints)."""
    cnt = 0
    for a in range(res):
        u = (a + 0.5) / res
        for b in range(res):
            v = (b + 0.5) / res
            if maxgap_pts(grid_config(u, v, n1, n2)) > TH:
                cnt += 1
    return cnt / (res * res)

print("=" * 96)
print("PART 1 -- step-independence of mu_{1/7}(grid(n1,n2; d1,d2)), measured")
print("=" * 96)
print(f"  m_P = {MP:.6f}")
for (n1, n2) in [(2, 11), (7, 2), (5, 3)]:
    vals = []
    for (d1, d2) in [(11, 2), (13, 2), (23, 3), (17, 5), (29, 4), (37, 8), (101, 13), (211, 30)]:
        if math.gcd(d1, d2) != 1: continue
        m = mu_grid_1d(d1, d2, n1, n2, 20000)
        vals.append((d1, d2, m))
    lo = min(v[2] for v in vals); hi = max(v[2] for v in vals)
    print(f"  grid {n1}x{n2}: mu over steps = [{lo:.4f}, {hi:.4f}]  spread {hi-lo:.4f}")
    for d1, d2, m in vals[:4]:
        print(f"      (d1,d2)=({d1},{d2}): mu={m:.4f}")

print()
print("=" * 96)
print("PART 2 -- the 2-torus table mu2(n1,n2) (step-independent limit), n2<=6")
print("=" * 96)
print(f"  {'n1 x n2':>9} {'N=n1*n2':>8} {'mu2':>8}  vs m_P")
tab = {}
sel = {2: [7, 9, 11, 14, 18, 24, 31, 40], 3: [5, 7, 9, 12, 16, 21, 26],
       4: [4, 5, 7, 9, 12, 16, 20], 5: [3, 4, 6, 8, 11, 14, 16], 6: [3, 4, 6, 8, 10, 13]}
for n2, n1s in sel.items():
    for n1 in n1s:
        N = n1 * n2
        if N < 13 or N > 80: continue
        m = mu2_torus(n1, n2, 240)
        tab[(n1, n2)] = m
        flag = "OK" if m >= MP else "*** below m_P ***"
        print(f"  {n1:4d} x {n2:1d} {N:8d} {m:8.4f}  {flag}")
# find per-n2 crossing
print("\n  per-n2 crossing (largest n1 with mu2 >= m_P):")
for n2 in range(2, 7):
    ok = [n1 for (a, b), m in tab.items() if b == n2 and m >= MP for n1 in [a]]
    if ok:
        print(f"    n2={n2}: n1 up to {max(ok)} (N up to {max(ok)*n2})")

print()
print("=" * 96)
print("PART 3 -- GAP-coverability of spread/record families")
print("=" * 96)
def best_grid_cover(E, maxN=80):
    """smallest-N rank-2 grid cover a + {i d1 + j d2} found by direct search."""
    E = sorted(E); best = None
    D = E[-1] - E[0]
    E0 = [e - E[0] for e in E]
    for d1 in range(1, D + 1):
        for d2 in range(1, d1 + 1):
            if math.gcd(d1, d2) > 1: continue
            # greedy: can E0 be written as i*d1 + j*d2, 0<=i<n1, 0<=j<n2 with small n1,n2?
            for n2 in range(1, 7):
                need_i = set()
                ok = True
                for e in E0:
                    found = False
                    for j in range(n2):
                        r = e - j * d2
                        if r >= 0 and r % d1 == 0:
                            need_i.add(r // d1); found = True; break
                    if not found: ok = False; break
                if ok and need_i:
                    n1 = max(need_i) + 1
                    N = n1 * n2
                    if N <= maxN and (best is None or N < best[0]):
                        best = (N, n1, n2, d1, d2)
    return best

fams = {
    "monad record {2..22e}u{11,13}": [2,4,6,8,10,11,12,13,14,16,18,20,22],
    "parity diam 80: 8*{0..10}u{39,41}": sorted([8*i for i in range(11)]+[39,41]),
    "interlace diam 100: 10*{0..10}u{49,51}": sorted([10*i for i in range(11)]+[49,51]),
    "stretch {0,2..12,17,28}": [0]+list(range(2,13))+[17,28],
    "two-block {0..5}u{70..76}": list(range(6))+list(range(70,77)),
    "3-adic cascade (mac-mini EU-min)": [0,30,36,45,50,54,60,63,70,72,81,90,108],
}
for nm, E in fams.items():
    bc = best_grid_cover(E)
    if bc:
        N, n1, n2, d1, d2 = bc
        m2 = tab.get((n1, n2))
        m2txt = f"mu2({n1},{n2})={m2:.4f}" if m2 else f"mu2({n1},{n2}) not in table"
        clr = "clears m_P" if m2 and m2 >= MP else ("?" if not m2 else "BELOW")
        print(f"  {nm:40s} cover: {n1}x{n2} steps ({d1},{d2}) N={N}  {m2txt}  {clr}")
    else:
        print(f"  {nm:40s} NO small grid cover (<=80 pts, n2<=6) -- decorrelation lane")

print()
print("=" * 96)
print("PART 4 -- error anatomy: |mu(grid;d1,d2) - mu2| vs the 1/(d1 d2) pair-bound scale")
print("=" * 96)
for (n1, n2) in [(2, 11), (7, 2)]:
    m2 = mu2_torus(n1, n2, 500)
    worst = 0; rows = []
    for (d1, d2) in [(3, 2), (5, 2), (7, 3), (11, 2), (17, 5), (29, 4), (101, 13)]:
        if math.gcd(d1, d2) != 1: continue
        m = mu_grid_1d(d1, d2, n1, n2, 20000)
        err = abs(m - m2); sc = err * d1 * d2
        rows.append((d1, d2, m, err, sc))
        worst = max(worst, sc)
    print(f"  grid {n1}x{n2} (mu2={m2:.4f}):")
    for d1, d2, m, err, sc in rows:
        print(f"    ({d1:3d},{d2:2d}): mu={m:.4f} err={err:.4f}  err*(d1*d2)={sc:.2f}")
    print(f"    => err ~ C/(d1*d2) with C <= {worst:.1f} (empirical; Fourier-provable shape)")
print()
print("DONE.")
