"""
lrc14_pz_general_integrator_opus_S148.py   (opus-2026-07-08-S148, HYP-5327)

GENERAL exact Paley-Zygmund covering-floor integrator (generalizes mac-mini-S57's
block-only lrc14_block_PZ_exact to ARBITRARY integer families E).

W(x) = sum_i (g_i(x) - 1/7)_+, g_i = circular gaps of {frac(e x): e in E};
mu_{1/7}(E) = P(W>0) >= E[W]^2/E[W^2] = PZ (THM-660).  This computes E[W], E[W^2], PZ as
EXACT rationals for any E, via Farey-cell piecewise-linear integration:
  - phase order of {frac(ex)} + each floor(ex) are constant between breakpoints
    x = p/d, d in {e_i} u {e_i - e_j} (all elements and positive differences), p integer;
  - on each cell frac(ex) = e x - floor(e x_mid) is linear; each circular gap g = a+bx;
  - split at g = 1/7 crossings; W is linear per sub-cell; integrate W, W^2 exactly.

USE: (a) validate vs mac-mini block values; (b) exhaustive compact-shape min PZ for
k=12,13 (the [compact check] half klein did only for k=11); (c) PZ vs diameter -> the
decorrelation-tail threshold D_k where PZ >= bar is guaranteed.
"""
from fractions import Fraction as F
from math import floor, gcd
from itertools import combinations

TH = F(1, 7)

# exact honest bars (MISTAKE-123): bar_k = m_P + 1 - min_{|P|=13-k} meas(G_P)
# k=11/12/13 (from THM-660 / mac-mini exact):
BAR = {8: F(1702763, 2522520), 9: F(35456, 63063), 10: F(114041, 252252),
       11: F(83549, 252252), 12: F(50285, 252252), 13: F(14249, 252252)}


def pz_exact(E):
    """Exact (E[W], E[W^2], PZ) for integer family E (list). Translation-invariant; we use
       the raw e-values (phases frac(ex))."""
    E = sorted(E)
    k = len(E)
    e0 = E[0]
    # normalize to start at 0 (translation of speeds shifts all phases equally -> same gaps)
    Es = [e - e0 for e in E]  # e0-based; frac((e-e0)x) is a rigid rotation of frac(ex)-frac(e0 x)...
    # NOTE: gaps of {frac(e x)} = gaps of {frac((e-e0)x)} (rigid rotation by frac(e0 x)),
    # so W is identical. Use Es (starts at 0) -> denominators are the differences.
    dens = set()
    for e in Es:
        if e > 0:
            dens.add(e)
    for i in range(k):
        for j in range(i + 1, k):
            dens.add(Es[j] - Es[i])
    bps = set([F(0), F(1)])
    for d in dens:
        for p in range(0, d + 1):
            bps.add(F(p, d))
    bps = sorted(bps)
    EW = F(0); EW2 = F(0)
    for c in range(len(bps) - 1):
        x0, x1 = bps[c], bps[c + 1]
        xm = (x0 + x1) / 2
        lin = []
        for e in Es:
            ce = floor(e * xm)
            lin.append((F(-ce), F(e)))   # frac(ex) = e x - ce = a + b x
        order = sorted(range(k), key=lambda j: lin[j][0] + lin[j][1] * xm)
        sp = [lin[j] for j in order]
        gaps = []
        for i in range(k - 1):
            gaps.append((sp[i + 1][0] - sp[i][0], sp[i + 1][1] - sp[i][1]))
        gaps.append((F(1) + sp[0][0] - sp[k - 1][0], sp[0][1] - sp[k - 1][1]))
        subs = set([x0, x1])
        for (a, b) in gaps:
            if b != 0:
                xs = (TH - a) / b
                if x0 < xs < x1:
                    subs.add(xs)
        subs = sorted(subs)
        for s in range(len(subs) - 1):
            u0, u1 = subs[s], subs[s + 1]
            um = (u0 + u1) / 2
            A = F(0); B = F(0)
            for (a, b) in gaps:
                if a + b * um > TH:
                    A += (a - TH); B += b
            EW += A * (u1 - u0) + B * (u1 * u1 - u0 * u0) / 2
            EW2 += (A * A * (u1 - u0) + A * B * (u1 * u1 - u0 * u0)
                    + B * B * (u1**3 - u0**3) / 3)
    PZ = EW * EW / EW2 if EW2 != 0 else F(0)
    return EW, EW2, PZ


def additive_energy_reduced(E):
    """R2 = sum_{d != 0} r_d^2, r_d = #{(i,j): e_i - e_j = d} (ordered)."""
    from collections import Counter
    c = Counter()
    for i in range(len(E)):
        for j in range(len(E)):
            if i != j:
                c[E[i] - E[j]] += 1
    return sum(v * v for v in c.values())


def main():
    print("=" * 92)
    print("(0) VALIDATION: general integrator vs mac-mini block-only exact values")
    print("=" * 92)
    macmini = {11: F(3400663, 9797402), 12: F(82210303, 266892098), 13: F(221828403, 815409784)}
    for k in (11, 12, 13):
        EW, EW2, PZ = pz_exact(list(range(k)))
        ok = (PZ == macmini[k])
        print(f"  k={k} block: PZ = {PZ} = {float(PZ):.6f}  vs mac-mini {float(macmini[k]):.6f}"
              f"  {'MATCH' if ok else '*** DIFF ***'}   bar {float(BAR[k]):.5f}"
              f"  margin {float(PZ-BAR[k]):+.5f}")

    print()
    print("=" * 92)
    print("(1) PZ vs DIAMETER (the decorrelation tail): min PZ over shapes at each diameter band")
    print("=" * 92)
    import random
    for k in (11, 12, 13):
        print(f"  k={k} (bar {float(BAR[k]):.4f}):")
        for D in (k - 1, 15, 25, 40, 70, 120):
            rng = random.Random(1000 * k + D)
            mn = None; arg = None
            trials = 400 if D > k - 1 else 1
            for _ in range(trials):
                if D == k - 1:
                    E = list(range(k))
                else:
                    S = sorted(rng.sample(range(1, D + 1), k - 1))
                    E = [0] + S
                    if E[-1] != D:  # ensure diameter exactly D sometimes
                        E[-1] = D
                    E = sorted(set(E))
                    if len(E) != k:
                        continue
                _, _, PZ = pz_exact(E)
                if mn is None or PZ < mn:
                    mn = PZ; arg = E
            tag = "block" if D == k - 1 else f"min/{trials}"
            print(f"    diam {D:4d}: min PZ = {float(mn):.4f} ({tag})  "
                  f"{'>= bar' if mn >= BAR[k] else '*** < bar ***'}  argmin {arg}")

    print()
    print("=" * 92)
    print("(2) ADDITIVE-ENERGY ORDERING: corr(R2, PZ) sign; does higher energy => lower PZ?")
    print("=" * 92)
    for k in (11, 13):
        rng = random.Random(7 * k)
        pts = []
        for _ in range(120):
            D = rng.randint(k - 1, 40)
            E = sorted(set([0] + rng.sample(range(1, D + 1), k - 1)))
            if len(E) != k:
                continue
            _, _, PZ = pz_exact(E)
            R2 = additive_energy_reduced(E)
            pts.append((R2, float(PZ)))
        # correlation
        n = len(pts)
        mr = sum(p[0] for p in pts) / n; mp = sum(p[1] for p in pts) / n
        cov = sum((p[0] - mr) * (p[1] - mp) for p in pts) / n
        sr = (sum((p[0] - mr)**2 for p in pts) / n)**.5
        spz = (sum((p[1] - mp)**2 for p in pts) / n)**.5
        corr = cov / (sr * spz) if sr * spz else 0
        blkR2 = additive_energy_reduced(list(range(k)))
        maxR2 = max(p[0] for p in pts)
        print(f"  k={k}: corr(R2, PZ) = {corr:+.3f} over {n} shapes; block R2 = {blkR2} "
              f"(= max? {blkR2 >= maxR2}); higher energy -> lower PZ: {corr < -0.5}")


if __name__ == "__main__":
    main()
