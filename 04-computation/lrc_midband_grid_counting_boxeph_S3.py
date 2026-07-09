#!/usr/bin/env python3
"""
HYP-5760: mid-band grid-counting survival (klein-S208 attack axis (b)).
boxeph-2026-07-09-S3.  Integer arithmetic + exact Fractions at the end.

Setup: covering 13-set S = P~ u L~, V = max S, cluster co-offsets
E = {V - u : u in L~} (0 in E, spread <= 5V/14 after the adaptive split),
killers P~ = S \\ L~ (slow + rider + MID-BAND members), |P~| = 13 - k~ <= 5.

Pipeline (all exact):
  1. Good = {j in [0,V) : 7*maxCircGap(e*j mod V) > V}  (strict good period).
  2. For j in Good: phi_j = center of the widest gap; tau_j = (j + phi_j)/V.
  3. TRUE kill per killer m: {j in Good : ||m tau_j|| < 1/14}.
  4. COUNTING BOUND per killer: V/7 + gcd(m, V)  (fixed-phi subgroup count;
     the true kill is measured against it -- the phi_j coupling is the gap).
  5. Survivors realize lonely tau: exact ||v tau_j|| >= 1/14 for ALL 13 v.

The lemma's overlap-agnostic survival floor:
  survivors >= |Good| - Sum_m kill_m >= V[mu_grid - |P~|/7 - Sum gcd(m,V)/V - c]
positive by the ladder margins when all gcd(m,V) small. Zero-survivor
instances should be exactly the large-gcd (harmonic/coarse) family.
"""
from fractions import Fraction as F
from math import gcd

ONE14 = F(1, 14)


def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))


def is_primitive(S):
    g = 0
    for v in S:
        g = gcd(g, v)
    return g == 1


def good_periods(E, V):
    """REALIZABLE good j: widest gap admits a drift-valid phi (monad-S2
    phi-interval, verified per-j by exact full-cluster clearance at
    tau = (j + phi*)/V). Returns {j: (num, den) of tau}."""
    s = max(E)
    out = {}
    for j in range(V):
        teeth = sorted(set((e * j) % V for e in E))
        n = len(teeth)
        best, ustart = -1, 0
        for i in range(n):
            a = teeth[i]
            b = teeth[(i + 1) % n] + (V if i == n - 1 else 0)
            if b - a > best:
                best, ustart = b - a, a
        if 7 * best <= V:
            continue
        # monad phi-interval: a=ustart/V, g=best/V, r=s/V ->
        # phi in ((a+1/14)/(1-r), (a+g-1/14)/(1+r)); take center, verify.
        a, g, r = F(ustart, V), F(best, V), F(s, V)
        lo = (a + ONE14) / (1 - r)
        hi = (a + g - ONE14) / (1 + r)
        if lo >= hi:
            continue
        phi = (lo + hi) / 2
        tau = (j + phi) / V
        # exact cluster + observer clearance
        ok = True
        for e in E:
            v = V - e
            x = (v * tau) % 1
            if min(x, 1 - x) < ONE14:
                ok = False
                break
        if ok:
            out[j] = tau
    return out


def analyze(name, P_killers, E, V, verbose_survivor=True):
    S = sorted(P_killers + [V - e for e in E])
    assert len(S) == 13, f"{name}: {len(S)} speeds"
    cov, prim = is_covering(S), is_primitive(S)
    k = len(E)
    s = max(E)
    mid_lo, mid_hi = V / 14, 9 * V / 14
    zones = {v: ('P' if v <= 13 else 'T' if v <= V / 14 else
                 'M' if v < 9 * V / 14 else 'L') for v in S}
    killers = sorted(P_killers)
    good = good_periods(E, V)
    nGood = len(good)
    print(f"\n=== {name} ===")
    print(f"  S={S}")
    print(f"  V={V} k~={k} spread(E)={s} (5V/14={5*V/14:.0f})  "
          f"covering={cov} primitive={prim}")
    print(f"  killers P~={killers} zones={[zones[m] for m in killers]}  "
          f"mid-band=({mid_lo:.0f},{mid_hi:.0f})")
    print(f"  |Good|={nGood} ({nGood/V:.3f} of V; need > sum kills)")
    # per-killer true kill among Good under gap-centered phi
    kill_sets = {}
    for m in killers:
        ks = set()
        for j, tau in good.items():
            x = (m * tau) % 1
            if min(x, 1 - x) < ONE14:
                ks.add(j)
        kill_sets[m] = ks
        g = gcd(m, V)
        bound = V / 7 + g
        print(f"    killer {m:6d} (zone {zones[m]}, gcd={g:4d}): "
              f"true kill among Good = {len(ks):5d}  "
              f"[counting bound V/7+g = {bound:7.1f}]  "
              f"{'<= bound OK' if len(ks) <= bound else 'EXCEEDS fixed-phi bound (phi-coupling)'}")
    union_kill = set().union(*kill_sets.values()) if kill_sets else set()
    survivors = [j for j in good if j not in union_kill]
    floor = nGood - sum(len(k_) for k_ in kill_sets.values())
    print(f"  union kill={len(union_kill)}  survivors={len(survivors)}  "
          f"(overlap-agnostic floor |Good|-Sum kills = {floor})")
    if survivors:
        # realize + verify one survivor exactly (ALL 13 speeds)
        j = survivors[len(survivors) // 2]
        tau = good[j]
        clr = min(min((v * tau) % 1, 1 - (v * tau) % 1) for v in S)
        ok = clr >= ONE14
        print(f"  REALIZED: j={j} tau={tau}  min clearance={clr} "
              f"= {float(clr):.6f}  {'>= 1/14 LONELY OK' if ok else 'FAIL'}")
        return len(survivors), ok
    else:
        gs = [gcd(m, V) for m in killers]
        print(f"  ZERO SURVIVORS -- gcds {gs}; max gcd/V = {max(gs)/V:.3f} "
              f"({'LARGE-gcd/coarse family as predicted' if max(gs) > V/20 else 'SMALL-gcd zero-survivor: LEMMA COUNTEREXAMPLE, investigate!'})")
        return 0, None


def main():
    print("HYP-5760: mid-band grid-counting survival -- exact pipeline")

    # --- bank 1: k~=10 cluster, 3 killers incl 2 mid-band, generic V ---
    V = 842                                     # 2*421
    E = [0, 11, 37, 68, 105, 133, 160, 191, 224, 260]  # k=10, spread 260 < 5V/14=300
    # covering duty: V covers 2; need 3..14 from cluster or killers
    # V-e residues: pick killers to patch. P~: {12, 13} P-zone + mid m
    # find mid-band m and adjust E residues for covering
    P_ = [12, 13, 350]                           # 350 in (60.1, 541.3) mid-band
    S = sorted(P_ + [V - e for e in E])
    # patch covering by scanning small changes to the mid member
    if not is_covering(S):
        for m in range(300, 500):
            S2 = sorted([12, 13, m] + [V - e for e in E])
            if is_covering(S2) and is_primitive(S2):
                P_ = [12, 13, m]
                break
    analyze("k~=10, 2P+1M generic", P_, E, V)

    # --- bank 2: k~=10, 3 mid-band killers, near-resonant m/V ~ 1/2, 1/3 ---
    V = 840 + 2                                 # keep 842
    for trial in [[V // 2 - 1, V // 3 + 1, V // 5 + 2]]:
        M3 = trial
        E2 = [0, 11, 37, 68, 105, 133, 160, 191, 224, 260]
        S2 = sorted(M3 + [V - e for e in E2])
        if not (is_covering(S2) and is_primitive(S2)):
            # scan neighborhood for covering
            import itertools as it
            found = False
            for d1, d2, d3 in it.product(range(0, 8), repeat=3):
                M3b = [V // 2 - 1 + d1, V // 3 + 1 + d2, V // 5 + 2 + d3]
                if len(set(M3b)) < 3:
                    continue
                S2 = sorted(M3b + [V - e for e in E2])
                if len(set(S2)) == 13 and is_covering(S2) and is_primitive(S2):
                    M3 = M3b
                    found = True
                    break
            if not found:
                print("\n(bank 2: no covering variant found in scan)")
                continue
        analyze("k~=10, 3M near-resonant", M3, E2, V)

    # --- bank 3: HARMONIC killers (large gcd) -- expect zero/low survivors ---
    V = 840                                     # highly divisible
    E3 = [0, 11, 37, 68, 105, 133, 160, 191, 224, 260]
    M3 = [420, 280, 168]                        # V/2, V/3, V/5 exact harmonics
    S3 = sorted(M3 + [V - e for e in E3])
    if is_covering(S3):
        analyze("k~=10, 3M EXACT harmonics", M3, E3, V)
    else:
        print("\n(bank 3: exact-harmonic set not covering -- adjusting)")
        M3 = [420, 280, 170]
        S3 = sorted(M3 + [V - e for e in E3])
        if is_covering(S3):
            analyze("k~=10, harmonics adj", M3, E3, V)

    # --- bank 4: k~=8, the full 5-killer budget (hardest counting case) ---
    V = 1006                                    # 2*503
    E4 = [0, 13, 41, 78, 120, 167, 219, 276]    # k=8, spread 276 < 5V/14=359
    # 5 killers: 2 slow + 3 mid-band; scan for covering
    import itertools as it
    base = [11, 13]
    found = None
    for m1 in range(120, 400, 7):
        for m2 in range(m1 + 50, 500, 11):
            for m3 in range(m2 + 40, 600, 13):
                P5 = base + [m1, m2, m3]
                S4 = sorted(P5 + [V - e for e in E4])
                if len(set(S4)) == 13 and is_covering(S4) and is_primitive(S4):
                    found = P5
                    break
            if found: break
        if found: break
    if found:
        analyze("k~=8, 2P+3M (5-killer budget)", found, E4, V)
    else:
        print("\n(bank 4: no covering 5-killer found in scan)")

    # --- bank 5: adversarial -- mid member tuned AGAINST the good comb ---
    # pick m = multiple of V/gcd structure aligned with cluster resonances
    V = 910                                     # 2*5*7*13 -- rich divisors
    E5 = [0, 13, 41, 78, 120, 167, 219, 276]
    M5 = [455, 364, 182]                        # V/2, 2V/5, V/5: all large gcd
    S5 = sorted(M5 + [V - e for e in E5])
    if len(set(S5)) == 13 and is_covering(S5):
        analyze("k~=8, 3M large-gcd adversarial", M5, E5, V)


if __name__ == '__main__':
    main()
