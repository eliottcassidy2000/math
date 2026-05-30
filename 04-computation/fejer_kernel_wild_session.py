#!/usr/bin/env python3
"""
fejer_kernel_wild_session.py

Loose but reproducible Fejer-kernel probes for circulant tournaments on Z_p.

The target is the gap left by the older S64 thread:

    additive energy / Fejer spectrum  ->  spatial cycle localization

This script keeps the experiments small enough to rerun quickly.  It checks
the exact Fejer/autocorrelation identities, then compares them with pinned
cycle profiles J_k(0,v), Hamiltonian path counts, and single-pair mutations.
"""

from __future__ import annotations

import cmath
import math
from collections import defaultdict
from itertools import combinations


def valid_sets(p: int):
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]
    for bits in range(1 << m):
        yield tuple(pairs[i][0] if bits & (1 << i) else pairs[i][1] for i in range(m))


def interval_set(p: int):
    return tuple(range(1, (p - 1) // 2 + 1))


def paley_set(p: int):
    return tuple(sorted(a for a in range(1, p) if pow(a, (p - 1) // 2, p) == 1))


def adjacency(p: int, S):
    S = set(S)
    A = [[0] * p for _ in range(p)]
    for i in range(p):
        for s in S:
            A[i][(i + s) % p] = 1
    return A


def eigenvalues(p: int, S):
    omega = cmath.exp(2j * cmath.pi / p)
    return [sum(omega ** (k * s) for s in S) for k in range(p)]


def q_spectrum(p: int, S):
    lam = eigenvalues(p, S)
    return [abs(lam[k]) ** 2 for k in range(1, p)]


def autocorrelation(p: int, S):
    S = set(S)
    return [sum(1 for a in S if (a + d) % p in S) for d in range(p)]


def additive_energy(p: int, S):
    r = autocorrelation(p, S)
    return sum(x * x for x in r)


def fejer_values(p: int):
    m = (p - 1) // 2
    vals = []
    for k in range(1, p):
        vals.append((math.sin(math.pi * m * k / p) / math.sin(math.pi * k / p)) ** 2)
    return vals


def entropy(vals):
    total = sum(vals)
    probs = [v / total for v in vals if v > 0]
    return -sum(q * math.log(q) for q in probs)


def concentration_metrics(p: int, S):
    q = q_spectrum(p, S)
    total = sum(q)
    ipr = sum(x * x for x in q) / (total * total)
    top = max(q) / total
    pr = 1.0 / ipr
    return {
        "top": top,
        "ipr": ipr,
        "pr": pr,
        "entropy": entropy(q),
        "energy": additive_energy(p, S),
    }


def hamiltonian_paths(p: int, S):
    A = adjacency(p, S)
    full = (1 << p) - 1
    dp = [[0] * p for _ in range(1 << p)]
    for v in range(p):
        dp[1 << v][v] = 1
    for mask in range(1 << p):
        for v in range(p):
            cur = dp[mask][v]
            if not cur:
                continue
            row = A[v]
            rem = full ^ mask
            w = 0
            while rem:
                if rem & 1 and row[w]:
                    dp[mask | (1 << w)][w] += cur
                rem >>= 1
                w += 1
    return sum(dp[full])


def pinned_cycles_profile(p: int, S, k: int):
    """J_k(0,v): directed simple k-cycles through both 0 and v."""
    A = adjacency(p, S)
    profile = [0] * p
    vertices = list(range(1, p))

    for others in combinations(vertices, k - 1):
        verts = (0,) + others
        idx = {v: i for i, v in enumerate(verts)}
        n = len(verts)
        dp = {(1, 0): 1}
        for _ in range(1, k):
            ndp = {}
            for (mask, pos), count in dp.items():
                u = verts[pos]
                for nxt_pos, w in enumerate(verts):
                    if mask & (1 << nxt_pos):
                        continue
                    if A[u][w]:
                        key = (mask | (1 << nxt_pos), nxt_pos)
                        ndp[key] = ndp.get(key, 0) + count
            dp = ndp

        full = (1 << n) - 1
        cycle_count = 0
        for (mask, pos), count in dp.items():
            if mask == full and A[verts[pos]][0]:
                cycle_count += count
        if cycle_count:
            for v in verts:
                profile[v] += cycle_count

    return profile


def profile_stats(profile):
    vals = profile[1:]
    mean = sum(vals) / len(vals)
    var = sum((x - mean) ** 2 for x in vals) / len(vals)
    return mean, var, max(vals) - min(vals)


def count_3cycles_through_zero_by_oriented_wedge(p: int, S):
    """A sanity expression for J_3: count cyclic oriented wedges explicitly."""
    S = set(S)
    out0 = S
    in0 = {(-s) % p for s in S}
    out_of = lambda v: {(v + s) % p for s in S}
    profile = [0] * p
    for v in range(1, p):
        if v in out0:
            profile[v] = len((out_of(v) & in0) - {0, v})
        else:
            # Cyclic orientation must be 0 <- v <- w <- 0.
            inv = {(v - s) % p for s in S}
            profile[v] = len((inv & out0) - {0, v})
    profile[0] = sum(profile[1:]) // 2
    return profile


def pearson(xs, ys):
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return float("nan")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)


def single_flip_neighbors(p: int, S):
    S = set(S)
    for s in range(1, (p - 1) // 2 + 1):
        T = set(S)
        if s in T:
            T.remove(s)
            T.add((-s) % p)
            label = f"{s}->{p-s}"
        else:
            T.remove((-s) % p)
            T.add(s)
            label = f"{p-s}->{s}"
        yield label, tuple(sorted(T))


def summarize_set(p: int, name: str, S, include_h: bool = True):
    metrics = concentration_metrics(p, S)
    j3 = pinned_cycles_profile(p, S, 3)
    j3_wedge = count_3cycles_through_zero_by_oriented_wedge(p, S)
    mean3, var3, range3 = profile_stats(j3)
    line = (
        f"{name:10s} top={metrics['top']:.6f} ipr={metrics['ipr']:.6f} "
        f"PR={metrics['pr']:.3f} entropy={metrics['entropy']:.3f} "
        f"E={metrics['energy']:4d} J3var={var3:8.3f} J3range={range3:5d}"
    )
    if include_h:
        line += f" H={hamiltonian_paths(p, S)}"
    print(line)
    if j3 != j3_wedge:
        print(f"  WARNING: wedge formula mismatch: exact={j3}, wedge={j3_wedge}")


def part_identity_checks():
    print("=" * 78)
    print("I. FEJER = FOURIER MAGNITUDE = AUTOCORRELATION TRANSFORM")
    print("=" * 78)
    for p in [7, 11, 13, 17, 19, 23, 29, 31]:
        S = interval_set(p)
        q = q_spectrum(p, S)
        f = fejer_values(p)
        ac = autocorrelation(p, S)
        omega = cmath.exp(2j * cmath.pi / p)
        ac_hat = []
        for k in range(1, p):
            ac_hat.append(sum(ac[d] * omega ** (k * d) for d in range(p)).real)
        err_fejer = max(abs(a - b) for a, b in zip(q, f))
        err_auto = max(abs(a - b) for a, b in zip(q, ac_hat))
        m = (p - 1) // 2
        print(
            f"p={p:2d} m={m:2d} max|Q-Fejer|={err_fejer:.2e} "
            f"max|Q-hat(auto)|={err_auto:.2e} "
            f"top={max(q)/sum(q):.6f}"
        )


def part_landscape():
    print("\n" + "=" * 78)
    print("II. SMALL PRIME LANDSCAPE: SPECTRAL CONCENTRATION VS CYCLE LOCALIZATION VS H")
    print("=" * 78)
    for p in [7, 11, 13]:
        rows = []
        for S in valid_sets(p):
            metrics = concentration_metrics(p, S)
            j3 = pinned_cycles_profile(p, S, 3)
            _, var3, range3 = profile_stats(j3)
            H = hamiltonian_paths(p, S)
            rows.append((metrics["ipr"], metrics["energy"], var3, range3, H, S))
        print(f"\np={p}: {len(rows)} circulant tournaments")
        print(f"  corr(IPR,H)      = {pearson([r[0] for r in rows], [r[4] for r in rows]): .4f}")
        print(f"  corr(J3var,H)    = {pearson([r[2] for r in rows], [r[4] for r in rows]): .4f}")
        print(f"  corr(IPR,J3var)  = {pearson([r[0] for r in rows], [r[2] for r in rows]): .4f}")
        rows_h = sorted(rows, key=lambda r: (-r[4], -r[0], r[5]))
        rows_i = sorted(rows, key=lambda r: (-r[0], -r[4], r[5]))
        print("  top by H:")
        for rank, r in enumerate(rows_h[:5], 1):
            print(
                f"    {rank:2d}. H={r[4]:10d} ipr={r[0]:.6f} "
                f"J3var={r[2]:8.3f} S={r[5]}"
            )
        print("  top by IPR:")
        for rank, r in enumerate(rows_i[:5], 1):
            print(
                f"    {rank:2d}. ipr={r[0]:.6f} H={r[4]:10d} "
                f"J3var={r[2]:8.3f} S={r[5]}"
            )


def part_j3_energy_law():
    print("\n" + "=" * 78)
    print("III. J_3 VARIANCE IS THE ADDITIVE-ENERGY AXIS")
    print("=" * 78)
    print("Fit over the full orientation cube: Var(J_3) = a*E(S) + b.")
    print("Empirically a = 1/(p-1) and b = -(p^2 - 2p + 5)/16 exactly.")
    for p in [7, 11, 13, 17, 19, 23]:
        rows = []
        for S in valid_sets(p):
            energy = additive_energy(p, S)
            j3var = profile_stats(pinned_cycles_profile(p, S, 3))[1]
            rows.append((energy, j3var))
        n = len(rows)
        mean_e = sum(e for e, _ in rows) / n
        mean_j = sum(j for _, j in rows) / n
        var_e = sum((e - mean_e) ** 2 for e, _ in rows)
        cov = sum((e - mean_e) * (j - mean_j) for e, j in rows)
        a = cov / var_e
        b = mean_j - a * mean_e
        b_formula = -(p * p - 2 * p + 5) / 16
        max_err = max(abs((a * e + b) - j) for e, j in rows)
        slope_err = abs(a - 1 / (p - 1))
        intercept_err = abs(b - b_formula)
        print(
            f"p={p:2d} cube={n:4d} a={a:.12f} "
            f"|a-1/(p-1)|={slope_err:.2e} b={b:.6f} "
            f"|b-b_p|={intercept_err:.2e} maxerr={max_err:.2e}"
        )


def part_profiles():
    print("\n" + "=" * 78)
    print("IV. PINNED CYCLE PROFILES: FEJER LOCALITY SURVIVES EXCLUSION CORRECTIONS")
    print("=" * 78)
    for p in [7, 11, 13]:
        print(f"\np={p}")
        candidates = [("Interval", interval_set(p))]
        if p % 4 == 3:
            candidates.append(("Paley", paley_set(p)))
        for name, S in candidates:
            print(f"  {name} S={S}")
            for k in range(3, min(p, 9) + 1, 2):
                profile = pinned_cycles_profile(p, S, k)
                mean, var, spread = profile_stats(profile)
                half = profile[1 : (p + 1) // 2 + 1]
                print(
                    f"    J_{k}: half-profile={half} mean={mean:.2f} "
                    f"var={var:.2f} spread={spread}"
                )


def part_mutations():
    print("\n" + "=" * 78)
    print("V. SINGLE-FLIP MUTATIONS OFF THE INTERVAL")
    print("=" * 78)
    for p in [11, 13, 17]:
        base = interval_set(p)
        base_metrics = concentration_metrics(p, base)
        base_j3 = profile_stats(pinned_cycles_profile(p, base, 3))[1]
        base_H = hamiltonian_paths(p, base) if p <= 13 else None
        rows = []
        for label, S in single_flip_neighbors(p, base):
            m = concentration_metrics(p, S)
            j3var = profile_stats(pinned_cycles_profile(p, S, 3))[1]
            H = hamiltonian_paths(p, S) if p <= 13 else None
            rows.append((m["ipr"] - base_metrics["ipr"], j3var - base_j3, None if H is None else H - base_H, label, S))
        print(f"\np={p}, interval={base}")
        print(f"  base: ", end="")
        summarize_set(p, "Interval", base, include_h=p <= 13)
        rows.sort(key=lambda r: r[0])
        print("  flips sorted by IPR loss:")
        for dipr, dj3, dH, label, S in rows:
            hpart = "" if dH is None else f" dH={dH:9d}"
            print(f"    {label:7s} dIPR={dipr:+.6f} dJ3var={dj3:+8.3f}{hpart} S={S}")


def part_tangents():
    print("\n" + "=" * 78)
    print("VI. ODD TANGENTS WORTH NOT FORGETTING")
    print("=" * 78)
    print("1. Autocorrelation is the honest Fejer object: Q is its Fourier transform.")
    print("   J_3 is not raw autocorrelation; it is an oriented wedge fold of it.")
    print("2. Interval maximizes IPR and J_3 variance together at p=13, but p=7,11")
    print("   show the old anti/locality regime: localization helps packing before it")
    print("   helps total H.")
    print("3. Single flips usually destroy the Fejer peak and the J_3 gradient in the")
    print("   same direction. This is a finite-difference version of rearrangement.")
    print("4. The next proof object may be a pinned-exclusion expansion: start with")
    print("   spectral pinned walk counts, then subtract repeated-vertex diagrams.")


def main():
    part_identity_checks()
    part_landscape()
    part_j3_energy_law()
    part_profiles()
    part_mutations()
    part_tangents()


if __name__ == "__main__":
    main()
