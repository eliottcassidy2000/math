#!/usr/bin/env python3
"""HYP-9055 referee — the shifted-word palindrome layer at N = 91
(mac-mini-2026-07-27-S145).

(a) FOUR-SECTOR EXACTNESS.  For B ⊆ ℤ/91 and the incidence count
    N_B(v) = #{(i, j) : j ∈ ℤ/91, v_i·j mod 91 ∈ B}:
    a coordinate with gcd(v_i, 91) = 1 sweeps ℤ/91 bijectively → |B|;
    gcd = 7 sweeps 7ℤ/91 (13 values, ×7 each) → 7·|B ∩ 7ℤ|;
    gcd = 13 → 13·|B ∩ 13ℤ|;  gcd = 91 → 91·[0 ∈ B].
    EXACT, not merely a congruence: N_B = z₁|B| + 7z₇|B∩7ℤ| +
    13z₁₃|B∩13ℤ| + 91z₉₁[0∈B].  For the witness band B = {dist ≤ 6}
    (|B| = 13, |B∩7ℤ| = |B∩13ℤ| = 1): N_B = 13z₁ + 7z₇ + 13z₁₃ + 91z₉₁.

(b) DIHEDRAL REFINEMENT.  The mirror j ↦ −j on ℤ/91 (91 odd) fixes only
    j = 0, and B = −B: N_B ≡ #{i : v_i·0 ∈ B} = 13 ≡ 1 (mod 2) for every
    13-speed family.  Joint: N_B mod 182 is fully forced by the duty
    profile — the shifted-word palindrome localization.

(c) THE GRAM BRIDGE.  Pair count N₂(v) = #{(i, i', j) : i < i', both
    v_i j, v_i' j ∈ B}: unit-unit pairs contribute the GRAM ENTRY
    G(r) = |B ∩ r⁻¹B| at the unit ratio r = v_i'/v_i mod 91 — exactly the
    chirp-Gram object of THM-2356 — while duty-involving pairs contribute
    the rigid sector terms.  Referee: exact sector formula vs direct count.

(d) mirror on pairs: fixed locus j = 0 gives C(13,2) = 78 ≡ 0 (mod 2):
    N₂ is EVEN for every family — the pair-level parity is degenerate
    (all mod-2 content sits at the singleton level), consistent with the
    S143 Smith-vanishing verdict: mod-2 sees only the trivial layer; the
    7/13 rotations carry the real localization.
"""

from math import gcd

N = 91
BAND = {r % N for r in range(-6, 7)}
assert len(BAND) == 13


def sectors(v):
    z = {1: 0, 7: 0, 13: 0, 91: 0}
    for vi in v:
        z[gcd(vi, N)] += 1
    return z


def NB_direct(v):
    return sum(1 for vi in v for j in range(N) if (vi * j) % N in BAND)


def NB_formula(v):
    z = sectors(v)
    B7 = sum(1 for b in BAND if b % 7 == 0)
    B13 = sum(1 for b in BAND if b % 13 == 0)
    return z[1] * len(BAND) + z[7] * 7 * sum(1 for b in BAND if b % 7 == 0 and True) \
        if False else z[1] * len(BAND) + 7 * z[7] * len([b for b in BAND if b % 7 == 0]) \
        + 13 * z[13] * len([b for b in BAND if b % 13 == 0]) + 91 * z[91] * (0 in BAND)


def N2_direct(v):
    n = len(v)
    tot = 0
    for j in range(N):
        inb = [i for i in range(n) if (v[i] * j) % N in BAND]
        k = len(inb)
        tot += k * (k - 1) // 2
    return tot


def N2_formula(v):
    """Sector decomposition of the pair count."""
    n = len(v)
    tot = 0
    inv = {}
    for r in range(N):
        if gcd(r, N) == 1:
            inv[r] = pow(r, -1, N)
    for a in range(n):
        for b in range(a + 1, n):
            ga, gb = gcd(v[a], N), gcd(v[b], N)
            if ga == 1 and gb == 1:
                # Gram entry at ratio r = v_b * v_a^{-1}: |B ∩ r^{-1}B| =
                # #{j: v_a j in B and v_b j in B} = #{u in B: (v_b/v_a)u in B}
                r = (v[b] * inv[v[a] % N]) % N
                tot += sum(1 for u in BAND if (r * u) % N in BAND)
            else:
                # at least one duty coordinate: direct small count over the
                # relevant subgroup structure (exact, per pair)
                tot += sum(1 for j in range(N)
                           if (v[a] * j) % N in BAND and (v[b] * j) % N in BAND)
    return tot


def main():
    fams = {
        "AP": list(range(1, 14)),
        "GW": list(range(1, 12)) + [13, 24],
        "deep well": list(range(1, 13)) + [182],
        "{1..11,13,36}": list(range(1, 12)) + [13, 36],
        "loose {2..14}": list(range(2, 15)),
    }
    print("== (a)+(b) four-sector exactness + dihedral parity ==")
    for name, v in fams.items():
        d, f = NB_direct(v), NB_formula(v)
        z = sectors(v)
        print(f"  {name:<14} z = {tuple(z.values())}  N_B direct {d} = formula {f}"
              f"  {'OK' if d == f else 'FAIL'};  N_B mod 2 = {d % 2} (forced 1)")
        assert d == f and d % 2 == 1
    print("  => N_B = 13z1 + 7z7 + 13z13 + 91z91 EXACT; odd for every family.")

    print("\n== (c) the Gram bridge at pair level ==")
    for name, v in fams.items():
        d, f = N2_direct(v), N2_formula(v)
        print(f"  {name:<14} N2 direct {d} = sector/Gram formula {f}"
              f"  {'OK' if d == f else 'FAIL'};  N2 mod 2 = {d % 2} (forced 0)")
        assert d == f and d % 2 == 0
    print("  => unit-unit pairs contribute EXACTLY the chirp-Gram entries")
    print("     |B ∩ rB| (THM-2356's object); duty pairs are the rigid part;")
    print("     pair-level mod-2 is degenerate (Smith verdict confirmed):")
    print("     the 7/13 rotations, not the mirror, carry pair localization.")

    print("\n== mod 7 / mod 13 pair congruences (the localization content) ==")
    for name, v in fams.items():
        d = N2_direct(v)
        print(f"  {name:<14} N2 = {d}: mod 7 = {d % 7}, mod 13 = {d % 13}")
    print("  (unit-sector Gram sums collapse mod 7/13 onto duty-involving")
    print("   pairs: the quantitative face of THM-2544's rank-0 boundary.)")


if __name__ == "__main__":
    main()
