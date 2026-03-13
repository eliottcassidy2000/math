#!/usr/bin/env python3
"""
additive_energy_disjointness_proof.py -- Toward proving HYP-480

KEY IDEA: The additive energy E(S) = sum r(s)^2 measures how "additively
structured" the connection set S is. Higher E(S) means vertices have
more similar neighborhoods, which creates more disjoint cycles.

For circulant tournaments:
- S = {1,...,m} (interval): E = m(2m^2+1)/3 ~ (2/3)m^3
- S = QR_p: E = m^2 + (p-1)*((m-1)/2)^2 ~ m^3/4

Ratio: E(Int)/E(QR) -> (2/3)/(1/4) = 8/3 ... no, that gives a different answer.

Let me compute more carefully.

E(Int) = sum_{s=-(m-1)}^{m-1} (m - |s|)^2
       = m^2 + 2*sum_{k=1}^{m-1} (m-k)^2
       = m^2 + 2*sum_{j=1}^{m-1} j^2
       = m^2 + 2*(m-1)*m*(2m-1)/6
       = m^2 + m*(m-1)*(2m-1)/3
       = (3m^2 + 2m^3 - 3m^2 + m)/3
       = m(2m^2+1)/3

E(QR) = m^2 + (p-1)*lambda^2  where lambda = (p-3)/4 = (m-1)/2
       = m^2 + 2m * ((m-1)/2)^2
       = m^2 + m*(m-1)^2/2
       = m*(2m + (m-1)^2)/2
       = m*(m^2+1)/2

Ratio = E(Int)/E(QR) = [m(2m^2+1)/3] / [m(m^2+1)/2] = 2(2m^2+1)/(3(m^2+1))
       -> 4/3 as m -> infinity

This matches the computed data!

Now: does higher additive energy IMPLY more disjoint cycle pairs?

THEOREM ATTEMPT: For circulant tournaments on Z_p, the number of
disjoint 3-cycle pairs is a MONOTONE FUNCTION of the additive energy
of the connection set.

Proof sketch:
- A 3-cycle on vertices {a, b, c} exists iff both {b-a, c-b, a-c} all
  have representatives in S (one from each complementary pair)
- Two 3-cycles {a,b,c} and {d,e,f} are disjoint iff their vertex sets
  don't overlap
- The probability (in a random-S sense) that two specific vertex-disjoint
  3-element subsets both form 3-cycles is proportional to E(S)
  (because shared neighborhoods amplify cycle probabilities)

Author: kind-pasteur-2026-03-12-S58
"""

import math


def additive_energy_formula():
    """Verify the additive energy formulas."""
    print("=" * 70)
    print("ADDITIVE ENERGY FORMULAS VERIFICATION")
    print("=" * 70)

    for p in [7, 11, 19, 23, 29, 31, 43, 47]:
        m = (p - 1) // 2

        # E(Interval) = m(2m^2+1)/3
        E_int = m * (2*m*m + 1) // 3

        # E(QR) = m(m^2+1)/2
        E_qr = m * (m*m + 1) // 2

        ratio = E_int / E_qr
        limit = 4/3

        print(f"  p={p:3d}, m={m:3d}: E(Int)={E_int:>10d}, E(QR)={E_qr:>10d}, "
              f"ratio={ratio:.6f} (limit 4/3={limit:.6f})")


def disjoint_3cycle_formula():
    """Compute alpha_2(3,3) for both tournaments analytically."""
    print(f"\n{'='*70}")
    print("DISJOINT 3-CYCLE PAIRS: ANALYTICAL APPROACH")
    print("=" * 70)

    print("""
  For a circulant tournament on Z_p with connection set S:
  - A 3-cycle on {a, a+d1, a+d1+d2} exists iff d1, d2, -(d1+d2) each
    map to the correct element of their complementary pair
  - c_3 depends only on |S| = m (all regular tournaments have
    c_3 = C(p,3) - p*C(m,2))
  - But WHICH triples are 3-cycles depends on S

  For the NUMBER of disjoint 3-cycle pairs:
  alpha_2(3,3) = #{(C1, C2) : C1, C2 are 3-cycles, V(C1) cap V(C2) = empty} / 2

  By circulant symmetry, we can fix one cycle and count disjoint others.
""")

    for p in [7, 11]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        # Build adjacency
        for name, S in [("Paley", S_qr), ("Interval", S_int)]:
            S_set = set(S)
            # Count 3-cycles through vertex 0
            cycles_through_0 = []
            for i in range(1, p):
                for j in range(i+1, p):
                    # Check if {0, i, j} is a 3-cycle
                    # 3-cycle: 0->i->j->0 or 0->j->i->0
                    d01, d02 = i, j
                    d10, d20 = p - i, p - j
                    d12, d21 = (j - i) % p, (i - j) % p
                    if d01 in S_set and d12 in S_set and d20 in S_set:
                        cycles_through_0.append(frozenset({0, i, j}))
                    elif d02 in S_set and d21 in S_set and d10 in S_set:
                        cycles_through_0.append(frozenset({0, i, j}))

            c3_through_0 = len(cycles_through_0)
            c3 = p * c3_through_0 // 3

            # Count disjoint pairs through vertex 0
            # A pair (C1, C2) with 0 in C1 but not in C2
            disjoint_from_0 = 0
            for cyc in cycles_through_0:
                # Count 3-cycles disjoint from cyc
                other_verts = [v for v in range(p) if v not in cyc]
                for i_idx in range(len(other_verts)):
                    for j_idx in range(i_idx+1, len(other_verts)):
                        a, b = other_verts[i_idx], other_verts[j_idx]
                        for c in other_verts[j_idx+1:]:
                            triple = (a, b, c)
                            # Check all 3 for 3-cycle
                            ab, bc, ca = (b-a)%p, (c-b)%p, (a-c)%p
                            ba, cb, ac = (a-b)%p, (b-c)%p, (c-a)%p
                            if (ab in S_set and bc in S_set and ca in S_set) or \
                               (ac in S_set and cb in S_set and ba in S_set):
                                disjoint_from_0 += 1

            # Each disjoint pair is counted: C1 through 0 times C2 not through 0
            # But we need to normalize...
            # alpha_2(3,3) = #{unordered pairs of disjoint 3-cycles}
            # By symmetry: each pair is counted p times (once per vertex of C1)
            # times 1 (C1 is a specific cycle through v0)
            # Actually: for each of the p*c3_through_0/3 = c3 total cycles C1,
            # the number of disjoint C2 is the same by symmetry.
            # So alpha_2(3,3) = c3 * (avg disjoint C2 per C1) / 2
            avg_disjoint = disjoint_from_0 / c3_through_0 if c3_through_0 > 0 else 0
            alpha2_33 = c3 * avg_disjoint // 2

            print(f"\n  p={p}, {name}:")
            print(f"    c3 = {c3}, cycles through 0 = {c3_through_0}")
            print(f"    disjoint_from_0 (total) = {disjoint_from_0}")
            print(f"    avg disjoint per cycle = {avg_disjoint:.2f}")
            print(f"    alpha_2(3,3) = {alpha2_33}")


def asymptotic_advantage():
    """Compute the asymptotic alpha_1 and alpha_2+ scaling."""
    print(f"\n{'='*70}")
    print("ASYMPTOTIC SCALING OF ISING DECOMPOSITION")
    print("=" * 70)

    # From computed data:
    data = {
        7: {"D1": 42, "D2plus": 28, "H_P": 189, "H_I": 175},
        11: {"D1": 5544, "D2plus": 3476, "H_P": 95095, "H_I": 93027},
        19: {"D1": 9043330440, "D2plus": 20560408288,
             "H_P": 1172695746915, "H_I": 1184212824763},
    }

    print(f"\n  {'p':>3s}  {'D1 (Paley adv)':>16s}  {'D2+ (Int adv)':>16s}  {'D1/D2+':>8s}  {'D1/H_P':>10s}  {'D2+/H_P':>10s}")
    for p in [7, 11, 19]:
        d = data[p]
        ratio = d["D1"] / d["D2plus"] if d["D2plus"] > 0 else float('inf')
        d1_frac = d["D1"] / d["H_P"]
        d2_frac = d["D2plus"] / d["H_P"]
        print(f"  {p:3d}  {d['D1']:>16,d}  {d['D2plus']:>16,d}  {ratio:>8.4f}  {d1_frac:>10.6f}  {d2_frac:>10.6f}")

    print("""
  KEY PATTERN:
  - D1/D2+ ratio: 1.500 (p=7), 1.595 (p=11), 0.440 (p=19)
  - The ratio INCREASES from p=7 to p=11, then DROPS sharply at p=19
  - This is the signature of a PHASE TRANSITION

  MECHANISM:
  - D1 grows as alpha_1 grows (roughly exponentially in p)
  - D2+ grows FASTER because the number of disjoint cycle packings
    grows super-exponentially (combinatorial explosion of independent sets)
  - At p <= 11: Paley's cycle count advantage dominates
  - At p >= 19: Interval's packing advantage dominates

  SCALING ARGUMENT:
  - alpha_1 ~ (p-1)! / 2^{p-1} (dominated by c_{p-2} and c_p)
  - Paley advantage in alpha_1: multiplicative factor ~1.1 (from quasi-random excess)
  - alpha_2 ~ alpha_1^2 * (disjointness rate)
  - Interval advantage in disjointness rate: additive energy factor ~4/3
  - So: D1 ~ 0.1 * alpha_1, D2+ ~ 0.33 * alpha_1^2 * (disjointness rate)
  - D2+ / D1 ~ 3.3 * alpha_1 * (disjointness rate) / 1
  - This GROWS with p, eventually exceeding 1
""")


def main():
    additive_energy_formula()
    disjoint_3cycle_formula()
    asymptotic_advantage()


if __name__ == '__main__':
    main()
