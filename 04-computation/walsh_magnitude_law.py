#!/usr/bin/env python3
"""
walsh_magnitude_law.py -- Universal magnitude law for Walsh coefficients

From the p=7 and p=11 data, degree-2 Walsh coefficients of H take only
TWO magnitudes (one per resonance class q):
  p=7:  |h_hat_H| = 3.5 (q=3 only)
  p=11: |h_hat_H| = 272.25 (q=3) and 8.25 (q=5)

Questions:
1. Are these magnitudes expressible in terms of p, m, q?
2. What's the relationship between degree-2 and degree-4 magnitudes?
3. Is there a universal formula?

From hhat_magnitude_formula.py: |h_hat_H(q=3)| ~ 9p^2/4, |h_hat_H(q=5)| ~ 3p/4
  p=7:  3.5 = 9*49/4 ??? = 110.25 NO.
  Actually 3.5 is the CORRECT value from simple cycles.

Let me also look at alpha_j magnitudes.

Author: kind-pasteur-2026-03-12-S60
"""


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def resonance_level(a, b, p):
    for q in range(1, p, 2):
        if (q * a - b) % p == 0 or (q * a + b) % p == 0:
            return q
        if q > 1 and ((a - q * b) % p == 0 or (a + q * b) % p == 0):
            return q
    return p


def main():
    print("=" * 70)
    print("WALSH MAGNITUDE LAW")
    print("=" * 70)

    # Collect all degree-2 data
    data_p7 = {
        'H': {(0,1): 3.5, (0,2): 3.5, (1,2): 3.5},
        'a1': {(0,1): 5.25, (0,2): 5.25, (1,2): 5.25},
        'a2': {(0,1): -1.75, (0,2): -1.75, (1,2): -1.75},
    }

    data_p11 = {
        'H': {(0,1): -8.25, (0,2): 272.25, (0,3): 272.25, (0,4): 8.25,
              (1,2): -272.25, (1,3): -8.25, (1,4): -272.25,
              (2,3): 8.25, (2,4): 8.25, (3,4): 272.25},
        'a1': {(0,1): -35.75, (0,2): 13.75, (0,3): 13.75, (0,4): 35.75,
               (1,2): -13.75, (1,3): -35.75, (1,4): -13.75,
               (2,3): 35.75, (2,4): 35.75, (3,4): 13.75},
        'a2': {(0,1): 47.4375, (0,2): 29.5625, (0,3): 29.5625, (0,4): -47.4375,
               (1,2): -29.5625, (1,3): 47.4375, (1,4): -29.5625,
               (2,3): -47.4375, (2,4): -47.4375, (3,4): 29.5625},
        'a3': {(0,1): -15.8125, (0,2): 15.8125, (0,3): 15.8125, (0,4): 15.8125,
               (1,2): -15.8125, (1,3): -15.8125, (1,4): -15.8125,
               (2,3): 15.8125, (2,4): 15.8125, (3,4): 15.8125},
    }

    # Degree-4 data at p=11
    deg4_p11 = {
        'H': 118.25,  # uniform magnitude, sign = chi(prod)
        'a1': 321.75,
        'a2': 19.9375,  # sign = -chi(prod) (anti)
        'a3': 55.6875,  # sign = -chi(prod) (anti)
    }

    for p, data in [(7, data_p7), (11, data_p11)]:
        m = (p - 1) // 2
        print(f"\n{'='*60}")
        print(f"p={p}, m={m}")
        print("=" * 60)

        # Group by resonance level
        pairs = [(a, b) for a in range(m) for b in range(a+1, m)]
        by_q = {}
        for a, b in pairs:
            q = resonance_level(a + 1, b + 1, p)
            if q not in by_q:
                by_q[q] = []
            by_q[q].append((a, b))

        print(f"\n  Resonance classes: {dict((q, len(ps)) for q, ps in sorted(by_q.items()))}")

        for name in sorted(data.keys()):
            print(f"\n  {name}:")
            for q in sorted(by_q.keys()):
                mags = [abs(data[name][pair]) for pair in by_q[q]]
                mag = mags[0]
                uniform = all(abs(m - mag) < 0.001 for m in mags)
                print(f"    q={q}: |h_hat| = {mag:.4f} "
                      f"({'UNIFORM' if uniform else f'VARIES: {mags}'})")

        # Express magnitudes as fractions
        print(f"\n  Magnitudes as fractions of p:")
        for name in sorted(data.keys()):
            for q in sorted(by_q.keys()):
                mag = abs(data[name][by_q[q][0]])
                # Try various forms
                for denom in [1, 2, 4, 8, 16, 32]:
                    val = mag * denom
                    if abs(val - round(val)) < 0.01:
                        print(f"    {name} q={q}: {mag:.4f} = {round(val)}/{denom}")
                        break

    # At p=11, degree-4 magnitudes
    print(f"\n\n{'='*60}")
    print(f"DEGREE-4 MAGNITUDES at p=11")
    print("=" * 60)
    for name in ['H', 'a1', 'a2', 'a3']:
        mag = deg4_p11[name]
        print(f"  {name}: |h_hat_4| = {mag:.4f}")
        for denom in [1, 2, 4, 8, 16, 32]:
            val = mag * denom
            if abs(val - round(val)) < 0.01:
                print(f"    = {round(val)}/{denom}")
                break

    # Relationships between degree-2 and degree-4
    print(f"\n  Relationships:")
    for name in ['H', 'a1', 'a2', 'a3']:
        d2_q3 = abs(data_p11[name][by_q[3][0]]) if name in data_p11 else 0
        d2_q5 = abs(data_p11[name][by_q[5][0]]) if name in data_p11 else 0
        d4 = deg4_p11[name]
        if d2_q3 > 0 and d2_q5 > 0:
            print(f"  {name}: d2_q3={d2_q3:.4f}, d2_q5={d2_q5:.4f}, d4={d4:.4f}")
            print(f"    d4/d2_q3 = {d4/d2_q3:.6f}")
            print(f"    d4/d2_q5 = {d4/d2_q5:.6f}")
            print(f"    d2_q3/d2_q5 = {d2_q3/d2_q5:.6f}")
            print(f"    d2_q3 * d2_q5 = {d2_q3 * d2_q5:.4f}")
            print(f"    d2_q3 + d2_q5 = {d2_q3 + d2_q5:.4f}")
            print(f"    d2_q3 - d2_q5 = {d2_q3 - d2_q5:.4f}")

    # Verify H OCF: h_H = 2*h_a1 + 4*h_a2 + 8*h_a3
    print(f"\n  OCF verification at degree-4:")
    h_H_d4_recon = 2 * deg4_p11['a1'] + 4 * deg4_p11['a2'] + 8 * deg4_p11['a3']
    print(f"    2*{deg4_p11['a1']} + 4*{deg4_p11['a2']} + 8*{deg4_p11['a3']} = {h_H_d4_recon:.4f}")
    print(f"    Actual h_H_d4 = {deg4_p11['H']}")
    # Note: sign matters!
    h_H_d4_signed = 2 * (-321.75) + 4 * (19.9375) + 8 * (55.6875)
    print(f"    Signed: 2*(-321.75) + 4*(19.9375) + 8*(55.6875) = {h_H_d4_signed:.4f}")
    h_H_d4_alt = 2 * (321.75) + 4 * (-19.9375) + 8 * (-55.6875)
    print(f"    Alt: 2*(321.75) + 4*(-19.9375) + 8*(-55.6875) = {h_H_d4_alt:.4f}")

    # The right version: for chi(prod)=+1:
    # a1: +321.75, a2: -19.9375, a3: -55.6875
    h_d4_pos = 2*321.75 + 4*(-19.9375) + 8*(-55.6875)
    print(f"    chi(prod)=+1: 2(321.75)+4(-19.9375)+8(-55.6875) = {h_d4_pos:.4f}")
    # Expect +118.25

    h_d4_neg = 2*(-321.75) + 4*(19.9375) + 8*(55.6875)
    print(f"    chi(prod)=-1: 2(-321.75)+4(19.9375)+8(55.6875) = {h_d4_neg:.4f}")
    # Expect -118.25

    print("\nDONE.")


if __name__ == '__main__':
    main()
