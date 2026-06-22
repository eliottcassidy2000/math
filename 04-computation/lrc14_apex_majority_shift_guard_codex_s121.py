#!/usr/bin/env python3
"""S121: exact shift-sieve guardrail for the LRC14 apex-majority branch.

We study rows

    S = 14 Q union R,       |Q| >= 7, |R| <= 6, 14 not divides r in R.

After the below-14 LRC theorem gives a strict 1/14-safe point `u` for Q,
the multiples stay safe at all lifts `t=(u+k)/14`.  The question is whether
some lift also avoids all residual speeds R.

The clean finite lemma is true if R has at most one half-step speed
(`gcd(r,14)=7`).  It is false in that exact form for two half-step speeds:
the pair {7,161} covers all 14 shifts at u=2/49.  This script audits the
single-speed shift facts and prints the counterexample, using exact rationals.

Tournament Analysis note:
  vertices are proof carriers, not runners.  The observable is dependency
  in the proof DAG; the switch is whether a carrier preserves the LRC predicate
  through the quotient `t=(u+k)/14`.  The tie Hamiltonian path is the proof
  order printed at the end.
"""

from fractions import Fraction as F
from math import gcd


THRESH = F(1, 14)


def dist1(x):
    r = x % 1
    return min(r, 1 - r)


def forbidden_shifts(speed, u):
    return tuple(k for k in range(14) if dist1(F(speed, 14) * (u + k)) < THRESH)


def parity_profile(shifts):
    return (sum(1 for k in shifts if k % 2 == 0), sum(1 for k in shifts if k % 2 == 1))


def audit_single_speed_facts(max_speed=420, denom=392):
    failures = []
    halfstep_profiles = set()
    ordinary_profiles = set()
    for speed in range(1, max_speed + 1):
        if speed % 14 == 0:
            continue
        h = gcd(speed, 14)
        for a in range(denom):
            u = F(a, denom)
            fs = forbidden_shifts(speed, u)
            pp = parity_profile(fs)
            if h == 7:
                halfstep_profiles.add((len(fs), pp))
                if len(fs) not in (0, 7) or (pp not in ((0, 0), (7, 0), (0, 7))):
                    failures.append(("halfstep", speed, u, fs, pp))
            else:
                ordinary_profiles.add((len(fs), pp))
                if len(fs) > 2 or pp[0] > 1 or pp[1] > 1:
                    failures.append(("ordinary", speed, u, fs, pp))
    return failures, halfstep_profiles, ordinary_profiles


def find_two_halfstep_cover(max_speed=420, denom=392):
    halfsteps = [s for s in range(1, max_speed + 1) if gcd(s, 14) == 7]
    target = set(range(14))
    for a in range(denom):
        u = F(a, denom)
        data = [(s, set(forbidden_shifts(s, u))) for s in halfsteps]
        data = [(s, fs) for s, fs in data if fs]
        for i, (s1, f1) in enumerate(data):
            for s2, f2 in data[i + 1 :]:
                if f1 | f2 == target:
                    return u, s1, tuple(sorted(f1)), s2, tuple(sorted(f2))
    return None


def verify_named_counterexample():
    u = F(2, 49)
    f7 = set(forbidden_shifts(7, u))
    f161 = set(forbidden_shifts(161, u))
    return u, tuple(sorted(f7)), tuple(sorted(f161)), tuple(sorted(f7 | f161))


def proof_carrier_tournament():
    carriers = [
        "below14_scaled_Q_margin",
        "fourteen_lifts_t=(u+k)/14",
        "ordinary_residual_shift_sieve",
        "single_halfstep_parity_guard",
        "generic_boundary_perturbation",
        "two_halfstep_phase_collision",
    ]
    scores = {name: len(carriers) - i - 1 for i, name in enumerate(carriers)}
    edges = []
    for i, a in enumerate(carriers):
        for b in carriers[i + 1 :]:
            edges.append((a, b))
    return carriers, scores, edges


def main():
    failures, half_profiles, ordinary_profiles = audit_single_speed_facts()
    print("LRC14 apex-majority shift guard (S121)")
    print()
    print("1. Single-speed exact shift facts")
    print(f"  ordinary profiles (gcd(speed,14)!=7): {sorted(ordinary_profiles)}")
    print(f"  half-step profiles (gcd(speed,14)=7): {sorted(half_profiles)}")
    print(f"  failures: {len(failures)}")
    if failures:
        print(f"  first failure: {failures[0]}")
    print()
    print("2. Consequence")
    print("  If |R|<=6 and at most one r in R has gcd(r,14)=7,")
    print("  then the union of forbidden shift classes has size at most 12:")
    print("    no half-step: six ordinary speeds forbid at most 6*2=12 shifts;")
    print("    one half-step: it forbids one parity class, and each ordinary")
    print("      speed contributes at most one new shift in the opposite parity,")
    print("      so at most 7+5=12 shifts.")
    print()
    print("3. Sharp guardrail: two half-step speeds can cover all shifts")
    found = find_two_halfstep_cover()
    print(f"  first grid cover found: {found}")
    u, f7, f161, union = verify_named_counterexample()
    print(f"  named exact cover at u={u}:")
    print(f"    speed 7 forbids   {f7}")
    print(f"    speed 161 forbids {f161}")
    print(f"    union             {union}")
    print()
    print("4. Tournament Analysis fingerprint")
    carriers, scores, edges = proof_carrier_tournament()
    print(f"  carriers={carriers}")
    print(f"  score_histogram={sorted(scores.values(), reverse=True)}")
    print(f"  directed_edges={len(edges)}")
    print("  Hamiltonian path:")
    print("    " + " > ".join(carriers))
    print()
    print("Readout: THM-570 closes the apex-majority branch except for the")
    print("two-or-more half-step residual collision, which is now the named")
    print("HYP-2911 phase-alignment target.")


if __name__ == "__main__":
    main()
