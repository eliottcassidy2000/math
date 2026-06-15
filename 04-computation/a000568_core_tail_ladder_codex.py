#!/usr/bin/env python3
"""
a000568_core_tail_ladder_codex.py

codex-2026-06-15-S12

Exact Burnside kernel extraction for A000568 after peeling the 1-tail.

Main objects:
  a(n) = sum over odd partitions lambda of n of 2^e(lambda)/z(lambda)

Split lambda = mu U 1^r, where mu has odd parts >= 3.
If |mu| = m and len(mu) = t, then
  e(lambda) = e(mu) + C(r,2) + r*t
  z(lambda) = z(mu) * r!

So
  a(n) = sum_{m,t} B[m,t] * 2^(C(n-m,2) + (n-m)t) / (n-m)!
with the exact core kernel
  B[m,t] = sum_{mu of size m, length t, odd parts >= 3} 2^e(mu) / z(mu).

It is cleaner integrally to define
  K[m,t] = m! * B[m,t] in Z_{>=0}.
Then
  n! a(n) = sum_{m,t} C(n,m) * K[m,t] * 2^(C(n-m,2) + (n-m)t).

This script:
  1. extracts K[m,t] exactly up to N;
  2. verifies the peeled formula against exact A000568 values up to n=20 and n=100;
  3. proves the active-state count at n=100 is exactly 834;
  4. verifies the next exact rung: peeling the 3-tail only needs one extra
     statistic c3 = number of remaining core parts divisible by 3.
"""

from collections import defaultdict
from fractions import Fraction
from math import comb, factorial, gcd


KNOWN = {
    1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
    9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168,
    13: 48542114686912, 14: 28401423719122304,
    15: 31021002160355166848,
    16: 63530415842308265100288,
    17: 244912778438520759443245824,
    18: 1783398846284777975419600287232,
    19: 24605641171260376770598003978281472,
    20: 645022068557873570931850526424042500096,
}

A100 = int(
    "1344234018830050151564780020895938096109808495427980734364623804072132907"
    "7533435494973955760104942457014606847008800075709883275222927949216976884"
    "8160000183618389033370455493448802541206980041054338793531485803213029091"
    "7174300999972198689579994381905422764066856105336835197886517573825587271"
    "6702792840904699583716750233247839360259050389291641888200682131023724211"
    "0626960831810984531870265547103956616258777873789708283632486931366819647"
    "4929122503512274081026319710955325388761318437959332825429272057097901052"
    "0733060528211442308711294856501378381685272845240406471872857412369411560"
    "8585993907002984756261082295132716438445495829476650486973791358206317897"
    "4602405695293172119994272761616735173554923428812677207497307717758567534"
    "0845558980759849689025904834723337323468530804943547000492785208349598988"
    "2804409727645055875197533759975341292069914508357536418467769620429805827"
    "3684086730624992672148802739877560069622154816882030932308413550477079496"
    "0906780269086232079344472937045734792643656334466440727093581478557227141"
    "8471746716666530123463267035947453857529101707034349312840042784124464088"
    "3018248530830082017669464820549904592243765418224667842784663553704671477"
    "4484251724089893346228928918107332790375472638258677940846887096466735079"
    "9476066499445464008983796843642335324911602613073272432978620214660148832"
    "6026841632941277184"
)


def odd_partition_count(n):
    ways = [0] * (n + 1)
    ways[0] = 1
    for part in range(1, n + 1, 2):
        for s in range(part, n + 1):
            ways[s] += ways[s - part]
    return ways[n]


def support_count_formula(n):
    total = 1
    for t in range(1, n // 3 + 1):
        total += (n - 3 * t) // 2 + 1
    return total


def build_kernels(max_n):
    facts = [factorial(i) for i in range(max_n + 1)]
    B = defaultdict(Fraction)   # (m,t) -> exact rational kernel
    B3 = defaultdict(Fraction)  # (m,t,c3) -> exact 3-free rational kernel

    def rec(max_part, mass, length, e_val, z_val, groups, c3_count):
        if mass:
            term = Fraction(1 << e_val, z_val)
            B[(mass, length)] += term
            if all(part != 3 for part, _ in groups):
                B3[(mass, length, c3_count)] += term

        start = min(max_part, max_n - mass)
        if start % 2 == 0:
            start -= 1
        for part in range(start, 2, -2):
            cross = 0
            for old_part, old_mult in groups:
                cross += old_mult * gcd(part, old_part)
            max_mult = (max_n - mass) // part
            part_pow = 1
            fact_mult = 1
            for mult in range(1, max_mult + 1):
                part_pow *= part
                fact_mult *= mult
                delta_e = (
                    mult * ((part - 1) // 2)
                    + (mult * (mult - 1) // 2) * part
                    + mult * cross
                )
                rec(
                    part - 2,
                    mass + mult * part,
                    length + mult,
                    e_val + delta_e,
                    z_val * part_pow * fact_mult,
                    groups + ((part, mult),),
                    c3_count + (mult if part % 3 == 0 else 0),
                )

    rec(max_n if max_n % 2 else max_n - 1, 0, 0, 0, 1, tuple(), 0)
    B[(0, 0)] = Fraction(1, 1)
    B3[(0, 0, 0)] = Fraction(1, 1)

    K = {}
    for (m, t), value in B.items():
        scaled = facts[m] * value
        assert scaled.denominator == 1
        K[(m, t)] = scaled.numerator

    K3 = {}
    for (m, t, c3), value in B3.items():
        scaled = facts[m] * value
        assert scaled.denominator == 1
        K3[(m, t, c3)] = scaled.numerator

    return B, B3, K, K3, facts


def a_from_K(n, K, facts):
    total = 0
    for (m, t), kval in K.items():
        if m > n:
            continue
        r = n - m
        total += comb(n, m) * kval * (1 << (r * (r - 1) // 2 + r * t))
    assert total % facts[n] == 0
    return total // facts[n]


def reconstruct_K_from_K3(max_n, K3, facts):
    recon = defaultdict(int)
    for total_mass in range(max_n + 1):
        for total_len in range(total_mass + 1):
            acc = 0
            max_s = min(total_len, total_mass // 3)
            for s in range(max_s + 1):
                base_mass = total_mass - 3 * s
                base_len = total_len - s
                if base_mass < 0 or base_len < 0:
                    continue
                pref = facts[total_mass] // (facts[base_mass] * (3 ** s) * factorial(s))
                for (m, t, c3), kval in K3.items():
                    if m != base_mass or t != base_len:
                        continue
                    delta = s + 3 * s * (s - 1) // 2 + s * (t + 2 * c3)
                    acc += pref * kval * (1 << delta)
            if acc:
                recon[(total_mass, total_len)] = acc
    return recon


def print_small_kernel(K, facts, upto=15):
    print("Small core kernel rows B[m,t] = K[m,t]/m!:")
    for m in range(upto + 1):
        entries = []
        for (mm, t), kval in sorted(K.items()):
            if mm != m:
                continue
            entries.append(f"t={t}:{kval}/{facts[m]}")
        if entries:
            print(f"  m={m:2d} -> " + ", ".join(entries))


def main():
    max_n = 100
    B, B3, K, K3, facts = build_kernels(max_n)

    print("=" * 72)
    print("A000568 CORE TAIL LADDER")
    print("=" * 72)

    odd_parts_100 = odd_partition_count(100)
    support_100 = support_count_formula(100)
    print(f"odd partitions of 100             : {odd_parts_100}")
    print(f"active (m,t) core states at n=100 : {support_100}")
    print(f"outer collapse factor             : {odd_parts_100 / support_100:.3f}x")
    print()

    print("Verifying unit-tail core formula against exact A000568:")
    ok = True
    for n in range(1, 21):
        val = a_from_K(n, K, facts)
        match = (val == KNOWN[n])
        ok &= match
        mark = "OK" if match else f"FAIL got {val}"
        print(f"  a({n:2d}) = {val} [{mark}]")
    print(f"  a(100) exact via K = {a_from_K(100, K, facts)}")
    print(f"  a(100) matches stored exact value: {a_from_K(100, K, facts) == A100}")
    print()

    print_small_kernel(K, facts)
    print()

    print("Integrality checks:")
    print(f"  active K states  : {len(K)}")
    print(f"  active K3 states : {len(K3)}")
    print(f"  all m! * B[m,t] integral  : {all((facts[m] * val).denominator == 1 for (m, _), val in B.items())}")
    print(
        "  all m! * B3[m,t,c3] integral : "
        f"{all((facts[m] * val).denominator == 1 for (m, _, _), val in B3.items())}"
    )
    print()

    print("Verifying the 3-tail refinement K <- C3 exactly through mass 100:")
    recon = reconstruct_K_from_K3(max_n, K3, facts)
    mismatch = [
        key for key in set(K) | set(recon)
        if K.get(key, 0) != recon.get(key, 0)
    ]
    print(f"  mismatches: {len(mismatch)}")
    if mismatch:
        print(f"  first mismatch: {mismatch[0]} -> K={K.get(mismatch[0], 0)} "
              f"recon={recon.get(mismatch[0], 0)}")
    else:
        print("  exact match on all (m,t) states.")
    print()

    print("Support sanity:")
    support_from_K = len([1 for (m, t) in K if m <= 100])
    print(f"  support(K up to 100) = {support_from_K}")
    print(f"  support formula      = {support_100}")
    print(f"  active 3-free states = {len(K3)}")
    print()

    print("DONE.")


if __name__ == "__main__":
    main()
