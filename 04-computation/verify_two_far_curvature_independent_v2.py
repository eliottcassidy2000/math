#!/usr/bin/env python3
"""INDEPENDENT adversarial verification of THM-548 / HYP-2679 two-far curvature claims.

Written from scratch (own exact-Fraction p0 and I_B) to try to REFUTE:

  (a) consecutive-pair curvature I_B(u,u+1) for B = (0,1,2,3,4,5,6,7) SATURATES to a
      finite limit ~0.0139 as u -> infinity (does NOT diverge).
      Compute for u = 20, 80, 320, 1280, 5120.

  (b) |I_B(u,v) - Phi_2(B)| * resdist(u,v) is BOUNDED (~0.01),
      resdist = min over (m,n) != (0,0), |m|,|n| <= 7, of |m*u + n*v|.
      Compute for several pairs incl. resonant (consecutive) and independent.

DEFINITION of p0 (re-derived from canon, NOT imported):
  7 sectors on the circle.  Runner e sits at sector floor(7*e*t) mod 7 at time t.
  The 6 INNER sectors are {1,2,3,4,5,6} (sector 0 excluded).
  p0(S) = measure of t in [0,1) such that EVERY inner sector {1..6} is occupied
          by at least one runner e in S.  (i.e. "missed = 0 sectors".)
  Exact: walls of runner e are at t = j/(7e), j=0..7e.  Pick a common refinement,
         evaluate the cover mask on each open subinterval via its midpoint.

p_t(B) = measure of t where B misses EXACTLY t of the 6 inner sectors.
Phi_2(B) = (2*p2(B) - p1(B)) / 49.

I_B(u,v) = p0(B u {u,v}) - p0(B u {u}) - p0(B u {v}) + p0(B).
"""

from __future__ import annotations

from fractions import Fraction as Fr
from functools import lru_cache, reduce
from math import gcd


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def sector_of(e: int, num: int, den: int) -> int:
    """Sector index floor(7*e*t) mod 7, where t = num/den exactly.

    floor(7*e*num/den) mod 7.
    """
    return (7 * e * num // den) % 7


@lru_cache(maxsize=None)
def missed_profile(row: tuple[int, ...]) -> tuple[Fr, ...]:
    """Return (p_0, p_1, ..., p_6) where p_t = measure of t missing exactly t inner sectors.

    Independent re-implementation. Walls live at t = j/(7e). We refine using the
    common denominator D = 7 * L where L = lcm of nonzero speeds; every wall t=j/(7e)
    is an integer multiple of (D/(7e)) = L/e in numerator units over D. We evaluate
    each cell's cover mask at its midpoint, using exact integer arithmetic.
    """
    nonzero = [e for e in row if e != 0]
    if not nonzero:
        # no runners -> all 6 inner sectors missed everywhere
        out = [Fr(0)] * 7
        out[6] = Fr(1)
        return tuple(out)

    L = reduce(lcm, nonzero)
    D = 7 * L  # common denominator for time t -> t = (integer)/D

    # collect breakpoints (numerators over D)
    bps = {0, D}
    for e in nonzero:
        # walls at t = j/(7e) -> numerator = j * (D/(7e)) = j * (L/e)
        step = L // e  # exact since e | L
        x = 0
        for _ in range(7 * e + 1):
            bps.add(x)
            x += step
    bp = sorted(bps)

    nums = [0] * 7  # accumulate cell lengths (in units of 1/D) by #missed inner sectors
    # midpoint of cell [lo,hi] over D is (lo+hi)/(2D); evaluate sector via num=lo+hi, den=2*D
    den_mid = 2 * D
    for lo, hi in zip(bp, bp[1:]):
        if hi <= lo:
            continue
        mid_num = lo + hi  # so t_mid = mid_num / (2D)
        mask = 0
        for e in nonzero:
            s = sector_of(e, mid_num, den_mid)
            mask |= 1 << s
        # inner sectors are 1..6
        inner_covered = bin(mask & 0b1111110).count("1")
        missed = 6 - inner_covered
        nums[missed] += hi - lo

    return tuple(Fr(n, D) for n in nums)


def p0(row: tuple[int, ...]) -> Fr:
    return missed_profile(tuple(sorted(set(row))))[0]


def p_t(row: tuple[int, ...], t: int) -> Fr:
    return missed_profile(tuple(sorted(set(row))))[t]


def Phi_2(B: tuple[int, ...]) -> Fr:
    prof = missed_profile(tuple(sorted(set(B))))
    p1, p2 = prof[1], prof[2]
    return (2 * p2 - p1) / 49


def I_B(B: tuple[int, ...], u: int, v: int) -> Fr:
    Bs = tuple(sorted(set(B)))
    return (
        p0(Bs + (u, v))
        - p0(Bs + (u,))
        - p0(Bs + (v,))
        + p0(Bs)
    )


def resdist(u: int, v: int, lim: int = 7) -> int:
    best = None
    for m in range(-lim, lim + 1):
        for n in range(-lim, lim + 1):
            if m == 0 and n == 0:
                continue
            val = abs(m * u + n * v)
            if best is None or val < best:
                best = val
    return best


def fmt(q: Fr) -> str:
    return f"{q} = {float(q):.9f}"


def cross_check_p0_against_canon() -> None:
    """Sanity: reproduce known exact p0 values from the canon output."""
    print("== cross-check: reproduce canon-reported exact p0 values ==")
    cases = [
        ((0, 4, 6, 8, 10, 12, 14), Fr(67, 1470), "core HYP2675"),
        ((0, 4, 6, 8, 10, 12, 14, 15), Fr(13, 120), "core+15"),
        ((0, 4, 6, 8, 10, 12, 14, 16), Fr(1609, 5880), "core+16"),
        ((0, 4, 6, 8, 10, 12, 14, 15, 16), Fr(321, 980), "full HYP2675 leader"),
        ((0, 3, 6, 9, 12, 14), Fr(0), "k8 core"),
        ((0, 2, 4, 6, 8, 10, 12), Fr(31, 210), "even AP core"),
    ]
    ok = True
    for row, expected, label in cases:
        got = p0(row)
        match = got == expected
        ok = ok and match
        print(f"  {'OK ' if match else 'FAIL'} {label:24s} p0={got}  expected={expected}")
    # curvature cross-check
    iv = I_B((0, 4, 6, 8, 10, 12, 14), 15, 16)
    print(f"  I_B(15,16) over HYP2675 core = {iv}  expected=-13/1470 -> "
          f"{'OK' if iv == Fr(-13, 1470) else 'FAIL'}")
    print(f"  all base p0 cross-checks pass: {ok}")
    print()


def claim_a() -> None:
    print("== CLAIM (a): consecutive-pair curvature I_B(u,u+1), B=(0,1,2,3,4,5,6,7) ==")
    print("  Does it SATURATE to ~0.0139 (bounded) or DIVERGE as u grows?")
    B = (0, 1, 2, 3, 4, 5, 6, 7)
    print(f"  Phi_2(B) = {fmt(Phi_2(B))}")
    print("    u    | I_B(u,u+1) exact            | decimal       | I_B - Phi_2")
    prev = None
    vals = []
    for u in (20, 80, 320, 1280, 5120):
        iv = I_B(B, u, u + 1)
        vals.append(iv)
        dev = iv - Phi_2(B)
        print(f"  {u:5d}  | {str(iv):27s} | {float(iv):.9f} | {fmt(dev)}")
        prev = iv
    # is it converging / bounded?
    mx = max(abs(v) for v in vals)
    print(f"  max |I_B(u,u+1)| over tested u = {fmt(mx)}")
    print(f"  claimed saturation value ~0.0139.  Last value decimal = {float(vals[-1]):.9f}")
    print(f"  sequence of decimals = {[round(float(v),9) for v in vals]}")
    print()


def claim_b() -> None:
    print("== CLAIM (b): |I_B(u,v) - Phi_2(B)| * resdist(u,v) bounded (~0.01)? ==")
    B = (0, 1, 2, 3, 4, 5, 6, 7)
    print(f"  base B={B}, Phi_2(B)={fmt(Phi_2(B))}")
    # mix of resonant (consecutive => resdist=1), small-resonance, and independent pairs
    pairs = [
        # consecutive (resdist=1)
        (20, 21), (80, 81), (320, 321), (1280, 1281), (5120, 5121),
        # near-resonant 2u ~ v (m=2,n=-1 => resdist small) e.g. v=2u
        (20, 40), (80, 160), (320, 640),
        # v = 2u+1
        (20, 41), (80, 161),
        # 3u-v small: v=3u
        (20, 60), (80, 240),
        # "independent"-ish large coprime gaps (resdist should be larger)
        (20, 33), (101, 211), (113, 251), (127, 311),
        (211, 499), (307, 701), (401, 907),
    ]
    print("    (u,v)        | resdist | I_B(u,v)-Phi_2 exact        | |dev|        | product |dev|*resdist")
    max_prod = Fr(0)
    max_pair = None
    for u, v in pairs:
        dev = I_B(B, u, v) - Phi_2(B)
        rd = resdist(u, v)
        prod = abs(dev) * rd
        if prod > max_prod:
            max_prod = prod
            max_pair = (u, v)
        print(f"  ({u:5d},{v:5d}) | {rd:7d} | {str(dev):27s} | {float(abs(dev)):.8f} | "
              f"{fmt(prod)}")
    print(f"  MAX product = {fmt(max_prod)} at pair {max_pair}")
    print()


def claim_b_resonance_sweep() -> None:
    """Adversarial: hunt for resonances that could blow up the product.

    Try to maximize |dev| for fixed small resdist by sweeping u with v chosen to
    create a target small relation, and also try other bases B.
    """
    print("== CLAIM (b) ADVERSARIAL: maximize product over resonant families ==")
    bases = [
        (0, 1, 2, 3, 4, 5, 6, 7),
        (0, 1, 2, 4, 8),
        (0, 2, 4, 6, 8, 10, 12, 14),  # dilated
        (0, 1, 3, 5, 7, 9, 11, 13),
        (0,),
    ]
    for B in bases:
        Bs = tuple(sorted(set(B)))
        ph2 = Phi_2(Bs)
        worst = Fr(0)
        worst_info = None
        # sweep consecutive (resdist=1) deep
        for u in (50, 200, 800, 3200):
            v = u + 1
            dev = I_B(Bs, u, v) - ph2
            rd = resdist(u, v)
            prod = abs(dev) * rd
            if prod > worst:
                worst = prod
                worst_info = (u, v, rd, dev)
        # sweep v=2u (m=2,n=-1) resdist = min(...) likely small
        for u in (50, 200, 800):
            v = 2 * u
            dev = I_B(Bs, u, v) - ph2
            rd = resdist(u, v)
            prod = abs(dev) * rd
            if prod > worst:
                worst = prod
                worst_info = (u, v, rd, dev)
        if worst_info is None:
            print(f"  B={Bs}: Phi_2={float(ph2):.6f}  (no nonzero product)")
            continue
        u, v, rd, dev = worst_info
        print(f"  B={Bs}: Phi_2={float(ph2):.6f}  worst product={fmt(worst)} "
              f"at (u,v)=({u},{v}) resdist={rd} dev={float(dev):.6f}")
    print()


def claim_b_offresonance_growth() -> None:
    """CRUX test of claim (b): does the product grow with resdist for independent pairs?

    The claimed form is |I_B - Phi_2| <= C / resdist (i.e. product <= C).
    For genuinely dissociated (u,v), |I_B - Phi_2| should -> 0, but the claim is
    that resdist*|dev| stays bounded. Test: scale up a fixed RATIO pair so resdist
    grows linearly while keeping the additive geometry fixed, and watch the product.
    """
    print("== CLAIM (b) CRUX: off-resonance, does product grow with resdist? ==")
    B = (0, 1, 2, 3, 4, 5, 6, 7)
    ph2 = Phi_2(B)
    # Use coprime well-separated ratio pairs of increasing size. For a fixed
    # 'irrational-like' ratio v/u, resdist grows ~ linearly in u (no small relation
    # with |m|,|n|<=7), so if |dev| ~ 1/resdist the product is constant; if |dev|
    # decays SLOWER than 1/resdist the product GROWS (refutes bound).
    print("  fixed-ratio scaling family v ~ phi*u (golden, hardest to approximate):")
    print("    (u,v)          | resdist | |dev|        | product")
    fam = [(13, 21), (34, 55), (89, 144), (233, 377), (610, 987), (1597, 2584)]
    for u, v in fam:
        dev = I_B(B, u, v) - ph2
        rd = resdist(u, v)
        prod = abs(dev) * rd
        print(f"  ({u:5d},{v:6d}) | {rd:7d} | {float(abs(dev)):.8f} | {fmt(prod)}")
    print()
    print("  prime-ish coprime pairs with growing min-relation (v ~ 2.3 u):")
    print("    (u,v)          | resdist | |dev|        | product")
    fam2 = [(101, 233), (211, 487), (401, 929), (809, 1873), (1601, 3697)]
    for u, v in fam2:
        dev = I_B(B, u, v) - ph2
        rd = resdist(u, v)
        prod = abs(dev) * rd
        print(f"  ({u:5d},{v:6d}) | {rd:7d} | {float(abs(dev)):.8f} | {fmt(prod)}")
    print()


if __name__ == "__main__":
    cross_check_p0_against_canon()
    claim_a()
    claim_b()
    claim_b_resonance_sweep()
    claim_b_offresonance_growth()
