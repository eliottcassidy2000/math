#!/usr/bin/env python3
"""
LRC(14) squarefree divisor profiles, Markov-Hurwitz, and K_n crossings.

The user's proposed bridge has three pieces:

  1. squarefree divisor profiles;
  2. the Markov-Hurwitz surface
       w^2 + x^2 + y^2 + z^2 = wxyz;
  3. the complete-graph crossing formula
       cr(K_n) = floor(n/2) floor((n-1)/2) floor((n-2)/2) floor((n-3)/2) / 4.

For a four-block tuple q=(w,x,y,z), the Hill/Zarankiewicz product q0 q1 q2 q3
is four times the crossing count.  The Markov-Hurwitz equation is exactly the
energy-critical surface where that four-block crossing product equals the
quadratic diagonal energy sum q_i^2.  This script tests which part of that
analogy is useful for the LRC14 modular/coimage thread.

Tournament Analysis declaration:
  vertices are proof quotients, not runners:
    squarefree_210_crossing_profile, lrc14_prime7_coimage_seam,
    markov_hurwitz_mutation, mod30_address, raw_crossing_number,
    raw_hill_blocks, raw_runner_vertices.
  Pairwise observable: (preserves LRC14 predicate, squarefree transfer value,
    recurrence legitimacy, crossing/Hurwitz explanatory value).
  Tie path: the listed order.

Assumption challenge:
  I considered complete-graph vertices, crossings, Hill four-block factors,
  Markov-Hurwitz coordinates, prime masks, coimage classes, and proof
  obligations.  The useful quotient is not raw crossings or raw runners.  It is
  the squarefree radical of the four-block product, because for n=14 the Hill
  product is 1260 = 6*210 and its radical is exactly the mod-210 address from
  HYP-2625/HYP-2626.
"""
from __future__ import annotations

import math
from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from operator import mul

PRIMES = (2, 3, 5, 7)


def section(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def prod(vals: tuple[int, ...] | list[int]) -> int:
    return reduce(mul, vals, 1)


def factor(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def fmt_factor(n: int) -> str:
    if n == 1:
        return "1"
    pieces = []
    for p, e in factor(n).items():
        pieces.append(str(p) if e == 1 else f"{p}^{e}")
    return "*".join(pieces)


def radical(n: int) -> int:
    r = 1
    for p in factor(n):
        r *= p
    return r


def mask(n: int) -> int:
    m = 0
    for i, p in enumerate(PRIMES):
        if n % p == 0:
            m |= 1 << i
    return m


def mask_name(m: int) -> str:
    vals = [str(p) for i, p in enumerate(PRIMES) if m & (1 << i)]
    return "{" + ",".join(vals) + "}" if vals else "{}"


def tuple_masks(q: tuple[int, ...]) -> str:
    return "[" + ", ".join(mask_name(mask(x)) for x in q) + "]"


def hill_tuple(n: int) -> tuple[int, int, int, int]:
    return tuple(
        sorted(
            (
                n // 2,
                (n - 1) // 2,
                (n - 2) // 2,
                (n - 3) // 2,
            )
        )
    )


def hurwitz_defect(q: tuple[int, int, int, int]) -> int:
    return prod(q) - sum(x * x for x in q)


def hurwitz_pressure(q: tuple[int, int, int, int]) -> Fraction | None:
    p = prod(q)
    if p == 0:
        return None
    return Fraction(sum(x * x for x in q), p)


def hill_report() -> None:
    section("HILL FOUR-BLOCK CROSSING PRODUCTS AND SQUAREFREE RADICALS")
    print(
        "For q_n=(floor((n-3)/2),...,floor(n/2)), Hill's product P_n=prod(q_n) "
        "satisfies cr(K_n)=P_n/4.  Markov-Hurwitz pressure is sum(q_i^2)/P_n."
    )
    print(
        f"{'n':>2} {'q_n':>17} {'P_n=4cr':>9} {'cr':>8} {'rad(P_n)':>9} "
        f"{'P/rad':>7} {'pressure':>12} {'defect':>8} {'masks':>32}"
    )
    first_rad_210 = None
    first_product_1260 = None
    for n in range(5, 23):
        q = hill_tuple(n)
        p = prod(q)
        rad = radical(p)
        pressure = hurwitz_pressure(q)
        if first_rad_210 is None and rad % 210 == 0:
            first_rad_210 = n
        if first_product_1260 is None and p == 1260:
            first_product_1260 = n
        cr = Fraction(p, 4)
        print(
            f"{n:>2} {str(q):>17} {p:>9} {str(cr):>8} {rad:>9} "
            f"{p // rad:>7} {str(pressure):>12} {hurwitz_defect(q):>8} "
            f"{tuple_masks(q):>32}"
        )
    print()
    print(f"first n with rad(P_n) divisible by 210: {first_rad_210}")
    print(f"first n with P_n=1260: {first_product_1260}")
    q14 = hill_tuple(14)
    print(
        "K_14 readout: q_14="
        f"{q14}, P_14={prod(q14)}={fmt_factor(prod(q14))}, "
        f"cr(K_14)={Fraction(prod(q14), 4)}, rad(P_14)={radical(prod(q14))}."
    )
    print(
        "This is the compact arithmetic packet behind the mod hierarchy: "
        "two 6-blocks give the mod-6 skeleton, the 5-block adds mod 30, "
        "and the 7-block turns the radical into 210."
    )


def is_hurwitz_solution(q: tuple[int, int, int, int]) -> bool:
    return sum(x * x for x in q) == prod(q)


def generate_hurwitz(max_coord: int) -> list[tuple[int, int, int, int]]:
    """Generate the positive Markov-Hurwitz orbit from (2,2,2,2)."""
    seed = (2, 2, 2, 2)
    seen = {seed}
    q: deque[tuple[int, int, int, int]] = deque([seed])
    out = []
    while q:
        cur = q.popleft()
        out.append(cur)
        for i in range(4):
            rest = [cur[j] for j in range(4) if j != i]
            nxt = prod(rest) - cur[i]
            if nxt <= 0 or nxt > max_coord:
                continue
            new = tuple(sorted(rest + [nxt]))
            if new not in seen:
                if not is_hurwitz_solution(new):
                    raise AssertionError(new)
                seen.add(new)
                q.append(new)
    return sorted(out, key=lambda t: (sum(t), t))


def hurwitz_report() -> None:
    section("MARKOV-HURWITZ MUTATION TREE AND PRIME MASKS")
    sols = generate_hurwitz(100_000_000)
    print(f"generated positive solutions with max coordinate <= 1e8: {len(sols)}")
    print(
        "All generated solutions are even; dividing by 2 gives the normalized "
        "surface a^2+b^2+c^2+d^2 = 4abcd.  The mutation is q_i -> product(other three)-q_i."
    )
    print(
        f"{'solution':>30} {'sum':>9} {'rad(prod)':>12} {'pressure':>8} {'masks after /2':>34}"
    )
    for s in sols[:18]:
        half = tuple(x // 2 for x in s)
        print(
            f"{str(s):>30} {sum(s):>9} {radical(prod(s)):>12} "
            f"{str(hurwitz_pressure(s)):>8} {tuple_masks(half):>34}"
        )

    any_5 = [s for s in sols if any(x % 5 == 0 for x in s)]
    any_7 = [s for s in sols if any(x % 7 == 0 for x in s)]
    print()
    print(f"solutions with some coordinate divisible by 5: {len(any_5)}")
    if any_5:
        print(f"  first: {any_5[0]}")
    print(f"solutions with some coordinate divisible by 7: {len(any_7)}")
    if any_7:
        print(f"  first: {any_7[0]}")

    branch = []
    for s in sols:
        if s[0] == 2 and s[1] == 2:
            branch.append(s[2] // 2)
            branch.append(s[3] // 2)
    branch = sorted(set(branch))[:10]
    rec_ok = all(branch[i + 1] == 4 * branch[i] - branch[i - 1] for i in range(1, len(branch) - 1))
    print()
    print(f"normalized (1,1,u,v) branch values: {branch}")
    print(f"recurrence u_(j+1)=4u_j-u_(j-1) holds on shown branch: {rec_ok}")
    print(
        "Readout: Markov-Hurwitz gives a genuine integer recurrence, but its "
        "small-prime profile is not the LRC14 mod-210 carrier.  In this orbit, "
        "prime 5 is absent through the generated range, while the K_14 crossing "
        "tuple needs the 5-block and the 7-block simultaneously."
    )


def telescope_report() -> None:
    section("THM-523 TWO-SPEED CLASH TELESCOPE")
    p14 = prod(hill_tuple(14))
    half = Fraction(15, 36) - Fraction(2, 5) - Fraction(1, 70) - Fraction(1, 504)
    print("The known LRC14 low-measure champion has half-gap")
    print("  15/36 - 2/5 - 1/70 - 1/504")
    print(f"  = {half} = 1/(2*{p14}).")
    print()
    print("The cancellation exposes the same four blocks:")
    print("  15/36 - 2/5 = 5/12 - 2/5 = 1/60")
    print("  1/60 - 1/70 = 1/420 = 1/(5*6*14)")
    print("  1/420 - 1/504 = 1/(5*6*14) - 1/(6*6*14)")
    print("                    = 1/(5*6*6*14)")
    print("                    = 1/(2*7*6*6*5).")
    print(f"After tau <-> 1-tau symmetry this is {2 * half} = 1/{p14}.")
    print(
        "Thus the denominator 1260 is not only the Hill K_14 product and not only "
        "the HYP-2625 squarefree carrier; it is exactly the THM-523 local "
        "two-speed-clash denominator after symmetry."
    )


@dataclass(frozen=True)
class Route:
    name: str
    preserves_lrc: int
    squarefree_transfer: int
    recurrence_legitimacy: int
    crossing_value: int

    @property
    def score(self) -> tuple[int, int, int, int]:
        return (
            self.preserves_lrc,
            self.squarefree_transfer,
            self.recurrence_legitimacy,
            self.crossing_value,
        )


def tournament_report() -> None:
    section("TOURNAMENT ANALYSIS: WHICH QUOTIENT SHOULD CARRY THE ANALOGY?")
    routes = [
        Route("squarefree_210_crossing_profile", 4, 4, 3, 4),
        Route("lrc14_prime7_coimage_seam", 4, 3, 4, 3),
        Route("markov_hurwitz_mutation", 2, 2, 4, 3),
        Route("mod30_address", 2, 3, 3, 2),
        Route("raw_crossing_number", 1, 1, 1, 4),
        Route("raw_hill_blocks", 1, 2, 1, 3),
        Route("raw_runner_vertices", 1, 1, 1, 1),
    ]
    wins = Counter()
    cycles = 0
    for a in routes:
        wins[a.name] += 0
        for b in routes:
            if a == b:
                continue
            if a.score > b.score:
                wins[a.name] += 1
    for a in routes:
        for b in routes:
            for c in routes:
                if len({a, b, c}) != 3:
                    continue
                ab = a.score > b.score
                bc = b.score > c.score
                ca = c.score > a.score
                if ab and bc and ca:
                    cycles += 1
    ordered = sorted(routes, key=lambda r: r.score, reverse=True)
    print("Hamiltonian proof path:")
    print("  " + " > ".join(r.name for r in ordered))
    print(f"score histogram: {dict(sorted(Counter(wins.values()).items()))}")
    print(f"directed 3-cycles: {cycles // 3}")
    print(
        "Interpretation: the Markov-Hurwitz equation is the right recurrence "
        "metaphor, but not the right carrier.  The carrier is the squarefree "
        "radical of the Hill four-block product, because K_14 is the first Hill "
        "tuple whose product radical is divisible by 210."
    )


def main() -> None:
    print("LRC14 squarefree / Markov-Hurwitz / complete-graph crossing probe - codex S20")
    hill_report()
    hurwitz_report()
    telescope_report()
    tournament_report()
    section("SESSION READOUT")
    print(
        "The relation is not that the K_14 crossing tuple solves "
        "w^2+x^2+y^2+z^2=wxyz.  It does not: its Hurwitz pressure is 73/630, "
        "well below the Markov pressure 1.  The useful fact is subtler: Hill's "
        "four-block product for K_14 is 1260, with radical 210, exactly the "
        "squarefree address HYP-2625 needed for mod 6 -> mod 30 -> mod 210."
    )
    print(
        "So Markov-Hurwitz supplies the mutation/recurrence archetype, while "
        "the complete-graph crossing formula supplies the four-block squarefree "
        "profile.  For LRC14, the latter is the live object: (5,6,6,7) packages "
        "5, the repeated 6 skeleton, and the prime-7 coimage seam in one product."
    )
    print(
        "The additional telescope check strengthens the bridge: THM-523's "
        "15/36-2/5-1/70-1/504 half-gap is exactly 1/(2*1260), hence the "
        "symmetrized low-measure scale is 1/(7*6*6*5)."
    )


if __name__ == "__main__":
    main()
