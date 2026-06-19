#!/usr/bin/env python3
"""
LRC(14) Euler-copy squarefree profiles.

User reframe:
  Create copies of natural numbers so that for each n,

      sum_{d|n} copies(d) = n.

Mobius inversion forces copies(n)=phi(n).  This script uses that identity as a
squarefree prime-mask ledger for the existing HYP-2625/2626/2627 LRC14 route.

The key question is not whether this proves LRC(14).  It does not.  The key
question is whether the copy ledger explains why the raw Hill product

      P_14 = 5*6*6*7 = 1260

is the proof-relevant carrier, while the divided crossing value cr(K_14)=315
loses information.

Tournament Analysis declaration:
  vertices are quotient/proof carriers, not runners:
    euler_copy_squarefree_profile, raw_hill_product, mod210_address,
    prime7_coimage_seam, markov_hurwitz_energy_surface, crossing_value,
    raw_divisor_lattice, raw_runner_vertices.
  Pairwise observable: preservation of the LRC14 predicate, dyadic gate,
    squarefree transfer, coimage compatibility, and recurrence value.
  Tie Hamiltonian path: the score order printed below.

Assumption challenge:
  I considered runners, divisors, individual divisor copies, squarefree masks,
  Hill blocks, the crossing number, Markov-Hurwitz coordinates, mod-7 coimage
  classes, and proof obligations.  The useful quotient is not the raw divisor
  lattice or the raw runners.  It is the Euler-copy mass grouped by squarefree
  masks, because this quotient preserves the mod-210 address while showing
  exactly where cr(K_14)=315 destroys the dyadic gate.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from operator import mul

PRIMES = (2, 3, 5, 7)
FULL_MASK = (1 << len(PRIMES)) - 1


def section(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def prod(vals) -> int:
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


def phi(n: int) -> int:
    result = n
    for p in factor(n):
        result = result // p * (p - 1)
    return result


def mobius(n: int) -> int:
    fac = factor(n)
    if any(e > 1 for e in fac.values()):
        return 0
    return -1 if len(fac) % 2 else 1


def divisors(n: int) -> list[int]:
    items = list(factor(n).items())
    out: list[int] = []

    def rec(i: int, cur: int) -> None:
        if i == len(items):
            out.append(cur)
            return
        p, e = items[i]
        x = 1
        for _ in range(e + 1):
            rec(i + 1, cur * x)
            x *= p

    rec(0, 1)
    return sorted(out)


def radical(n: int) -> int:
    return prod(factor(n).keys())


def mask_of(n: int) -> int:
    out = 0
    for i, p in enumerate(PRIMES):
        if n % p == 0:
            out |= 1 << i
    return out


def mask_name(mask: int) -> str:
    vals = [str(p) for i, p in enumerate(PRIMES) if mask & (1 << i)]
    return "{" + ",".join(vals) + "}" if vals else "{}"


def copy_profile(n: int) -> dict[int, int]:
    """Group sum_{d|n} phi(d) by the tracked squarefree prime mask of d."""
    prof = {mask: 0 for mask in range(1 << len(PRIMES))}
    for d in divisors(n):
        prof[mask_of(d)] += phi(d)
    return prof


def formula_profile(n: int) -> dict[int, int]:
    """
    Closed profile formula.

    If a tracked prime p has exponent a in n, its included-mask mass is p^a-1.
    Untracked primes are neutral here: all their divisor-copy mass is summed
    into every compatible tracked mask, contributing q^a.
    """
    fac = factor(n)
    neutral = prod(p**e for p, e in fac.items() if p not in PRIMES)
    prof: dict[int, int] = {}
    for mask in range(1 << len(PRIMES)):
        val = neutral
        ok = True
        for i, p in enumerate(PRIMES):
            e = fac.get(p, 0)
            if mask & (1 << i):
                if e == 0:
                    ok = False
                    break
                val *= p**e - 1
        prof[mask] = val if ok else 0
    return prof


def profile_line(n: int, mask: int = FULL_MASK) -> str:
    prof = copy_profile(n)
    full = prof[mask]
    frac = Fraction(full, n)
    return (
        f"N={n:<6} factor={fmt_factor(n):<18} phi={phi(n):<5} "
        f"rad={radical(n):<5} copy{mask_name(mask)}={full:<6} "
        f"copy/N={frac}"
    )


def hill_tuple(n: int) -> tuple[int, int, int, int]:
    return tuple(sorted((n // 2, (n - 1) // 2, (n - 2) // 2, (n - 3) // 2)))


def hill_product(n: int) -> int:
    return prod(hill_tuple(n))


def verify_copy_rule(limit: int = 60) -> None:
    section("COPY RULE IS EULER TOTIENT")
    print("Let c(n) satisfy sum_{d|n} c(d)=n.  Mobius inversion gives")
    print("  c(n) = sum_{d|n} mu(d)*(n/d),")
    print("which is phi(n).  Exact check:")
    ok = True
    for n in range(1, limit + 1):
        inv = sum(mobius(d) * (n // d) for d in divisors(n))
        divsum = sum(phi(d) for d in divisors(n))
        if inv != phi(n) or divsum != n:
            ok = False
            print(f"failure at n={n}: inv={inv}, phi={phi(n)}, divsum={divsum}")
            break
    print(f"verified c(n)=phi(n) and sum phi(d)=n for n<= {limit}: {ok}")
    print("First copy counts:")
    print("  " + ", ".join(f"{n}:{phi(n)}" for n in range(1, 21)))


def prime_extension_report() -> None:
    section("SQUAREFREE COPY RECURRENCE")
    print("For squarefree S, grouped copy mass is w_S(M)=prod_{p in M}(p-1).")
    print("Adding a new prime p keeps old masks and adds a shifted layer")
    print("  w_{S*p}(M union {p}) = (p-1) w_S(M).")
    cur = 1
    cur_mask = 0
    for p in PRIMES:
        before = copy_profile(cur)
        after = copy_profile(cur * p)
        shifted = cur_mask | (1 << PRIMES.index(p))
        print(
            f"add p={p}: total {cur}->{cur*p}, "
            f"full {mask_name(shifted)} mass {after[shifted]} = "
            f"{p-1}*{before[cur_mask]}"
        )
        cur *= p
        cur_mask = shifted
    print()
    print("The mod-address chain is therefore:")
    for n, mask in [(6, 0b0011), (30, 0b0111), (210, 0b1111)]:
        print("  " + profile_line(n, mask))


def profile_table_report() -> None:
    section("PROFILE TABLE: RAW PRODUCT VS DIVIDED CROSSING VALUE")
    rows = [
        ("mod6 radical", 6, 0b0011),
        ("mod30 radical", 30, 0b0111),
        ("mod210 radical", 210, FULL_MASK),
        ("K13 Hill product", 900, FULL_MASK),
        ("K14 raw Hill product", 1260, FULL_MASK),
        ("K14 crossing value", 315, FULL_MASK),
        ("K15 Hill product", 1764, FULL_MASK),
    ]
    print(
        f"{'object':<22} {'N':>7} {'factor':>18} {'phi':>6} {'rad':>6} "
        f"{'full-mask copies':>16} {'copy/N':>10}"
    )
    for name, n, mask in rows:
        prof = copy_profile(n)
        full = prof[mask]
        print(
            f"{name:<22} {n:>7} {fmt_factor(n):>18} {phi(n):>6} "
            f"{radical(n):>6} {full:>16} {str(Fraction(full, n)):>10}"
        )
    print()
    print("K14 comparison:")
    print("  P_14=1260 has copy mass on {2,3,5,7}: 576.")
    print("  cr(K_14)=315 has copy mass on {2,3,5,7}: 0.")
    print("  Reason: dividing P_14 by 4 removes the entire tracked 2-adic gate.")
    print("  This is exactly the HYP-2627 warning: keep the raw product ledger first.")


def detailed_k14_profile() -> None:
    section("DETAILED K14 COPY PROFILE")
    n = 1260
    prof = copy_profile(n)
    formula = formula_profile(n)
    print(f"P_14={n}={fmt_factor(n)}")
    print(f"sum over all mask copy masses: {sum(prof.values())}")
    print(f"closed formula matches divisor enumeration: {prof == formula}")
    print(
        "For P_14=2^2*3^2*5*7, a mask containing p receives factor "
        "p^a-1 from that prime."
    )
    print(f"{'mask':>13} {'copy_mass':>10} {'copy/P_14':>12}")
    for mask in range(1 << len(PRIMES)):
        if prof[mask]:
            print(f"{mask_name(mask):>13} {prof[mask]:>10} {str(Fraction(prof[mask], n)):>12}")
    print()
    full_rad = copy_profile(210)[FULL_MASK]
    full_p14 = prof[FULL_MASK]
    print(
        f"The radical profile has full-mask mass phi(210)={full_rad}; "
        f"P_14 thickens it to {full_p14}."
    )
    print(
        "The multiplier is ((2^2-1)/(2-1))*((3^2-1)/(3-1)) = 3*4 = 12, "
        "coming from the repeated 6 blocks."
    )
    print(
        f"By contrast phi(P_14)={phi(n)}={Fraction(phi(n), phi(210))}*phi(210); "
        "totient itself sees only the p-fold exponent thickening, while the "
        "full squarefree copy atom sees all nonzero p-power divisor layers."
    )


def hill_scan_report() -> None:
    section("HILL PRODUCT SCAN")
    print(
        f"{'n':>2} {'q_n':>17} {'P_n':>8} {'rad(P)':>8} "
        f"{'copy210(P)':>12} {'copy/P':>9} {'copy210(cr)':>12}"
    )
    first_full = None
    first_raw_not_crossing = None
    for n in range(5, 33):
        q = hill_tuple(n)
        p = prod(q)
        prof_p = copy_profile(p)
        full_p = prof_p[FULL_MASK]
        cr = p // 4
        full_cr = copy_profile(cr)[FULL_MASK] if p % 4 == 0 else 0
        if first_full is None and full_p:
            first_full = n
        if first_raw_not_crossing is None and full_p and not full_cr:
            first_raw_not_crossing = n
        print(
            f"{n:>2} {str(q):>17} {p:>8} {radical(p):>8} "
            f"{full_p:>12} {str(Fraction(full_p, p)):>9} {full_cr:>12}"
        )
    print()
    print(f"first Hill row with nonzero full {{2,3,5,7}} copy mass: n={first_full}")
    print(
        "first row where the raw product has the full mask but the divided "
        f"crossing value loses it: n={first_raw_not_crossing}"
    )


def hurwitz_copy_report() -> None:
    section("MARKOV-HURWITZ IN COPY LANGUAGE")
    q = hill_tuple(14)
    product_copies = prod(sum(phi(d) for d in divisors(x)) for x in q)
    square_energy = sum(x * x for x in q)
    primitive_blocks = tuple(phi(x) for x in q)
    primitive_product = prod(primitive_blocks)
    primitive_square_energy = sum(x * x for x in primitive_blocks)
    print(f"K14 Hill tuple q={q}")
    print(
        "Since x=sum_{d|x} phi(d), the product wxyz is the number of four-block "
        "Euler-copy assignments."
    )
    print(f"  product copy assignments wxyz = {product_copies}")
    print(f"  Markov-Hurwitz diagonal energy sum q_i^2 = {square_energy}")
    print(f"  product minus diagonal energy = {product_copies - square_energy}")
    print()
    print(f"Primitive block counts phi(q_i)={primitive_blocks}")
    print(f"  primitive product = {primitive_product}")
    print(f"  primitive square energy = {primitive_square_energy}")
    print(
        "Readout: the totient-copy frame does not make q_14 Markov-Hurwitz. "
        "It explains the RHS product as a copy-assignment ledger, and it shows "
        "which squarefree assignment layers survive projection."
    )


@dataclass(frozen=True)
class Route:
    name: str
    lrc_predicate: int
    dyadic_gate: int
    squarefree_transfer: int
    coimage_compatibility: int
    recurrence_value: int

    @property
    def score(self) -> tuple[int, int, int, int, int]:
        return (
            self.lrc_predicate,
            self.dyadic_gate,
            self.squarefree_transfer,
            self.coimage_compatibility,
            self.recurrence_value,
        )


def tournament_report() -> None:
    section("TOURNAMENT ANALYSIS")
    routes = [
        Route("euler_copy_squarefree_profile", 5, 5, 5, 4, 5),
        Route("raw_hill_product", 4, 5, 4, 3, 3),
        Route("mod210_address", 4, 4, 5, 4, 4),
        Route("prime7_coimage_seam", 4, 2, 3, 5, 3),
        Route("markov_hurwitz_energy_surface", 2, 2, 2, 1, 5),
        Route("crossing_value", 2, 0, 2, 1, 2),
        Route("raw_divisor_lattice", 1, 3, 2, 1, 1),
        Route("raw_runner_vertices", 1, 1, 1, 1, 1),
    ]
    wins = Counter()
    for a in routes:
        wins[a.name] += 0
        for b in routes:
            if a != b and a.score > b.score:
                wins[a.name] += 1

    cycles = 0
    for a in routes:
        for b in routes:
            for c in routes:
                if len({a, b, c}) != 3:
                    continue
                if a.score > b.score and b.score > c.score and c.score > a.score:
                    cycles += 1
    ordered = sorted(routes, key=lambda r: r.score, reverse=True)
    print("Hamiltonian proof path:")
    print("  " + " > ".join(r.name for r in ordered))
    print(f"score histogram: {dict(sorted(Counter(wins.values()).items()))}")
    print(f"directed 3-cycles: {cycles // 3}")
    print(
        "The quotient preserves the proof-relevant mod-210 copy address.  It "
        "destroys exact divisor representatives and witness-time geometry, which "
        "must be recovered only after the finite wall/coimage ledgers are fixed."
    )


def main() -> None:
    print("LRC14 Euler-copy squarefree profile probe - codex S21")
    verify_copy_rule()
    prime_extension_report()
    profile_table_report()
    detailed_k14_profile()
    hill_scan_report()
    hurwitz_copy_report()
    tournament_report()
    section("SESSION READOUT")
    print(
        "The divisor-copy reframe is exactly Euler totient by Mobius inversion. "
        "Its useful LRC14 form is the squarefree copy profile "
        "sum_{d|N, mask(d)=M} phi(d)."
    )
    print(
        "For K14, raw P_14=1260 has full {2,3,5,7} copy mass 576, while "
        "cr(K_14)=315 has full copy mass 0.  The division by 4 deletes the "
        "dyadic gate, so HYP-2627's raw-product discipline is necessary."
    )
    print(
        "The recurrence hidden in the user's reframe is prime extension: adding "
        "p appends a shifted copy layer multiplied by p-1.  Thus "
        "{2,3}->{2,3,5}->{2,3,5,7} is a genuine Euler-copy recurrence, not just "
        "a list of moduli."
    )


if __name__ == "__main__":
    main()
