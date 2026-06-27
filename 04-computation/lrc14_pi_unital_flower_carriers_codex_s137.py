#!/usr/bin/env python3
"""S137: pi approximants, flower families, and unital carriers for LRC14.

The prompt adds three signals to the current LRC14 Farey/product/shell route:

* 22/7 is a good approximation to pi.
* cbrt(31) is also a good approximation to pi.
* turning by 1/pi radians before each petal gives 22 visible families.

This script keeps those signals as quotients with declared loss.  The 22-family
claim is real in the unit-radian residue quotient because 1/pi ~= 7/22.  It is
not a true full-circle closure, whose rotation number is 1/(2*pi^2).

The script also separates geometric unitals, algebraic unital/unit-preserving
maps, and unit groups.  The useful LRC14 connection is q=3 geometric unital
parameters: (q^3+1,q+1,1)=(28,4,1), so C=27 is one point short of the unital
point count and 31=q^3+(q+1) is the cubic package behind cbrt(31).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, getcontext
from fractions import Fraction
from itertools import combinations


getcontext().prec = 80

PI = Decimal(
    "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899"
)
N = 14
C = 2 * N - 1


def cbrt_decimal(n: int) -> Decimal:
    x = Decimal(str(n ** (1 / 3)))
    Ndec = Decimal(n)
    for _ in range(40):
        x = (2 * x + Ndec / (x * x)) / 3
    return +x


def dec(frac: Fraction) -> Decimal:
    return Decimal(frac.numerator) / Decimal(frac.denominator)


def continued_fraction(x: Decimal, terms: int = 12) -> list[int]:
    out: list[int] = []
    y = +x
    for _ in range(terms):
        a = int(y)
        out.append(a)
        y -= a
        if y == 0:
            break
        y = 1 / y
    return out


def best_fraction(x: Decimal, max_den: int) -> tuple[Fraction, Decimal]:
    best: tuple[Fraction, Decimal] | None = None
    for q in range(1, max_den + 1):
        center = int((x * q).to_integral_value(rounding="ROUND_HALF_UP"))
        for p in (center - 1, center, center + 1):
            f = Fraction(p, q)
            err = abs(dec(f) - x)
            if best is None or err < best[1]:
                best = (f, err)
    assert best is not None
    return best


def circular_distance_mod1(a: Decimal, b: Decimal) -> Decimal:
    d = abs(a - b) % 1
    return min(d, Decimal(1) - d)


def flower_family_stats() -> dict[str, object]:
    theta = 1 / PI
    rational = Fraction(7, 22)
    per_step_error = theta - dec(rational)
    deviations = []
    for k in range(22):
        actual = (theta * Decimal(k)) % 1
        model = dec(Fraction((7 * k) % 22, 22))
        deviations.append(circular_distance_mod1(actual, model))
    rho = 1 / (2 * PI * PI)
    return {
        "theta": theta,
        "radian_model": rational,
        "per_step_error": per_step_error,
        "twenty_two_step_drift_from_7_radians": Decimal(22) * theta - Decimal(7),
        "max_family_deviation_mod_1_radian": max(deviations),
        "circle_rotation_number": rho,
        "circle_best_den_le_22": best_fraction(rho, 22),
        "circle_best_den_le_100": best_fraction(rho, 100),
        "circle_best_den_le_300": best_fraction(rho, 300),
        "circle_best_den_le_1000": best_fraction(rho, 1000),
        "circle_model_from_22_7": Fraction(49, 968),
    }


def unital_parameters(q: int) -> dict[str, int]:
    v = q**3 + 1
    k = q + 1
    lam = 1
    b = v * (v - 1) * lam // (k * (k - 1))
    r = lam * (v - 1) // (k - 1)
    return {"q": q, "v": v, "k": k, "lambda": lam, "blocks": b, "rep": r}


@dataclass(frozen=True)
class Carrier:
    name: str
    role: str
    score: tuple[int, int, int, int, int, int]


def carriers() -> list[Carrier]:
    return [
        Carrier(
            "exact_M_Farey_branch",
            "keeps q and e=14p-q as theorem scale",
            (6, 5, 5, 2, 5, 5),
        ),
        Carrier(
            "marked_C27_shell_transfer",
            "unit/nonunit shell holes and doubles from HYP-2937",
            (5, 5, 5, 4, 4, 5),
        ),
        Carrier(
            "bigraded_relation_signature",
            "visible and hidden summand/multiplicand channels",
            (5, 5, 4, 4, 4, 5),
        ),
        Carrier(
            "Kpq_K33_incidence",
            "complete-bipartite owner packet for p>=3",
            (5, 4, 3, 4, 5, 4),
        ),
        Carrier(
            "geometric_unital_28_4_1",
            "q=3 pair-unique design with 28=2*14 points",
            (4, 4, 5, 5, 4, 5),
        ),
        Carrier(
            "algebraic_unital_maps",
            "identity/unit-preserving quotient discipline",
            (4, 4, 5, 3, 4, 5),
        ),
        Carrier(
            "octa_Clebsch_packet",
            "support-six current and folded covariance carriers",
            (4, 4, 3, 4, 5, 4),
        ),
        Carrier(
            "PZ_degree4_gateway",
            "second moment only after upgrade to labelled degree-4 region",
            (3, 3, 2, 3, 3, 4),
        ),
        Carrier(
            "flower_22_radian_quotient",
            "1/pi ~= 7/22 gives 22 near families modulo one radian",
            (2, 3, 3, 2, 2, 5),
        ),
        Carrier(
            "cuberoot31_cubic_shell",
            "pi^3 ~= 31 = 27+4, a q=3 unital package mnemonic",
            (2, 2, 3, 2, 2, 5),
        ),
    ]


def role_score(c: Carrier) -> tuple[int, int, int, int, int, int]:
    return c.score


def tournament_fingerprint(items: list[Carrier]) -> dict[str, object]:
    wins = {c.name: 0 for c in items}
    c3 = 0
    for a, b in combinations(items, 2):
        if role_score(a) >= role_score(b):
            wins[a.name] += 1
        else:
            wins[b.name] += 1
    for a, b, c in combinations(items, 3):
        ab = role_score(a) >= role_score(b)
        bc = role_score(b) >= role_score(c)
        ca = role_score(c) >= role_score(a)
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            c3 += 1
    order = sorted(items, key=role_score, reverse=True)
    return {
        "score_hist": dict(sorted(Counter(wins.values()).items())),
        "directed_3cycles": c3,
        "hamiltonian_order": [c.name for c in order],
    }


def print_decimal(x: Decimal, digits: int = 18) -> str:
    quant = Decimal(10) ** -digits
    return str(x.quantize(quant))


def main() -> None:
    c31 = cbrt_decimal(31)
    approx_rows = [
        ("22/7", dec(Fraction(22, 7)), "continued-fraction/radian-family"),
        ("cuberoot(31)", c31, "cubic C27+4 / unital mnemonic"),
        ("355/113", dec(Fraction(355, 113)), "control: strong pi convergent"),
    ]

    print("S137 LRC14 PI / UNITAL / FLOWER CARRIER SYNTHESIS")
    print("=" * 78)
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    petals, radian residues, full-circle phases, pi approximants,")
    print("    q=3 unital points/blocks, algebraic unit maps, C27 shells,")
    print("    Kpq owners, support-six packets, and proof obligations.")
    print("  chosen quotient:")
    print("    proof carriers with exact notes on what the pi/unital quotients forget.")
    print("  preserved predicate:")
    print("    whether the carrier keeps LRC14 branch, unit identity, pair coverage,")
    print("    or state-lift owner data.")
    print("  destroyed predicate:")
    print("    raw safe-time geometry; exact M/Farey branch remains attached.")

    print("\n[1] Pi approximants")
    print(f"  pi = {print_decimal(PI, 30)}")
    print(f"  pi^3 - 31 = {print_decimal(PI**3 - Decimal(31), 30)}")
    print("  name             value                 signed error          abs error       role")
    for name, value, role in approx_rows:
        err = value - PI
        print(
            f"  {name:<15} {print_decimal(value, 18):>20} "
            f"{print_decimal(err, 18):>21} {print_decimal(abs(err), 18):>16}  {role}"
        )
    print("  readout:")
    print("    cuberoot(31) is numerically closer to pi than 22/7 in absolute error,")
    print("    but 22/7 is the quotient that explains the 22-family flower signal.")

    print("\n[2] Flower from turning 1/pi radians")
    fstats = flower_family_stats()
    print(f"  theta=1/pi = {print_decimal(fstats['theta'], 30)} radians")
    print("  radian-residue quotient:")
    print("    1/pi ~= 7/22, so the model step on Z/22 is +7.")
    print("    gcd(7,22)=1, hence all 22 radian-residue families are visited.")
    print(
        "    after 22 petals, actual angle differs from 7 radians by "
        f"{print_decimal(fstats['twenty_two_step_drift_from_7_radians'], 24)}."
    )
    print(
        "    max deviation from the 7/22 family grid over the first 22 petals is "
        f"{print_decimal(fstats['max_family_deviation_mod_1_radian'], 24)}."
    )
    print("  full-circle quotient guardrail:")
    print(
        "    rotation number is theta/(2*pi)=1/(2*pi^2)="
        f"{print_decimal(fstats['circle_rotation_number'], 30)}."
    )
    for key in (
        "circle_best_den_le_22",
        "circle_best_den_le_100",
        "circle_best_den_le_300",
        "circle_best_den_le_1000",
    ):
        frac, err = fstats[key]
        print(f"    best {key.removeprefix('circle_best_den_'):>7}: {frac} err={print_decimal(err, 24)}")
    print(f"    substituting 22/7 gives full-circle model {fstats['circle_model_from_22_7']}.")
    print("  readout:")
    print("    the 22 families are a real radian-address quotient, not a true")
    print("    circle-closing period.  This is analogous to C27 shell transfer:")
    print("    a useful quotient is allowed only after its forgotten geometry is named.")

    print("\n[3] Geometric and algebraic unital distinctions")
    for q in (2, 3, 4):
        p = unital_parameters(q)
        special = ""
        if q == 3:
            special = "  <-- v=28=2*14, v-1=27=C, k=4, b=63=7*9"
        print(
            f"  q={q}: 2-({p['v']},{p['k']},{p['lambda']}) "
            f"blocks={p['blocks']} rep={p['rep']}{special}"
        )
    print("  q=3 package:")
    print("    q^3+1=28 equals the 2*14 antipodal/apex point count.")
    print("    q^3=27 equals the C=2*14-1 shell in HYP-2937.")
    print("    q+1=4 is the block size, and q^3+(q+1)=31.")
    print("    Thus cuberoot(31) is a cubic-shell mnemonic: 31=27+4.")
    print("  algebraic unital guardrail:")
    print("    geometric unitals, unit groups, and unital maps are different.")
    print("    For LRC14, the useful algebraic rule is unit preservation:")
    print("    quotient maps should send the identity/floor object to the")
    print("    identity carrier rather than erasing AP/GW or unit-visible shells.")

    print("\n[4] Carrier Tournament Analysis")
    cs = carriers()
    fp = tournament_fingerprint(cs)
    print("  vertices: proof carriers, not runners.")
    print("  observable:")
    print("    (branch retention, typed visibility, unit preservation, pair coverage,")
    print("     state-lift fit, scalar-decoy resistance).")
    print("  switch/gauge: lexicographically larger role score wins.")
    for c in sorted(cs, key=role_score, reverse=True):
        print(f"    {c.name:<30} score={c.score} role={c.role}")
    print(
        f"  fingerprint: score_hist={fp['score_hist']} "
        f"c3={fp['directed_3cycles']}"
    )
    print("  Hamiltonian carrier order:")
    print("    " + " > ".join(fp["hamiltonian_order"]))

    print("\n[5] Proof readout")
    print("  Add the pi/unital ideas as conservative carriers, not proof scalars.")
    print("  The flower's 22-family quotient is useful because 1/pi ~= 7/22")
    print("  and +7 is a unit step on Z/22; it should be compared to unit-visible")
    print("  C27 shell holes, not to a literal full-circle period.")
    print("  The cuberoot(31) signal points to the q=3 unital package")
    print("  31=27+4, linking C=27 with a block size 4 pair-design carrier.")
    print("  The next theorem target is unchanged but better typed:")
    print("    after exact M/Farey branch and C27 transfer are attached,")
    print("    try a q=3 unital/pair-unique packet as a finite pair-coverage")
    print("    interface for the p>=3 K33 owner branch or HYP-2908 state lift.")
    print("  Guardrail:")
    print("    none of 22/7, cuberoot(31), or the word 'unital' replaces q,")
    print("    exact M(S), or the C27/K33 branch labels.")


if __name__ == "__main__":
    main()
