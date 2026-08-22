#!/usr/bin/env python3
"""Exact n<=101 slope screen and four explicit positive Pell rays.

This is a standard-library continuation of the finite computation accompanying
THM-3376.  The congruence screen is deliberately weaker than a local-global
classification: a survivor has only been checked modulo every prime p<=199.
In particular, the seven screen survivors outside ``RAY_SLOPES`` are not
claimed to be globally insoluble or everywhere locally soluble.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from math import gcd
from pathlib import Path


PRIME_LIMIT = 199
EXPECTED_CANDIDATE_COUNT = 528
EXPECTED_SURVIVORS = (
    (14, 23),
    (26, 29),
    (26, 47),
    (38, 47),
    (38, 53),
    (50, 53),
    (50, 71),
    (62, 95),
    (74, 95),
    (74, 101),
    (98, 101),
)
EXPECTED_OBSTRUCTION_COUNTS = {
    3: 395,
    5: 21,
    7: 44,
    11: 11,
    13: 7,
    17: 8,
    23: 2,
    29: 5,
    31: 4,
    37: 3,
    41: 2,
    47: 1,
    53: 2,
    59: 3,
    61: 2,
    71: 2,
    79: 1,
    83: 2,
    89: 2,
}

# Each row gives an integral point on W^2-D*u^2=C and a positive norm-one
# unit P+Q*sqrt(D).  All entries are independently checked below.
RAYS = {
    (14, 23): {
        "D": 44965,
        "C": 676,
        "W": 16715091502370452792679,
        "u": 78826357654129385469,
        "P": 1062835917426709745982924462536293407897885896493741491943228169,
        "Q": 5012206133771469409285758349474043581500293093464114075203972,
    },
    (26, 29): {
        "D": 522261,
        "C": 2692,
        "W": 5059,
        "u": 7,
        "P": 454592181623521601375,
        "Q": 629039857366305528,
    },
    (26, 47): {
        "D": 364485,
        "C": 2116,
        "W": 1576376215602550757032908052542957159292644701,
        "u": 2611079190314729909962215655566812890419809,
        "P": 24107434053023296515322082382772446705687987537902965213862639312622666158519,
        "Q": 39931089269621661093589841956801775898467910099449466227642200335186151476,
    },
    (98, 101): {
        "D": 95940405,
        "C": 38404,
        "W": 220788520986775462919614657060547064697173954066965079893,
        "u": 22541131703851573370870656292410483656470170855181183,
        "P": 105189902320645145211266433397204047145119454659522482848985261748433399672601631439105713727025690111,
        "Q": 10739233323941537083546286758788229602860630133788513945467635000056068525534702555933010597208288,
    },
}
RAY_SLOPES = tuple(RAYS)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primes_up_to(limit: int) -> tuple[int, ...]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * (
                (limit - p * p) // p + 1
            )
    return tuple(p for p in range(2, limit + 1) if sieve[p])


def candidates() -> tuple[tuple[int, int], ...]:
    """Primitive parity-correct slopes n/2<m<n with odd n and even m."""
    return tuple(
        (m, n)
        for n in range(3, 102, 2)
        for m in range(n // 2 + 1, n)
        if m % 2 == 0 and gcd(m, n) == 1
    )


def conic_coefficients(m: int, n: int) -> tuple[int, int]:
    """Return A,B in 3*W^2=A*u^2+B."""
    return n * n * (4 * m * m - n * n), 4 * (
        2 * m * m + 2 * m * n - n * n
    )


def soluble_mod(m: int, n: int, prime: int) -> bool:
    """Exhaust the exact residue equation 3*w^2=A*u^2+B (mod prime)."""
    a_coefficient, constant = conic_coefficients(m, n)
    left = {(3 * w * w) % prime for w in range(prime)}
    right = {
        (a_coefficient * u * u + constant) % prime for u in range(prime)
    }
    return not left.isdisjoint(right)


def first_obstruction(m: int, n: int, primes: tuple[int, ...]) -> int | None:
    return next((p for p in primes if not soluble_mod(m, n, p)), None)


def advance(w: int, u: int, pell_d: int, p: int, q: int) -> tuple[int, int]:
    """Multiply W+u*sqrt(D) by the norm-one unit P+Q*sqrt(D)."""
    return p * w + pell_d * q * u, q * w + p * u


def cube_member(m: int, n: int, w: int, u: int) -> dict[str, int]:
    """Compile one Pell point into x^3+y^3=(2r+1)^2+2."""
    d = n * n * u * u + 2
    e = u * w
    a = m * n * n * u**3 + (2 * m + n) * u
    q_value = m * m * n * n * u**4 + 2 * m * (m + n) * u * u + 1
    require((d + e) % 2 == 0 and (d - e) % 2 == 0, "cube parity")
    require((a - 1) % 2 == 0, "depth parity")
    x = (d + e) // 2
    y = (d - e) // 2
    r = (a - 1) // 2
    return {
        "d": d,
        "e": e,
        "a": a,
        "q": q_value,
        "x": x,
        "y": y,
        "r": r,
        "value": a * a + 2,
    }


def verify_ray(slope: tuple[int, int], iterations: int = 6) -> dict[str, int]:
    m, n = slope
    data = RAYS[slope]
    pell_d = data["D"]
    pell_c = data["C"]
    w = data["W"]
    u = data["u"]
    p = data["P"]
    q = data["Q"]

    raw_d, raw_c = conic_coefficients(m, n)
    require(raw_d % 3 == 0 and raw_c % 3 == 0, f"{slope}: division by 3")
    require(raw_d // 3 == pell_d and raw_c // 3 == pell_c, f"{slope}: Pell data")
    require(0 < pell_d < n**4 and pell_c > 0, f"{slope}: positive conic")
    require(w * w - pell_d * u * u == pell_c, f"{slope}: seed norm")
    require(p * p - pell_d * q * q == 1, f"{slope}: unit norm")
    require(p > 1 and q > 0 and w > 0 and u > 0, f"{slope}: positive data")
    require(p % 2 == 1 and q % 2 == 0 and pell_d % 2 == 1, f"{slope}: unit parity")
    require(w % 2 == 1 and u % 2 == 1, f"{slope}: seed parity")

    # Put L=n^2*u-W and H=n^2*W-D*u.  The recurrence obeys
    # L'=P*L+Q*H and H'=P*H+D*Q*L.  Thus L,H>0 is an invariant cone;
    # L>0 gives 0<e=uW<n^2*u^2<d and therefore x>y>0 forever.
    lower_defect = n * n * u - w
    upper_defect = n * n * w - pell_d * u
    require(lower_defect > 0 and upper_defect > 0, f"{slope}: seed chamber")

    first: dict[str, int] | None = None
    previous_u = 0
    for step in range(iterations):
        require(w * w - pell_d * u * u == pell_c, f"{slope}: orbit norm {step}")
        require(w % 2 == 1 and u % 2 == 1, f"{slope}: orbit parity {step}")
        require(u > previous_u, f"{slope}: strict orbit growth {step}")
        lower_defect = n * n * u - w
        upper_defect = n * n * w - pell_d * u
        require(lower_defect > 0 and upper_defect > 0, f"{slope}: chamber {step}")

        row = cube_member(m, n, w, u)
        require(row["a"] ** 2 + 2 == row["d"] * row["q"], f"{slope}: factor identity")
        require(4 * row["q"] - row["d"] ** 2 == 3 * row["e"] ** 2, f"{slope}: square identity")
        require(row["x"] > row["y"] > 0, f"{slope}: positive cubes")
        require(row["x"] ** 3 + row["y"] ** 3 == row["value"], f"{slope}: cube identity")
        require(row["value"] == (2 * row["r"] + 1) ** 2 + 2, f"{slope}: target identity")
        if first is None:
            first = row

        next_w, next_u = advance(w, u, pell_d, p, q)
        next_lower = n * n * next_u - next_w
        next_upper = n * n * next_w - pell_d * next_u
        require(next_lower == p * lower_defect + q * upper_defect, f"{slope}: L recurrence")
        require(next_upper == p * upper_defect + pell_d * q * lower_defect, f"{slope}: H recurrence")
        require(next_w * next_w - pell_d * next_u * next_u == pell_c, f"{slope}: norm recurrence")
        previous_u, w, u = u, next_w, next_u

    require(first is not None, f"{slope}: first member")
    return first


def main() -> None:
    primes = primes_up_to(PRIME_LIMIT)
    require(len(primes) == 46 and primes[0] == 2 and primes[-1] == 199, "prime table")
    slopes = candidates()
    require(len(slopes) == EXPECTED_CANDIDATE_COUNT, "candidate count")

    obstruction_by_slope = {
        slope: first_obstruction(*slope, primes) for slope in slopes
    }
    survivors = tuple(
        slope for slope in slopes if obstruction_by_slope[slope] is None
    )
    counts = Counter(
        obstruction for obstruction in obstruction_by_slope.values() if obstruction is not None
    )
    require(survivors == EXPECTED_SURVIVORS, "screen survivors")
    require(dict(sorted(counts.items())) == EXPECTED_OBSTRUCTION_COUNTS, "obstruction counts")
    require(sum(counts.values()) == 517, "excluded count")
    require(all(slope in survivors for slope in RAY_SLOPES), "ray survives screen")

    first_members = {slope: verify_ray(slope) for slope in RAY_SLOPES}
    require(len({(row["x"], row["y"]) for row in first_members.values()}) == 4, "distinct first members")

    atlas_text = ";".join(
        f"{m}/{n}:{obstruction_by_slope[(m, n)] or 0}" for m, n in slopes
    )
    ray_text = ";".join(
        ":".join(
            str(value)
            for value in (
                m,
                n,
                RAYS[(m, n)]["D"],
                RAYS[(m, n)]["C"],
                RAYS[(m, n)]["W"],
                RAYS[(m, n)]["u"],
                RAYS[(m, n)]["P"],
                RAYS[(m, n)]["Q"],
                first_members[(m, n)]["x"],
                first_members[(m, n)]["y"],
                first_members[(m, n)]["r"],
                first_members[(m, n)]["value"],
            )
        )
        for m, n in RAY_SLOPES
    )
    source = Path(__file__).read_bytes().replace(b"\r\n", b"\n")

    print("BERGGREN POSITIVE-CUBE SLOPE SCREEN THROUGH 101")
    print("status=FINITE-EXACT+VERIFIED-EXACT")
    print("universe=primitive_slopes;3<=n<=101;n_odd;m_even;n/2<m<n;gcd(m,n)=1")
    print("prime_screen=all_primes_p<=199;prime_count=46")
    print("candidate_count=528;screen_excluded=517;screen_survivor_count=11")
    print("screen_survivors=" + repr(list(survivors)).replace(" ", ""))
    print(
        "first_obstruction_counts="
        + ",".join(f"p{prime}:{count}" for prime, count in sorted(counts.items()))
    )
    print(f"atlas_semantic_sha256={sha256(atlas_text.encode('ascii')).hexdigest()}")
    print("explicit_positive_ray_count=4;ray_slopes=" + repr(list(RAY_SLOPES)).replace(" ", ""))
    print(
        "recurrence_certificate=norm_multiplier_1;odd_parity_preserved;"
        "L'=P*L+Q*H;H'=P*H+D*Q*L;L,H_positive;u_strictly_grows"
    )
    for m, n in RAY_SLOPES:
        row = first_members[(m, n)]
        data = RAYS[(m, n)]
        print(
            f"ray_{m}_{n}=W^2-{data['D']}u^2={data['C']};"
            f"seed=({data['W']},{data['u']});unit=({data['P']},{data['Q']})"
        )
        print(
            f"first_member_{m}_{n}=x:{row['x']};y:{row['y']};r:{row['r']};"
            f"x^3+y^3=(2r+1)^2+2:{row['value']}"
        )
    print(f"ray_semantic_sha256={sha256(ray_text.encode('ascii')).hexdigest()}")
    print(
        "boundary=the_11_slopes_only_pass_the_displayed_prime_screen;"
        "no_global_insolubility_or_everywhere_local_solubility_claim_for_the_other_7;"
        "explicit_seeds_not_claimed_minimal"
    )
    print(f"source_lf_sha256={sha256(source).hexdigest()}")


if __name__ == "__main__":
    main()
