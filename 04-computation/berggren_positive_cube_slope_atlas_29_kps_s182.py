#!/usr/bin/env python3
"""Finite exact slope atlas and second positive Pell ray for THM-3376."""

from __future__ import annotations

from hashlib import sha256
from math import gcd
from pathlib import Path


OBSTRUCTION_PRIMES = (3, 5, 7, 11, 13, 23)
EXPECTED_GROUPS = {
    3: ((2, 3), (4, 5), (6, 7), (8, 9), (6, 11), (10, 11),
        (8, 13), (12, 13), (8, 15), (14, 15), (10, 17),
        (12, 17), (16, 17), (12, 19), (14, 19), (18, 19),
        (16, 21), (20, 21), (12, 23), (16, 23), (18, 23),
        (22, 23), (14, 25), (18, 25), (24, 25), (14, 27),
        (16, 27), (20, 27), (22, 27), (26, 27), (16, 29),
        (18, 29), (22, 29), (24, 29), (28, 29)),
    5: ((4, 7), (14, 17)),
    7: ((10, 13), (20, 23), (16, 25)),
    11: ((8, 11), (20, 29)),
    13: ((10, 19), (16, 19)),
    23: ((22, 25),),
}
SURVIVORS = ((14, 23), (26, 29))

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
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def candidates() -> tuple[tuple[int, int], ...]:
    return tuple(
        (m, n)
        for n in range(3, 30, 2)
        for m in range(n // 2 + 1, n)
        if m % 2 == 0 and gcd(m, n) == 1
    )


def conic_coefficients(m: int, n: int) -> tuple[int, int]:
    return n * n * (4 * m * m - n * n), 4 * (2 * m * m + 2 * m * n - n * n)


def soluble_mod(m: int, n: int, prime: int) -> bool:
    a, b = conic_coefficients(m, n)
    return any(
        (3 * w * w - a * u * u - b) % prime == 0
        for u in range(prime)
        for w in range(prime)
    )


def first_obstruction(m: int, n: int) -> int | None:
    return next((p for p in OBSTRUCTION_PRIMES if not soluble_mod(m, n, p)), None)


def advance(w: int, u: int, d: int, p: int, q: int) -> tuple[int, int]:
    return p * w + d * q * u, q * w + p * u


def member(m: int, n: int, w: int, u: int) -> dict[str, int]:
    d = n * n * u * u + 2
    e = u * w
    a = m * n * n * u**3 + (2 * m + n) * u
    q = m * m * n * n * u**4 + 2 * m * (m + n) * u * u + 1
    require((d + e) % 2 == 0 and (d - e) % 2 == 0, "cube parity")
    x = (d + e) // 2
    y = (d - e) // 2
    require((a - 1) % 2 == 0, "depth parity")
    r = (a - 1) // 2
    return {"d": d, "e": e, "a": a, "q": q, "x": x, "y": y, "r": r}


def main() -> None:
    slopes = candidates()
    require(len(slopes) == 47, "candidate count")
    groups: dict[int, list[tuple[int, int]]] = {p: [] for p in OBSTRUCTION_PRIMES}
    survivors: list[tuple[int, int]] = []
    for m, n in slopes:
        obstruction = first_obstruction(m, n)
        if obstruction is None:
            survivors.append((m, n))
        else:
            groups[obstruction].append((m, n))
    require(tuple(survivors) == SURVIVORS, "survivor slopes")
    require(
        {p: tuple(group) for p, group in groups.items()} == EXPECTED_GROUPS,
        "obstruction atlas",
    )

    ray_summaries = []
    for slope in SURVIVORS:
        m, n = slope
        data = RAYS[slope]
        d_conic, c3 = conic_coefficients(m, n)
        require(d_conic % 3 == 0 and c3 % 3 == 0, "Pell division by three")
        require(d_conic // 3 == data["D"] and c3 // 3 == data["C"], "Pell data")
        require(data["W"] ** 2 - data["D"] * data["u"] ** 2 == data["C"],
                "seed norm")
        require(data["P"] ** 2 - data["D"] * data["Q"] ** 2 == 1, "unit norm")
        require(data["W"] % 2 == data["u"] % 2 == data["P"] % 2 == 1,
                "odd seed/unit scalar")
        require(data["Q"] % 2 == 0, "even unit radical")
        require(data["W"] < n * n * data["u"], "positive ratio gate")
        w, u = data["W"], data["u"]
        members = []
        for _ in range(6):
            row = member(m, n, w, u)
            require(w * w - data["D"] * u * u == data["C"], "orbit norm")
            require(row["a"] ** 2 + 2 == row["d"] * row["q"], "factor identity")
            require(4 * row["q"] - row["d"] ** 2 == 3 * row["e"] ** 2,
                    "square identity")
            require(row["x"] > row["y"] > 0, "positive chamber")
            require(row["x"] ** 3 + row["y"] ** 3 == row["a"] ** 2 + 2,
                    "cube identity")
            members.append(row)
            w, u = advance(w, u, data["D"], data["P"], data["Q"])
        require(all(members[i]["r"] < members[i + 1]["r"] for i in range(5)),
                "strict depth growth")
        ray_summaries.append((slope, members[0]))

    first_1423 = ray_summaries[0][1]
    require(first_1423["x"] == 2302290678332599157988107845344339684655411,
            "14/23 first x")
    require(first_1423["y"] == 984700897345246967397057762882950552473960,
            "14/23 first y")
    require(first_1423["r"] == 1813711014853445365774251238977913323291598894943013194651502886,
            "14/23 first depth")

    # The two limiting chamber ratios are different, so these are distinct rays.
    require(44965 * 841**2 != 522261 * 529**2, "distinct limiting ratios")

    atlas_text = ";".join(
        f"{m}/{n}:{first_obstruction(m,n) or 0}" for m, n in slopes
    )
    ray_text = ";".join(
        f"{m}/{n}:{row['r']}:{row['x']}:{row['y']}"
        for (m, n), row in ray_summaries
    )
    source = Path(__file__).read_bytes().replace(b"\r\n", b"\n")
    print("BERGGREN POSITIVE-CUBE SLOPE ATLAS THROUGH 29")
    print("status=PROVED+FINITE-EXACT+VERIFIED-EXACT")
    print("universe=primitive_slopes;3<=n<=29;n_odd;m_even;n/2<m<n;gcd(m,n)=1")
    print("candidate_count=47;modularly_excluded=45;survivors=[(14,23),(26,29)]")
    print("obstruction_counts=" + ",".join(f"p{p}:{len(groups[p])}" for p in OBSTRUCTION_PRIMES))
    print(f"atlas_semantic_sha256={sha256(atlas_text.encode('ascii')).hexdigest()}")
    print(
        "ray_14_23=Pell(W^2-44965u^2=676);"
        f"seed=({RAYS[(14,23)]['W']},{RAYS[(14,23)]['u']});"
        f"first_r={first_1423['r']};first_x={first_1423['x']};first_y={first_1423['y']}"
    )
    print("ray_26_29=Pell(W^2-522261u^2=2692);seed=(5059,7);continuation=THM3375")
    print("both_units=norm_one_with_odd_scalar_even_radical;both_positive_orbits=infinite")
    print("limiting_chamber_ratios=distinct")
    print(f"ray_semantic_sha256={sha256(ray_text.encode('ascii')).hexdigest()}")
    print(
        "boundary=finite_slope_atlas_only_through_n29;explicit_seeds_not_"
        "claimed_minimal;no_counting_asymptotic_or_other_frontier_transfer"
    )
    print(f"source_lf_sha256={sha256(source).hexdigest()}")


if __name__ == "__main__":
    main()
