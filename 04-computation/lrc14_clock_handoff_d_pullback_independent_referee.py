#!/usr/bin/env python3
"""Independent exact referee for the physical THM-2670 D-handoff.

This script calls the separately implemented interval engine in
``lrc14_thm2670_physical_reversed_two_edge_scout.py`` and exhausts all four
ordered guard-cospan sector pairs.  It verifies the ten-clock physical atlas,
the full labelled safe/safe atom census, the danger/danger base-support zero,
both mixed legs, and the absence of any composable pair of physical two-edge
clock triples.

The event chronology is

    E_(b,a;j,h) intersect D^(-1)E_(c,b;h,k),  D(x)={13x}.

All source pairs and state labels are retained.  The carrier is pulled back on
the exact denominator 13*T grid, and the delayed words are multiplied as
Q_(a,h)(z) Q_(b,k)(13z), z={13^6 x}, before the prefix integral is evaluated.
No product of marginal masses is used.
"""

from itertools import product

import lrc14_thm2670_physical_reversed_two_edge_scout as scout


S10 = (
    (3, 0, 3),
    (3, 0, 4),
    (3, 0, 5),
    (3, 0, 6),
    (3, 1, 0),
    (4, 0, 1),
    (4, 0, 2),
    (4, 0, 3),
    (4, 0, 4),
    (4, 6, 0),
)

EXPECTED_PHYSICAL = {
    ("safe", "safe"): (792, 1320, 1320, 1584, 2640,
                         1848, 1848, 1848, 1584, 2376),
    ("safe", "danger"): (264, 528, 528, 528, 264,
                           528, 528, 528, 264, 528),
    ("danger", "safe"): (660, 660, 660, 660, 396,
                           132, 132, 132, 132, 132),
    ("danger", "danger"): (),
}

EXPECTED_BASE_TOTAL = {
    ("safe", "safe"): 31944,
    ("safe", "danger"): 4752,
    ("danger", "safe"): 5016,
    ("danger", "danger"): 0,
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def clock_vector(result):
    return tuple(result["physical_by_clock"].get(clock, 0) for clock in S10)


def excluded_source_pairs(result, clock):
    support = {
        pair
        for (triple, pair), count in result["physical_by_clock_source"].items()
        if triple == clock and count
    }
    universe = set(product(range(1, scout.P), repeat=2))
    return tuple(sorted(universe - support))


def main():
    print("THM-2670 physical D-handoff independent referee -- exact")
    print(f"T={scout.T} Tprime=13*T={scout.TD} R=13^6={scout.SUCCESSOR_SPEED}")
    results = {}
    for sector_pair in EXPECTED_PHYSICAL:
        left, right = sector_pair
        result = scout.global_physical_census(
            left,
            right,
            verbose=False,
            audit_atom_pairs=(sector_pair == ("safe", "safe")),
        )
        results[sector_pair] = result
        totals = result["totals"]
        vector = clock_vector(result)
        expected = EXPECTED_PHYSICAL[sector_pair]
        if expected:
            require(
                tuple(clock for clock, count in sorted(
                    result["physical_by_clock"].items()
                ) if count) == S10,
                f"{sector_pair} physical clock atlas changed",
            )
            require(vector == expected,
                    f"{sector_pair} physical clock counts changed")
        else:
            require(not result["physical_by_clock"],
                    "danger/danger acquired a physical clock")
        require(
            totals.get("base_positive", 0) == EXPECTED_BASE_TOTAL[sector_pair],
            f"{sector_pair} base-positive total changed",
        )
        print(
            f"sector={left}->{right} atoms={result['atom_counts']} "
            f"formal_state_keys={totals['formal']} "
            f"base_positive={totals.get('base_positive',0)} "
            f"word_empty={totals.get('word_empty',0)} "
            f"skew_zero={totals.get('skew_zero',0)} "
            f"physical_state_keys={totals.get('physical',0)}"
        )
        print(f"  physical_by_clock={vector}")

    safe = results[("safe", "safe")]
    safe_totals = safe["totals"]
    require(safe_totals["physical"] == 17160,
            "safe/safe physical state-key total changed")
    require(safe_totals["physical_atom_pairs"] == 17160,
            "safe/safe no longer has one positive labelled atom pair per key")
    require(safe_totals["base_atom_pairs"] == 32736,
            "safe/safe labelled pre-word atom-pair total changed")
    require(safe_totals["word_compatible_atom_pairs"] == 17688,
            "safe/safe joint-word-compatible atom-pair total changed")
    require(safe_totals["atom_pair_skew_zero"] == 528,
            "safe/safe weighted atom-pair skew-zero total changed")
    require(safe_totals["word_empty"] == 14784,
            "safe/safe union-event word-empty total changed")
    require(len(safe["physical_by_source"]) == 144,
            "global safe atlas lost a source pair")

    for index, clock in enumerate(S10):
        missing = excluded_source_pairs(safe, clock)
        forbidden = 6 if index < 5 else 7
        expected_missing = tuple((source0, forbidden)
                                 for source0 in range(1, scout.P))
        require(missing == expected_missing,
                f"safe excluded-source column changed at {clock}")
        require(
            safe["physical_by_clock"][clock] // 132
            == EXPECTED_PHYSICAL[("safe", "safe")][index] // 132,
            "safe per-family state support division changed",
        )
        active_counts = tuple(
            count
            for (triple, _), count in safe["physical_by_clock_source"].items()
            if triple == clock
        )
        expected_family_count = (
            EXPECTED_PHYSICAL[("safe", "safe")][index] // 132
        )
        require(
            len(active_counts) == 132
            and set(active_counts) == {expected_family_count},
            f"safe source-family factorization changed at {clock}",
        )
        print(
            f"safe clock={clock} excluded_s1={forbidden} "
            f"state_support_per_active_source="
            f"{safe['physical_by_clock'][clock] // 132}"
        )

    danger = results[("danger", "danger")]
    require(danger["totals"]["formal"] == 27640,
            "danger/danger formal state-key universe changed")
    require(danger["totals"]["base_zero"] == 27640,
            "a danger/danger chain survived the base D-pullback")
    print("danger/danger: all 27640 formal state chains are base-zero")
    print(
        "safe/safe labelled atom pairs: "
        "base-positive=32736 word-compatible=17688 "
        "skew-zero=528 physical-positive=17160"
    )

    live_sector_pairs = (
        ("safe", "safe"),
        ("safe", "danger"),
        ("danger", "safe"),
    )
    for left_index, left_pair in enumerate(live_sector_pairs):
        for right_pair in live_sector_pairs[left_index + 1:]:
            require(
                results[left_pair]["physical_state_keys"].isdisjoint(
                    results[right_pair]["physical_state_keys"]
                ),
                f"sector-state supports overlap: {left_pair}, {right_pair}",
            )
    free_by_clock = tuple(
        sum(results[pair]["physical_by_clock"].get(clock, 0)
            for pair in live_sector_pairs)
        for clock in S10
    )
    expected_free_by_clock = (
        1716, 2508, 2508, 2772, 3300,
        2508, 2508, 2508, 1980, 3036,
    )
    require(free_by_clock == expected_free_by_clock,
            "guard-free disjoint-union clock counts changed")
    require(sum(free_by_clock) == 25344,
            "guard-free disjoint-union total changed")
    print(f"guard-free disjoint sector union by clock={free_by_clock}")
    print("guard-free total=25344; per-active-source states="
          "(13,19,19,21,25,19,19,19,15,23)")

    union_clocks = set()
    for result in results.values():
        union_clocks.update(result["physical_by_clock"])
    require(union_clocks == set(S10), "cospan-union clock atlas changed")
    composable = tuple(
        (first, second)
        for first in S10
        for second in S10
        if first[1:] == second[:2]
    )
    require(not composable, "a physical three-edge clock chain appeared")
    print("cospan union clock atlas=S10; composable two-edge pairs=0")
    print("verdict=PASS: exact physical chronology has positive two-edge handoffs but no three-edge chain")


if __name__ == "__main__":
    main()
