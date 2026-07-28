#!/usr/bin/env python3
"""Exact target-b coefficient census and BABA primitive-unit certificate.

This is the missing normalization companion for the intrinsic guard-cospan
debt carrier in ``lrc14_delayed_tail_scc_factor_diagnostic_20260728.py``.  It
changes only
the delayed target word in the complete THM-2623 successor-half-cell universe:

    target a (inherited): danger at 13^3, safe at 2*13^5;
    target b (this file): danger at 2*13^5, safe at 13^3.

The guard safe/danger sectors, five ordinary safeties, seven owner-clock
words, 162 THM-2584 rails, two tooth halves, two future half-digits, thirteen
future digits, and twelve nonzero roots are unchanged.  All 18,398,016
target-b subdigit coefficients are recomputed exactly.  Their gcd is 520, so
combining them with the proved target-a bank of content 26 has gcd 26.

At the four exact BABA states, the selected /26 coefficient vectors are then
recomputed from their actual rail, guard, half, future digit, half-digit, and
private root.  After private-root normalization, all four multiplication
determinants in F_13[z]/(Phi_7) are nonzero.  This certifies the primitive-unit
sidecar only; it does not provide transport between those units, a terminal-
word transition, endpoint transport, a row exclusion, or LRC(14).
"""

from concurrent.futures import ProcessPoolExecutor
from functools import reduce
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_successor_halfcell_carry_no_go_thm2623 as probe


core = probe.core
old = probe.old
P, Q7, T, N = probe.P, probe.Q7, probe.T, probe.N


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def merge_support(intervals):
    out = []
    for left, right in sorted(intervals):
        if out and left <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], right))
        else:
            out.append((left, right))
    return out


def target_b_sector_words(module):
    """Guard-safe/danger words with the two deep roles exchanged."""
    base = module.make_comb(module.C3, 182, -13, 13)
    words = []
    for guard_danger in (False, True):
        if guard_danger:
            word = module.intersect_comb(
                base, module.W[module.GUARD], 91, -13, 13
            )
        else:
            word = module.subtract_comb(
                base, module.W[module.GUARD], 91, -13, 13
            )
        for index in module.UNIT_IDX:
            word = module.subtract_comb(
                word, module.W[index], 182, -13, 13
            )
        word = module.subtract_comb(word, module.C2, 182, -13, 13)
        words.append(word)
    return tuple(words)


def target_b_danger_split_words(module):
    """Split broad guard danger into its narrow tooth and annulus."""
    base = module.make_comb(module.C3, 182, -13, 13)
    broad = module.intersect_comb(
        base, module.W[module.GUARD], 91, -13, 13
    )
    narrow = module.intersect_comb(
        base, module.W[module.GUARD], 182, -13, 13
    )
    annulus = module.subtract_comb(
        broad, module.W[module.GUARD], 182, -13, 13
    )
    result = []
    for word in (narrow, annulus):
        for index in module.UNIT_IDX:
            word = module.subtract_comb(
                word, module.W[index], 182, -13, 13
            )
        word = module.subtract_comb(word, module.C2, 182, -13, 13)
        result.append(word)
    return tuple(result)


def build_subprefixes(module, words):
    result = []
    for word in words:
        by_clock = []
        for ell in range(Q7):
            qell = module.subtract_comb(
                word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
            )
            pieces = tuple(
                module.make_prefix(old.sat.intersect_interval(
                    qell, j * T // N, (j + 1) * T // N
                )) for j in range(N)
            )
            require(
                sum(piece[2][-1] for piece in pieces)
                == sum(right - left for left, right in qell),
                "future half-digits stopped partitioning an owner word",
            )
            by_clock.append(pieces)
        result.append(tuple(by_clock))
    return tuple(result)


EXPECTED_WORD_CENSUS = (
    (188_634, 10_807_572_370_920),
    (86_950, 4_981_755_322_960),
)
EXPECTED_SPLIT_CENSUS = (
    (44_222, 2_533_638_252_600),
    (42_728, 2_448_117_070_360),
)
EXPECTED_SHARDS = (
    ((0, 41), 520, 787_088, 44_772, 4_656_288, 0),
    ((41, 82), 520, 720_690, 44_772, 4_656_288, 0),
    ((82, 123), 520, 720_690, 44_772, 4_656_288, 0),
    ((123, 162), 520, 758_384, 42_588, 4_429_152, 0),
)
EXPECTED_REFINED_SHARDS = (
    (520, 834_432, 6_984_432, 0),
    (520, 764_910, 6_984_432, 0),
    (520, 764_910, 6_984_432, 0),
    (520, 803_616, 6_643_728, 0),
)


def shard(bounds):
    start, stop = bounds
    module, _safe, _whole, _masses, rails, present, starts = (
        core.build_carrier_data()
    )
    safe_word = target_b_sector_words(module)[0]
    narrow_word, annulus_word = target_b_danger_split_words(module)
    prefixes = build_subprefixes(
        module, (safe_word, narrow_word, annulus_word)
    )
    caches = [[dict() for _ in range(Q7)] for _ in range(3)]
    broad_content = 0
    broad_positive = 0
    refined_content = 0
    refined_positive = 0
    checks = 0
    broad_values_tested = 0
    refined_values_tested = 0
    broad_nondivisible_26 = 0
    refined_nondivisible_26 = 0
    for j in range(start, stop):
        _source, _owner, _deep, pieces = rails[j]
        for h in range(P):
            for ell5 in range(Q7):
                overlap = old.intersect_weighted_union(
                    pieces, present[ell5, (-h) % P], starts[ell5, (-h) % P]
                )
                for root in range(1, P):
                    halves = (
                        old.intersect_weighted_comb(
                            overlap, module.C3, 182,
                            14 * root - 13, 14 * root,
                        ),
                        old.intersect_weighted_comb(
                            overlap, module.C3, 182,
                            14 * root, 14 * root + 13,
                        ),
                    )
                    for half in halves:
                        sector_values = tuple(
                            probe.delayed_all_subdigits(
                                half, prefixes[sector][ell5],
                                caches[sector][ell5],
                            ) for sector in range(3)
                        )
                        safe_values, narrow_values, annulus_values = sector_values
                        broad_values = tuple(
                            narrow + annulus
                            for narrow, annulus
                            in zip(narrow_values, annulus_values)
                        )
                        broad_values_tested += len(safe_values) + len(broad_values)
                        refined_values_tested += sum(map(len, sector_values))
                        for value in safe_values + broad_values:
                            if value:
                                broad_positive += 1
                                broad_content = gcd(broad_content, value)
                                broad_nondivisible_26 += int(value % 26 != 0)
                        for values in sector_values:
                            for value in values:
                                if value:
                                    refined_positive += 1
                                    refined_content = gcd(refined_content, value)
                                    refined_nondivisible_26 += int(value % 26 != 0)
                    checks += 1
    return (
        bounds, broad_content, broad_positive, checks,
        broad_values_tested, broad_nondivisible_26,
        refined_content, refined_positive,
        refined_values_tested, refined_nondivisible_26,
    )


# target, guard sector (0 safe/1 danger), rail (source,owner,deep),
# tooth half (0 left/1 right), future half-digit, future digit, private root.
MIXED_SLICES = (
    ("b", 1, (1, 1, 0), 1, 1, 0, 2),
    ("a", 0, (1, 0, 12), 1, 1, 6, 3),
    ("b", 1, (1, 1, 0), 0, 0, 12, 2),
    ("a", 0, (1, 0, 12), 0, 0, 6, 3),
)
EXPECTED_BROAD_MIXED_VECTORS = (
    (9, 0, 0, 0, 0, 1, 0),
    (0, 0, 0, 2, 0, 0, 0),
    (9, 0, 0, 0, 10, 10, 0),
    (0, 0, 0, 11, 0, 0, 0),
)
EXPECTED_BROAD_MIXED_DETERMINANTS = (12, 12, 10, 12)
EXPECTED_DEBT_MIXED_VECTORS = (
    (9, 0, 0, 0, 0, 12, 0),
    (0, 0, 0, 2, 0, 0, 0),
    (9, 0, 0, 0, 2, 4, 0),
    (0, 0, 0, 11, 0, 0, 0),
)
EXPECTED_DEBT_MIXED_DETERMINANTS = (12, 12, 8, 12)
EXPECTED_ANNULUS_B_VECTORS = (
    (0, 0, 0, 0, 0, 2, 0),
    (0, 0, 0, 0, 8, 6, 0),
)


def selected_mixed_rows():
    module, _safe, _whole, _masses, rails, present, starts = (
        core.build_carrier_data()
    )
    target_a_prefixes = probe.build_subprefixes(module)
    target_b_broad_prefixes = build_subprefixes(
        module, target_b_sector_words(module)
    )
    target_b_split_prefixes = build_subprefixes(
        module, target_b_danger_split_words(module)
    )
    caches = {
        label: [dict() for _ in range(Q7)]
        for label in ("a", "broad", "narrow", "annulus")
    }
    broad_rows = []
    debt_rows = []
    annulus_rows = []

    def compute_raw(by_clock, cache, pieces, edge, kappa, h, root):
        raw = []
        for ell5 in range(Q7):
            overlap = old.intersect_weighted_union(
                pieces, present[ell5, (-h) % P], starts[ell5, (-h) % P]
            )
            low = 14 * root - 13 if edge == 0 else 14 * root
            high = 14 * root if edge == 0 else 14 * root + 13
            half = old.intersect_weighted_comb(
                overlap, module.C3, 182, low, high
            )
            values = probe.delayed_all_subdigits(
                half, by_clock[ell5], cache[ell5]
            )
            raw.append(values[2 * h + kappa])
        return tuple(raw)

    def certify(raw, root):
        require(any(raw) and all(value % 26 == 0 for value in raw),
                "selected mixed row lost joint /26 integrality")
        vector = tuple((value // 26) % P for value in raw)
        scaled = tuple(value * pow(root, -1, P) % P for value in vector)
        reduced = tuple((scaled[i] - scaled[-1]) % P for i in range(Q7 - 1))
        determinant = old.sat.multiplication_determinant_7(reduced)
        return raw, vector, determinant

    for target, sector, rail_key, edge, kappa, h, root in MIXED_SLICES:
        matches = tuple(
            pieces for source, owner, deep, pieces in rails
            if (source, owner, deep) == rail_key
        )
        require(len(matches) == 1, "selected THM-2584 rail stopped being unique")
        pieces = matches[0]
        if target == "a":
            raw = compute_raw(
                target_a_prefixes[sector], caches["a"],
                pieces, edge, kappa, h, root,
            )
            row = certify(raw, root)
            broad_rows.append(row)
            debt_rows.append(row)
        else:
            broad_raw = compute_raw(
                target_b_broad_prefixes[sector], caches["broad"],
                pieces, edge, kappa, h, root,
            )
            narrow_raw = compute_raw(
                target_b_split_prefixes[0], caches["narrow"],
                pieces, edge, kappa, h, root,
            )
            annulus_raw = compute_raw(
                target_b_split_prefixes[1], caches["annulus"],
                pieces, edge, kappa, h, root,
            )
            require(tuple(narrow + annulus
                          for narrow, annulus
                          in zip(narrow_raw, annulus_raw)) == broad_raw,
                    "selected narrow/annulus rows stopped summing to broad")
            broad_rows.append(certify(broad_raw, root))
            debt_rows.append(certify(narrow_raw, root))
            annulus_rows.append(certify(annulus_raw, root))
    return tuple(broad_rows), tuple(debt_rows), tuple(annulus_rows)


def main():
    module, _safe, _whole, _masses, _rails, _present, _starts = (
        core.build_carrier_data()
    )
    words = target_b_sector_words(module)
    census = tuple(
        (len(word), sum(right - left for left, right in word))
        for word in words
    )
    require(census == EXPECTED_WORD_CENSUS,
            "target-b safe/danger word census changed")
    split_words = target_b_danger_split_words(module)
    split_census = tuple(
        (len(word), sum(right - left for left, right in word))
        for word in split_words
    )
    require(split_census == EXPECTED_SPLIT_CENSUS,
            "target-b narrow/annulus word census changed")
    require(merge_support(split_words[0] + split_words[1]) == words[1]
            and sum(row[1] for row in split_census) == census[1][1],
            "narrow/annulus words stopped partitioning broad danger")
    for sector, word in enumerate(words + split_words):
        owner_words = tuple(
            module.subtract_comb(
                word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
            ) for ell in range(Q7)
        )
        require(
            merge_support(piece for row in owner_words for piece in row) == word,
            f"target-b sector {sector} owner words stopped covering its base",
        )

    with ProcessPoolExecutor(max_workers=len(core.SHARDS)) as pool:
        results = tuple(pool.map(shard, core.SHARDS))
    require(tuple(row[:6] for row in results) == EXPECTED_SHARDS,
            "target-b broad shard census/content changed")
    require(tuple(row[6:] for row in results) == EXPECTED_REFINED_SHARDS,
            "target-b refined shard census/content changed")
    global_content = reduce(gcd, (row[1] for row in results), 0)
    positive = sum(row[2] for row in results)
    checks = sum(row[3] for row in results)
    values_tested = sum(row[4] for row in results)
    nondivisible_26 = sum(row[5] for row in results)
    refined_content = reduce(gcd, (row[6] for row in results), 0)
    refined_positive = sum(row[7] for row in results)
    refined_values_tested = sum(row[8] for row in results)
    refined_nondivisible_26 = sum(row[9] for row in results)
    require(
        (global_content, positive, checks, values_tested, nondivisible_26)
        == (520, 2_986_852, 176_904, 18_398_016, 0),
        "target-b aggregate coefficient census changed",
    )
    require((refined_content, refined_positive, refined_values_tested,
             refined_nondivisible_26)
            == (520, 3_167_868, 27_597_024, 0),
            "target-b refined coefficient census/content changed")
    require(gcd(26, global_content) == 26,
            "joint target-a/target-b primitive content changed")
    require(gcd(26, refined_content) == 26,
            "joint target-a/refined-target-b content changed")

    broad_rows, debt_rows, annulus_rows = selected_mixed_rows()
    broad_vectors = tuple(row[1] for row in broad_rows)
    broad_determinants = tuple(row[2] for row in broad_rows)
    debt_vectors = tuple(row[1] for row in debt_rows)
    debt_determinants = tuple(row[2] for row in debt_rows)
    annulus_vectors = tuple(row[1] for row in annulus_rows)
    require(broad_vectors == EXPECTED_BROAD_MIXED_VECTORS,
            "selected broad BABA /26 coefficient vectors changed")
    require(broad_determinants == EXPECTED_BROAD_MIXED_DETERMINANTS,
            "selected broad BABA unit determinants changed")
    require(debt_vectors == EXPECTED_DEBT_MIXED_VECTORS,
            "selected narrow-debt BABA /26 coefficient vectors changed")
    require(debt_determinants == EXPECTED_DEBT_MIXED_DETERMINANTS,
            "selected narrow-debt BABA unit determinants changed")
    require(annulus_vectors == EXPECTED_ANNULUS_B_VECTORS,
            "selected annulus B-row vectors changed")
    require(all(debt_determinants),
            "selected narrow-debt BABA row lost unitness")

    print("LRC14 target-b successor-half-cell content and BABA unit certificate")
    print("status=VERIFIED-EXACT COMPANION; canonical THM-2623 universe")
    print(f"target_b_word_census_safe_danger={census}")
    print(f"target_b_danger_split_census_narrow_annulus={split_census}")
    print(f"shards={results}")
    print(f"target_b_global_content={global_content} positives={positive} "
          f"fine_checks={checks} subdigit_values={values_tested} "
          f"nondivisible_by_26={nondivisible_26}")
    print(f"target_b_refined_content={refined_content} "
          f"refined_positives={refined_positive} "
          f"refined_subdigit_values={refined_values_tested} "
          f"refined_nondivisible_by_26={refined_nondivisible_26}")
    print("inherited_target_a_content=26 "
          f"joint_a_b_broad_content=gcd(26,{global_content})="
          f"{gcd(26, global_content)} "
          f"joint_a_b_refined_content=gcd(26,{refined_content})="
          f"{gcd(26, refined_content)}")
    print(f"mixed_slices={MIXED_SLICES}")
    print(f"mixed_broad_raw_rows={tuple(row[0] for row in broad_rows)}")
    print(f"mixed_debt_raw_rows={tuple(row[0] for row in debt_rows)}")
    print(f"mixed_annulus_b_raw_rows={tuple(row[0] for row in annulus_rows)}")
    print(f"mixed_broad_vectors_div26_mod13={broad_vectors}")
    print(f"mixed_debt_vectors_div26_mod13={debt_vectors}")
    print(f"mixed_annulus_b_vectors_div26_mod13={annulus_vectors}")
    print(f"mixed_broad_root_normalized_determinants={broad_determinants}")
    print(f"mixed_debt_root_normalized_determinants={debt_determinants}")
    print("verdict=all four narrow BABA debt-sector rows have nonzero unit "
          "determinant after joint /26 normalization")
    print("SCOPE: no unit transport, terminal-word transition, endpoint "
          "transport, row exclusion, or LRC14 conclusion")


if __name__ == "__main__":
    main()
