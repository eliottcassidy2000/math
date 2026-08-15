#!/usr/bin/env python3
"""Exact stopping audit for six-clock completions and dyadic samples.

This companion composes with the *current refined* THM-3366 k=2 and k=3
row keys.  It does not subtract a raw support census.  On each row surviving
the proved <=5 complement-clock screen it asks whether exactly six clocks
from 1..14 cover the unsupported open cells and separates:

* the literal six-body set itself as a completion;
* a different six-clock completion; and
* no six-clock completion from the pool.

It then audits the half sample u=1/2 and the paired quarter samples
u=1/4,3/4.  These compile to literal half-twist sheets on D and 2D,
respectively.  Only rank-support contradictions from proved THM-3415,
THM-3416, and THM-3434 are counted.  Aligned parity is deliberately not
inferred from the current (body,D) row key.
"""

from __future__ import annotations

from collections import Counter
from functools import reduce
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

K2_PATH = ROOT / "04-computation/lrc14_k2_refined_complement_clock_composition_kps_s175.py"
K3_PATH = ROOT / "04-computation/lrc14_k3_refined_complement_clock_composition_kps_s176.py"
BASE_PATH = ROOT / "04-computation/lrc14_k1_body_complement_clock_scan_kps_s172.py"
THM3366_PATH = ROOT / "01-canon/theorems/THM-3366-all-sector-complement-clock-completion.md"
THM3415_PATH = ROOT / "01-canon/theorems/THM-3415-zero-mode-cochain-global-rank-five-support.md"
THM3416_PATH = ROOT / "01-canon/theorems/THM-3416-zero-mode-cochain-global-rank-six-support.md"
THM3434_PATH = ROOT / "01-canon/theorems/THM-3434-seventeen-fibre-two-sided-mass-closure.md"

EXPECTED_DEPENDENCIES = {
    K2_PATH: "414b3777cc44f2e059d5bc4258ab555fc055de8c0f8b3ebfd99ce0c90c7da14c",
    K3_PATH: "27e4bff52705189bf8ff73db42d76d4e2fc94c44330d295f166f8d4217cb1804",
    BASE_PATH: "bdb2001cf22f7e92884e895b0095021e42e8f1febd9adbf779b250a2f6c53507",
    THM3366_PATH: "9d6a356cf7440f0eae271ec8163b53895df9ac70332be27a1b0fe62da5172e43",
    THM3415_PATH: "de8d2f615695070dc75cad38ad4a9c22308d8bf900cd6f9a69cd2003f815a14d",
    THM3416_PATH: "42a9309145de51d1bb6fca0b7c1945302ff37a63a3183e1dfed838c07118e8bf",
    THM3434_PATH: "52a5a5caed75ad48ab35ce287c15f6cece074c88dddda2b87329792a3df70af7",
}

POOL = tuple(range(1, 15))
RANK5_BASES = (8, 9, 10, 12)
RANK6_BASES = (8, 9, 10, 11, 12, 15, 23, 25)
ODD_RANK7_BASES = (9, 11, 13, 15, 23, 25, 29, 51)

EXPECTED_K2_PRE = (26_899_164_786, 200_141_092_521, 4_354, 2_966, 147)
EXPECTED_K2_POST = (4_056, 200_069_517_203)
EXPECTED_K2_CATEGORIES = {
    "body-only": (2_968, 199_075_988_191),
    "nontrivial-only": (382, 133_049_669),
    "none": (706, 860_479_343),
    "body+nontrivial": (0, 0),
}
EXPECTED_K2_SEMANTIC = "4325e0a02db45ca1c8e89383317141eb820dabcf35415a7cd532c50baa4ccdb5"

EXPECTED_K3_PRE = (398_241_574, 2_548_901_482, 1_904, 1_823, 107)
EXPECTED_K3_POST = (1_897, 2_548_893_834)
EXPECTED_K3_CATEGORIES = {
    "body-only": (1_823, 2_547_058_578),
    "nontrivial-only": (20, 313_120),
    "none": (54, 1_522_136),
    "body+nontrivial": (0, 0),
}
EXPECTED_K3_SEMANTIC = "9420bb276b8652f19550ba8c18f1db1b63b59b5246dd3a8312d120a99c66d415"


def require(condition: bool, payload) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("module", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def supported(modulus: int, bases: tuple[int, ...]) -> bool:
    return any(modulus % base == 0 for base in bases)


def first_support(modulus: int, bases: tuple[int, ...]) -> int | None:
    return next((base for base in bases if modulus % base == 0), None)


def arrangement_data(base):
    points = base.arrangement_points(POOL)
    atoms = tuple(zip(points, points[1:]))
    samples = points + tuple((left + right) / 2 for left, right in atoms)
    clock_masks = {
        clock: sum(
            (1 << index for index, value in enumerate(samples) if base.danger(clock, value)),
            0,
        )
        for clock in POOL
    }
    six = tuple(combinations(POOL, 6))
    six_unions = {
        clocks: reduce(int.__or__, (clock_masks[clock] for clock in clocks), 0)
        for clocks in six
    }
    return points, atoms, samples, clock_masks, six, six_unions


def k2_refined(c2, d6, base, points, atoms, solve):
    rows_by_divisor, _nonseptimal = d6.base.build_rows()
    queries = d6.base.engine_queries(rows_by_divisor)
    _raw, answers, engine_output = c2.run_engine_portable(d6.base, queries)
    before = c2.aggregate(d6, rows_by_divisor, answers)
    require(before["summary"] == EXPECTED_K2_PRE, ("k2 pre", before["summary"]))

    rowdata = {}
    for divisor, rows in rows_by_divisor.items():
        for body, ruler, support_count, _histogram, _values in rows:
            key = (body, divisor)
            if key not in before["rows"]:
                continue
            check_ruler, ranges = d6.base.support.safe_cell_ranges(body)
            require(check_ruler == ruler, ("k2 ruler", key, ruler, check_ruler))
            arcs = base.residue_arcs(divisor, ranges)
            require(
                sum(right - left for left, right in arcs) == support_count,
                ("k2 support", key),
            )
            gaps = base.unsupported_gaps(divisor, arcs)
            target = c2.target_mask(base, divisor, gaps, points, atoms)
            rowdata[key] = (ruler, support_count, arcs, gaps, target)

    old_terminals = {key for key, data in rowdata.items() if solve(data[-1]) is not None}
    require(len(old_terminals) == 298, ("k2 old rows", len(old_terminals)))
    require(
        sum(before["row_counts"][key] for key in old_terminals) == 71_575_318,
        "k2 old occurrences",
    )
    post = before["rows"] - old_terminals
    post_occurrences = sum(before["row_counts"][key] for key in post)
    require((len(post), post_occurrences) == EXPECTED_K2_POST, "k2 post")
    return before, rowdata, post, sha256(engine_output.encode()).hexdigest()


def k3_refined(c3, k3, base, points, atoms, solve):
    by_divisor, body_count, body_divisor_rows = k3.q7.build_rows()
    require((body_count, body_divisor_rows) == (3_003, 251_536), "k3 universe")
    before = c3.aggregate(k3, by_divisor)
    require(before["summary"] == EXPECTED_K3_PRE, ("k3 pre", before["summary"]))

    rowdata = {}
    for divisor, rows in by_divisor.items():
        for support_count, body, ruler, arcs in rows:
            key = (body, divisor)
            if key not in before["rows"]:
                continue
            gaps = base.unsupported_gaps(divisor, arcs)
            target = c3.target_mask(base, divisor, gaps, points, atoms)
            rowdata[key] = (ruler, support_count, arcs, gaps, target)

    old_terminals = {key for key, data in rowdata.items() if solve(data[-1]) is not None}
    require(len(old_terminals) == 7, ("k3 old rows", len(old_terminals)))
    require(
        sum(before["row_counts"][key] for key in old_terminals) == 7_648,
        "k3 old occurrences",
    )
    post = before["rows"] - old_terminals
    post_occurrences = sum(before["row_counts"][key] for key in post)
    require((len(post), post_occurrences) == EXPECTED_K3_POST, "k3 post")
    return before, rowdata, post


def classify_six(post, rowdata, row_counts, six, six_unions, expected, expected_semantic):
    category_rows = Counter()
    category_occurrences = Counter()
    coindex = Counter()
    examples = {}
    classifications = {}
    semantic = sha256()

    for body, divisor in sorted(post, key=lambda key: (key[1], key[0])):
        ruler, support_count, _arcs, _gaps, target = rowdata[(body, divisor)]
        body_covers = target & ~six_unions[body] == 0
        nontrivial = None
        for clocks in six:
            if clocks != body and target & ~six_unions[clocks] == 0:
                nontrivial = clocks
                break
        if body_covers and nontrivial is not None:
            category = "body+nontrivial"
        elif body_covers:
            category = "body-only"
        elif nontrivial is not None:
            category = "nontrivial-only"
        else:
            category = "none"
        key = (body, divisor)
        count = row_counts[key]
        classifications[key] = (category, nontrivial)
        category_rows[category] += 1
        category_occurrences[category] += count
        coindex[(category, ruler // divisor)] += 1
        examples.setdefault(category, (key, (ruler, support_count), nontrivial, count))
        semantic.update(
            f"{divisor}|{body}|{ruler}|{support_count}|{category}|{nontrivial}|{count}\n".encode()
        )

    observed = {
        category: (category_rows[category], category_occurrences[category])
        for category in expected
    }
    require(observed == expected, ("six categories", observed, expected))
    require(semantic.hexdigest() == expected_semantic, "six semantic")
    return observed, coindex, examples, classifications, semantic.hexdigest()


def first_seven_cover(key, rowdata, clock_masks):
    target = rowdata[key][-1]
    for clocks in combinations(POOL, 7):
        union = reduce(int.__or__, (clock_masks[clock] for clock in clocks), 0)
        if target & ~union == 0:
            return clocks
    return None


def complement_from_arcs(divisor: int, arcs) -> tuple[int, ...]:
    support = bytearray(divisor)
    for left, right in arcs:
        support[left:right] = b"\x01" * (right - left)
    return tuple(index for index, value in enumerate(support) if not value)


def half_block_size(modulus: int, residue: int) -> int:
    period = 2 * modulus
    count = 0
    for sheet in range(modulus):
        phase = residue * (2 * sheet + 1) % period
        if 7 * min(phase, period - phase) < modulus:
            count += 1
    return count


def maximum_transverse_half_block(modulus: int):
    best = (-1, None, None)
    for residue in range(1, 2 * modulus):
        if residue == modulus:
            continue
        size = half_block_size(modulus, residue)
        quotient_order = modulus // gcd(modulus, residue)
        if size > best[0]:
            best = (size, residue, quotient_order)
    return best


def support_gate_histogram(post, row_counts, multiplier, bases):
    rows = Counter()
    occurrences = Counter()
    for key in post:
        witness = first_support(multiplier * key[1], bases)
        rows[witness] += 1
        occurrences[witness] += row_counts[key]
    return tuple(sorted(rows.items(), key=lambda item: (item[0] is None, item[0] or 0))), tuple(
        sorted(occurrences.items(), key=lambda item: (item[0] is None, item[0] or 0))
    )


def main() -> None:
    for path, expected in EXPECTED_DEPENDENCIES.items():
        require(lf_sha256(path) == expected, ("dependency", path, lf_sha256(path), expected))

    c2 = load_module("refined_six_c2", K2_PATH)
    c3 = load_module("refined_six_c3", K3_PATH)
    d6 = load_module("refined_six_d6", c2.D6_PATH)
    k3 = load_module("refined_six_k3", c3.K3_PATH)
    base = load_module("refined_six_base", BASE_PATH)

    points, atoms, _samples, clock_masks, six, six_unions = arrangement_data(base)
    _solver_points, _solver_atoms, solve, exact = c2.build_solver(base)
    require(points == _solver_points and atoms == _solver_atoms, "arrangement identity")

    before2, rowdata2, post2, engine_hash = k2_refined(c2, d6, base, points, atoms, solve)
    observed2, coindex2, examples2, _types2, semantic2 = classify_six(
        post2,
        rowdata2,
        before2["row_counts"],
        six,
        six_unions,
        EXPECTED_K2_CATEGORIES,
        EXPECTED_K2_SEMANTIC,
    )

    before3, rowdata3, post3 = k3_refined(c3, k3, base, points, atoms, solve)
    observed3, coindex3, examples3, _types3, semantic3 = classify_six(
        post3,
        rowdata3,
        before3["row_counts"],
        six,
        six_unions,
        EXPECTED_K3_CATEGORIES,
        EXPECTED_K3_SEMANTIC,
    )

    expected_examples2 = {
        "body-only": (((1, 2, 3, 4, 8, 12), 336), (336, 160), None, 2_095),
        "nontrivial-only": (
            ((1, 2, 3, 4, 7, 12), 588),
            (1_176, 358),
            (1, 2, 6, 7, 8, 10),
            181,
        ),
        "none": (((1, 2, 4, 6, 9, 10), 1_260), (2_520, 646), None, 918),
    }
    expected_examples3 = {
        "body-only": (((1, 3, 4, 5, 6, 10), 840), (840, 366), None, 388),
        "nontrivial-only": (
            ((1, 2, 5, 7, 8, 10), 1_960),
            (3_920, 1_110),
            (1, 4, 5, 7, 8, 12),
            184,
        ),
        "none": (((2, 3, 5, 7, 10, 12), 2_940), (5_880, 1_574), None, 863),
    }
    require(
        {key: examples2[key] for key in expected_examples2} == expected_examples2,
        "k2 examples",
    )
    require(
        {key: examples3[key] for key in expected_examples3} == expected_examples3,
        "k3 examples",
    )

    k2_none_seven = first_seven_cover(expected_examples2["none"][0], rowdata2, clock_masks)
    k3_none_seven = first_seven_cover(expected_examples3["none"][0], rowdata3, clock_masks)
    require(k2_none_seven == (1, 2, 3, 5, 8, 9, 10), "k2 seven witness")
    require(k3_none_seven == (1, 5, 6, 7, 8, 10, 12), "k3 seven witness")

    k2_half_eligible = sorted(
        (key for key in post2 if not supported(key[1], RANK6_BASES)),
        key=lambda key: (key[1], key[0]),
    )
    require(
        k2_half_eligible
        == [
            ((1, 2, 3, 6, 7, 13), 3_822),
            ((2, 3, 6, 7, 13, 14), 3_822),
        ],
        ("k2 half eligible", k2_half_eligible),
    )
    k2_half_hostiles = []
    maximum_3822 = maximum_transverse_half_block(3_822)
    require(maximum_3822 == (1_274, 2_548, 3), ("3822 maximum", maximum_3822))
    for key in k2_half_eligible:
        complement = complement_from_arcs(key[1], rowdata2[key][2])
        require(len(complement) > maximum_3822[0], ("3822 hostile", key, len(complement)))
        k2_half_hostiles.append(
            (key, rowdata2[key][1], len(complement), maximum_3822, before2["row_counts"][key])
        )
    require([row[2] for row in k2_half_hostiles] == [1_530, 1_560], "3822 sizes")

    k2_quarter_eligible = [key for key in post2 if not supported(2 * key[1], RANK6_BASES)]
    k3_half_rank5_eligible = [key for key in post3 if not supported(key[1], RANK5_BASES)]
    k3_half_rank6_eligible = [key for key in post3 if not supported(key[1], RANK6_BASES)]
    k3_quarter_rank5_eligible = [
        key for key in post3 if not supported(2 * key[1], RANK5_BASES)
    ]
    k3_quarter_rank6_eligible = [
        key for key in post3 if not supported(2 * key[1], RANK6_BASES)
    ]
    k3_odd_rank7_eligible = [
        key
        for key in post3
        if key[1] % 2 and not supported(key[1], ODD_RANK7_BASES)
    ]
    require(not k2_quarter_eligible, "k2 quarter base-free")
    require(not k3_half_rank5_eligible, "k3 half rank5 base-free")
    require(not k3_half_rank6_eligible, "k3 half rank6 base-free")
    require(not k3_quarter_rank5_eligible, "k3 quarter rank5 base-free")
    require(not k3_quarter_rank6_eligible, "k3 quarter rank6 base-free")
    require(not k3_odd_rank7_eligible, "k3 odd rank7 base-free")

    k2_rank6_hist = support_gate_histogram(post2, before2["row_counts"], 1, RANK6_BASES)
    k2_quarter_rank6_hist = support_gate_histogram(post2, before2["row_counts"], 2, RANK6_BASES)
    k3_rank5_hist = support_gate_histogram(post3, before3["row_counts"], 1, RANK5_BASES)
    k3_rank6_hist = support_gate_histogram(post3, before3["row_counts"], 1, RANK6_BASES)
    k3_quarter_rank5_hist = support_gate_histogram(post3, before3["row_counts"], 2, RANK5_BASES)
    k3_quarter_rank6_hist = support_gate_histogram(post3, before3["row_counts"], 2, RANK6_BASES)

    # These are the smallest aligned shapes missed by both sufficient charts.
    require(not (all(a % 2 for a in (1, 4))), "k2 half hostile")
    require(any(a % 4 == 0 for a in (1, 4)), "k2 quarter hostile")
    require(not (all(a % 2 for a in (1, 2, 4))), "k3 half hostile")
    require(any(a % 4 == 0 for a in (1, 2, 4)), "k3 quarter hostile")

    print("LRC14 REFINED SIX-POOL AND DYADIC SAMPLE STOPPING AUDIT")
    print(
        "status=FINITE-EXACT unnumbered stopping package;exact current refined key composition;"
        "no new LRC14 row terminal"
    )
    print(
        "compiler=U_D six-pool cover plus hypothetical quotient cover gives at most 13 clocks;"
        "the cited theorem stops at 12;no forced clock collision is retained by (body,D)"
    )
    print(f"k2_post=(rows,occurrences)={EXPECTED_K2_POST};categories={observed2}")
    print(f"k2_category_coindex={tuple(sorted(coindex2.items()))}")
    print(f"k2_examples={expected_examples2};k2_none_first_seven={k2_none_seven}")
    print(f"k2_semantic_sha256={semantic2}")
    print(f"k3_post=(rows,occurrences)={EXPECTED_K3_POST};categories={observed3}")
    print(f"k3_category_coindex={tuple(sorted(coindex3.items()))}")
    print(f"k3_examples={expected_examples3};k3_none_first_seven={k3_none_seven}")
    print(f"k3_semantic_sha256={semantic3}")
    print(
        "dyadic_map=half u=1/2 maps supported r to sheet r of B_D;"
        "quarter u=1/4,3/4 maps supported r to both sheets 2r,2r+1 of B_(2D)"
    )
    print(f"k2_half_rank6_basefree={k2_half_hostiles}")
    print(
        "k2_dyadic_result=no one-aux hit because every arbitrary transverse D=3822 half block has "
        f"size at most {maximum_3822[0]};k2_quarter_rank6_basefree_rows=0"
    )
    print(
        "k3_dyadic_result=half_rank5_basefree_rows=0;half_rank6_basefree_rows=0;"
        "quarter_rank5_basefree_rows=0;quarter_rank6_basefree_rows=0;"
        "odd_half_rank7_basefree_rows=0"
    )
    print(f"k2_half_rank6_support_hist={k2_rank6_hist}")
    print(f"k2_quarter_rank6_support_hist={k2_quarter_rank6_hist}")
    print(f"k3_half_rank5_support_hist={k3_rank5_hist}")
    print(f"k3_half_rank6_support_hist={k3_rank6_hist}")
    print(f"k3_quarter_rank5_support_hist={k3_quarter_rank5_hist}")
    print(f"k3_quarter_rank6_support_hist={k3_quarter_rank6_hist}")
    print(
        "parity_loss=current keys omit aligned multiplier parities;"
        "half chart needs all odd;quarter chart needs none divisible by 4;"
        "minimal missed shapes k2=(1,4),k3=(1,2,4)"
    )
    print(
        "typed_boundary=the quarter compiler retains both oriented child sheets;"
        "a base-support projection alone would lose the dyadic coset"
    )
    print(f"k2_engine_output_sha256={engine_hash};solver_cache={exact.cache_info()}")
    print(
        "dependency_lf_sha256="
        + repr(tuple((str(path.relative_to(ROOT)), digest) for path, digest in EXPECTED_DEPENDENCIES.items()))
    )
    print("verdict=PASS;new_unconditional_terminals=0;new_parity_conditional_terminals=0")


if __name__ == "__main__":
    main()
