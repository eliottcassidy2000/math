#!/usr/bin/env python3
"""Validate the exploratory 47/103 Schoof stratifier on a lower-rank box.

The five known 5-rank-five specializations were used to *discover* the
47/103 bad-fibre label, so their 4-of-5 concentration is post-selected.  This
script applies the already frozen label to a disjoint response tier: every
reduced positive rational x/y with x,y <= B.  It computes each distinct
imaginary quadratic field once with PARI, then prints both specialization-
level and field-level rank histograms split by the bad-fibre label.

PARI ``quadclassunit`` is a GRH-assisted discovery path here.  The script
makes no negative class-group theorem and emits no challenge certificate.

Reproduction:

    python3 04-computation/classrank_schoof_bad_fibre_validation_20260824.py \
      --bound 200 --pari-gp /opt/homebrew/bin/gp
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from math import gcd
from pathlib import Path
import subprocess
import sys
from fractions import Fraction
from typing import DefaultDict, Dict, Iterable, Sequence, Set, Tuple


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from classrank_schoof_rank5_bad_fibre_scout_20260824 import (  # noqa: E402
    BAD_PRIMES,
    fundamental_discriminant_abs,
    label_at,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def is_bad(numerator: int, denominator: int) -> bool:
    return any(
        label_at(numerator, denominator, prime) != "good"
        for prime in BAD_PRIMES
    )


def gp_ranks(gp_path: str, discriminants: Sequence[int]) -> Dict[int, int]:
    literal = ",".join(map(str, discriminants))
    program = f"""
default(parisizemax,2000000000); allocatemem(128000000);
C=[{literal}];
for(i=1,#C,{{
  D=C[i];
  if(!isfundamental(-D),error("nonfundamental discriminant"));
  q=quadclassunit(-D);
  rk=sum(j=1,#q.cyc,valuation(q.cyc[j],5)>0);
  print("ROW|",D,"|",rk)
}});
print("DONE|",#C);
quit
"""
    process = subprocess.run(
        [gp_path, "-fq"],
        input=program,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    require(process.returncode == 0, f"PARI exited {process.returncode}")
    errors = [
        row for row in process.stdout.splitlines()
        if row.lstrip().startswith("***") and "Warning:" not in row
    ]
    require(not errors, "PARI errors:\n" + "\n".join(errors))
    ranks: Dict[int, int] = {}
    done_rows = []
    for row in process.stdout.splitlines():
        if row.startswith("ROW|"):
            _, raw_D, raw_rank = row.split("|")
            D = int(raw_D)
            require(D not in ranks, f"duplicate PARI row for D={D}")
            ranks[D] = int(raw_rank)
        elif row.startswith("DONE|"):
            done_rows.append(int(row.split("|", 1)[1]))
    require(done_rows == [len(discriminants)], "missing or bad PARI DONE row")
    require(len(ranks) == len(discriminants), "missing PARI rank rows")
    return ranks


def histogram(values: Iterable[int]) -> str:
    counts = Counter(values)
    return ",".join(f"{rank}:{counts[rank]}" for rank in sorted(counts))


def threshold_summary(values: Sequence[int]) -> str:
    return ",".join(
        f"ge{threshold}:{sum(rank >= threshold for rank in values)}"
        for threshold in (2, 3, 4, 5, 6)
    )


def rate_ratio(bad: Sequence[int], good: Sequence[int], threshold: int) -> str:
    bad_hits = sum(rank >= threshold for rank in bad)
    good_hits = sum(rank >= threshold for rank in good)
    if not bad or not good or good_hits == 0:
        return "undefined"
    ratio = Fraction(bad_hits * len(good), good_hits * len(bad))
    return f"{ratio}={float(ratio):.12f}"


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bound", type=int, default=200)
    parser.add_argument("--pari-gp", required=True)
    args = parser.parse_args(argv)
    require(args.bound >= 1, "bound must be positive")

    specializations = []
    labels_by_D: DefaultDict[int, Set[bool]] = defaultdict(set)
    multiplicity_by_D: Counter[int] = Counter()
    squareful_specializations = 0
    for numerator in range(1, args.bound + 1):
        for denominator in range(1, args.bound + 1):
            if gcd(numerator, denominator) != 1:
                continue
            D, square_factor = fundamental_discriminant_abs(
                numerator, denominator
            )
            bad = is_bad(numerator, denominator)
            specializations.append((D, bad))
            labels_by_D[D].add(bad)
            multiplicity_by_D[D] += 1
            squareful_specializations += square_factor != 1

    discriminants = tuple(sorted(labels_by_D))
    ranks = gp_ranks(args.pari_gp, discriminants)

    spec_bad = [ranks[D] for D, bad in specializations if bad]
    spec_good = [ranks[D] for D, bad in specializations if not bad]
    field_bad = [
        ranks[D] for D, labels in labels_by_D.items() if labels == {True}
    ]
    field_good = [
        ranks[D] for D, labels in labels_by_D.items() if labels == {False}
    ]
    field_mixed = [
        ranks[D] for D, labels in labels_by_D.items() if len(labels) > 1
    ]

    print("SCHOOF 47/103 BAD-FIBRE LOWER-TIER VALIDATION -- 2026-08-24")
    print("STATUS=FINITE-COMPUTED_GRH_ASSISTED_DISCOVERY")
    print(
        f"UNIVERSE=reduced_positive_rationals;1<=numerator,denominator<="
        f"{args.bound};bad_primes={BAD_PRIMES}"
    )
    print(
        f"COUNTS=specializations:{len(specializations)};"
        f"distinct_fields:{len(discriminants)};"
        f"duplicate_specializations:{len(specializations)-len(discriminants)};"
        f"squareful_specializations:{squareful_specializations};"
        f"max_field_multiplicity:{max(multiplicity_by_D.values(), default=0)}"
    )
    print(
        f"SPECIALIZATION_BAD=n:{len(spec_bad)};hist:{histogram(spec_bad)};"
        f"thresholds:{threshold_summary(spec_bad)}"
    )
    print(
        f"SPECIALIZATION_GOOD=n:{len(spec_good)};hist:{histogram(spec_good)};"
        f"thresholds:{threshold_summary(spec_good)}"
    )
    print(
        f"FIELD_BAD_ONLY=n:{len(field_bad)};hist:{histogram(field_bad)};"
        f"thresholds:{threshold_summary(field_bad)}"
    )
    print(
        f"FIELD_GOOD_ONLY=n:{len(field_good)};hist:{histogram(field_good)};"
        f"thresholds:{threshold_summary(field_good)}"
    )
    print(
        f"FIELD_MIXED=n:{len(field_mixed)};hist:{histogram(field_mixed)};"
        f"thresholds:{threshold_summary(field_mixed)}"
    )
    print(
        "FIELD_BAD_TO_GOOD_RATE_RATIO="
        f"ge3:{rate_ratio(field_bad, field_good, 3)};"
        f"ge4:{rate_ratio(field_bad, field_good, 4)}"
    )
    print(f"MAX_DISCOVERED_5_RANK={max(ranks.values(), default=-1)}")
    print(
        "DISTINCT_RANK4_FIELDS="
        f"{sum(rank >= 4 for rank in ranks.values())}"
    )
    print(f"PARI_FIELDS_CHECKED={len(ranks)}")
    print(
        "INTERPRETATION=lower_rank_out_of_sample_response_for_a_post_selected_"
        "residue_label;compare_rates_descriptively_and_retest_on_new_high_box"
    )
    print(
        "BOUNDARY=quadclassunit_full_ranks_are_GRH_assisted;no_rank_six_claim;"
        "no_negative_theorem;no_unadjusted_p_value"
    )


if __name__ == "__main__":
    main()
