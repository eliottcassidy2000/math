#!/usr/bin/env python3
"""Exact bounded audit for the n=12 nonextremal-core (sporadic) branch.

Scope is deliberately small: primitive eleven-speed cores in ``{1,...,13}``.
For every nonextremal core P, compute ``mu=M(P)`` exactly and enumerate every
top speed allowed by THM-759,

    w <= floor(max(P)/(13*mu-1)).

The exact M evaluator includes pair sums, pair differences, and single-runner
half-turn cusps.  THM-766's first-window cones are checked independently.

Tournament Analysis is not used as a decision layer here: pairwise
orientations of runners, cores, or tooth labels do not preserve the exact
predicate ``M(P union {w})=1/13``.  The faithful object is component-tooth
incidence with rational endpoints; THM-766 records its preserve/destroy audit.
This script is a bounded regression, not the uniform branch theorem.
"""

from collections import Counter
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd

from lrc14_certificates import M_exact


TARGET = F(1, 13)
CORE_FLOOR = F(1, 12)


def primitive(values: tuple[int, ...]) -> bool:
    return reduce(gcd, values) == 1


def cap_for(core: tuple[int, ...], mu: F) -> int:
    return max(core) // (13 * mu - 1)


def denominator_quantized_cap(b: int) -> int:
    """Universal sporadic cap from THM-668 + THM-759.

    If ``mu=m/q>1/12`` and the pair-sum ruler has ``q<=2b``, then
    ``mu>=1/12+1/(24b)`` and hence
    ``w<=floor(24b^2/(2b+13))``.
    """
    return (24 * b * b) // (2 * b + 13)


def tooth_labels(m: int, b: int, w: int) -> tuple[int, ...]:
    return tuple(
        k
        for k in range(1, 12)
        if (13 * k - 1) * m <= w and 12 * w <= (13 * k + 1) * b
    )


def main() -> None:
    cores = []
    for core in combinations(range(1, 14), 11):
        if not primitive(core):
            continue
        mu = M_exact(core)
        if mu > CORE_FLOOR:
            cores.append((core, mu))

    branch_cores = Counter()
    branch_candidates = Counter()
    narrow_label_hist = Counter()
    candidates = []
    max_cap_ratio = (F(0), None)
    coarse_ratio_candidates = 0
    denominator_cap_candidates = 0

    for core, mu in cores:
        m, b = min(core), max(core)
        cap = cap_for(core, mu)
        denominator_cap = denominator_quantized_cap(b)
        assert cap <= denominator_cap < 12 * b
        ratio = F(cap, b)
        if ratio > max_cap_ratio[0]:
            max_cap_ratio = (ratio, (core, mu, cap))

        branch = "wide" if b >= 12 * m else "narrow"
        branch_cores[branch] += 1
        coarse_ratio_candidates += sum(
            primitive(core + (w,)) for w in range(b + 1, 12 * b)
        )
        denominator_cap_candidates += sum(
            primitive(core + (w,)) for w in range(b + 1, denominator_cap + 1)
        )
        for w in range(b + 1, cap + 1):
            speeds = core + (w,)
            if not primitive(speeds):
                continue
            labels = tooth_labels(m, b, w) if branch == "narrow" else ()
            if branch == "narrow":
                narrow_label_hist[labels] += 1
            branch_candidates[branch] += 1
            candidates.append((speeds, M_exact(speeds)))

    tight = [speeds for speeds, margin in candidates if margin == TARGET]
    minimum = min(margin for _, margin in candidates)
    minimizers = [speeds for speeds, margin in candidates if margin == minimum]

    assert len(cores) == 77
    assert len(candidates) == 790
    assert coarse_ratio_candidates == 10813
    assert denominator_cap_candidates == 6897
    assert branch_cores == Counter({"wide": 65, "narrow": 12})
    assert branch_candidates == Counter({"wide": 750, "narrow": 40})
    assert narrow_label_hist == Counter({(): 40})
    assert not tight
    assert minimum == F(1, 12)
    assert len(minimizers) == 3
    assert max_cap_ratio[0] == F(11, 2)

    print("n=12 sporadic exact-cap audit (bounded core box {1,...,13})")
    print(f"nonextremal primitive 11-cores: {len(cores)}")
    print(f"core branches: wide={branch_cores['wide']} narrow={branch_cores['narrow']}")
    print(f"THM-759-capped primitive completions: {len(candidates)}")
    print(
        "uniform top-speed caps on this same core bank: "
        f"coarse w<12b gives {coarse_ratio_candidates}; "
        "pair-ruler quantization w<=floor(24b^2/(2b+13)) gives "
        f"{denominator_cap_candidates}"
    )
    print(
        "candidate branches: "
        f"wide={branch_candidates['wide']} narrow={branch_candidates['narrow']}"
    )
    print("narrow candidates with a THM-766 tooth label: 0/40")
    print("tight completions M=1/13: 0")
    print(f"minimum M in the bank: {minimum}; minimizers={minimizers}")
    ratio, (core, mu, cap) = max_cap_ratio
    print(
        "largest exact cap/max(core) ratio: "
        f"{ratio} at core={core}, M(core)={mu}, cap={cap}"
    )
    print("scope: bounded 77-core bank only; THM-763 gives uniform finiteness, not this enumeration")


if __name__ == "__main__":
    main()
