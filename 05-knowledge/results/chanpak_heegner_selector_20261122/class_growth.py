"""Factor the twelve certified endpoint heights and check the large remaining prime witness."""
from math import isqrt, prod
from pathlib import Path
from runpy import run_path

from cm_crossings import CHARTS, verify_literal
from e15_pair import p, strong

# Repeated prime factors are listed separately.
FACTORS = {
    33506: ((5, 6701), (3, 3, 3, 17, 73)),
    33512: ((23, 31, 47), (3, 11171)),
    885480: ((7, 7, 17, 1063), (941, 941)),
    884736744: ((11, 787, 102199), (5, 23, 31, 41, 6053)),
    147197952744: ((147197952743,), (5, 47, 59, 10616513)),
    262537412640768744: ((3677, 4787, 91433, 163129), (5, 103, 509781383768483)),
}


def is_prime(n):
    return n == 2 or (n > 2 and n % 2 == 1 and all(n % d for d in range(3, isqrt(n)+1, 2)))


def neighbor_scores(rows):
    scores = [r.bit_count() for r in rows]
    return sorted((scores[i], tuple(sorted(s for j, s in enumerate(scores) if r >> j & 1)))
                  for i, r in enumerate(rows))


def main():
    results = Path(__file__).resolve().parents[1]
    for prime in {p for pair in FACTORS.values() for factors in pair for p in factors}:
        assert is_prime(prime), prime
    for center, pair in FACTORS.items():
        assert [prod(factors) for factors in pair] == [center-1, center+1]
        print(center, "eventual degrees",
              [len(factors) - int(7 in factors) for factors in pair])

    # The F1 bridge formula in the targeted two-apex construction.
    b, e = 29, 12179923056399
    assert b*e+(1+p(b-1))*(1+p(e-1))-(b&e&1) == 509781383768483
    assert (b-1).bit_length()+(e-1).bit_length()+2 == 51
    incoming = run_path(str(results / "chanpak_heegner_absolute_targets_20261123/heegner-witness.py"))
    for prime in FACTORS[262537412640768744][0]:
        rows = incoming["CORES"][prime]
        n = len(rows)
        reverse = tuple(r ^ ((1 << n)-1) ^ (1 << i) for i, r in enumerate(rows))
        assert strong(rows) and neighbor_scores(rows) != neighbor_scores(reverse)
    # Lower and upper 58-vertex witnesses add a different one-component construction.
    chart = CHARTS[-1]
    verify_literal(chart, chart.first, chart.second, chart.cube_root**3 + 744)
    print("Distinct strong factor orientations at all four lower primes; upper large factor on 51 vertices")


if __name__ == "__main__":
    main()
