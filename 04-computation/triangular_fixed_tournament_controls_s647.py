#!/usr/bin/env python3
"""S647: triangular fixed controls, tournaments, and Hamiltonian paths.

S644 treated even perfect numbers as fixed controls in the triangular
pair-count carrier A=C(n,2).  This continuation moves the same carrier into
tournament language.

Every tournament on m vertices has exactly C(m,2) arcs.  If one Hamiltonian
path is fixed as a base order, then m-1 arcs are forced and the remaining
C(m-1,2) arcs form the off-path deformation fiber.

At the second even-perfect control, m=8:

    C(8,2) = 28 = 7 + 21.

Thus the permanent H-gap values 7 and 21 appear as the spine length and the
off-path triangular deformation budget of a perfect arc-count control, while
the exhaustive n=8 H-spectrum still says they are not Hamiltonian-path counts.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import factorial, isqrt
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "05-knowledge" / "results"
N8_SPECTRUM = RESULTS / "h_spectrum_n8_exhaustive_monad.out"


def choose2(n: int) -> int:
    return n * (n - 1) // 2


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    d = 5
    step = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += step
        step = 6 - step
    return True


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


def sigma(n: int) -> int:
    total = 1
    for p, e in factor(n).items():
        total *= (p ** (e + 1) - 1) // (p - 1)
    return total


def aliquot(n: int) -> int:
    return sigma(n) - n


def is_power_of_two(n: int) -> bool:
    return n > 0 and n & (n - 1) == 0


def log2_int(n: int) -> int:
    return n.bit_length() - 1


def perfect_control_rows(limit_power: int = 8) -> list[tuple[int, int, int, int, int, int, str]]:
    rows = []
    for p in range(1, limit_power + 1):
        m = 1 << p
        arcs = choose2(m)
        spine = m - 1
        fiber = choose2(m - 1)
        shell = 2 * m - 1
        mersenne = m - 1
        status = "perfect" if is_prime(p) and is_prime(mersenne) else "control-only"
        rows.append((p, m, arcs, spine, fiber, shell, status))
    return rows


def adjacency_from_fixed_path_bits(n: int, bits: int) -> list[int]:
    """Tournament adjacency bitsets with base path 0->1->...->n-1 forced."""
    adj = [0] * n
    for i in range(n - 1):
        adj[i] |= 1 << (i + 1)
    k = 0
    for i in range(n):
        for j in range(i + 2, n):
            if (bits >> k) & 1:
                adj[j] |= 1 << i
            else:
                adj[i] |= 1 << j
            k += 1
    return adj


def h_count(adj: list[int]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            ways = dp[mask][last]
            if not ways:
                continue
            avail = adj[last] & ~mask
            while avail:
                bit = avail & -avail
                nxt = bit.bit_length() - 1
                avail -= bit
                dp[mask | bit][nxt] += ways
    return sum(dp[-1])


def fixed_path_hist_by_enumeration(n: int) -> Counter[int]:
    free = choose2(n - 1)
    hist: Counter[int] = Counter()
    for bits in range(1 << free):
        hist[h_count(adjacency_from_fixed_path_bits(n, bits))] += 1
    return hist


def parse_n8_labeled_h_spectrum() -> Counter[int]:
    """Parse exhaustive labelled n=8 counts from monad S1."""
    hist: Counter[int] = Counter()
    pattern = re.compile(r"^\s*H=\s*(\d+):\s*(\d+)\s*$")
    with N8_SPECTRUM.open() as f:
        for line in f:
            match = pattern.match(line)
            if match:
                hist[int(match.group(1))] = int(match.group(2))
    if not hist:
        raise RuntimeError(f"no H histogram rows parsed from {N8_SPECTRUM}")
    return hist


def fixed_path_hist_from_labeled(n: int, labeled_hist: Counter[int]) -> Counter[int]:
    """Convert labelled tournament H-counts to a fixed-base-path fiber.

    Double-count pairs (T,P) where P is a Hamiltonian path of T.  For every
    fixed labelled path P0, symmetry gives the same number of tournaments with
    P0 present and H(T)=h.  Hence:

        fixed_count_h = h * labelled_count_h / n!
    """
    denom = factorial(n)
    out: Counter[int] = Counter()
    for h, count in labeled_hist.items():
        numerator = h * count
        if numerator % denom:
            raise RuntimeError(f"nonintegral fixed-path count for H={h}")
        out[h] = numerator // denom
    return out


def hist_summary(hist: Counter[int]) -> dict[str, object]:
    total = sum(hist.values())
    h_mass = sum(h * count for h, count in hist.items())
    mode_h, mode_count = max(hist.items(), key=lambda item: (item[1], item[0]))
    odd_gaps = [
        h for h in range(min(hist), max(hist) + 1, 2)
        if h not in hist
    ]
    return {
        "total": total,
        "distinct": len(hist),
        "min": min(hist),
        "max": max(hist),
        "mean": h_mass / total,
        "mode": (mode_h, mode_count),
        "odd_gaps": odd_gaps,
    }


@dataclass(frozen=True)
class Lens:
    name: str
    section: int
    deformation: int
    proof_use: int
    h_guardrail: int
    scalar_risk: int


LENSES = [
    Lens("fixed_hamiltonian_section", 5, 4, 5, 4, 1),
    Lens("off_path_triangular_fiber", 5, 5, 5, 5, 1),
    Lens("n8_perfect_control_28_equals_7_plus_21", 4, 5, 5, 5, 1),
    Lens("h_gap_role_mismatch_guardrail", 3, 5, 5, 5, 1),
    Lens("weighted_fixed_path_fiber_from_h_spectrum", 4, 4, 4, 5, 1),
    Lens("euclid_euler_arc_count_controls", 3, 4, 4, 3, 1),
    Lens("raw_7_21_numerology", 1, 2, 2, 1, 5),
]


def lens_score(lens: Lens) -> int:
    return (
        3 * lens.section
        + 3 * lens.deformation
        + 2 * lens.proof_use
        + 2 * lens.h_guardrail
        - 3 * lens.scalar_risk
    )


def beats(a: Lens, b: Lens) -> bool:
    sa = lens_score(a)
    sb = lens_score(b)
    if sa != sb:
        return sa > sb
    return a.name < b.name


def directed_triangles(lenses: list[Lens]) -> int:
    total = 0
    n = len(lenses)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                ij = beats(lenses[i], lenses[j])
                jk = beats(lenses[j], lenses[k])
                ki = beats(lenses[k], lenses[i])
                if ij and jk and ki:
                    total += 1
                ji = not ij
                kj = not jk
                ik = not ki
                if ji and kj and ik:
                    total += 1
    return total


def hamiltonian_paths_lens_tournament(lenses: list[Lens]) -> int:
    n = len(lenses)
    dp: dict[tuple[int, int], int] = {(1 << i, i): 1 for i in range(n)}
    for mask in range(1 << n):
        for last in range(n):
            ways = dp.get((mask, last), 0)
            if not ways:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if beats(lenses[last], lenses[nxt]):
                    key = (mask | (1 << nxt), nxt)
                    dp[key] = dp.get(key, 0) + ways
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def score_hist(lenses: list[Lens]) -> dict[int, int]:
    scores = {
        lens.name: sum(1 for other in lenses if lens is not other and beats(lens, other))
        for lens in lenses
    }
    return dict(sorted(Counter(scores.values()).items()))


def print_header() -> None:
    print("S647 triangular fixed controls and Hamiltonian paths")
    print("====================================================")
    print()
    print("Core identities")
    print("---------------")
    print("Tournament arc count: A_m = C(m,2)")
    print("Fixed Hamiltonian path spine: m-1 forced arcs")
    print("Off-path deformation fiber: C(m,2)-(m-1)=C(m-1,2) free arcs")
    print("Fixed-path fiber size: 2^C(m-1,2)")
    print()


def print_perfect_controls() -> None:
    print("Perfect triangular tournament controls")
    print("--------------------------------------")
    print("  p  m=2^p  arcs=C(m,2)  spine=m-1  offpath=C(m-1,2)  shell=2m-1  status")
    for p, m, arcs, spine, fiber, shell, status in perfect_control_rows():
        print(f"  {p:<2} {m:<6} {arcs:<12} {spine:<10} {fiber:<18} {shell:<10} {status}")
    print()
    print("The exact perfect controls occur when 2^p-1 is Mersenne prime.")
    print("At m=8, the perfect scalar splits as 28 = 7 + 21.")
    print("Those are dimensions of a path section and its deformation fiber,")
    print("not Hamiltonian-path counts.")
    print()


def print_fixed_path_spectra() -> None:
    print("Fixed-base-path deformation spectra")
    print("-----------------------------------")
    print("Rows n=3..7 are enumerated over the 2^C(n-1,2) fixed-path fiber.")
    print("Row n=8 is reconstructed from the full labelled n=8 H-spectrum by")
    print("fixed_count_h = h * labelled_count_h / 8! .")
    print()
    print("  n  free_bits  fiber_size  distinct_H  H_range    mean_H      mode       gaps_low")
    fixed_hists: dict[int, Counter[int]] = {}
    for n in range(3, 8):
        fixed_hists[n] = fixed_path_hist_by_enumeration(n)
    fixed_hists[8] = fixed_path_hist_from_labeled(8, parse_n8_labeled_h_spectrum())
    for n in range(3, 9):
        hist = fixed_hists[n]
        summary = hist_summary(hist)
        low_gaps = [h for h in summary["odd_gaps"] if h <= 65]
        mode_h, mode_count = summary["mode"]
        print(
            f"  {n:<2} {choose2(n-1):<10} {summary['total']:<11} "
            f"{summary['distinct']:<11} [{summary['min']},{summary['max']}]"
            f"{summary['mean']:>11.4f}  {mode_h}:{mode_count:<7} {low_gaps}"
        )
    print()
    hist8 = fixed_hists[8]
    summary8 = hist_summary(hist8)
    print("n=8 perfect-control details")
    print("---------------------------")
    print(f"  total fixed-path tournaments = {summary8['total']} = 2^21")
    print(f"  distinct H values = {summary8['distinct']}")
    print(f"  H range = [{summary8['min']}, {summary8['max']}]")
    print(f"  mean H in fixed-path fiber = {summary8['mean']:.4f}")
    print(f"  fixed-fiber count with H=7: {hist8.get(7, 0)}")
    print(f"  fixed-fiber count with H=21: {hist8.get(21, 0)}")
    print("  selected fixed-fiber counts:")
    for h in (1, 3, 5, 7, 21, 35, 39, 49, 63, 189, 661):
        print(f"    H={h:<3} count={hist8.get(h, 0)}")
    print()
    print("Double-counting check")
    print("---------------------")
    print("Every tournament with H(T)=h contributes h labelled Hamiltonian paths.")
    print("Dividing by 8! gives the exact count in one fixed-base-path fiber.")
    print("The reconstructed total is 2^21, matching the off-path triangular budget.")
    print()


def print_role_mismatch() -> None:
    print("Role mismatch: count versus dimension")
    print("-------------------------------------")
    print("H=7 and H=21 are forbidden as Hamiltonian-path counts.")
    print("But at the m=8 perfect tournament control:")
    print("  7  = the forced spine length of a Hamiltonian path")
    print("  21 = the number of off-path arcs after choosing that spine")
    print("  28 = the perfect total arc count")
    print()
    print("This suggests a guardrail for LRC/unit-distance transfers:")
    print("before treating a visible scalar as a target count, ask whether it is")
    print("actually a section length, a fiber dimension, a shell modulus, or a")
    print("side-channel mass.")
    print()


def print_tournament_analysis() -> None:
    print("Tournament Analysis over proof lenses")
    print("-------------------------------------")
    lenses = list(LENSES)
    scores = {
        lens.name: sum(1 for other in lenses if lens is not other and beats(lens, other))
        for lens in lenses
    }
    ranking = sorted(lenses, key=lambda lens: (-scores[lens.name], lens.name))
    print(f"vertices={len(lenses)}")
    print(f"score_hist={score_hist(lenses)}")
    print(f"directed_3cycles={directed_triangles(lenses)}")
    print(f"hamiltonian_paths={hamiltonian_paths_lens_tournament(lenses)}")
    print("tie Hamiltonian path:")
    for lens in ranking:
        print(f"  score={scores[lens.name]} {lens.name}")
    print()
    print("Pairwise observable: which lens best preserves the distinction between")
    print("Hamiltonian-path count, fixed section, off-path fiber, and scalar control.")
    print("Switch/gauge: proof transfer value minus scalar-numerology risk.")
    print("Vertices are proof lenses, not tournament vertices or LRC runners.")
    print()


def print_assumption_challenge() -> None:
    print("Assumption challenge")
    print("--------------------")
    print("Alternate vertices considered: tournament vertices, arcs, fixed path edges,")
    print("off-path cells, Hamiltonian paths, conflict-graph cycles, LRC shell residues,")
    print("unit-distance spine/bulk obligations, and proof lenses.")
    print("Chosen computation: off-path cells are the finite tiling variables, while")
    print("Tournament Analysis uses proof lenses.")
    print("Preserved information: path-section size, deformation-fiber size, exact")
    print("fixed-path H-spectrum through the perfect m=8 control, and the 7/21 gap")
    print("guardrail.")
    print("Destroyed information: tournament isomorphism class, endpoint distribution,")
    print("SCC decomposition, Omega side channels, and LRC carry/owner labels.")
    print("Challenged assumption: a scalar appearing in several problems has the same")
    print("role in each. Here 7 and 21 are forbidden H-counts but natural dimensions.")


def main() -> None:
    print_header()
    print_perfect_controls()
    print_fixed_path_spectra()
    print_role_mismatch()
    print_tournament_analysis()
    print_assumption_challenge()


if __name__ == "__main__":
    main()
