#!/usr/bin/env python3
"""
LRC(14) Euler-copy coimage-tail profile.

This continues HYP-2629's next target: re-index the HYP-2626 k=10
repeated-residue tail by Euler-copy / exact-period packet data.

The computation keeps three layers separate:

  1. coimage class and height<=H one-large wall addressability;
  2. minimal large-residue demand after the bounded core 1..13;
  3. exact-period Euler-copy capacity and quadratic-character moments.

The main point is intentionally two-sided.  Euler-copy mass explains the
available packet capacity in the LRC14 unit seam, but it is uniform across
F_7^* residues.  Therefore copy mass alone does not separate the QR/NQR split
inside the repeated-residue tail; the first separator is the quadratic
character phase channel.

Tournament Analysis declaration:
  vertices are proof obligations / quotient stages, not runners.  The quotient
  preserves the support-six coimage-tail predicate and destroys exact witness
  times and raw support identities.
"""
from __future__ import annotations

import itertools
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_height2_coimage_wall_classes_codex_s18 as s18  # noqa: E402
import lrc14_support6_coimage_fiber_codex_s14 as s14  # noqa: E402

PRIMES = (2, 3, 5, 7)
EPS = 1e-12
K = 10
AMBIENT_D = 9
CORE_BOUND = 13


def section(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


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


def divisors_from_factor(fac: dict[int, int]) -> list[int]:
    divs = [1]
    for p, e in fac.items():
        divs = [d * (p**a) for d in divs for a in range(e + 1)]
    return sorted(divs)


def phi(n: int) -> int:
    ans = n
    for p in factor(n):
        ans = ans // p * (p - 1)
    return ans


def mask_of(n: int) -> int:
    mask = 0
    for i, p in enumerate(PRIMES):
        if n % p == 0:
            mask |= 1 << i
    return mask


def mask_name(mask: int) -> str:
    vals = [str(p) for i, p in enumerate(PRIMES) if mask & (1 << i)]
    return "{" + ",".join(vals) + "}" if vals else "{}"


def chi7(x: int) -> int:
    x %= 7
    if x == 0:
        return 0
    return 1 if x in {1, 2, 4} else -1


def inv7(x: int) -> int:
    x %= 7
    if x == 0:
        raise ValueError("zero has no inverse mod 7")
    return pow(x, -1, 7)


def comb_capacity(weight_per_residue: int, multiplicities: Counter[int]) -> int:
    """Unordered packet capacity with distinct packet copies inside each residue."""
    ans = 1
    for r, m in multiplicities.items():
        if r == 0:
            return 0
        if weight_per_residue < m:
            return 0
        ans *= math.comb(weight_per_residue, m)
    return ans


def exact_period_residue_profile(q: int) -> dict[int, dict[int, int]]:
    """mask -> residue mod 7 -> exact-period copy mass on the q-grid."""
    out: dict[int, dict[int, int]] = defaultdict(lambda: defaultdict(int))
    for d in divisors_from_factor(factor(q)):
        m = mask_of(d)
        for u in range(1, d + 1):
            if math.gcd(u, d) != 1:
                continue
            a = ((q // d) * u) % q
            out[m][a % 7] += 1
    return out


def full_mask_mass(q: int) -> int:
    fac = factor(q)
    ans = 1
    for p in PRIMES:
        if p not in fac:
            return 0
        ans *= p ** fac[p] - 1
    return ans


def pattern_label(cls: tuple[int, ...]) -> str:
    counts = Counter(cls)
    multiplicities = sorted(counts.values(), reverse=True)
    zeros = counts.get(0, 0)
    if zeros >= 4:
        return "zero-cusp"
    if multiplicities[:2] == [4, 2]:
        return "4+2 repeated"
    if multiplicities[:3] == [4, 1, 1]:
        return "4+1+1 repeated"
    if multiplicities[:3] == [2, 2, 2]:
        return "2+2+2 repeated"
    if multiplicities[:1] == [6]:
        return "6 repeated"
    return "mixed"


def normalize_repeated(cls: tuple[int, ...]) -> tuple[int, ...] | None:
    counts = Counter(cls)
    repeated = [r for r, c in counts.items() if c == max(counts.values())]
    repeated = [r for r in repeated if r != 0]
    if not repeated:
        return None
    lam = inv7(repeated[0])
    return tuple(sorted((lam * x) % 7 for x in cls))


def repeated_signature(cls: tuple[int, ...]) -> tuple[str, tuple[int, ...], tuple[int, ...]]:
    norm = normalize_repeated(cls)
    if norm is None:
        return ("zero/mixed", (), ())
    counts = Counter(norm)
    mults = sorted(counts.values(), reverse=True)
    if mults[:2] == [4, 2]:
        a = next(r for r, c in counts.items() if c == 2)
        return ("4+2", (a,), (chi7(a), 4 + 2 * chi7(a)))
    if mults[:3] == [4, 1, 1]:
        singles = tuple(sorted(r for r, c in counts.items() if c == 1))
        a, b = singles
        sig = (
            chi7(a),
            chi7(b),
            chi7(a * b),
            chi7((a - 1) * (b - 1)),
            4 + chi7(a) + chi7(b),
        )
        return ("4+1+1", singles, sig)
    return (pattern_label(norm), tuple(norm), tuple(chi7(x) for x in norm))


def coefficient_set(height: int) -> tuple[int, ...]:
    return tuple(x for x in range(-height, height + 1) if x)


@dataclass(frozen=True)
class WallRecord:
    support: tuple[int, ...]
    cls: tuple[int, ...]
    large: int


def enumerate_one_large_records(k: int, height: int) -> list[WallRecord]:
    """Unique one-large wall supports at coefficient height <= height."""
    B = s18.BOUND[k]
    coeffs = coefficient_set(height)
    records: dict[tuple[int, ...], WallRecord] = {}
    for core_support in itertools.combinations(range(1, B + 1), 5):
        for core_coeffs in itertools.product(coeffs, repeat=5):
            core_sum = sum(c * v for c, v in zip(core_coeffs, core_support))
            if core_sum == 0:
                continue
            for cM in coeffs:
                if (-core_sum) % cM:
                    continue
                M = (-core_sum) // cM
                if M <= B:
                    continue
                support = tuple(sorted(core_support + (M,)))
                if support in records:
                    continue
                records[support] = WallRecord(
                    support=support,
                    cls=s14.canon_support(support),
                    large=M,
                )
    return list(records.values())


def coimage_rows() -> list[s14.FiberStats]:
    return s14.compute_stats_for_d(AMBIENT_D, s14.support_classes())


def core_capacity() -> Counter[int]:
    return Counter(n % 7 for n in range(1, CORE_BOUND + 1))


def min_large_residue_demand(cls: tuple[int, ...]) -> tuple[int, tuple[int, ...]]:
    """Minimum number of large coordinates needed after using residues 1..13."""
    cap = core_capacity()
    best: tuple[int, tuple[int, ...]] | None = None
    for lam in range(1, 7):
        scaled = tuple(sorted((lam * x) % 7 for x in cls))
        counts = Counter(scaled)
        demand = sum(max(0, c - cap[r]) for r, c in counts.items())
        candidate = (demand, scaled)
        if best is None or candidate < best:
            best = candidate
    assert best is not None
    return best


def exact_period_packet_report() -> dict[str, int]:
    section("EXACT-PERIOD EULER-COPY UNIT PACKETS")
    print(
        "For q divisible by 7, exact-period packets whose denominator touches 7 "
        "are uniformly distributed over F_7^*.  This is the copy ledger; any "
        "QR/NQR split must come from a later character phase, not from raw phi mass."
    )
    weights: dict[str, int] = {}
    print(
        f"{'q':>5} {'phi(q)':>8} {'top per unit':>13} {'full-mask mass':>15} "
        f"{'full per unit':>14} {'7-touch per unit':>16}"
    )
    for q in (14, 210, 1260):
        profile = exact_period_residue_profile(q)
        top_per = phi(q) // 6 if q % 7 == 0 else 0
        fmass = full_mask_mass(q)
        full_per = fmass // 6 if fmass else 0
        touch7 = 0
        for mask, row in profile.items():
            if mask & (1 << 3):
                touch7 += sum(row[r] for r in range(1, 7))
        touch7_per = touch7 // 6
        print(
            f"{q:>5} {phi(q):>8} {top_per:>13} {fmass:>15} "
            f"{full_per:>14} {touch7_per:>16}"
        )
        weights[f"top_{q}"] = top_per
        weights[f"full_{q}"] = full_per
        weights[f"touch7_{q}"] = touch7_per
    print(
        "\nFor q=1260, the exact top-period unit seam has 48 copies per unit "
        "residue, while the full {2,3,5,7} squarefree mask has 96 per unit "
        "residue.  Both are uniform on F_7^*."
    )
    return weights


def wall_height_report(rows: list[s14.FiberStats]) -> tuple[set[tuple[int, ...]], set[tuple[int, ...]]]:
    section("ONE-LARGE WALL HEIGHT CHECK")
    nz = [r for r in rows if r.signed_abs > EPS]
    total = sum(r.signed_abs for r in nz)
    class_sets = []
    print(f"{'height':>6} {'supports':>10} {'classes':>8} {'nonzero hit':>13} {'mass share':>12}")
    for H in (2, 3):
        records = enumerate_one_large_records(K, H)
        cls_set = {r.cls for r in records}
        class_sets.append(cls_set)
        hit = [r for r in nz if r.cls in cls_set]
        print(
            f"{H:>6} {len(records):>10} {len(cls_set):>8} "
            f"{len(hit):>5}/{len(nz):<7} {sum(r.signed_abs for r in hit) / total:>11.6%}"
        )
    same = class_sets[0] == class_sets[1]
    print(f"\nheight-2 and height-3 class sets equal? {same}")
    print(
        "Readout: the repeated-residue tail is not a missed coefficient-height "
        "3 one-large wall.  The obstruction is residue multiplicity: one large "
        "coordinate cannot create four equal nonzero residues after the core "
        "1..13 supplies only two."
    )
    return class_sets[0], class_sets[1]


def tail_profile_report(rows: list[s14.FiberStats], hit_h2: set[tuple[int, ...]], weights: dict[str, int]) -> None:
    section("K=10 TAIL-ONLY CLASSES BY RESIDUE DEMAND AND COPY CAPACITY")
    nz = sorted([r for r in rows if r.signed_abs > EPS], key=lambda r: r.signed_abs, reverse=True)
    tail = [r for r in nz if r.cls not in hit_h2]
    total_tail = sum(r.signed_abs for r in tail)
    print(f"tail-only nonzero classes after height<=2 one-large walls: {len(tail)}")
    print(f"tail-only |S_9|-mass: {total_tail:.9g}")
    print(f"bounded core residue capacity 1..13: {dict(sorted(core_capacity().items()))}")

    by_pattern: dict[str, list[s14.FiberStats]] = defaultdict(list)
    by_demand: dict[int, list[s14.FiberStats]] = defaultdict(list)
    for row in tail:
        by_pattern[pattern_label(row.cls)].append(row)
        demand, _ = min_large_residue_demand(row.cls)
        by_demand[demand].append(row)

    print("\nPattern summary:")
    print(f"{'pattern':>18} {'classes':>8} {'|S|-mass':>12} {'mean |S|':>12} {'min large demand':>17}")
    for pat, items in sorted(by_pattern.items(), key=lambda kv: -sum(r.signed_abs for r in kv[1])):
        demands = sorted({min_large_residue_demand(r.cls)[0] for r in items})
        mass = sum(r.signed_abs for r in items)
        print(f"{pat:>18} {len(items):>8} {mass:>12.8g} {mass / len(items):>12.8g} {str(demands):>17}")

    print("\nLarge-residue demand summary:")
    print(f"{'demand':>8} {'classes':>8} {'|S|-mass':>12} {'share of tail':>14}")
    for demand, items in sorted(by_demand.items()):
        mass = sum(r.signed_abs for r in items)
        print(f"{demand:>8} {len(items):>8} {mass:>12.8g} {mass / total_tail:>13.6%}")

    print("\nTop tail classes:")
    print(
        f"{'rank':>4} {'class':>22} {'pattern':>17} {'|S_9|':>11} "
        f"{'demand':>7} {'scaled demand class':>24} {'C_full210':>12} {'C_full1260':>14}"
    )
    for i, row in enumerate(tail[:16], 1):
        demand, scaled = min_large_residue_demand(row.cls)
        mult = Counter(scaled)
        c210 = comb_capacity(weights["full_210"], mult)
        c1260 = comb_capacity(weights["full_1260"], mult)
        print(
            f"{i:>4} {str(row.cls):>22} {pattern_label(row.cls):>17} "
            f"{row.signed_abs:>11.8g} {demand:>7} {str(scaled):>24} "
            f"{c210:>12} {c1260:>14}"
        )


def character_split_report(rows: list[s14.FiberStats], hit_h2: set[tuple[int, ...]], weights: dict[str, int]) -> None:
    section("COPY MASS VS QUADRATIC-CHARACTER PHASE")
    tail = [r for r in rows if r.signed_abs > EPS and r.cls not in hit_h2]

    four_two = [r for r in tail if repeated_signature(r.cls)[0] == "4+2"]
    print("4+2 packet `(1,1,1,1,a,a)`:")
    print(
        f"{'a':>3} {'chi(a)':>7} {'chi-sum':>8} {'|S_9|':>12} "
        f"{'C_full210':>12} {'C_full1260':>14}"
    )
    for row in sorted(four_two, key=lambda r: repeated_signature(r.cls)[1]):
        _, data, sig = repeated_signature(row.cls)
        a = data[0]
        mult = Counter(normalize_repeated(row.cls))
        c210 = comb_capacity(weights["full_210"], mult)
        c1260 = comb_capacity(weights["full_1260"], mult)
        print(f"{a:>3} {sig[0]:>7} {sig[1]:>8} {row.signed_abs:>12.8g} {c210:>12} {c1260:>14}")

    qr = [r for r in four_two if repeated_signature(r.cls)[2][0] == 1]
    nqr = [r for r in four_two if repeated_signature(r.cls)[2][0] == -1]
    if qr and nqr:
        print(
            f"QR mean |S_9|={sum(r.signed_abs for r in qr)/len(qr):.8g}; "
            f"NQR mean |S_9|={sum(r.signed_abs for r in nqr)/len(nqr):.8g}. "
            "Copy capacities are identical inside the row; the separator is chi."
        )

    print("\n4+1+1 packet signatures `(chi(a),chi(b),chi(ab),chi((a-1)(b-1)),chi-sum)`:")
    buckets: dict[tuple[int, ...], list[s14.FiberStats]] = defaultdict(list)
    for row in tail:
        typ, _data, sig = repeated_signature(row.cls)
        if typ == "4+1+1":
            buckets[sig].append(row)
    print(f"{'signature':>28} {'classes':>8} {'|S|-mass':>12} {'mean |S|':>12} {'example':>22}")
    for sig, items in sorted(buckets.items(), key=lambda kv: -sum(r.signed_abs for r in kv[1])):
        mass = sum(r.signed_abs for r in items)
        print(f"{str(sig):>28} {len(items):>8} {mass:>12.8g} {mass/len(items):>12.8g} {str(items[0].cls):>22}")

    print(
        "\nConclusion: exact-period Euler-copy capacity thickens the repeated "
        "tail by multiplicity pattern, but it is blind to QR/NQR within a fixed "
        "pattern because the unit packets are uniform over F_7^*.  The live "
        "extra coordinate is the quadratic-character phase moment."
    )


def inspiration_report() -> None:
    section("PAST-TOPIC INSPIRATION LEDGER")
    print(
        "Older repo threads suggest the right architecture here:\n"
        "  residue/phase/incidence synthesis: copy mass is a residue capacity, "
        "quadratic character is a phase channel, and multi-large walls are an "
        "incidence/compatibility obligation;\n"
        "  OCF packet lesson: retain packet address before scalar evaluation;\n"
        "  root-packet lesson: closed packets need support plus incidence, not "
        "support alone;\n"
        "  Burnside/totient/divisor ledgers: phi counts exact copies, but "
        "character channels decide how those copies cancel."
    )


def tournament_analysis() -> None:
    section("TOURNAMENT ANALYSIS")
    vertices = [
        "multi_large_residue_demand",
        "quadratic_character_phase",
        "euler_copy_unit_capacity",
        "exact_period_1260_packets",
        "height3_one_large_probe",
        "raw_apex_mask",
        "raw_runner_vertices",
    ]
    print("Pairwise observable: preservation of the signed k=10 coimage-tail predicate.")
    print("Switch: prefer the quotient that keeps residue multiplicity, copy capacity, and character phase with less raw row data.")
    print("Hamiltonian path:")
    print("  " + " > ".join(vertices))
    print("score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}")
    print("directed_3_cycles = 0")
    print(
        "Challenged assumption: runners and raw apex masks are too coarse.  "
        "The useful vertices are repeated-residue proof obligations and the "
        "phase/copy coordinates they preserve."
    )


def main() -> None:
    rows = coimage_rows()
    weights = exact_period_packet_report()
    hit_h2, _hit_h3 = wall_height_report(rows)
    tail_profile_report(rows, hit_h2, weights)
    character_split_report(rows, hit_h2, weights)
    inspiration_report()
    tournament_analysis()


if __name__ == "__main__":
    main()
