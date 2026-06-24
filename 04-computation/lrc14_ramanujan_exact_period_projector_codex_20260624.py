#!/usr/bin/env python3
"""HYP-2979: Ramanujan exact-period projector packets for LRC14.

This script is the retained-packet companion to HYP-2978's quotient guardrail.
HYP-2978 shows that scalar divisor/Ramanujan profiles collide across LRC14
proof routes.  Here we test a more active exact-period packet:

  * primitive phase witness profiles on q-th unit residues;
  * q=14 endpoint-sum traces via c_14(r+s);
  * Ramanujan-shell energies of the danger count over Z/qZ.

Arithmetic seed:
  c_q(n) = sum_{(a,q)=1} exp(2*pi*i*a*n/q)
         = sum_{d | gcd(q,n)} d*mu(q/d).

Projector seed:
  For f on Z/qZ, the primitive-frequency shell energy is
      E_q(f) = sum_{a,b mod q} f(a) f(b) c_q(a-b).
  This is an integer exact-period packet.  It keeps more than a scalar qdiv
  count, but less than endpoint-owner geometry.

Tournament Analysis:
  Vertices are proof carriers / exact-period modes rather than runners.
  The pairwise observable is retained LRC predicate payload: primitive weak
  witness status, strict-open status, boundary equality labels, endpoint
  ownership, and compatibility with harmonic/packet duals.  Ties follow the
  declared Hamiltonian path PROOF_TIE_PATH.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from math import gcd, isqrt
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]
N = 14
AP = tuple(range(1, 14))
THRESHOLD = Fraction(1, 14)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s146 = load_module(
    "h2979_haar_baire_boundary",
    REPO / "04-computation" / "lrc14_haar_baire_taut_boundary_s146.py",
)


def section(title: str) -> None:
    print("\n" + "=" * 96)
    print(title)
    print("=" * 96)


def fmt_frac(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def divisors(n: int) -> list[int]:
    out: list[int] = []
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
    return sorted(out)


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


def mobius(n: int) -> int:
    fs = factor(n)
    if any(e > 1 for e in fs.values()):
        return 0
    return -1 if len(fs) % 2 else 1


def phi(n: int) -> int:
    out = n
    for p in factor(n):
        out = out // p * (p - 1)
    return out


def ramanujan_sum(q: int, n: int) -> int:
    return sum(d * mobius(q // d) for d in divisors(gcd(q, n)))


def circular_dist_num(num: int, den: int) -> int:
    r = num % den
    return min(r, den - r)


def slack_num(v: int, a: int, q: int) -> int:
    """Integer threshold slack: nonnegative iff ||v*a/q|| >= 1/14."""
    return N * circular_dist_num(v * a, q) - q


def q_threshold(speeds: tuple[int, ...], cap: int = 240) -> int:
    for q in range(2, cap + 1):
        if all(v % q for v in speeds):
            return q
    return cap + 1


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def replace_one(drop: int, add: int) -> tuple[int, ...]:
    return tuple(sorted((set(AP) - {drop}) | {add}))


def replace_many(drops: tuple[int, ...], adds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((set(AP) - set(drops)) | set(adds)))


@dataclass(frozen=True)
class Row:
    name: str
    route: str
    speeds: tuple[int, ...]


def named_rows() -> list[Row]:
    return [
        Row("AP", "BOUNDARY-AP-GW", AP),
        Row("GW 12->24", "BOUNDARY-AP-GW", replace_one(12, 24)),
        Row("residue liar 12->26", "Q-WITNESS", replace_one(12, 26)),
        Row("near 12->36", "K33-STATE-LIFT", replace_one(12, 36)),
        Row("petal 10->20", "BOUNDARY-PETAL", replace_one(10, 20)),
        Row("petal 13->26", "BOUNDARY-PETAL", replace_one(13, 26)),
        Row("P10+GW", "BOUNDARY-PETAL", replace_many((10, 12), (20, 24))),
        Row("P10+K33", "K33-STATE-LIFT", replace_many((10, 12), (20, 36))),
        Row("covering 12->84", "COVERING-MOMENT", replace_one(12, 84)),
        Row("covering 12->168", "COVERING-MOMENT", replace_one(12, 168)),
        Row("covering 6->98", "COVERING-MOMENT", replace_one(6, 98)),
    ]


def ap_neighborhood_rows(single_limit: int = 180, two_limit: int = 36) -> list[Row]:
    rows: dict[tuple[int, ...], Row] = {row.speeds: row for row in named_rows()}
    for drop in AP:
        for add in range(14, single_limit + 1):
            speeds = replace_one(drop, add)
            if len(speeds) == 13 and primitive(speeds):
                route = "BOUNDARY-AP-GW" if speeds == replace_one(12, 24) else "ONE-SWAP"
                rows.setdefault(speeds, Row(f"single {drop}->{add}", route, speeds))
    for d1, d2 in combinations(AP, 2):
        for a1, a2 in combinations(range(14, two_limit + 1), 2):
            speeds = replace_many((d1, d2), (a1, a2))
            if len(speeds) == 13 and primitive(speeds):
                rows.setdefault(speeds, Row(f"two ({d1},{d2})->({a1},{a2})", "TWO-SWAP", speeds))
    return list(rows.values())


def primitive_phase_packet(speeds: tuple[int, ...], q: int) -> Counter[tuple[int, int]]:
    """Counts primitive phases by (danger_count, boundary_count)."""
    states: Counter[tuple[int, int]] = Counter()
    for a in range(1, q):
        if gcd(a, q) != 1:
            continue
        slacks = [slack_num(v, a, q) for v in speeds]
        danger = sum(1 for s in slacks if s < 0)
        boundary = sum(1 for s in slacks if s == 0)
        states[(danger, boundary)] += 1
    return states


def phase_summary(speeds: tuple[int, ...], q: int) -> tuple[int, int, int, tuple[tuple[tuple[int, int], int], ...]]:
    packet = primitive_phase_packet(speeds, q)
    weak = sum(c for (danger, _), c in packet.items() if danger == 0)
    strict = packet.get((0, 0), 0)
    boundary_only = weak - strict
    return weak, strict, boundary_only, tuple(sorted(packet.items()))


def first_witness_q(speeds: tuple[int, ...], qmax: int = 42, strict: bool = False) -> int | None:
    for q in range(2, qmax + 1):
        weak, strict_count, _, _ = phase_summary(speeds, q)
        if strict:
            if strict_count:
                return q
        elif weak:
            return q
    return None


def danger_counts_on_grid(speeds: tuple[int, ...], q: int) -> list[int]:
    counts: list[int] = []
    for a in range(q):
        counts.append(sum(1 for v in speeds if slack_num(v, a, q) < 0))
    return counts


def weak_safe_indicator_on_grid(speeds: tuple[int, ...], q: int) -> list[int]:
    vals: list[int] = []
    for a in range(q):
        vals.append(1 if all(slack_num(v, a, q) >= 0 for v in speeds) else 0)
    return vals


def ramanujan_energy(values: list[int], q: int) -> int:
    total = 0
    for a, fa in enumerate(values):
        if not fa:
            continue
        for b, fb in enumerate(values):
            if fb:
                total += fa * fb * ramanujan_sum(q, a - b)
    return total


def endpoint_sum_trace(speeds: tuple[int, ...], q: int = 14) -> Counter[int]:
    vals: Counter[int] = Counter()
    for r, s in combinations(speeds, 2):
        vals[ramanujan_sum(q, r + s)] += 1
    return vals


def endpoint_diff_trace(speeds: tuple[int, ...], q: int = 14) -> Counter[int]:
    vals: Counter[int] = Counter()
    for r, s in combinations(speeds, 2):
        vals[ramanujan_sum(q, r - s)] += 1
    return vals


def named_audit() -> None:
    section("NAMED LRC14 ROWS: PRIMITIVE PHASE AND RAMANUJAN PROJECTOR PACKETS")
    print(
        "row                  route              qdiv firstW firstS safe_mu      "
        "q14_states                    E14(N) E27(N) E41(N) E14(weak)"
    )
    for row in named_rows():
        comps = s146.safe_open_components(row.speeds)
        safe_mu = s146.interval_measure(comps)
        fw = first_witness_q(row.speeds, strict=False)
        fs = first_witness_q(row.speeds, strict=True)
        _, _, _, states14 = phase_summary(row.speeds, 14)
        e14 = ramanujan_energy(danger_counts_on_grid(row.speeds, 14), 14)
        e27 = ramanujan_energy(danger_counts_on_grid(row.speeds, 27), 27)
        e41 = ramanujan_energy(danger_counts_on_grid(row.speeds, 41), 41)
        s14 = ramanujan_energy(weak_safe_indicator_on_grid(row.speeds, 14), 14)
        print(
            f"{row.name:<20} {row.route:<18} {q_threshold(row.speeds):>4} "
            f"{str(fw):>6} {str(fs):>6} {fmt_frac(safe_mu):>11} "
            f"{str(states14):<29} {e14:>6} {e27:>6} {e41:>6} {s14:>9}"
        )

    print("\nq=14 endpoint c_14(r+s) traces on unordered pairs:")
    for row in named_rows():
        print(f"  {row.name:<20} {dict(sorted(endpoint_sum_trace(row.speeds).items()))}")

    print("\nq=14 difference c_14(r-r') traces on unordered pairs:")
    for row in named_rows():
        print(f"  {row.name:<20} {dict(sorted(endpoint_diff_trace(row.speeds).items()))}")

    print("\nReadout:")
    print("  q=14 primitive phases are weak boundary witnesses for AP/GW and many q14-front rows.")
    print("  The (danger,boundary) phase packet is sharper than a c_14 speed multiset,")
    print("  while Ramanujan energies keep exact-period variation of the whole danger count.")


def bank_audit() -> None:
    section("AP-NEIGHBORHOOD STRESS: PRIMITIVE EXACT-PERIOD WITNESS PROFILE")
    rows = ap_neighborhood_rows()
    route_counts = Counter(row.route for row in rows)
    first_weak = Counter()
    first_strict = Counter()
    no_weak: list[Row] = []
    no_strict: list[Row] = []
    qdiv_hist = Counter()
    q14_state_collisions: defaultdict[tuple, Counter[str]] = defaultdict(Counter)

    for row in rows:
        qd = q_threshold(row.speeds, 80)
        qdiv_hist[qd] += 1
        fw = first_witness_q(row.speeds, strict=False)
        fs = first_witness_q(row.speeds, strict=True)
        first_weak[fw] += 1
        first_strict[fs] += 1
        if fw is None:
            no_weak.append(row)
        if fs is None:
            no_strict.append(row)
        q14_state_collisions[phase_summary(row.speeds, 14)[3]][row.route] += 1

    print(f"rows audited={len(rows)}")
    print(f"route_counts={dict(sorted(route_counts.items()))}")
    print(f"qdiv_hist_first_80={dict(sorted(qdiv_hist.items(), key=lambda kv: str(kv[0])))}")
    print(f"first_weak_q_hist={dict(sorted(first_weak.items(), key=lambda kv: str(kv[0])))}")
    print(f"first_strict_q_hist={dict(sorted(first_strict.items(), key=lambda kv: str(kv[0])))}")
    print(f"no_weak_q<=42={len(no_weak)}")
    print(f"no_strict_q<=42={len(no_strict)}")
    if no_weak[:10]:
        print("first no-weak examples:", [r.name for r in no_weak[:10]])
    if no_strict[:10]:
        print("first no-strict examples:", [r.name for r in no_strict[:10]])

    mixed = []
    for sig, routes in q14_state_collisions.items():
        if len(routes) > 1:
            mixed.append((sum(routes.values()), sig, dict(routes)))
    mixed.sort(reverse=True, key=lambda x: x[0])
    print("\nq=14 primitive phase packet route-mixing signatures:", len(mixed))
    for size, sig, routes in mixed[:8]:
        print(f"  size={size:<5} routes={routes} sig={sig}")

    print("\nStress readout:")
    print("  Weak primitive witnesses are abundant and reproduce the twist-ladder signal.")
    print("  Strict primitive witnesses are rarer; boundary-only primitive packets must hand")
    print("  off to endpoint owners, C27/K33 labels, or Toeplitz/Fejer duals.")


PROOF_VERTICES = (
    "labelled_packet_sheaf",
    "toeplitz_fejer_dual",
    "spectral_shadow",
    "ramanujan_danger_energy",
    "primitive_phase_packet",
    "c14_endpoint_sum_trace",
    "carmichael_autocorrelation",
    "raw_qdiv",
    "raw_runner_residues",
)

PROOF_SCORES: dict[str, tuple[int, ...]] = {
    # weak witness, strict-open, boundary equality, harmonic compatibility, anti-scalar guard
    "labelled_packet_sheaf": (9, 9, 9, 9, 9),
    "toeplitz_fejer_dual": (7, 9, 4, 9, 8),
    "spectral_shadow": (7, 9, 3, 8, 8),
    "ramanujan_danger_energy": (6, 6, 5, 8, 7),
    "primitive_phase_packet": (8, 5, 7, 6, 7),
    "c14_endpoint_sum_trace": (5, 3, 9, 5, 6),
    "carmichael_autocorrelation": (5, 5, 5, 7, 6),
    "raw_qdiv": (6, 2, 2, 2, 3),
    "raw_runner_residues": (3, 2, 3, 1, 2),
}

PROOF_TIE_PATH = (
    "labelled_packet_sheaf",
    "toeplitz_fejer_dual",
    "spectral_shadow",
    "ramanujan_danger_energy",
    "primitive_phase_packet",
    "carmichael_autocorrelation",
    "c14_endpoint_sum_trace",
    "raw_qdiv",
    "raw_runner_residues",
)


def compare(a: str, b: str) -> int:
    sa = sum(x > y for x, y in zip(PROOF_SCORES[a], PROOF_SCORES[b]))
    sb = sum(x < y for x, y in zip(PROOF_SCORES[a], PROOF_SCORES[b]))
    if sa != sb:
        return 1 if sa > sb else -1
    return 1 if PROOF_TIE_PATH.index(a) < PROOF_TIE_PATH.index(b) else -1


def tournament_analysis() -> None:
    section("TOURNAMENT ANALYSIS ON EXACT-PERIOD PROOF CARRIERS")
    vertices = list(PROOF_VERTICES)
    edges = {(a, b): compare(a, b) for a in vertices for b in vertices if a != b}
    scores = Counter()
    for a in vertices:
        scores[sum(1 for b in vertices if a != b and edges[(a, b)] > 0)] += 1

    c3 = 0
    for a, b, c in combinations(vertices, 3):
        if edges[(a, b)] > 0 and edges[(b, c)] > 0 and edges[(c, a)] > 0:
            c3 += 1
        if edges[(a, c)] > 0 and edges[(c, b)] > 0 and edges[(b, a)] > 0:
            c3 += 1

    hp = 0
    for perm in permutations(vertices):
        if all(edges[(perm[i], perm[i + 1])] > 0 for i in range(len(perm) - 1)):
            hp += 1

    ordered = sorted(vertices, key=lambda v: -sum(1 for b in vertices if v != b and edges[(v, b)] > 0))
    print("vertices:")
    for v in ordered:
        print(f"  {v}: scores={PROOF_SCORES[v]}")
    print(f"score_hist={dict(sorted(scores.items()))}")
    print(f"directed_3_cycles={c3}")
    print(f"hamiltonian_paths={hp}")
    print("Hamiltonian path:")
    print("  " + " > ".join(ordered))
    print("\nAssumption challenge:")
    print("  considered vertices: runners, residues, primitive phases, q-denominators,")
    print("    Ramanujan modes, endpoint owner pairs, autocorrelation shifts, and proof obligations.")
    print("  chosen vertices: proof carriers / exact-period modes.")
    print("  preserved predicate: weak or strict LRC witness status, boundary equality,")
    print("    and compatibility with harmonic or labelled-packet discharge.")
    print("  destroyed data: raw runner names and endpoint ownership unless explicitly reattached.")


def theorem_target() -> None:
    section("HYP-2979 THEOREM TARGET")
    print("Primitive exact-period packet theorem target:")
    print()
    print("  After current Moon-core reductions, every LRC14 residual either")
    print("    (1) has a primitive weak or strict q-phase packet certifying the row;")
    print("    (2) has a negative Toeplitz/Fejer or positive spectral-shadow dual;")
    print("    (3) has an AP/GW c_14 endpoint-sum boundary packet;")
    print("    (4) has a Ramanujan danger-energy defect that forces a labelled handoff; or")
    print("    (5) carries the K33/HYP-2908/THM-572 state-lift debt.")
    print()
    print("What changed in this pass:")
    print("  Ramanujan sums should not be used only as speed-profile scalars.")
    print("  The useful object is the exact-period kernel on functions of phase:")
    print("      E_q(f)=sum_{a,b} f(a) f(b)c_q(a-b).")
    print("  That kernel can sit between the q-ladder and the full Toeplitz dual,")
    print("  while HYP-2978 tells us when the quotient has forgotten too much.")


def main() -> None:
    section("RAMANUJAN IDENTITY CHECK")
    failures = []
    for q in range(2, 50):
        for n in range(0, 120):
            lhs = ramanujan_sum(q, n)
            g = gcd(q, n)
            rhs = mobius(q // g) * phi(q) // phi(q // g)
            if lhs != rhs:
                failures.append((q, n, lhs, rhs))
                break
    print(f"checked q<50, n<120; failures={failures[:3]}")
    print("c_14 gcd-class table:")
    for g in (1, 2, 7, 14):
        print(f"  gcd(14,n)={g:<2} c_14={ramanujan_sum(14, g)}")
    named_audit()
    bank_audit()
    tournament_analysis()
    theorem_target()


if __name__ == "__main__":
    main()
