#!/usr/bin/env python3
"""Pisano quotient diagnostics for the LRC n=14 band ladder.

This extends THM-491/492 and HYP-2438 on the exact one-stranger family
    S(r) = 7*{1,...,12} union {r}.

Main targets:
  * isolate the shell-27 obstruction as a 9-class antipodal Pisano quotient;
  * measure whether the known S(r) residual is exactly the missing quotient
    class plus the 13-clock block;
  * compare band-1 shells, the fibered Q27/Q41 lattices, and B' gates;
  * record Tournament Analysis over proof gates rather than over runners.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from functools import reduce
from itertools import combinations, permutations
from math import gcd

N = 14
TARGET = F(1, N)
CORE = [7 * k for k in range(1, 13)]
MAX_R = 13 * 84


def pisano(m: int) -> int:
    if m == 1:
        return 1
    a, b = 0, 1
    for i in range(6 * m * m + 1):
        a, b = b, (a + b) % m
        if a == 0 and b == 1:
            return i + 1
    raise RuntimeError(f"pisano search exhausted for m={m}")


def ord_mod(a: int, m: int) -> int:
    if gcd(a, m) != 1:
        raise ValueError("multiplicative order requires a unit")
    x = 1
    for k in range(1, 10 * m + 1):
        x = (x * a) % m
        if x == 1:
            return k
    raise RuntimeError(f"order search exhausted for a={a}, m={m}")


def antipodal_class(r: int, m: int) -> int | None:
    r %= m
    if r == 0 or gcd(r, m) != 1:
        return None
    return min(r, (-r) % m)


def unit_classes_mod_pm(m: int) -> list[int]:
    return sorted(
        {antipodal_class(r, m) for r in range(1, m) if gcd(r, m) == 1}
    )


def orbit_by_two_mod_pm(m: int, start: int = 1) -> list[int]:
    out: list[int] = []
    r = start % m
    while True:
        c = antipodal_class(r, m)
        assert c is not None
        if c in out:
            return out
        out.append(c)
        r = (2 * r) % m


def class_coverage(S: list[int], m: int = 27) -> list[int]:
    return sorted({c for v in S if (c := antipodal_class(v, m)) is not None})


def inv_class(c: int, m: int = 27) -> int:
    for u in range(1, m):
        if gcd(u, m) == 1 and (c * u) % m in (1, m - 1):
            out = antipodal_class(u, m)
            assert out is not None
            return out
    raise RuntimeError(f"no inverse antipodal class for {c} mod {m}")


def witness_at_q(S: list[int], q: int) -> int | None:
    """Strict band witness at q: all va avoid +-floor(q/14)."""
    band = q // N
    prof = [v % q for v in S]
    for a in range(1, q):
        if gcd(a, q) != 1:
            continue
        if all(min((p * a) % q, q - (p * a) % q) > band for p in prof):
            return a
    return None


def min_witness_modulus(S: list[int], qmax: int = 600) -> int | None:
    for q in range(2, qmax + 1):
        if witness_at_q(S, q) is not None:
            return q
    return None


def first_witness_in(S: list[int], Q: list[int]) -> tuple[int, int] | None:
    for q in Q:
        a = witness_at_q(S, q)
        if a is not None:
            return q, a
    return None


def q_lattice(max_m: int) -> list[int]:
    return sorted({d * m for d in (1, 2, 7, 14) for m in range(1, max_m + 1)} - {1})


def safe_components(U: list[int]) -> list[tuple[F, F]]:
    arcs: list[tuple[F, F]] = []
    for u in U:
        for k in range(0, u + 1):
            a = F(N * k - 1, N * u)
            b = F(N * k + 1, N * u)
            arcs.append((max(a, F(0)), min(b, F(1))))
    arcs.sort()
    merged: list[tuple[F, F]] = []
    for a, b in arcs:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    return [
        (merged[i][1], merged[i + 1][0])
        for i in range(len(merged) - 1)
        if merged[i + 1][0] > merged[i][1]
    ]


def bprime_margin_for_runner(S: list[int], v: int) -> F:
    comps = safe_components([u for u in S if u != v])
    if not comps:
        return F(-1)
    threshold = F(2, N * v)
    return max((b - a) - threshold for a, b in comps)


def bprime_mult14(S: list[int]) -> bool:
    for v in [u for u in S if u % N == 0]:
        if bprime_margin_for_runner(S, v) > 0:
            return True
    return False


def bprime_any(S: list[int], preferred: int | None = None) -> tuple[bool, int | None, F]:
    """Return a B' certificate, trying a preferred runner first when supplied."""
    checked: set[int] = set()
    if preferred is not None:
        margin = bprime_margin_for_runner(S, preferred)
        checked.add(preferred)
        if margin > 0:
            return True, preferred, margin
    best_v: int | None = None
    best_margin = F(-1)
    for v in S:
        if v in checked:
            continue
        margin = bprime_margin_for_runner(S, v)
        if margin > 0:
            return True, v, margin
        if margin > best_margin:
            best_v, best_margin = v, margin
    return False, best_v, best_margin


@dataclass(frozen=True)
class FamilyRow:
    r: int
    min_q: int
    min_a: int
    q27: tuple[int, int] | None
    q41: tuple[int, int] | None
    bmult14: bool
    bany: tuple[bool, int | None, F]

    @property
    def S(self) -> list[int]:
        return sorted(CORE + [self.r])


def valid_family_rows() -> list[FamilyRow]:
    Q27 = q_lattice(27)
    Q41 = q_lattice(41)
    rows: list[FamilyRow] = []
    for r in range(1, MAX_R + 1):
        if r % 7 == 0:
            continue
        S = sorted(CORE + [r])
        if len(set(S)) != 13 or reduce(gcd, S) != 1:
            continue
        mq = min_witness_modulus(S)
        assert mq is not None, f"no witness q<=600 for r={r}"
        ma = witness_at_q(S, mq)
        assert ma is not None
        bmult = bprime_mult14(S)
        rows.append(
            FamilyRow(
                r=r,
                min_q=mq,
                min_a=ma,
                q27=first_witness_in(S, Q27),
                q41=first_witness_in(S, Q41),
                bmult14=bmult,
                bany=bprime_any(S, preferred=r),
            )
        )
    return rows


def fraction_short(x: F) -> str:
    return f"{x.numerator}/{x.denominator}"


def gate_tournament(rows: list[FamilyRow]) -> list[str]:
    old_evader_set = {
        row.r for row in rows if row.min_q > 27 and not row.bmult14
    }
    shell_blockers = {row.r for row in rows if row.min_q > 27}
    gates = [
        ("Q41_fiber_lattice", lambda row: row.q41 is not None, 1),
        ("band2_min_q_le_41", lambda row: row.min_q <= 41, 2),
        ("Q27_fiber_lattice", lambda row: row.q27 is not None, 3),
        ("Bprime_any_runner", lambda row: row.bany[0], 4),
        ("band1_plain_q_le_27", lambda row: row.min_q <= 27, 5),
        ("Bprime_mult14_only", lambda row: row.bmult14, 6),
        ("not_missing_shell27_class", lambda row: row.r % 27 not in {0, 10, 17}, 7),
        ("not_13_clock_block", lambda row: row.r % 13 != 0, 8),
    ]

    metrics: dict[str, tuple[int, int, int, int]] = {}
    covered: dict[str, set[int]] = {}
    for name, pred, tie_rank in gates:
        s = {row.r for row in rows if pred(row)}
        covered[name] = s
        metrics[name] = (
            len(s),
            len(s & old_evader_set),
            len(s & shell_blockers),
            -tie_rank,
        )

    names = [name for name, _, _ in gates]
    arcs: dict[tuple[str, str], str] = {}
    flips = 0
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            winner = a if metrics[a] >= metrics[b] else b
            arcs[(a, b)] = winner
            if winner != a:
                flips += 1

    outdeg = {name: 0 for name in names}
    for (a, b), winner in arcs.items():
        outdeg[winner] += 1
    score_hist = Counter(outdeg.values())

    cycles = 0
    for a_i in range(len(names)):
        for b_i in range(a_i + 1, len(names)):
            for c_i in range(b_i + 1, len(names)):
                tri = [names[a_i], names[b_i], names[c_i]]
                wins = {x: 0 for x in tri}
                for x, y in ((tri[0], tri[1]), (tri[0], tri[2]), (tri[1], tri[2])):
                    winner = arcs[(x, y)] if (x, y) in arcs else arcs[(y, x)]
                    wins[winner] += 1
                if sorted(wins.values()) == [1, 1, 1]:
                    cycles += 1

    def beats(a: str, b: str) -> bool:
        if (a, b) in arcs:
            return arcs[(a, b)] == a
        return arcs[(b, a)] == a

    hp_count = 0
    for path in permutations(names):
        if all(beats(path[i], path[i + 1]) for i in range(len(path) - 1)):
            hp_count += 1

    lines = []
    lines.append("== Tournament Analysis: proof gates, not runners ==")
    lines.append(
        "  vertices: Q41, band2, Q27, B'(any), band1, B'(mult14), shell27-residue, 13-clock"
    )
    lines.append(
        "  observable tuple per gate: (family coverage, old-evader coverage, shell-blocker coverage, -tie-rank)"
    )
    lines.append(
        "  switch: A -> B iff A has lexicographically larger tuple; tie Hamiltonian path is declaration order"
    )
    for name in names:
        lines.append(f"    {name}: metric={metrics[name]}, outdeg={outdeg[name]}")
    lines.append(f"  score histogram: {dict(sorted(score_hist.items()))}")
    lines.append(f"  directed 3-cycles: {cycles}; Hamiltonian paths: {hp_count}; edge flips vs declaration path: {flips}/28")
    lines.append(
        "  assumption challenge: alternate vertices considered = runners, gaps, residues, shell unit classes,"
    )
    lines.append(
        "    Pisano quotient classes, band-rung events, B' targets, and proof obligations."
    )
    lines.append(
        "    This quotient preserves exact S(r) coverage and the residue mechanism; it destroys arbitrary"
    )
    lines.append(
        "    multi-stranger geometry, endpoint-owner data, and nonlocal pressure dependencies."
    )
    return lines


def main() -> None:
    print("== Pisano / shell-27 quotient ==")
    pis = {m: pisano(m) for m in (3, 9, 27, 81)}
    print(
        "  Pisano tower pi(3),pi(9),pi(27),pi(81) = "
        + ",".join(str(pis[m]) for m in (3, 9, 27, 81))
    )
    print(
        f"  ratios: pi(9)/pi(3)={pis[9] // pis[3]}, "
        f"pi(27)/pi(9)={pis[27] // pis[9]}, pi(27)/pi(3)={pis[27] // pis[3]}"
    )
    print(f"  ord_27(2)={ord_mod(2, 27)}; quotient by +- has orbit length {len(orbit_by_two_mod_pm(27))}")
    classes = unit_classes_mod_pm(27)
    orbit = orbit_by_two_mod_pm(27)
    core_classes = class_coverage(CORE, 27)
    missing = sorted(set(classes) - set(core_classes))
    print(f"  (Z/27)^*/+- classes ({len(classes)}): {classes}")
    print(f"  multiply-by-2 orbit on classes: {orbit}")
    print(f"  core 7*{{1..12}} covers {len(core_classes)}/{len(classes)} unit classes: {core_classes}")
    print(f"  missing shell-27 class: {missing}; inverse multiplier class: {[inv_class(c, 27) for c in missing]}")
    core_q27_a = witness_at_q(CORE, 27)
    print(f"  q=27 witness for the 12-runner core alone: a={core_q27_a}")
    print()

    rows = valid_family_rows()
    assert len(rows) == 936
    min_hist = Counter(row.min_q for row in rows)
    shell_blockers = [row for row in rows if row.min_q > 27]
    old_evaders = [row for row in rows if row.min_q > 27 and not row.bmult14]
    assert [row.r for row in old_evaders] == [611, 702, 793, 962, 1053]

    print("== Exact S(r)=7*{1..12} union {r} family ==")
    print(f"  valid primitive rows r<=1092: {len(rows)}")
    print(f"  all loose by exact band witness q<=600: {sum(row.min_q is not None for row in rows)}/{len(rows)}")
    print(f"  max minimal witness modulus: {max(row.min_q for row in rows)}")
    print(f"  minimal-q histogram: {dict(sorted(min_hist.items()))}")
    print(f"  rows blocking all plain q<=27 shells: {len(shell_blockers)}")
    print(f"  old evaders of [plain q<=27 shells union B'(mult14)]: {len(old_evaders)}")
    print()

    residue_counter = Counter((row.r % 13, row.r % 27) for row in shell_blockers)
    print("== Residue signature of the shell-27 residual ==")
    print(f"  shell-blocker residues (r mod 13, r mod 27): {dict(sorted(residue_counter.items()))}")
    print("  old evader details:")
    for row in old_evaders:
        q27 = row.q27 if row.q27 is not None else ("none", "none")
        q41 = row.q41 if row.q41 is not None else ("none", "none")
        bany_ok, bany_v, bany_margin = row.bany
        cls = antipodal_class(row.r, 27)
        print(
            f"    r={row.r}: min=(q={row.min_q}, a={row.min_a}); "
            f"r mod 13={row.r % 13}, r mod 27={row.r % 27}, class27={cls}; "
            f"Q27={q27}; Q41={q41}; B'any={bany_ok} at v={bany_v} "
            f"margin={fraction_short(bany_margin)}"
        )
    print(
        "  interpretation: the 7-core misses exactly class +-10 mod 27; the five old evaders"
    )
    print(
        "    have r mod 27 in {0,+-10} and r mod 13=0, so they either plug the missing"
    )
    print(
        "    shell-27 unit class or collapse to zero while also blocking the generic 13-clock."
    )
    print()

    q27_hits = sum(row.q27 is not None for row in rows)
    q41_hits = sum(row.q41 is not None for row in rows)
    bany_hits = sum(row.bany[0] for row in rows)
    bany_stranger_hits = sum(row.bany[0] and row.bany[1] == row.r for row in rows)
    q27_misses = [row.r for row in rows if row.q27 is None]
    q41_misses = [row.r for row in rows if row.q41 is None]
    print("== Fibered lattice and B' gate comparison ==")
    print(f"  Q27={{d*m: d|14, m<=27}} size={len(q_lattice(27))}, family coverage={q27_hits}/{len(rows)}")
    print(f"  Q41={{d*m: d|14, m<=41}} size={len(q_lattice(41))}, family coverage={q41_hits}/{len(rows)}")
    print(f"  Q27 misses: {q27_misses}")
    print(f"  Q41 misses: {q41_misses}")
    print(f"  B'(any runner) coverage over family: {bany_hits}/{len(rows)}")
    print(
        f"  B'(any runner) certificates: stranger target {bany_stranger_hits}, "
        f"core-runner target {bany_hits - bany_stranger_hits}"
    )
    print(
        "  exact closure on this family: "
        f"Q27 covers all rows = {q27_hits == len(rows)}; "
        f"Q41 covers all rows = {q41_hits == len(rows)}; "
        f"B'(any) covers all rows = {bany_hits == len(rows)}"
    )
    print()

    print("== Tiny two-stranger stress probe ==")
    hard_rs = [row.r for row in rows if row.min_q > 27]
    core11 = [7 * k for k in range(1, 12)]
    pair_records = []
    for r, s in combinations(hard_rs, 2):
        S = sorted(core11 + [r, s])
        if len(set(S)) != 13 or reduce(gcd, S) != 1:
            continue
        qhit = first_witness_in(S, q_lattice(27))
        bany = bprime_any(S)
        pair_records.append((r, s, qhit, bany[0]))
    qhit_hist = Counter(qhit[0] if qhit is not None else None for _, _, qhit, _ in pair_records)
    pair_misses = [(r, s) for r, s, qhit, bany_ok in pair_records if qhit is None and not bany_ok]
    print(f"  paired hard single-stranger residues over 7*{{1..11}}: {hard_rs}")
    print(f"  tested primitive pairs: {len(pair_records)}")
    print(f"  first Q27 witness-q histogram: {dict(sorted(qhit_hist.items(), key=lambda kv: str(kv[0])))}")
    print(f"  Q27 plus B'(any) misses: {pair_misses}")
    print(
        "  lesson: simply pairing the one-stranger blockers is too easy; deleting 84 opens"
    )
    print(
        "    the q=12 divisor clock. A hard multi-stranger search must keep low-clock coverage"
    )
    print("    while spending shell-27/Pisano resources.")
    print()

    for line in gate_tournament(rows):
        print(line)

    print()
    print("== Proof-shaped takeaway ==")
    print("  HYP-2438's Q lattice is exact on the whole one-stranger family already at Q27;")
    print("  the two rows with minimal plain witness q=41 are nevertheless caught at q=91.")
    print("  Q41 is therefore redundant for S(r), though it remains the visible band-2 stress test.")
    print("  The Pisano quotient explains the residual rather than merely listing it: shell 27")
    print("  leaves one antipodal unit class, +-10, and the old evaders occupy that class (or 0)")
    print("  while also killing the 13-clock.  Multi-stranger configurations remain the live gap.")


if __name__ == "__main__":
    main()
