#!/usr/bin/env python3
"""Exact THM-4006 physical Rule 30 cross-scale/off-ray observer census.

The finite theorem universe is the canonical highest-shell realization
``n = 2**m + t`` with ``0 <= t < 2**m``.  Every record keeps the exact gap
branch, a finite phase-owner portrait, a larger block of odd normalized
units, the ordinary center-carry state, a finite projective scale history,
and routed off-ray data.  Equal proposed states are grouped before comparing
the physical center consumer and the next selected bank.

Nothing here proves an all-scale quotient, a bounded gap, or a Rule 30 prize.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys


sys.stdout.reconfigure(newline="\n")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEPENDENCY = ROOT / "04-computation" / "rule30_orbit_signalizer_gap_thm3511.py"
DEPENDENCY_SHA256 = "2ce110f0b8e9c71c3d298aaf07e8e6c02b70d33e5671bc763f3f3b490caa5445"

MAX_SCALE = 9
OWNER_DEPTH = 4
SELECTED_RAYS = (0, 1, 2)
ODD_BLOCK_LENGTH = 4
ODD_PRECISION = 4
PROJECTIVE_HISTORY = 3
PROJECTIVE_PRECISION = 4


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def load_dependency():
    require(hashlib.sha256(DEPENDENCY.read_bytes()).hexdigest() == DEPENDENCY_SHA256,
            "THM-3511 dependency hash")
    spec = importlib.util.spec_from_file_location("rule30_thm3511", DEPENDENCY)
    require(spec is not None and spec.loader is not None, "dependency loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rows_through(module, stop: int) -> list[int]:
    rows = [1]
    row = 1
    # Exact Python integers: support of R_t is at most bit 2t.
    for _ in range(stop):
        row = row ^ ((row << 1) | (row << 2))
        rows.append(row)
    return rows


def unit(rows: list[int], v: tuple[int, ...], scale: int, phase: int) -> int:
    q = 1 << scale
    numerator = rows[phase + q] - rows[phase]
    require(numerator % (1 << v[scale]) == 0, ("unit integrality", scale, phase))
    answer = numerator >> v[scale]
    require(answer & 1, ("unit odd", scale, phase))
    return answer


def owner_portrait(module, owner: bytearray, depth: int) -> tuple[int, ...]:
    return tuple(module.apply_word(owner, ray, depth) for ray in range(1 << depth))


def selected_bank(module, owner: bytearray, depth: int = OWNER_DEPTH) -> tuple[int, ...]:
    return tuple(module.apply_word(owner, ray, depth) for ray in SELECTED_RAYS)


def center_carry_state(module, rows: list[int], v: tuple[int, ...], n: int):
    alpha = module.nu2(n)
    va = v[alpha]
    center = (rows[n] >> n) & 1
    if va > n:
        require(center == 0, ("pre-visibility center", n, va))
        return ("pre",), center, ()

    level = n - va
    modulus = 1 << (level + 1)
    terms = []
    target_xor = 0
    lower_sum = 0
    for scale in range(n.bit_length()):
        if not ((n >> scale) & 1):
            continue
        phase = n & ((1 << scale) - 1)
        raw_unit = unit(rows, v, scale, phase)
        residue = ((1 << (v[scale] - va)) * raw_unit) % modulus
        target = (residue >> level) & 1
        lower = residue & ((1 << level) - 1)
        target_xor ^= target
        lower_sum += lower
        terms.append((scale, target, lower))
    carry_parity = (lower_sum >> level) & 1 if level else 0
    require(center == (target_xor ^ carry_parity), ("ordinary carry", n))
    # The state keeps the direct target parity and the ordinary carry parity,
    # not merely their XOR.  The term list is retained as a hostile audit but
    # is not part of the proposed finite observer.
    return ("visible", target_xor, carry_parity), center, tuple(terms)


def projective_history(rows: list[int], v: tuple[int, ...], gaps: tuple[int, ...],
                       scale: int, phase: int):
    entries: list[object] = [None] * max(0, PROJECTIVE_HISTORY - scale - 1)
    exact: list[object] = [None] * max(0, PROJECTIVE_HISTORY - scale - 1)
    first = max(0, scale - PROJECTIVE_HISTORY + 1)
    for level in range(first, scale + 1):
        d = gaps[level]
        precision = PROJECTIVE_PRECISION + d
        modulus = 1 << precision
        u0 = unit(rows, v, level, phase)
        u1 = unit(rows, v, level, phase + (1 << level))
        g = (-u1 * pow(u0, -1, modulus)) % modulus
        require(((g - 1) & ((1 << d) - 1)) == 0,
                ("projective gap congruence", level, phase))
        z = ((1 - g) >> d) % (1 << PROJECTIVE_PRECISION)
        entries.append((d, g % (1 << PROJECTIVE_PRECISION), z))
        g_num, g_den = -u1, u0
        common = math.gcd(abs(g_num), g_den)
        g_num //= common
        g_den //= common
        z_num, z_den = u0 + u1, (1 << d) * u0
        common = math.gcd(abs(z_num), z_den)
        z_num //= common
        z_den //= common
        exact.append((d, g_num, g_den, z_num, z_den))
    require(len(entries) == PROJECTIVE_HISTORY, "projective history length")
    require(len(exact) == PROJECTIVE_HISTORY, "exact projective history length")
    return tuple(entries), tuple(exact)


def routed_chains(module, owner: bytearray, gap: int,
                  extension: int) -> tuple[tuple[int, int, int], ...]:
    # Depth is just large enough to carry the target depth after stripping the
    # exact fixed zero prefix.  These are physical, state-dependent queries.
    depth = OWNER_DEPTH + gap
    modulus = 1 << depth
    chains = []
    require(0 <= extension < (1 << gap), "phase-owned extension")
    for ray in SELECTED_RAYS:
        # The marked origin has extension zero.  A general physical phase is
        # based at its actual tail prefix a, as in THM-3511 (29); silently
        # reusing the zero ray is false already at (m,t)=(2,2).
        x = extension | (ray << gap)
        first = module.apply_word(owner, x, depth)
        second = module.apply_word(owner, first, depth)
        require(0 <= first < modulus and 0 <= second < modulus, "routed range")
        chains.append((x, first, second))
    return tuple(chains)


@dataclass(frozen=True)
class Record:
    n: int
    scale: int
    phase: int
    gap: int
    owner_bank: tuple[int, ...]
    owner_portrait: tuple[int, ...]
    odd_block: tuple[int, ...]
    exact_odd_block: tuple[int, ...]
    carry_state: tuple[object, ...]
    projective: tuple[object, ...]
    exact_projective: tuple[object, ...]
    chains: tuple[tuple[int, int, int], ...]
    center: int
    next_bank: tuple[int, ...]
    exact_owner_sha256: str
    exact_owner: bytes
    carry_terms: tuple[tuple[int, int, int], ...]

    def target(self) -> tuple[object, ...]:
        return self.center, self.next_bank


def build_records(module) -> tuple[list[Record], tuple[int, ...], tuple[int, ...]]:
    words, raw_gaps, _, _ = module.signalizer_words()
    v = tuple(module.EXPECTED_V_0_TO_23)
    gaps = tuple(raw_gaps)
    max_q = 1 << MAX_SCALE
    # Larger block and projective ratios need up to phase + 4q.
    rows = rows_through(module, 5 * max_q)
    records = []
    modulus = 1 << ODD_PRECISION

    for scale in range(MAX_SCALE + 1):
        q = 1 << scale
        gap = gaps[scale]
        power_word = bytearray((module.A,)) * q
        doubled_word = bytearray((module.A,)) * (2 * q)
        for phase in range(q):
            n = q + phase
            prefix = rows[phase] & ((1 << v[scale]) - 1)
            owner = module.section_at_prefix(power_word, prefix, v[scale])
            require(module.activity(owner) == 1, ("active phase owner", scale, phase))

            extension = (rows[phase] >> v[scale]) & ((1 << gap) - 1)
            induced_next = module.section_at_prefix(owner + owner, extension, gap)
            long_prefix = rows[phase] & ((1 << (v[scale] + gap)) - 1)
            direct_next = module.section_at_prefix(
                doubled_word, long_prefix, v[scale] + gap,
            )
            require(induced_next == direct_next,
                    ("physical phase transition", scale, phase))

            chains = routed_chains(module, owner, gap, extension)
            recovered = tuple(
                (chain[2] >> gap) & ((1 << OWNER_DEPTH) - 1)
                for chain in chains
            )
            next_bank = selected_bank(module, induced_next)
            require(recovered == next_bank,
                    ("routed next bank", scale, phase))

            carry_state, center, carry_terms = center_carry_state(module, rows, v, n)
            exact_odd_block = tuple(
                unit(rows, v, scale, phase + j * q)
                for j in range(ODD_BLOCK_LENGTH)
            )
            odd_block = tuple(value % modulus for value in exact_odd_block)
            finite_projective, exact_projective = projective_history(
                rows, v, gaps, scale, phase,
            )
            records.append(
                Record(
                    n=n,
                    scale=scale,
                    phase=phase,
                    gap=gap,
                    owner_bank=selected_bank(module, owner),
                    owner_portrait=owner_portrait(module, owner, OWNER_DEPTH),
                    odd_block=odd_block,
                    exact_odd_block=exact_odd_block,
                    carry_state=carry_state,
                    projective=finite_projective,
                    exact_projective=exact_projective,
                    chains=chains,
                    center=center,
                    next_bank=next_bank,
                    exact_owner_sha256=hashlib.sha256(owner).hexdigest(),
                    exact_owner=bytes(owner),
                    carry_terms=carry_terms,
                )
            )

    require(tuple(item.n for item in records) == tuple(range(1, (1 << (MAX_SCALE + 1)))),
            "canonical chronological universe")
    return records, v, gaps


def state(record: Record, level: str) -> tuple[object, ...]:
    pieces: list[object] = [record.gap]
    if level in {"bank", "portrait", "odd", "carry", "projective", "base", "one", "two", "exact"}:
        pieces.append(record.owner_bank)
    if level in {"portrait", "odd", "carry", "projective", "base", "one", "two", "exact"}:
        pieces.append(record.owner_portrait)
    if level in {"odd", "carry", "projective", "base", "one", "two"}:
        pieces.append(record.odd_block)
    if level in {"carry", "projective", "base", "one", "two"}:
        pieces.append(record.carry_state)
    if level in {"projective", "base", "one", "two"}:
        pieces.append(record.projective)
    if level in {"base", "one", "two"}:
        pieces.append(record.chains[0])
    if level in {"one", "two"}:
        pieces.append(record.chains[1])
    if level == "two":
        pieces.append(record.chains[2])
    if level == "exact":
        pieces.extend(
            (
                record.exact_owner,
                record.exact_odd_block,
                (record.carry_state, record.carry_terms),
                record.exact_projective,
            )
        )
    return tuple(pieces)


def quotient(records: list[Record], level: str):
    fibres: dict[tuple[object, ...], list[Record]] = defaultdict(list)
    for record in records:
        fibres[state(record, level)].append(record)
    mismatches = []
    for key, rows in fibres.items():
        targets: dict[tuple[object, ...], Record] = {}
        for row in rows:
            targets.setdefault(row.target(), row)
        if len(targets) > 1:
            representatives = sorted(targets.values(), key=lambda item: item.n)
            mismatches.append((representatives[1].n, representatives[0].n, key,
                               representatives[0], representatives[1]))
    mismatches.sort(key=lambda item: (item[0], item[1], repr(item[2])))
    return fibres, mismatches


def mismatch_scope(fibres):
    same = cross = 0
    first_cross = None
    for key, rows in fibres.items():
        has_same = has_cross = False
        for index, left in enumerate(rows):
            for right in rows[index + 1:]:
                if left.target() == right.target():
                    continue
                if left.scale == right.scale:
                    has_same = True
                else:
                    has_cross = True
                    candidate = (max(left.n, right.n), min(left.n, right.n), key, left, right)
                    if first_cross is None or candidate[:2] < first_cross[:2]:
                        first_cross = candidate
        same += has_same
        cross += has_cross
    return same, cross, first_cross


def compact_record(record: Record) -> tuple[object, ...]:
    return (
        record.n, record.scale, record.phase, record.gap,
        record.owner_bank, record.odd_block, record.carry_state,
        record.projective, record.chains, record.target(),
        record.exact_owner_sha256,
    )


def main() -> None:
    module = load_dependency()
    records, v, gaps = build_records(module)

    # THM-3824 same-scale positive control: exact free defects are strictly
    # increasing on every canonical phase interval used by the scout.
    strict_checks = 0
    same_scale_ranges = []
    rows = rows_through(module, 5 * (1 << MAX_SCALE))
    for scale in range(MAX_SCALE + 1):
        q = 1 << scale
        defects = tuple(
            unit(rows, v, scale, phase + q) - unit(rows, v, scale, phase)
            for phase in range(q)
        )
        require(all(value > 0 for value in defects), ("positive defect", scale))
        require(all(a < b for a, b in zip(defects, defects[1:])),
                ("strict defect", scale))
        strict_checks += len(defects) + max(0, len(defects) - 1)
        same_scale_ranges.append(
            (
                scale,
                len(defects),
                defects[0].bit_length(),
                defects[-1].bit_length(),
                defects[0] % 65536,
                defects[-1] % 65536,
            )
        )

    levels = ("branch", "bank", "portrait", "odd", "carry", "projective",
              "base", "one", "two", "exact")
    census = []
    first_hostiles = []
    hostile_by_level = {}
    cross_by_level = {}
    for level in levels:
        fibres, mismatches = quotient(records, level)
        same_mismatch_fibres, cross_mismatch_fibres, first_cross = mismatch_scope(fibres)
        collision_fibres = sum(len(rows) > 1 for rows in fibres.values())
        census.append(
            (
                level, len(fibres), collision_fibres, len(mismatches),
                same_mismatch_fibres, cross_mismatch_fibres,
            )
        )
        if mismatches:
            _, _, key, left, right = mismatches[0]
            hostile_by_level[level] = (key, left, right)
            first_hostiles.append(
                (level, key, compact_record(left), compact_record(right))
            )
        if first_cross is not None:
            _, _, key, left, right = first_cross
            cross_by_level[level] = (key, left, right)
            first_hostiles.append(
                (level + "_FIRST_CROSS_SCALE", key,
                 compact_record(left), compact_record(right))
            )

    # Two routed off-ray chains must always determine the next bank by the
    # inherited exact section formula.  The center is independently fixed by
    # carry_state.  This is a bounded positive candidate, not a universal
    # finite observer theorem.
    _, final_mismatches = quotient(records, "two")
    require(not final_mismatches, "full bounded candidate target congruence")

    # Cheapest physical warning against transporting the marked-origin zero
    # ray to every phase.  The n=6 record is (m,t)=(2,2).
    phase_hostile = records[5]
    require((phase_hostile.scale, phase_hostile.phase) == (2, 2),
            "phase-owned-base hostile address")
    naive_chains = routed_chains(
        module, bytearray(phase_hostile.exact_owner), phase_hostile.gap, 0,
    )
    naive_bank = tuple(
        (chain[2] >> phase_hostile.gap) & ((1 << OWNER_DEPTH) - 1)
        for chain in naive_chains
    )
    require(naive_bank != phase_hostile.next_bank,
            "zero-ray transport hostile")
    require(
        tuple(
            (chain[2] >> phase_hostile.gap) & ((1 << OWNER_DEPTH) - 1)
            for chain in phase_hostile.chains
        ) == phase_hostile.next_bank,
        "phase-owned-base repair",
    )

    projective_key, projective_left, projective_right = hostile_by_level["projective"]
    carry_cross_key, carry_cross_left, carry_cross_right = cross_by_level["carry"]
    require(projective_left.scale == projective_right.scale == MAX_SCALE,
            "terminal finite mismatch scale")
    require(carry_cross_left.scale != carry_cross_right.scale,
            "cross-scale finite mismatch")

    semantic = (
        DEPENDENCY_SHA256, MAX_SCALE, OWNER_DEPTH, SELECTED_RAYS,
        ODD_BLOCK_LENGTH, ODD_PRECISION, PROJECTIVE_HISTORY,
        PROJECTIVE_PRECISION, tuple(v[: MAX_SCALE + 2]),
        tuple(gaps[: MAX_SCALE + 1]), tuple(same_scale_ranges),
        tuple(census), tuple(first_hostiles),
        tuple(compact_record(record) for record in records),
    )
    semantic_sha256 = hashlib.sha256(
        json.dumps(semantic, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("RULE30_CROSS_SCALE_FIRST_MISMATCH_20260824")
    print("status=THM-4006_FINITE-EXACT;physical_highest-shell_chronology;no_all-scale_claim;no_rule30_prize")
    print(f"dependency_sha256={DEPENDENCY_SHA256}")
    print(
        f"universe=n=1..{records[-1].n};scales=0..{MAX_SCALE};"
        f"records={len(records)};phase_rule=0<=t<2^m;n=2^m+t"
    )
    print(
        f"owner_depth={OWNER_DEPTH};selected_rays={SELECTED_RAYS};"
        f"odd_block_length={ODD_BLOCK_LENGTH};odd_precision={ODD_PRECISION};"
        f"projective_history={PROJECTIVE_HISTORY};"
        f"projective_precision={PROJECTIVE_PRECISION}"
    )
    print(f"v={tuple(v[: MAX_SCALE + 2])};gaps={tuple(gaps[: MAX_SCALE + 1])}")
    print("same_scale_defect_ranges=" + repr(tuple(same_scale_ranges)).replace(" ", ""))
    print(f"same_scale_strict_checks={strict_checks}")
    print("quotient_census=" + repr(tuple(census)).replace(" ", ""))
    print(
        "finite_first_mismatch_after_projective="
        + repr(
            (
                projective_key,
                compact_record(projective_left),
                compact_record(projective_right),
            )
        ).replace(" ", "")
    )
    print(
        "finite_first_cross_scale_mismatch_after_carry="
        + repr(
            (
                carry_cross_key,
                compact_record(carry_cross_left),
                compact_record(carry_cross_right),
            )
        ).replace(" ", "")
    )
    print(
        "phase_owned_base_hostile="
        + repr(
            (
                (phase_hostile.n, phase_hostile.scale, phase_hostile.phase),
                naive_bank,
                phase_hostile.next_bank,
                phase_hostile.chains,
            )
        ).replace(" ", "")
    )
    print(
        "projective_mismatch_phase_base_repair="
        + repr(
            (
                projective_left.chains[0], projective_right.chains[0],
                projective_left.next_bank, projective_right.next_bank,
            )
        ).replace(" ", "")
    )
    print(f"semantic_sha256={semantic_sha256}")
    print(
        "interpretation=finite physical cross-scale quotient;branch+owner+odd-block+carry+"
        "projective-history+phase-base-chain are tested before routed off-ray chains;"
        "two-chain candidate has no target mismatch only in declared universe"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
