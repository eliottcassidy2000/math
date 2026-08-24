#!/usr/bin/env python3
"""Exact physical Rule 30 depth-15 observer census for THM-4013.

The phase owner is advanced through the exact conjugacy
    s_(m,t) = r_(m,t)^(-1) s_m r_(m,t)
on the required finite tree.  This avoids rebuilding a 2**m-letter section
at every phase while retaining the actual tail-prefix route.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


sys.stdout.reconfigure(newline="\n")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEPENDENCY = ROOT / "04-computation" / "rule30_orbit_signalizer_gap_thm3511.py"
DEPENDENCY_SHA256 = "2ce110f0b8e9c71c3d298aaf07e8e6c02b70d33e5671bc763f3f3b490caa5445"

MAX_SCALE = 15
OWNER_DEPTH = 4
ODD_BLOCK_LENGTH = 4
ARITHMETIC_PRECISION = 4
MAX_PROJECTIVE_HISTORY = 6
LOW_WIDTH = 64


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


def exact_rows(stop: int) -> list[int]:
    rows = [1]
    row = 1
    for _ in range(stop):
        row = row ^ ((row << 1) | (row << 2))
        rows.append(row)
    return rows


def low_rows(stop: int) -> list[int]:
    mask = (1 << LOW_WIDTH) - 1
    rows = [1]
    row = 1
    for _ in range(stop):
        row = (row ^ ((row << 1) | (row << 2))) & mask
        rows.append(row)
    return rows


def unit_exact(rows: list[int], v: tuple[int, ...], scale: int, phase: int) -> int:
    numerator = rows[phase + (1 << scale)] - rows[phase]
    require(numerator % (1 << v[scale]) == 0,
            ("exact unit integrality", scale, phase))
    answer = numerator >> v[scale]
    require(answer & 1, ("exact unit odd", scale, phase))
    return answer


def unit_mod(rows: list[int], v: tuple[int, ...], scale: int, phase: int,
             precision: int) -> int:
    require(v[scale] + precision <= LOW_WIDTH, "low-width unit precision")
    numerator = (rows[phase + (1 << scale)] - rows[phase]) % (1 << LOW_WIDTH)
    require(numerator % (1 << v[scale]) == 0,
            ("modular unit integrality", scale, phase))
    answer = (numerator >> v[scale]) % (1 << precision)
    require(answer & 1, ("modular unit odd", scale, phase))
    return answer


def center_carry(rows: list[int], v: tuple[int, ...], n: int):
    alpha = (n & -n).bit_length() - 1
    va = v[alpha]
    center = (rows[n] >> n) & 1
    if va > n:
        require(center == 0, ("pre-visible center", n))
        return ("pre",), center, 0

    level = n - va
    modulus = 1 << (level + 1)
    target_xor = 0
    lower_sum = 0
    term_digest = hashlib.sha256()
    for scale in range(n.bit_length()):
        if not ((n >> scale) & 1):
            continue
        phase = n & ((1 << scale) - 1)
        residue = (
            (1 << (v[scale] - va)) * unit_exact(rows, v, scale, phase)
        ) % modulus
        target = (residue >> level) & 1
        lower = residue & ((1 << level) - 1)
        target_xor ^= target
        lower_sum += lower
        lower_bytes = lower.to_bytes((level + 7) // 8, "little") if level else b""
        term_digest.update(scale.to_bytes(2, "little"))
        term_digest.update(bytes((target,)))
        term_digest.update(len(lower_bytes).to_bytes(4, "little"))
        term_digest.update(lower_bytes)
    carry = (lower_sum >> level) & 1 if level else 0
    require(center == (target_xor ^ carry), ("ordinary carry", n))
    return ("visible", target_xor, carry), center, int.from_bytes(
        term_digest.digest()[:8], "big"
    )


def projective_history(rows: list[int], v: tuple[int, ...], gaps: tuple[int, ...],
                       scale: int, phase: int) -> tuple[object, ...]:
    entries: list[object] = [None] * max(0, MAX_PROJECTIVE_HISTORY - scale - 1)
    first = max(0, scale - MAX_PROJECTIVE_HISTORY + 1)
    for level in range(first, scale + 1):
        d = gaps[level]
        precision = ARITHMETIC_PRECISION + d
        modulus = 1 << precision
        u0 = unit_mod(rows, v, level, phase, precision)
        u1 = unit_mod(rows, v, level, phase + (1 << level), precision)
        g = (-u1 * pow(u0, -1, modulus)) % modulus
        require((g - 1) % (1 << d) == 0,
                ("projective valuation", level, phase))
        z = ((1 - g) >> d) % (1 << ARITHMETIC_PRECISION)
        entries.append((d, g % (1 << ARITHMETIC_PRECISION), z))
    require(len(entries) == MAX_PROJECTIVE_HISTORY, "projective history length")
    return tuple(entries)


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    """Right action product: apply left, then right."""
    return tuple(right[left[value]] for value in range(len(left)))


def invert(permutation: tuple[int, ...]) -> tuple[int, ...]:
    out = [0] * len(permutation)
    for source, target in enumerate(permutation):
        out[target] = source
    return tuple(out)


@dataclass(frozen=True)
class Record:
    n: int
    scale: int
    phase: int
    gap: int
    extension: int
    owner_bank: tuple[int, ...]
    owner_portrait: tuple[int, ...]
    odd_block: tuple[int, ...]
    carry_state: tuple[object, ...]
    carry_digest: int
    projective: tuple[object, ...]
    chains: tuple[tuple[int, int, int], ...]
    center: int
    next_bank: tuple[int, ...]

    @property
    def target(self):
        return self.center, self.next_bank


def build_records(module):
    words, raw_gaps, _, _ = module.signalizer_words()
    v = tuple(module.EXPECTED_V_0_TO_23)
    gaps = tuple(raw_gaps)
    maximum_n = (1 << (MAX_SCALE + 1)) - 1
    maximum_q = 1 << MAX_SCALE
    rows = exact_rows(maximum_n)
    low = low_rows(5 * maximum_q)
    require(all((rows[index] % (1 << LOW_WIDTH)) == low[index]
                for index in range(maximum_n + 1)), "exact/low row agreement")

    records: list[Record] = []
    by_scale_phase: dict[tuple[int, int], Record] = {}
    owner_controls = conjugacy_controls = route_controls = 0

    for scale in range(MAX_SCALE + 1):
        q = 1 << scale
        gap = gaps[scale]
        depth = OWNER_DEPTH + gap
        size = 1 << depth
        actions = module.generator_actions(depth)
        marked = module.permutation(words[scale], actions, {})
        transport = tuple(range(size))
        transport_inverse = transport

        for phase in range(q):
            # r^-1 s r in right-action convention: apply r^-1, then s, then r.
            owner = tuple(
                transport[marked[transport_inverse[value]]]
                for value in range(size)
            )
            require(sorted(owner) == list(range(size)),
                    ("owner permutation", scale, phase))
            conjugacy_controls += size

            tail = (low[phase] >> v[scale]) % size
            first_tail = (low[phase + q] >> v[scale]) % size
            second_tail = (low[phase + 2 * q] >> v[scale]) % size
            require(owner[tail] == first_tail and owner[owner[tail]] == second_tail,
                    ("physical owner tail", scale, phase))
            owner_controls += 2

            extension = tail & ((1 << gap) - 1)
            chains = []
            recovered = []
            for ray in (0, 1, 2):
                x = extension | (ray << gap)
                first_image = owner[x]
                second_image = owner[first_image]
                require(second_image & ((1 << gap) - 1) == extension,
                        ("fixed physical base", scale, phase, ray))
                chains.append((x, first_image, second_image))
                recovered.append((second_image >> gap) & 15)
                route_controls += 1

            n = q + phase
            carry_state, center, carry_digest = center_carry(rows, v, n)
            odd_block = tuple(
                unit_mod(low, v, scale, phase + j * q, ARITHMETIC_PRECISION)
                for j in range(ODD_BLOCK_LENGTH)
            )
            record = Record(
                n=n,
                scale=scale,
                phase=phase,
                gap=gap,
                extension=extension,
                owner_bank=tuple(owner[ray] & 15 for ray in (0, 1, 2)),
                owner_portrait=tuple(owner[ray] & 15 for ray in range(16)),
                odd_block=odd_block,
                carry_state=carry_state,
                carry_digest=carry_digest,
                projective=projective_history(low, v, gaps, scale, phase),
                chains=tuple(chains),
                center=center,
                next_bank=tuple(recovered),
            )
            records.append(record)
            by_scale_phase[(scale, phase)] = record

            # r_(t+1)=r_t * (A|_(u_m^(A^t))).
            prefix = rows[phase] & ((1 << v[scale]) - 1)
            generator_word = module.section_at_prefix(
                bytearray((module.A,)), prefix, v[scale],
            )
            require(len(generator_word) == 1, "one-state transport increment")
            generator = actions[generator_word[0]]
            generator_inverse = invert(generator)
            transport = compose(transport, generator)
            transport_inverse = compose(generator_inverse, transport_inverse)
            require(all(transport_inverse[transport[x]] == x for x in range(size)),
                    ("transport inverse", scale, phase))

    require(tuple(record.n for record in records) == tuple(range(1, maximum_n + 1)),
            "chronological universe")

    transition_controls = 0
    for record in records:
        if record.scale == MAX_SCALE:
            continue
        direct_next = by_scale_phase[(record.scale + 1, record.phase)].owner_bank
        require(record.next_bank == direct_next,
                ("cross-scale direct owner bank", record.scale, record.phase))
        transition_controls += 1

    return (
        records, rows, low, v, gaps,
        (owner_controls, conjugacy_controls, route_controls, transition_controls),
    )


def base_state(record: Record, history: int = 3):
    require(0 <= history <= MAX_PROJECTIVE_HISTORY, "history range")
    projective = record.projective[-history:] if history else ()
    return (
        record.gap,
        record.owner_bank,
        record.owner_portrait,
        record.odd_block,
        record.carry_state,
        projective,
    )


def state(record: Record, label: str):
    if label.startswith("history"):
        return base_state(record, int(label.removeprefix("history")))
    answer = list(base_state(record, 3))
    if label in {"extension", "base_first", "base_second", "base_chain",
                 "base_plus_c1_first", "base_plus_c2_first", "base_plus_both_first",
                 "base_plus_c1", "base_plus_c2", "base_plus_both"}:
        answer.append(record.extension)
    if label in {"base_first", "base_chain", "base_plus_c1_first",
                 "base_plus_c2_first", "base_plus_both_first", "base_plus_c1",
                 "base_plus_c2", "base_plus_both"}:
        answer.append(record.chains[0][1])
    if label in {"base_second", "base_chain", "base_plus_c1_first",
                 "base_plus_c2_first", "base_plus_both_first", "base_plus_c1",
                 "base_plus_c2", "base_plus_both"}:
        answer.append(record.chains[0][2])
    if label in {"base_plus_c1_first", "base_plus_both_first"}:
        answer.append(record.chains[1][1])
    if label in {"base_plus_c2_first", "base_plus_both_first"}:
        answer.append(record.chains[2][1])
    if label in {"base_plus_c1", "base_plus_both"}:
        answer.append(record.chains[1])
    if label in {"base_plus_c2", "base_plus_both"}:
        answer.append(record.chains[2])
    return tuple(answer)


def census(records: list[Record], label: str):
    fibres: dict[tuple[object, ...], list[Record]] = defaultdict(list)
    for record in records:
        fibres[state(record, label)].append(record)

    mismatches = same = cross = 0
    first = first_cross = None
    for key, rows in fibres.items():
        targets = {row.target for row in rows}
        if len(targets) <= 1:
            continue
        mismatches += 1
        scale_targets: dict[int, set[tuple[object, ...]]] = defaultdict(set)
        for row in rows:
            scale_targets[row.scale].add(row.target)
        if any(len(values) > 1 for values in scale_targets.values()):
            same += 1
        scales = sorted(scale_targets)
        has_cross = any(
            any(left != right for left in scale_targets[a] for right in scale_targets[b])
            for index, a in enumerate(scales) for b in scales[index + 1:]
        )
        cross += has_cross

        seen: dict[tuple[object, ...], Record] = {}
        seen_by_scale: dict[int, dict[tuple[object, ...], Record]] = defaultdict(dict)
        for row in rows:
            for target, old in seen.items():
                if target != row.target:
                    candidate = (row.n, old.n, key, old, row)
                    if first is None or candidate[:2] < first[:2]:
                        first = candidate
            for old_scale, old_targets in seen_by_scale.items():
                if old_scale == row.scale:
                    continue
                for target, old in old_targets.items():
                    if target != row.target:
                        candidate = (row.n, old.n, key, old, row)
                        if first_cross is None or candidate[:2] < first_cross[:2]:
                            first_cross = candidate
            seen.setdefault(row.target, row)
            seen_by_scale[row.scale].setdefault(row.target, row)

    nontrivial = sum(len(rows) > 1 for rows in fibres.values())
    return (len(fibres), nontrivial, mismatches, same, cross), first, first_cross


def conditional_census(records: list[Record], first_stage_only: bool):
    fibres: dict[tuple[object, ...], list[Record]] = defaultdict(list)
    for record in records:
        fibres[state(record, "base_chain")].append(record)

    labels = (
        "base_enough", "either_one", "c1_only", "c2_only", "both_needed",
        "pair_unresolved",
    )
    fibre_counts = {label: 0 for label in labels}
    record_counts = {label: 0 for label in labels}
    controls = {}

    def observation(record: Record, index: int):
        return record.chains[index][1] if first_stage_only else record.chains[index]

    def determines(rows: list[Record], indices: tuple[int, ...]) -> bool:
        quotient: dict[tuple[object, ...], set[tuple[object, ...]]] = defaultdict(set)
        for row in rows:
            quotient[tuple(observation(row, index) for index in indices)].add(row.target)
        return all(len(targets) == 1 for targets in quotient.values())

    for rows in fibres.values():
        if len({row.target for row in rows}) == 1:
            label = "base_enough"
        else:
            c1 = determines(rows, (1,))
            c2 = determines(rows, (2,))
            pair = determines(rows, (1, 2))
            if not pair:
                label = "pair_unresolved"
            elif c1 and c2:
                label = "either_one"
            elif c1:
                label = "c1_only"
            elif c2:
                label = "c2_only"
            else:
                label = "both_needed"
            controls.setdefault(label, tuple(compact(row) for row in rows[:3]))
        fibre_counts[label] += 1
        record_counts[label] += len(rows)

    return (
        tuple((label, fibre_counts[label], record_counts[label]) for label in labels),
        tuple((label, controls.get(label)) for label in labels if label != "base_enough"),
    )


def compact(record: Record):
    return (
        record.n, record.scale, record.phase, record.gap, record.extension,
        record.owner_bank, record.odd_block, record.carry_state,
        record.projective[-3:], record.chains, record.target,
    )


def main() -> None:
    module = load_dependency()
    records, rows, low, v, gaps, controls = build_records(module)

    # Same-scale control: verify the exact top-bit inequalities used by
    # THM-3824 through every time needed for m<=MAX_SCALE, then check its
    # positive comparison polynomial for each physical scale.
    row = 1
    top_bit_checks = 0
    top_stop = 3 * (1 << MAX_SCALE)
    for time in range(1, top_stop + 1):
        row = row ^ ((row << 1) | (row << 2))
        require((3 * (1 << (2 * time))) < 2 * row < (4 * (1 << (2 * time))),
                ("THM-3824 top-bit bounds", time))
        top_bit_checks += 1
    scale_inequalities = []
    for scale in range(MAX_SCALE + 1):
        q = 1 << scale
        c = 1 << (2 * q)
        margin = 4 * c * c - 13 * c + 4
        require(margin > 0, ("strict-defect margin", scale))
        scale_inequalities.append((scale, q, margin.bit_length(), margin % 65536))

    labels = tuple(f"history{k}" for k in range(MAX_PROJECTIVE_HISTORY + 1)) + (
        "extension", "base_first", "base_second", "base_chain",
        "base_plus_c1_first", "base_plus_c2_first", "base_plus_both_first",
        "base_plus_c1", "base_plus_c2", "base_plus_both",
    )
    censuses = []
    witnesses = {}
    cross_witnesses = {}
    for label in labels:
        counts, first, first_cross = census(records, label)
        censuses.append((label,) + counts)
        if first is not None:
            witnesses[label] = (first[2], compact(first[3]), compact(first[4]))
        if first_cross is not None:
            cross_witnesses[label] = (
                first_cross[2], compact(first_cross[3]), compact(first_cross[4])
            )

    # Identify the first declared phase-owned refinement with no mismatch.
    phase_labels = (
        "extension", "base_first", "base_second", "base_chain",
        "base_plus_c1_first", "base_plus_c2_first", "base_plus_both_first",
        "base_plus_c1", "base_plus_c2", "base_plus_both",
    )
    no_mismatch = [label for label in phase_labels
                   if next(row[3] for row in censuses if row[0] == label) == 0]
    require(no_mismatch, "phase-owned no-mismatch candidate")
    first_phase_repair = no_mismatch[0]

    require("history3" in witnesses, "three-level projective hostile")
    conditional_first = conditional_census(records, True)
    conditional_full = conditional_census(records, False)

    semantic = (
        DEPENDENCY_SHA256, MAX_SCALE, OWNER_DEPTH, ODD_BLOCK_LENGTH,
        ARITHMETIC_PRECISION, MAX_PROJECTIVE_HISTORY,
        tuple(v[: MAX_SCALE + 2]), tuple(gaps[: MAX_SCALE + 1]), controls,
        top_bit_checks, tuple(scale_inequalities), tuple(censuses),
        first_phase_repair, conditional_first, conditional_full,
        tuple((label, witnesses.get(label)) for label in labels),
        tuple((label, cross_witnesses.get(label)) for label in labels),
        tuple(compact(record) for record in records),
    )
    semantic_sha256 = hashlib.sha256(
        json.dumps(semantic, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("RULE30_DEPTH15_HISTORY_ADAPTIVE_ROUTE_THM4013")
    print("status=FINITE-EXACT_CANONICAL;no_all-scale_claim;no_rule30_prize")
    print(f"dependency_sha256={DEPENDENCY_SHA256}")
    print(
        f"universe=n=1..{records[-1].n};scales=0..{MAX_SCALE};records={len(records)};"
        "phase_rule=0<=t<2^m;n=2^m+t"
    )
    print(
        f"owner_depth={OWNER_DEPTH};odd_block={ODD_BLOCK_LENGTH};"
        f"arithmetic_precision={ARITHMETIC_PRECISION};"
        f"maximum_projective_history={MAX_PROJECTIVE_HISTORY}"
    )
    print(f"v={tuple(v[:MAX_SCALE+2])};gaps={tuple(gaps[:MAX_SCALE+1])}")
    print(f"owner_conjugacy_route_transition_controls={controls}")
    print(f"same_scale_top_bit_checks={top_bit_checks}")
    print("same_scale_margin_controls=" + repr(tuple(scale_inequalities)).replace(" ", ""))
    print("censuses=" + repr(tuple(censuses)).replace(" ", ""))
    print(f"first_phase_owned_no_mismatch_refinement={first_phase_repair}")
    print("conditional_first_address_census=" + repr(conditional_first[0]).replace(" ", ""))
    print("conditional_full_chain_census=" + repr(conditional_full[0]).replace(" ", ""))
    first_controls = dict(conditional_first[1])
    full_controls = dict(conditional_full[1])
    print(
        "first_address_pair_unresolved_control="
        + repr(tuple(row[0] for row in first_controls["pair_unresolved"])).replace(" ", "")
    )
    print(
        "conditional_full_direction_controls="
        + repr(
            (
                tuple(row[0] for row in full_controls["c1_only"]),
                tuple(row[0] for row in full_controls["c2_only"]),
            )
        ).replace(" ", "")
    )
    print("history3_first_mismatch=" + repr(witnesses["history3"]).replace(" ", ""))
    print(
        "history3_first_cross_scale_mismatch="
        + (
            repr(cross_witnesses["history3"]).replace(" ", "")
            if "history3" in cross_witnesses else "NONE_IN_DECLARED_UNIVERSE"
        )
    )
    if first_phase_repair in witnesses:
        print("repair_first_mismatch=" + repr(witnesses[first_phase_repair]).replace(" ", ""))
    else:
        print("repair_first_mismatch=NONE_IN_DECLARED_UNIVERSE")
    print(f"semantic_sha256={semantic_sha256}")
    print(
        "interpretation=three-level_projective_history_is_tested_after_carry;"
        "phase-owned refinements are compared componentwise;"
        "first no-mismatch refinement is finite-universe only"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
