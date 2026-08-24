#!/usr/bin/env python3
"""No-import exact audit of the THM-4006 Rule 30 cross-scale observer census.

The Mealy action, signalizer recursion, physical phase routing, Rule 30 rows,
ordinary carry, projective ratios, quotient fibres, and chronological first-
mismatch tests are rebuilt here.  No producer module is imported.
"""

from collections import defaultdict
from dataclasses import dataclass
from hashlib import sha256
import json
import math
import sys


sys.stdout.reconfigure(newline="\n")

A, B, C = 0, 1, 2
MAX_SCALE = 9
PORTRAIT_DEPTH = 4
RAYS = (0, 1, 2)
BLOCK = 4
MOD_BITS = 4
HISTORY = 3
EXPECTED_SEMANTIC_SHA256 = "47c79b0e6ab72cc6519b32602608b77ddb9da7b9043c4579b97d2e159cc2231b"
GATES = 0

# (output bit, next section state), read low bit first.
MEALY = {
    A: ((0, A), (1, B)),
    B: ((1, C), (0, B)),
    C: ((1, A), (0, B)),
}


def gate(condition, label):
    global GATES
    GATES += 1
    if not bool(condition):
        raise RuntimeError(label)


def v2(value):
    gate(value != 0, "v2 nonzero")
    return (value & -value).bit_length() - 1


def root_section(word, input_bit):
    """Section of a right-action product, derived from the Mealy table."""
    current = input_bit
    section = []
    for state in word:
        output, child = MEALY[state][current]
        section.append(child)
        current = output
    return tuple(section), current


def section_at_integer(word, prefix, length):
    answer = tuple(word)
    for depth in range(length):
        answer, _ = root_section(answer, (prefix >> depth) & 1)
    return answer


def active(word):
    return sum(state != A for state in word) & 1


def act_state(state, value, width):
    output_value = 0
    current_state = state
    for depth in range(width):
        bit = (value >> depth) & 1
        output, current_state = MEALY[current_state][bit]
        output_value |= output << depth
    return output_value


def act_word(word, value, width):
    answer = value
    for state in word:
        answer = act_state(state, answer, width)
    return answer


def signalizer_step(word):
    gate(active(word) == 1, "active signalizer input")
    left, left_root = root_section(word, 0)
    right, right_root = root_section(word, 1)
    gate((left_root, right_root) == (1, 0), "active child routing")
    candidate = left + right
    gap = 1
    while active(candidate) == 0:
        candidate, root = root_section(candidate, 0)
        gate(root == 0, "fixed prefix descent")
        gap += 1
        gate(gap < 64, "finite signalizer probe")
    return candidate, gap


def signalizer_prefix():
    word = (B,)
    valuation = [1]
    gaps = []
    words = []
    for scale in range(MAX_SCALE + 1):
        gate(len(word) == 1 << scale, f"signalizer length {scale}")
        gate(active(word) == 1, f"signalizer activity {scale}")
        words.append(word)
        word, gap = signalizer_step(word)
        gaps.append(gap)
        valuation.append(valuation[-1] + gap)
    return tuple(words), tuple(valuation), tuple(gaps)


def rule30_rows(stop):
    row = 1
    rows = [row]
    for _ in range(stop):
        row = row ^ ((row << 1) | (row << 2))
        rows.append(row)
    return rows


def normalized_unit(rows, valuations, scale, phase):
    q = 1 << scale
    numerator = rows[phase + q] - rows[phase]
    divisor = 1 << valuations[scale]
    gate(numerator % divisor == 0, ("unit integral", scale, phase))
    answer = numerator // divisor
    gate(answer & 1, ("unit odd", scale, phase))
    return answer


def image_bank(owner):
    return tuple(act_word(owner, ray, PORTRAIT_DEPTH) for ray in RAYS)


def portrait(owner):
    return tuple(act_word(owner, ray, PORTRAIT_DEPTH)
                 for ray in range(1 << PORTRAIT_DEPTH))


def route(owner, gap, extension):
    gate(0 <= extension < 1 << gap, "extension range")
    depth = PORTRAIT_DEPTH + gap
    mask = (1 << depth) - 1
    result = []
    for ray in RAYS:
        source = extension + (ray << gap)
        first = act_word(owner, source, depth)
        second = act_word(owner, first, depth)
        gate(0 <= first <= mask and 0 <= second <= mask, "route range")
        result.append((source, first, second))
    return tuple(result)


def bank_from_routes(routes, gap):
    return tuple((second >> gap) & 15 for _source, _first, second in routes)


def carry_observer(rows, valuations, n):
    alpha = v2(n)
    base_v = valuations[alpha]
    center = (rows[n] >> n) & 1
    if base_v > n:
        gate(center == 0, ("pre-visible center", n))
        return ("pre",), center, ()

    level = n - base_v
    modulus = 1 << (level + 1)
    target_xor = 0
    low_sum = 0
    terms = []
    for scale in range(n.bit_length()):
        if not ((n >> scale) & 1):
            continue
        phase = n % (1 << scale) if scale else 0
        odd = normalized_unit(rows, valuations, scale, phase)
        residue = ((1 << (valuations[scale] - base_v)) * odd) % modulus
        target = (residue >> level) & 1
        low = residue & ((1 << level) - 1)
        target_xor ^= target
        low_sum += low
        terms.append((scale, target, low))
    carry = ((low_sum >> level) & 1) if level else 0
    gate(center == (target_xor ^ carry), ("ordinary carry", n))
    return ("visible", target_xor, carry), center, tuple(terms)


def scale_history(rows, valuations, gaps, scale, phase):
    finite = [None] * max(0, HISTORY - scale - 1)
    exact = [None] * max(0, HISTORY - scale - 1)
    for level in range(max(0, scale - HISTORY + 1), scale + 1):
        gap = gaps[level]
        modulus = 1 << (MOD_BITS + gap)
        first = normalized_unit(rows, valuations, level, phase)
        second = normalized_unit(rows, valuations, level, phase + (1 << level))
        ratio = (-second * pow(first, -1, modulus)) % modulus
        gate((ratio - 1) % (1 << gap) == 0,
             ("ratio gap", level, phase))
        quotient = ((1 - ratio) >> gap) % (1 << MOD_BITS)
        finite.append((gap, ratio % 16, quotient))

        g_num, g_den = -second, first
        common = math.gcd(abs(g_num), g_den)
        g_num, g_den = g_num // common, g_den // common
        z_num, z_den = first + second, (1 << gap) * first
        common = math.gcd(abs(z_num), z_den)
        z_num, z_den = z_num // common, z_den // common
        exact.append((gap, g_num, g_den, z_num, z_den))
    gate(len(finite) == HISTORY and len(exact) == HISTORY,
         "history length")
    return tuple(finite), tuple(exact)


@dataclass(frozen=True)
class Snapshot:
    n: int
    scale: int
    phase: int
    gap: int
    extension: int
    bank: tuple
    full_portrait: tuple
    odd_shadow: tuple
    exact_odds: tuple
    carry: tuple
    carry_terms: tuple
    projective: tuple
    exact_projective: tuple
    chains: tuple
    center: int
    next_bank: tuple
    owner: tuple
    owner_digest: str

    @property
    def target(self):
        return self.center, self.next_bank


def build_universe():
    words, valuations, gaps = signalizer_prefix()
    gate(valuations == (1, 3, 4, 6, 7, 9, 15, 16, 24, 25, 27),
         "valuation prefix")
    gate(gaps == (2, 1, 2, 1, 2, 6, 1, 8, 1, 2), "gap prefix")
    gate(tuple(len(word) for word in words) == tuple(1 << m for m in range(10)),
         "signalizer word lengths")

    qmax = 1 << MAX_SCALE
    rows = rule30_rows(5 * qmax)
    records = []
    for scale in range(MAX_SCALE + 1):
        q = 1 << scale
        gap = gaps[scale]
        power = (A,) * q
        doubled = (A,) * (2 * q)
        for phase in range(q):
            n = q + phase
            prefix = rows[phase] % (1 << valuations[scale])
            owner = section_at_integer(power, prefix, valuations[scale])
            gate(active(owner) == 1, ("physical owner active", scale, phase))

            extension = (rows[phase] >> valuations[scale]) % (1 << gap)
            induced = section_at_integer(owner + owner, extension, gap)
            long_prefix = rows[phase] % (1 << (valuations[scale] + gap))
            direct = section_at_integer(
                doubled, long_prefix, valuations[scale] + gap
            )
            gate(induced == direct, ("phase transition", scale, phase))

            chains = route(owner, gap, extension)
            next_bank = image_bank(induced)
            gate(bank_from_routes(chains, gap) == next_bank,
                 ("route recovers bank", scale, phase))

            carry, center, terms = carry_observer(rows, valuations, n)
            exact_odds = tuple(
                normalized_unit(rows, valuations, scale, phase + j * q)
                for j in range(BLOCK)
            )
            odd_shadow = tuple(value % 16 for value in exact_odds)
            finite_projective, exact_projective = scale_history(
                rows, valuations, gaps, scale, phase
            )
            owner_bytes = bytes(owner)
            records.append(Snapshot(
                n=n,
                scale=scale,
                phase=phase,
                gap=gap,
                extension=extension,
                bank=image_bank(owner),
                full_portrait=portrait(owner),
                odd_shadow=odd_shadow,
                exact_odds=exact_odds,
                carry=carry,
                carry_terms=terms,
                projective=finite_projective,
                exact_projective=exact_projective,
                chains=chains,
                center=center,
                next_bank=next_bank,
                owner=owner,
                owner_digest=sha256(owner_bytes).hexdigest(),
            ))

    gate(tuple(record.n for record in records) == tuple(range(1, 1024)),
         "highest-shell chronology is n=1..1023 exactly once")
    return records, rows, valuations, gaps


def observer(record, level):
    state = [record.gap]
    if level in {"bank", "portrait", "odd", "carry", "projective",
                 "base", "one", "two", "exact"}:
        state.append(record.bank)
    if level in {"portrait", "odd", "carry", "projective",
                 "base", "one", "two", "exact"}:
        state.append(record.full_portrait)
    if level in {"odd", "carry", "projective", "base", "one", "two"}:
        state.append(record.odd_shadow)
    if level in {"carry", "projective", "base", "one", "two"}:
        state.append(record.carry)
    if level in {"projective", "base", "one", "two"}:
        state.append(record.projective)
    if level in {"base", "one", "two"}:
        state.append(record.chains[0])
    if level in {"one", "two"}:
        state.append(record.chains[1])
    if level == "two":
        state.append(record.chains[2])
    if level == "exact":
        state.extend((record.owner, record.exact_odds,
                      (record.carry, record.carry_terms),
                      record.exact_projective))
    return tuple(state)


def analyze(records, level):
    fibres = defaultdict(list)
    for record in records:
        fibres[observer(record, level)].append(record)

    collision_fibres = [rows for rows in fibres.values() if len(rows) > 1]
    mismatch_fibres = []
    same_count = 0
    cross_count = 0
    mismatch_pairs = []
    cross_pairs = []
    for key, rows in fibres.items():
        has_mismatch = has_same = has_cross = False
        for i, left in enumerate(rows):
            for right in rows[i + 1:]:
                if left.target == right.target:
                    continue
                has_mismatch = True
                pair = (right.n, left.n, repr(key), left, right)
                mismatch_pairs.append(pair)
                if left.scale == right.scale:
                    has_same = True
                else:
                    has_cross = True
                    cross_pairs.append(pair)
        if has_mismatch:
            mismatch_fibres.append(rows)
        same_count += has_same
        cross_count += has_cross

    first_repeat = None
    for rows in collision_fibres:
        ordered = sorted(rows, key=lambda item: item.n)
        candidate = (ordered[1].n, ordered[0].n, ordered[0], ordered[1])
        if first_repeat is None or candidate[:2] < first_repeat[:2]:
            first_repeat = candidate
    first_mismatch = min(mismatch_pairs, key=lambda item: item[:3]) if mismatch_pairs else None
    first_cross = min(cross_pairs, key=lambda item: item[:3]) if cross_pairs else None
    return {
        "fibres": fibres,
        "collision_fibres": collision_fibres,
        "mismatch_fibres": mismatch_fibres,
        "same": same_count,
        "cross": cross_count,
        "first_repeat": first_repeat,
        "first_mismatch": first_mismatch,
        "first_cross": first_cross,
    }


def compact(record):
    return (
        record.n, record.scale, record.phase, record.gap, record.bank,
        record.odd_shadow, record.carry, record.projective,
        record.chains, record.target, record.owner_digest,
    )


def integer_digest(value):
    return sha256(str(value).encode("ascii")).hexdigest()


def main():
    records, rows, valuations, gaps = build_universe()

    # Independent THM-3824 same-scale control.
    ranges = []
    strict_control_gates = 0
    defect_by_scale = {}
    for scale in range(MAX_SCALE + 1):
        q = 1 << scale
        defects = tuple(
            normalized_unit(rows, valuations, scale, phase + q)
            - normalized_unit(rows, valuations, scale, phase)
            for phase in range(q)
        )
        for index, value in enumerate(defects):
            gate(value > 0, ("positive exact defect", scale, index))
            strict_control_gates += 1
        for index in range(len(defects) - 1):
            gate(defects[index] < defects[index + 1],
                 ("strict exact defect", scale, index))
            strict_control_gates += 1
        ranges.append((scale, len(defects), defects[0].bit_length(),
                       defects[-1].bit_length(), defects[0] % 65536,
                       defects[-1] % 65536))
        defect_by_scale[scale] = defects
    gate(strict_control_gates == 2036, "THM3824 strict-control gate count")
    expected_ranges = (
        (0, 1, 3, 3, 6, 6), (1, 2, 6, 8, 44, 196),
        (2, 4, 13, 19, 6378, 58722), (3, 8, 27, 41, 12416, 54368),
        (4, 16, 58, 88, 34118, 61694), (5, 32, 120, 182, 22038, 23446),
        (6, 64, 242, 368, 42356, 51636), (7, 128, 497, 751, 58834, 15834),
        (8, 256, 1001, 1511, 26148, 38836),
        (9, 512, 2024, 3046, 25710, 54918),
    )
    gate(tuple(ranges) == expected_ranges, "same-scale defect ranges")

    levels = ("branch", "bank", "portrait", "odd", "carry", "projective",
              "base", "one", "two", "exact")
    analyses = {level: analyze(records, level) for level in levels}
    census = tuple(
        (level, len(analyses[level]["fibres"]),
         len(analyses[level]["collision_fibres"]),
         len(analyses[level]["mismatch_fibres"]),
         analyses[level]["same"], analyses[level]["cross"])
        for level in levels
    )
    expected_census = (
        ("branch", 4, 4, 4, 4, 2), ("bank", 246, 185, 185, 176, 58),
        ("portrait", 395, 229, 226, 224, 15),
        ("odd", 753, 190, 183, 181, 2),
        ("carry", 870, 128, 109, 107, 2),
        ("projective", 1019, 4, 1, 1, 0),
        ("base", 1022, 1, 0, 0, 0),
        ("one", 1023, 0, 0, 0, 0),
        ("two", 1023, 0, 0, 0, 0),
        ("exact", 1023, 0, 0, 0, 0),
    )
    gate(census == expected_census, "complete refinement census")

    projective = analyses["projective"]
    gate(projective["first_mismatch"] is not None, "projective mismatch exists")
    later, earlier, _key, left, right = projective["first_mismatch"]
    gate((earlier, later) == (943, 951), "first projective target mismatch chronology")
    gate((left.scale, left.phase, right.scale, right.phase) == (9, 431, 9, 439),
         "projective mismatch addresses")
    gate(observer(left, "projective") == observer(right, "projective"),
         "projective observer collision")
    gate(left.target != right.target, "projective target mismatch")
    gate(left.gap == right.gap == 2, "projective mismatch gap")
    gate(left.bank == right.bank == (15, 6, 1), "projective mismatch bank")
    gate(left.full_portrait == right.full_portrait
         == (15, 6, 1, 8, 3, 10, 13, 4, 7, 14, 9, 0, 11, 2, 5, 12),
         "projective mismatch portrait")
    gate(left.odd_shadow == right.odd_shadow == (5, 7, 5, 15),
         "projective mismatch odd shadow")
    gate(left.carry == right.carry == ("visible", 1, 0),
         "projective mismatch carry")
    gate(left.projective == right.projective
         == ((8, 1, 5), (1, 3, 15), (2, 5, 11)),
         "projective mismatch history")
    gate((left.center, right.center) == (1, 1), "projective centers agree")
    gate((left.next_bank, right.next_bank) == ((7, 10, 13), (11, 14, 1)),
         "projective next banks differ")
    gate((left.chains[0], right.chains[0]) == ((1, 54, 29), (1, 22, 45)),
         "phase base chain repairs mismatch")
    gate(observer(left, "exact") != observer(right, "exact"),
         "exact faithful states remain distinct")
    gate(defect_by_scale[9][431] < defect_by_scale[9][439],
         "THM3824 exact defect separates mismatch phases")

    carry = analyses["carry"]
    gate(carry["first_cross"] is not None, "carry cross-scale mismatch exists")
    later, earlier, _key, cross_left, cross_right = carry["first_cross"]
    gate((earlier, later) == (20, 574), "first carry cross-scale mismatch")
    gate((cross_left.scale, cross_left.phase, cross_right.scale, cross_right.phase)
         == (4, 4, 9, 62), "carry cross-scale addresses")
    gate(observer(cross_left, "carry") == observer(cross_right, "carry"),
         "carry observer collision")
    gate(cross_left.target != cross_right.target, "carry target mismatch")
    gate(cross_left.projective != cross_right.projective,
         "projective history repairs first cross-scale mismatch")
    gate(projective["cross"] == 0, "no cross-scale mismatch after projective")

    # Chronological distinction between a state repeat and a target mismatch.
    projective_collisions = tuple(
        tuple(record.n for record in sorted(fibre, key=lambda item: item.n))
        for fibre in sorted(projective["collision_fibres"],
                            key=lambda rows_: tuple(r.n for r in rows_))
    )
    projective_collision_details = tuple(
        tuple((record.n, record.scale, record.phase, record.target)
              for record in sorted(fibre, key=lambda item: item.n))
        for fibre in sorted(projective["collision_fibres"],
                            key=lambda rows_: tuple(r.n for r in rows_))
    )
    first_repeat = projective["first_repeat"]
    gate(first_repeat is not None, "projective observer has a repeated state")
    first_repeat_pair = (first_repeat[1], first_repeat[0])
    gate(first_repeat_pair == (128, 132), "literal first projective state collision")
    gate(tuple(len({entry[3] for entry in fibre})
               for fibre in projective_collision_details) == (1, 1, 1, 2),
         "first three projective collisions are harmless; fourth mismatches")

    # Exact audit of the phase-owned base correction.
    nonzero_extensions = []
    naive_mismatches = []
    accidental_naive_equal = []
    for record in records:
        naive_routes = route(record.owner, record.gap, 0)
        naive_bank = bank_from_routes(naive_routes, record.gap)
        if record.extension:
            nonzero_extensions.append(record)
            if naive_bank == record.next_bank:
                accidental_naive_equal.append(record.n)
        if naive_bank != record.next_bank:
            naive_mismatches.append((record, naive_bank))
    gate(nonzero_extensions and nonzero_extensions[0].n == 6,
         "n=6 is first nonzero physical phase extension")
    gate(naive_mismatches and naive_mismatches[0][0].n == 6,
         "n=6 is first zero-base next-bank failure")
    gate(len(nonzero_extensions) == 684 and len(naive_mismatches) == 676,
         "zero-base census")
    gate(len(accidental_naive_equal) == 8,
         "eight nonzero phases have accidental selected-bank equality")
    six, six_naive = naive_mismatches[0]
    gate((six.scale, six.phase, six.extension) == (2, 2, 1),
         "n=6 phase-tail address")
    gate(six_naive == (11, 12, 9), "n=6 naive zero-base bank")
    gate(six.next_bank == (15, 8, 5), "n=6 physical next bank")
    gate(six.chains == ((1, 16, 61), (5, 4, 33), (9, 24, 21)),
         "n=6 physical routed chains")
    gate(bank_from_routes(six.chains, six.gap) == six.next_bank,
         "n=6 routed recovery")

    carry_one_records = [record for record in records
                         if len(record.carry) == 3 and record.carry[2] == 1]
    gate(carry_one_records and carry_one_records[0].n == 6,
         "n=6 is independently first ordinary-carry correction")

    # Hash every exact record through an independent streaming serialization.
    record_stream = sha256()
    for record in records:
        exact_payload = (
            record.n, record.scale, record.phase, record.gap, record.extension,
            record.bank, record.full_portrait, record.exact_odds,
            record.carry, record.carry_terms, record.exact_projective,
            record.chains, record.target, record.owner_digest,
        )
        record_stream.update(
            (json.dumps(exact_payload, separators=(",", ":")) + "\n").encode("ascii")
        )

    semantic_record = {
        "universe": [1, 1023],
        "valuations": list(valuations),
        "gaps": list(gaps),
        "census": [list(row) for row in census],
        "ranges": [list(row) for row in ranges],
        "projective_collisions": [list(row) for row in projective_collisions],
        "projective_collision_details": [
            [[entry[0], entry[1], entry[2], list(entry[3])]
             for entry in fibre]
            for fibre in projective_collision_details
        ],
        "projective_first_repeat": list(first_repeat_pair),
        "projective_first_mismatch": [943, 951],
        "carry_first_cross": [20, 574],
        "nonzero_extension_count": len(nonzero_extensions),
        "naive_bank_mismatch_count": len(naive_mismatches),
        "accidental_naive_equal": accidental_naive_equal,
        "first_zero_base_failure": [6, 2, 2, 1],
        "first_carry_one": carry_one_records[0].n,
        "record_stream_sha256": record_stream.hexdigest(),
        "defect_431": [defect_by_scale[9][431].bit_length(),
                       defect_by_scale[9][431] % 65536,
                       integer_digest(defect_by_scale[9][431])],
        "defect_439": [defect_by_scale[9][439].bit_length(),
                       defect_by_scale[9][439] % 65536,
                       integer_digest(defect_by_scale[9][439])],
    }
    semantic_digest = sha256(
        json.dumps(semantic_record, sort_keys=True,
                   separators=(",", ":")).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        gate(semantic_digest == EXPECTED_SEMANTIC_SHA256, "semantic freeze")

    print("RULE30_CROSS_SCALE_FIRST_MISMATCH_INDEPENDENT_AUDIT_20260824")
    print("status=THM-4006_INDEPENDENT_AUDIT_PASS;NO_ALL_SCALE_CLAIM;NO_RULE30_PRIZE")
    print("universe=n_1_through_1023_exactly_once;highest_shell_n=2^m+t;0<=m<=9;0<=t<2^m")
    print(f"valuations={valuations};gaps={gaps}")
    print("same_scale_defect_ranges=" + repr(tuple(ranges)).replace(" ", ""))
    print(f"same_scale_strict_checks={strict_control_gates}")
    print("quotient_census=" + repr(census).replace(" ", ""))
    print(f"projective_collision_fibres={projective_collisions}")
    print("projective_collision_details=" + repr(projective_collision_details).replace(" ", ""))
    print(f"projective_first_state_repeat={first_repeat_pair}")
    print("projective_first_target_mismatch=(943,951);same_scale=(9,431)/(9,439)")
    print("projective_mismatch_target=((1,(7,10,13)),(1,(11,14,1)))")
    print("carry_first_cross_scale_target_mismatch=(20,574);projective_repairs_it=True")
    print(f"zero_base_counts=(nonzero_extensions={len(nonzero_extensions)},naive_bank_mismatches={len(naive_mismatches)})")
    print(f"zero_base_accidental_selected_bank_equal={tuple(accidental_naive_equal)}")
    print("zero_base_first_failure=(n=6,m=2,t=2,extension=1);naive=(11,12,9);physical=(15,8,5)")
    print("zero_base_scope=routed_chain_coordinates_and_recovered_next_bank_only;projective_shadow_carry_and_center_unchanged")
    print("n6_coincidence=also_first_ordinary_carry_one_but_not_the_same_correction")
    print("n6_physical_chains=((1,16,61),(5,4,33),(9,24,21))")
    print("record_stream_sha256=" + record_stream.hexdigest())
    print("semantic_sha256=" + semantic_digest)
    print(f"gates={GATES}")
    print("RESULT=PASS;FIRST_TARGET_MISMATCH=943/951;FINITE_BOUNDARY=1023")


if __name__ == "__main__":
    main()
