#!/usr/bin/env python3
"""Exact companion for THM-2884.

The script couples the endpoint E3/complement bit to the XOR of the two
outer semantic danger bits.  It reconstructs the 449-sheet THM-2835 horn,
adds the fully bare outer word, verifies a four-state V4-shaped diagonal,
and finds the 20 physical cells on which I,J0,J3,J11,J7 retain every
other native factor.  It then evaluates the actual two-origin endpoint
coefficient over both finite-field embeddings.

The positive result is a lawful correlated Boolean carrier and a
canonical stepped-only q3 section.  The negative result is that q11 and
q7 are origin-even on this carrier, so their signed endpoint coefficients
still cancel.  No row exclusion or LRC(14) conclusion is asserted.
"""

from __future__ import annotations

from bisect import bisect_right
from hashlib import sha256
from itertools import permutations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
sys.path.insert(0, str(COMP))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINNED = {
    COMP / "lrc14_q11_semantic_word_horn_thm2835.py":
        "207dd94f086338ae1e80b7d7196f99bf41e795893d13b6d48e4e7d516af03523",
    RESULTS / "lrc14_q11_semantic_word_horn_thm2835.out":
        "1ebe0cbaf7d4ef13defed0bdb5b37df1364880acdbfc6139b243ab9df65f6bf6",
    COMP / "lrc14_event_twisted_all_q_coefficient_carry_lift_thm2882.py":
        "3ed346e0c631b34bd61f0c4d27d7f161e8d35b70decfb95f5207c5f57893d005",
    RESULTS / "lrc14_event_twisted_all_q_coefficient_carry_lift_thm2882.out":
        "0faa0a24f6ba8b6c88b6bbfc4f225e38667097b1a937d977741453499884901c",
}
for path, expected in PINNED.items():
    require(lf_hash(path) == expected, f"pinned dependency changed: {path.name}")


import lrc14_q11_semantic_word_horn_thm2835 as horn
import lrc14_event_twisted_all_q_coefficient_carry_lift_thm2882 as carry


allocation = horn.allocation
word_atlas = horn.atlas
endpoint_atlas = carry.atlas
endpoint_base = carry.endpoint_base
endpoint = carry.endpoint
P = 13
DEPTH = P**5
WORDS = ("Q0", "QA", "QB", "QAB")
WORD_BITS = {
    "Q0": (0, 0),
    "QA": (1, 0),
    "QB": (0, 1),
    "QAB": (1, 1),
}
H = {
    (e, u, v)
    for e in range(2)
    for u in range(2)
    for v in range(2)
    if e ^ u ^ v == 0
}


def build_empty_word(base_module):
    """Build guard-good, units/owner/outer-slots-all-safe Q_empty."""
    pattern = dict(word_atlas.PAT_QAB)
    pattern[7] = "out"
    pattern[8] = "out"
    intervals = [(0, word_atlas.T_DEN)]
    for index, mode in pattern.items():
        if mode == "gout":
            intervals = base_module.subtract_comb(
                intervals, word_atlas.W[index], 91, -13, 13
            )
    for _speed, index in sorted(
        (word_atlas.W[index], index)
        for index, mode in pattern.items()
        if mode == "out"
    ):
        intervals = base_module.subtract_comb(
            intervals, word_atlas.W[index], 182, -13, 13
        )
    return tuple(intervals)


def flag_on_labels(interval, offset, period, pieces, starts, labels):
    """Whole-cylinder flags on a specified ancestry list."""
    flags = []
    for ancestry in labels:
        low = interval[0] + offset + ancestry * period
        high = interval[1] + offset + ancestry * period
        midpoint = (low + high) // 2
        index = bisect_right(starts, midpoint // DEPTH) - 1
        midpoint_hit = (
            index >= 0
            and pieces[index][0] * DEPTH <= midpoint
            < pieces[index][1] * DEPTH
        )
        whole_hit = (
            midpoint_hit
            and pieces[index][0] * DEPTH <= low
            and high <= pieces[index][1] * DEPTH
        )
        require(
            midpoint_hit == whole_hit,
            "midpoint and whole-cylinder membership differ",
        )
        flags.append(bool(whole_hit))
    return tuple(flags)


def xor(left, right):
    return tuple(a ^ b for a, b in zip(left, right))


def permutation_power(permutation, exponent):
    result = tuple(range(len(permutation)))
    base = permutation
    while exponent:
        if exponent & 1:
            result = tuple(base[result[index]] for index in range(len(base)))
        base = tuple(base[base[index]] for index in range(len(base)))
        exponent //= 2
    return result


def main() -> None:
    (
        _module, full_module, _details, e3, clocks, _q_pairs
    ) = allocation.build_geometry()
    period = full_module.T
    unit = period // P
    interval = allocation.ATOM_INTERVAL
    target = tuple(value + allocation.physical.SHIFT for value in interval)
    target_atom = ((*target, 1),)
    old = allocation.physical.relative.lift.m.core.old
    origins = (endpoint_atlas.ORIGIN, endpoint_atlas.STEPPED)

    require(
        period == word_atlas.T_DEN == endpoint_base.T
        and old.T == period
        and old.rail.DEPTH == DEPTH,
        "shared physical scale changed",
    )

    # Reconstruct the three inherited nonempty words and the missing fully
    # bare word Q0.  All four are honest Boolean combinations of the same
    # nine row factors.
    patterns = {
        "QA": word_atlas.PAT_QA,
        "QB": word_atlas.PAT_QB,
        "QAB": word_atlas.PAT_QAB,
    }
    pieces = {
        word: tuple(old.base.build_set(pattern, old.base.ZELL))
        for word, pattern in patterns.items()
    }
    pieces["Q0"] = build_empty_word(old.base)
    starts = {
        word: tuple(left for left, _right in pieces[word])
        for word in WORDS
    }
    require(
        tuple(WORD_BITS[word] for word in ("Q0", "QA", "QB", "QAB"))
        == ((0, 0), (1, 0), (0, 1), (1, 1))
        and len(pieces["Q0"]) == 132262
        and sum(right - left for left, right in pieces["Q0"])
        == 45425355597700,
        "four-word Boolean atlas changed",
    )

    source_qb = horn.word_flags(
        interval, 0, period, pieces["QB"], starts["QB"]
    )
    target_q11_qa = horn.word_flags(
        interval,
        allocation.physical.SHIFT + 11 * unit,
        period,
        pieces["QA"],
        starts["QA"],
    )
    labels = tuple(
        ancestry
        for ancestry in range(DEPTH)
        if source_qb[ancestry] and target_q11_qa[ancestry]
    )
    require(
        len(labels) == 449
        and labels[0] == 59306
        and labels[-1] == 311961
        and {ancestry % 169 for ancestry in labels} == {156},
        "THM-2835 horn labels changed",
    )

    # Determine the complete word, not merely the two outer bits, on every
    # horn label at every q.  The omitted q=3 state is exactly Q0.
    source_word_flags = {
        word: flag_on_labels(
            interval, 0, period, pieces[word], starts[word], labels
        )
        for word in WORDS
    }
    target_word_flags = {
        (word, q): flag_on_labels(
            interval,
            allocation.physical.SHIFT + q * unit,
            period,
            pieces[word],
            starts[word],
            labels,
        )
        for word in WORDS
        for q in range(P)
    }
    source_words = tuple(
        tuple(
            word for word in WORDS
            if source_word_flags[word][index]
        )
        for index in range(len(labels))
    )
    target_words_by_q = tuple(
        tuple(
            tuple(
                word for word in WORDS
                if target_word_flags[word, q][index]
            )
            for index in range(len(labels))
        )
        for q in range(P)
    )
    expected_word_by_q = (
        "QB", "Q0", "Q0", "Q0", "Q0", "Q0", "QA",
        "QAB", "QA", "QA", "QA", "QA", "QA",
    )
    require(
        set(source_words) == {("QB",)}
        and all(
            set(target_words_by_q[q]) == {(expected_word_by_q[q],)}
            for q in range(P)
        ),
        "all-q four-word horn atlas changed",
    )

    q3_empty_global = horn.word_flags(
        interval,
        allocation.physical.SHIFT + 3 * unit,
        period,
        pieces["Q0"],
        starts["Q0"],
    )
    q7_qab_plus_one = flag_on_labels(
        interval,
        allocation.physical.SHIFT + 7 * unit,
        period,
        pieces["QAB"],
        starts["QAB"],
        tuple((ancestry + 1) % DEPTH for ancestry in labels),
    )
    require(
        sum(q3_empty_global) == 66104
        and all(q3_empty_global[ancestry] for ancestry in labels)
        and all(q7_qab_plus_one),
        "Q0 completion or one-borrow q7 filler changed",
    )

    # The endpoint E3 bit is total on every selected atom at both origins.
    carriers = {"I": ((*interval, 1),)}
    carriers.update({
        q: allocation.physical.overlap.shift_weighted(
            target_atom, q * unit
        )
        for q in range(P)
    })
    origin_blocks = {}
    macro_truth = {}
    origin_truth_piece = {}
    for origin_index, address in enumerate(origins):
        ell = endpoint_base.REPS[address]
        e3_block = tuple(endpoint.build_set(endpoint_base.PAT_E3, ell))
        not_e3_block = carry.complement(e3_block, endpoint_base.T)
        origin_blocks[origin_index] = {
            1: e3_block,
            0: not_e3_block,
        }
        for key, carrier in carriers.items():
            restricted = {
                bit: carry.restricted_piece(carrier, block)
                for bit, block in origin_blocks[origin_index].items()
            }
            containing = tuple(
                bit for bit, piece in restricted.items()
                if piece == carrier
            )
            require(
                len(containing) == 1
                and not restricted[1 ^ containing[0]],
                "endpoint E3/complement bit ceased to be total",
            )
            macro_truth[origin_index, key] = containing[0]
            for bit, piece in restricted.items():
                origin_truth_piece[origin_index, key, bit] = piece

    zero_truth = tuple(macro_truth[0, q] for q in range(P))
    stepped_truth = tuple(macro_truth[1, q] for q in range(P))
    require(
        macro_truth[0, "I"] == macro_truth[1, "I"] == 1
        and tuple(q for q, bit in enumerate(zero_truth) if bit)
        == (0, 3, 11)
        and tuple(q for q, bit in enumerate(stepped_truth) if bit)
        == (0, 11),
        "two-origin macro truth atlas changed",
    )

    source_state = {
        origin_index:
            (macro_truth[origin_index, "I"], *WORD_BITS["QB"])
        for origin_index in range(2)
    }
    states = {
        (origin_index, q):
            (
                macro_truth[origin_index, q],
                *WORD_BITS[expected_word_by_q[q]],
            )
        for origin_index in range(2)
        for q in range(P)
    }
    expected_zero_states = (
        (1, 0, 1), (0, 0, 0), (0, 0, 0), (1, 0, 0),
        (0, 0, 0), (0, 0, 0), (0, 1, 0), (0, 1, 1),
        (0, 1, 0), (0, 1, 0), (0, 1, 0), (1, 1, 0),
        (0, 1, 0),
    )
    expected_stepped_states = list(expected_zero_states)
    expected_stepped_states[3] = (0, 0, 0)
    expected_stepped_states = tuple(expected_stepped_states)
    h_support = {
        origin_index: tuple(
            q for q in range(P) if states[origin_index, q] in H
        )
        for origin_index in range(2)
    }
    require(
        H == {(0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0)}
        and source_state == {0: (1, 0, 1), 1: (1, 0, 1)}
        and tuple(states[0, q] for q in range(P)) == expected_zero_states
        and tuple(states[1, q] for q in range(P))
        == expected_stepped_states
        and h_support[0] == (0, 1, 2, 4, 5, 7, 11)
        and h_support[1] == (0, 1, 2, 3, 4, 5, 7, 11),
        "macro-semantic diagonal support changed",
    )

    # The q3/q11/q7 target triangle plus the fixed source runs through all
    # four elements of H.  The edge labels are flat under XOR.
    h0 = states[1, 3]
    h_qa = states[1, 11]
    h_qab = states[1, 7]
    h_qb = source_state[1]
    edge_3_11 = xor(h0, h_qa)
    edge_11_7 = xor(h_qa, h_qab)
    edge_3_7 = xor(h0, h_qab)
    require(
        {h0, h_qa, h_qab, h_qb} == H
        and (h0, h_qa, h_qab, h_qb)
        == ((0, 0, 0), (1, 1, 0), (0, 1, 1), (1, 0, 1))
        and (edge_3_11, edge_11_7, edge_3_7)
        == (h_qa, h_qb, h_qab)
        and xor(edge_3_11, edge_11_7) == edge_3_7,
        "V4-shaped horn edge law changed",
    )

    # On the four seam checkpoints, the affine parity required by the
    # positive THM-2882 selector is the nonlinear indicator of H\\{0}.
    # This is a useful geometric description, but it is only seam-local:
    # q=1,2,4,5 give exact counterexamples to an all-q extension.
    seam = (0, 3, 11, 7)
    desired_selector_parity = tuple(1 ^ int(q == 3) for q in range(P))
    seam_nonzero_indicator = tuple(
        int(states[1, q] != (0, 0, 0)) for q in seam
    )
    all_q_nonzero_indicator = tuple(
        int(
            states[1, q] in H
            and states[1, q] != (0, 0, 0)
        )
        for q in range(P)
    )
    require(
        seam_nonzero_indicator
        == tuple(desired_selector_parity[q] for q in seam)
        == (1, 0, 1, 1)
        and tuple(
            q for q in range(P)
            if all_q_nonzero_indicator[q]
            != desired_selector_parity[q]
        ) == (1, 2, 4, 5, 6, 8, 9, 10, 12)
        and len(H - {(0, 0, 0)}) == 3,
        "nonzero-V4 seam selector boundary changed",
    )

    # Reconstruct the native six-factor cell bank efficiently.  The E3
    # factor is the only q7 failure on the final 20-cell intersection.
    physical_intervals = {
        "I": interval,
        "J": target,
        3: horn.circular_shift(target, 3 * unit, period),
        7: horn.circular_shift(target, 7 * unit, period),
        11: horn.circular_shift(target, 11 * unit, period),
    }
    e3_bit = {
        key: allocation.contained(value, e3)
        for key, value in physical_intervals.items()
    }
    clock_bit = {
        (key, clock): allocation.contained(value, clocks[clock])
        for key, value in physical_intervals.items()
        for clock in range(7)
    }
    q1_bit = {
        (key, s): horn.safe_comb_contains(
            value, full_module, full_module.W[1], 182,
            -14 * s - 13, -14 * s + 13,
        )
        for key, value in physical_intervals.items()
        for s in allocation.COMMON_S
    }
    q2_bit = {
        (key, t): horn.safe_comb_contains(
            value, full_module, full_module.W[2], 182,
            -14 * t - 13, -14 * t + 13,
        )
        for key, value in physical_intervals.items()
        for t in allocation.COMMON_T
    }
    c2_bit = {
        (key, s): horn.safe_comb_contains(
            value, full_module, full_module.C2, 182,
            14 * s - 13, 14 * s + 13,
        )
        for key, value in physical_intervals.items()
        for s in allocation.COMMON_S
    }
    c3_bit = {
        (key, t): horn.safe_comb_contains(
            value, full_module, full_module.C3, 182,
            14 * t - 13, 14 * t + 13,
        )
        for key, value in physical_intervals.items()
        for t in allocation.COMMON_T
    }

    def other_five(keys, s, t, clock):
        return all(
            clock_bit[key, clock]
            and q1_bit[key, s]
            and q2_bit[key, t]
            and c2_bit[key, s]
            and c3_bit[key, t]
            for key in keys
        )

    cells = tuple(
        (s, t, clock)
        for s in allocation.COMMON_S
        for t in allocation.COMMON_T
        for clock in range(7)
    )
    full_q3 = {
        cell for cell in cells
        if all(e3_bit[key] for key in ("I", "J", 3))
        and other_five(("I", "J", 3), *cell)
    }
    full_q11 = {
        cell for cell in cells
        if all(e3_bit[key] for key in ("I", "J", 11))
        and other_five(("I", "J", 11), *cell)
    }
    q7_e3_only = {
        cell for cell in cells
        if e3_bit["I"] and e3_bit["J"] and not e3_bit[7]
        and other_five(("I", "J", 7), *cell)
    }
    diagonal_cells = full_q3 & full_q11 & q7_e3_only
    expected_diagonal_cells = {
        (s, t, 1)
        for s in (0, 3, 8, 9, 12)
        for t in (5, 6, 9, 10)
    }
    require(
        e3_bit == {"I": True, "J": True, 3: True, 7: False, 11: True}
        and len(full_q3) == 56
        and len(full_q11) == 63
        and len(full_q3 & full_q11) == 42
        and len(q7_e3_only) == 49
        and len(full_q11 & q7_e3_only) == 35
        and diagonal_cells == expected_diagonal_cells,
        "20-cell diagonal physical packet changed",
    )

    # The diagonal selector chooses E3 when the outer parity is one and
    # complement when it is zero.  This is one origin-independent rule.
    # Its pointwise signed endpoint coefficient is supported only at q3.
    parity_by_q = tuple(
        WORD_BITS[word][0] ^ WORD_BITS[word][1]
        for word in expected_word_by_q
    )
    diagonal_full = {
        (origin_index, q):
            states[origin_index, q] in H
        for origin_index in range(2)
        for q in range(P)
    }
    require(
        all(
            diagonal_full[origin_index, q]
            == (
                parity_by_q[q] == macro_truth[origin_index, q]
            )
            for origin_index in range(2)
            for q in range(P)
        ),
        "diagonal selector is not equality of the two Boolean bits",
    )

    actual_samples = tuple(
        sample
        for value in endpoint_atlas.UNIT_SECTIONS
        for sample in (value, value + 1)
    )
    field_rows = []
    for embedding in endpoint.MODS:
        prime, _root = embedding
        source_value = endpoint_atlas.COMMON_SOURCE[prime]
        signed_supports = []
        q3_currents = []
        for sample in actual_samples:
            frequency = -(12 + 26 * sample)
            signed = []
            for q in range(P):
                selected_bit = parity_by_q[q]
                values = tuple(
                    allocation.endpoint_sum(
                        origin_truth_piece[
                            origin_index, q, selected_bit
                        ],
                        frequency,
                        embedding,
                    )
                    for origin_index in range(2)
                )
                base_value = allocation.endpoint_sum(
                    carriers[q], frequency, embedding
                )
                require(
                    base_value != 0
                    and all(
                        values[origin_index]
                        == (
                            base_value
                            if diagonal_full[origin_index, q]
                            else 0
                        )
                        for origin_index in range(2)
                    ),
                    "diagonal endpoint restriction/value changed",
                )
                signed.append((values[0] - values[1]) % prime)
            signed_supports.append(
                tuple(q for q, value in enumerate(signed) if value)
            )
            require(
                signed[3]
                == -allocation.endpoint_sum(
                    carriers[3], frequency, embedding
                ) % prime,
                "stepped-only q3 orientation changed",
            )
            q3_current = source_value * signed[3] % prime
            require(q3_current != 0, "full q3 current vanished")
            q3_currents.append(q3_current)
        require(
            tuple(signed_supports) == ((3,),) * len(actual_samples),
            "diagonal signed support is no longer delta_q3",
        )
        field_rows.append((
            prime,
            len(actual_samples),
            tuple(sorted(set(signed_supports))),
            len(set(q3_currents)),
        ))

    # A local V4 labelling does not become a C13 action: a generator of a
    # C13 action on four points would be a permutation pi with pi^13=1.
    # Exhausting S4 leaves only the identity.
    identity = tuple(range(4))
    c13_actions = tuple(
        permutation
        for permutation in permutations(range(4))
        if permutation_power(permutation, P) == identity
    )
    require(
        c13_actions == (identity,),
        "a nontrivial C13 action on the four diagonal states appeared",
    )

    print(
        "scale="
        f"(p={P},depth={DEPTH},T={period},U={unit}); "
        f"horn_labels={len(labels)}; range=({labels[0]},{labels[-1]}); "
        "residue_mod169=156"
    )
    print(
        "four_words="
        "Q0=(0,0),QA=(1,0),QB=(0,1),QAB=(1,1); "
        f"Q0_intervals={len(pieces['Q0'])}; "
        f"Q0_measure_numerator={sum(r-l for l,r in pieces['Q0'])}"
    )
    print(
        f"horn_word_by_q={expected_word_by_q}; "
        "source_word=QB; q7_zero_and_one_borrow=QAB; "
        f"q3_Q0_global={sum(q3_empty_global)}; q3_Q0_horn=449/449"
    )
    print(
        f"zero_macro_E3={tuple(q for q,b in enumerate(zero_truth) if b)}; "
        f"stepped_macro_E3={tuple(q for q,b in enumerate(stepped_truth) if b)}"
    )
    print(
        f"H=ker(e_XOR_u_XOR_v)={tuple(sorted(H))}; "
        f"source={source_state[1]}; q3_step={h0}; "
        f"q11={h_qa}; q7={h_qab}; "
        f"edges=(3to11={edge_3_11},11to7={edge_11_7},"
        f"3to7={edge_3_7})"
    )
    print(
        f"H_support_zero={h_support[0]}; "
        f"H_support_stepped={h_support[1]}; "
        "origin_parity=delta_q3"
    )
    print(
        f"seam_q={seam}; H_nonzero_indicator={seam_nonzero_indicator}; "
        "desired_positive_selector_parity=(1,0,1,1); "
        "H_nonzero_is_size3_not_subgroup_coset_or_character_fibre; "
        "all_q_counterexamples=(1,2,4,5,6,8,9,10,12)"
    )
    print(
        f"cells=(q3={len(full_q3)},q11={len(full_q11)},"
        f"q3capq11={len(full_q3 & full_q11)},"
        f"q7_E3_only={len(q7_e3_only)},"
        f"q11capq7_E3_only={len(full_q11 & q7_e3_only)},"
        f"diagonal={len(diagonal_cells)}); "
        f"diagonal_cells={tuple(sorted(diagonal_cells))}"
    )
    print(
        f"complete_stepped_packets=20*449={20 * 449}; "
        f"two_q7_fillers={2 * 20 * 449}; "
        "intervals=(I,J0,J3,J11,J7); "
        "other_native_factors=(clock,q1,q2,c2,c3)"
    )
    print(
        f"field_rows={tuple(field_rows)}; "
        "signed_diagonal_support=q3_only; "
        "q3_orientation=zero_minus_stepped=-base_nonzero; "
        "q11_q7=origin_even_cancel"
    )
    print(
        f"C13_actions_on_H={len(c13_actions)}_identity_only; "
        "V4_edge_triangle=flat_but_not_global_C13_action; "
        "carry_boundary=q7_a_and_a+1_have_same_QAB_H_state"
    )
    print(
        "CONCLUSION: the correlated physical selector "
        "H=ker(E3_bit XOR outer_slot7 XOR outer_slot8) repairs the "
        "THM-2835 macro/word mismatch.  On every one of the 449 horn "
        "sheets, the stepped-origin q3 target is the fully bare state "
        "000, while q11, q7, and the fixed source are the other three "
        "V4 states.  Twenty native cells retain I,J0,J3,J11,J7 and all "
        "five remaining factors, giving 8980 complete stepped packets.  "
        "The same origin-independent rule is a canonical origin-odd q3 "
        "selector, but q11 and q7 remain origin-even and cancel in the "
        "signed coefficient.  The V4 edge labels do not extend to a "
        "C13 action and remain blind to the q7 ancestry carry.  Thus this "
        "is a physical carrier/clutch candidate, not a row exclusion or "
        "an LRC(14) proof."
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
