#!/usr/bin/env python3
"""Exact 449-sheet macro/semantic diagonal and V4-boundary audit.

This scratch companion reconstructs the THM-2835 449-sheet horn and the
previously omitted all-out semantic word Qempty.  It then separates:

* the algebraic subgroup H = ker(e XOR u XOR v);
* its sheetwise realization on the four horn states; and
* the stronger, still-unproved assertion that these incidences define a
  physical V4 action on the complete current packet.

It proves no row exclusion or LRC(14).
"""

from __future__ import annotations

from hashlib import sha256
from itertools import product
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
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
    COMP / "lrc14_169_twist_variance_opus_20260726.py":
        "c0aa9c55c3e7d38dc28b4705e58c776a979f17d2874d1b762f96b97d2364e5e9",
    COMP / "lrc14_literal_fixed_sheet_allocation_thm2806.py":
        "311d0d85500f0c65945ebe5913f09d34a16293119c942b42eeaa854fbf85f71e",
}
for path, expected in PINNED.items():
    require(lf_hash(path) == expected, f"pinned dependency changed: {path.name}")


import lrc14_q11_semantic_word_horn_thm2835 as horn


F2 = (0, 1)
P = 13
DEPTH = P**5


def add(left, right):
    return tuple(x ^ y for x, y in zip(left, right))


def chi(state):
    e3, slot7, slot8 = state
    return e3 ^ slot7 ^ slot8


def macro_truth(q, origin):
    support = (0, 3, 11) if origin == 0 else (0, 11)
    return int(q in support)


def build_all_out(base, pattern, ell):
    """Build a semantic cell with no positive `in` anchor.

    `base.build_set` starts from the least-speed positive comb and therefore
    intentionally does not accept an all-out pattern.  Here the ambient
    circle is the positive starting object, after which the same guard and
    ordinary-comb subtractions are performed.
    """
    require(
        set(pattern.values()) <= {"out", "gout"}
        and "gout" in pattern.values(),
        "all-out reconstruction received a positive or guardless pattern",
    )
    intervals = [(0, base.T_DEN)]
    for index, mode in pattern.items():
        if mode == "gout":
            intervals = base.subtract_comb(
                intervals,
                base.W[index],
                91,
                -13 - 7 * ell[index],
                13 - 7 * ell[index],
            )
    for _speed, index in sorted(
        (base.W[index], index)
        for index, mode in pattern.items()
        if mode == "out"
    ):
        intervals = base.subtract_comb(
            intervals,
            base.W[index],
            182,
            -13 - 14 * ell[index],
            13 - 14 * ell[index],
        )
    return tuple(intervals)


def main() -> None:
    _module, full, _details, _e3, _clocks, _q_pairs = (
        horn.allocation.build_geometry()
    )
    period = full.T
    unit = period // P
    atom = horn.allocation.ATOM_INTERVAL
    target_shift = horn.allocation.physical.SHIFT
    old = horn.allocation.physical.relative.lift.m.core.old
    base = old.base
    require(
        period == horn.atlas.T_DEN == old.T == base.T_DEN
        and old.rail.DEPTH == DEPTH,
        "shared physical/ancestry scale changed",
    )

    patterns = {
        "QA": dict(horn.atlas.PAT_QA),
        "QB": dict(horn.atlas.PAT_QB),
        "QAB": dict(horn.atlas.PAT_QAB),
    }
    qempty_pattern = {
        0: "gout",
        1: "out",
        2: "out",
        3: "out",
        4: "out",
        5: "out",
        6: "out",
        7: "out",
        8: "out",
    }
    require(
        tuple(
            (patterns[word][7], patterns[word][8])
            for word in ("QB", "QA", "QAB")
        )
        == (("out", "in"), ("in", "out"), ("in", "in")),
        "outer semantic bit order changed",
    )

    pieces = {
        word: tuple(base.build_set(pattern, base.ZELL))
        for word, pattern in patterns.items()
    }
    pieces["Qempty"] = build_all_out(base, qempty_pattern, base.ZELL)
    starts = {
        word: tuple(left for left, _right in intervals)
        for word, intervals in pieces.items()
    }

    cache = {}

    def flags(word, q=None):
        key = word, q
        if key not in cache:
            offset = 0 if q is None else target_shift + q * unit
            cache[key] = horn.word_flags(
                atom,
                offset,
                period,
                pieces[word],
                starts[word],
            )
        return cache[key]

    source_qb = flags("QB")
    q11_qa = flags("QA", 11)
    labels = tuple(
        ancestry
        for ancestry in range(DEPTH)
        if source_qb[ancestry] and q11_qa[ancestry]
    )
    require(
        len(labels) == 449
        and labels[0] == 59306
        and labels[-1] == 311961
        and {ancestry % 169 for ancestry in labels} == {156},
        "pinned 449-sheet cospan changed",
    )

    q7_qab_plus_one_hits = sum(
        flags("QAB", 7)[(ancestry + 1) % DEPTH]
        for ancestry in labels
    )
    qempty_global_counts = tuple(sum(flags("Qempty", q)) for q in range(P))
    qempty_label_hits = tuple(
        sum(flags("Qempty", q)[ancestry] for ancestry in labels)
        for q in range(1, 6)
    )
    require(
        q7_qab_plus_one_hits == 449
        and qempty_label_hits == (449, 449, 449, 449, 449),
        "449-sheet QAB/Qempty completion changed",
    )

    outer = {
        "Qempty": (0, 0),
        "QB": (0, 1),
        "QA": (1, 0),
        "QAB": (1, 1),
    }
    semantic = {
        "Qempty": (0, *outer["Qempty"]),
        "QB": (1, *outer["QB"]),
        "QA": (1, *outer["QA"]),
        "QAB": (0, *outer["QAB"]),
    }
    ambient = tuple(product(F2, repeat=3))
    diagonal = tuple(state for state in ambient if chi(state) == 0)
    require(
        diagonal == ((0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0))
        and set(semantic.values()) == set(diagonal)
        and semantic["QB"] == add(semantic["QA"], semantic["QAB"])
        and semantic["QA"] == add(semantic["QB"], semantic["QAB"])
        and semantic["QAB"] == add(semantic["QB"], semantic["QA"])
        and all(
            add(left, right) in diagonal
            for left in diagonal for right in diagonal
        ),
        "macro/semantic diagonal stopped being V4",
    )
    nonzero_horn = {
        semantic["QB"], semantic["QA"], semantic["QAB"]
    }
    require(
        nonzero_horn == set(diagonal) - {semantic["Qempty"]}
        and len(nonzero_horn) == 3
        and len(nonzero_horn) not in (1, 2, 4, 8),
        "three horn words stopped being the punctured diagonal",
    )

    # At q1,q2,q4,q5, Qempty is diagonal at both origins.  At q3 the
    # zero origin uses E3 and lies in the opposite coset, while the stepped
    # origin uses complement and supplies the missing identity of H.
    qempty_states = {
        (q, origin): (
            macro_truth(q, origin),
            *outer["Qempty"],
        )
        for q in range(1, 6)
        for origin in F2
    }
    qempty_chi_pairs = tuple(
        (
            chi(qempty_states[q, 0]),
            chi(qempty_states[q, 1]),
        )
        for q in range(1, 6)
    )
    require(
        qempty_chi_pairs
        == ((0, 0), (0, 0), (1, 0), (0, 0), (0, 0))
        and qempty_states[3, 1] == semantic["Qempty"]
        and qempty_states[3, 0] == (1, 0, 0)
        and qempty_states[3, 0] not in diagonal,
        "q3 stopped being the unique Qempty diagonal-coset split",
    )

    # On the four seam checkpoints, the character chi is the exact
    # origin-difference bit delta_q3.  The desired THM-2882 positive
    # selector parity is its affine complement, so the diagonal explains
    # the already-lawful q3 asymmetry but not the missing q11/q7 one.
    checkpoint_word = {
        0: "QB",
        3: "Qempty",
        11: "QA",
        7: "QAB",
    }
    checkpoint_states = {
        q: tuple(
            (macro_truth(q, origin), *outer[word])
            for origin in F2
        )
        for q, word in checkpoint_word.items()
    }
    checkpoint_chi = {
        q: tuple(chi(state) for state in states)
        for q, states in checkpoint_states.items()
    }
    checkpoint_defect_parity = {
        q: pair[0] ^ pair[1]
        for q, pair in checkpoint_chi.items()
    }
    desired_selector_parity = {
        q: 1 ^ int(q == 3)
        for q in checkpoint_word
    }
    diagonal_signed_amplitude = {
        q: (1 ^ pair[0]) - (1 ^ pair[1])
        for q, pair in checkpoint_chi.items()
    }
    require(
        checkpoint_chi
        == {0: (0, 0), 3: (1, 0), 11: (0, 0), 7: (0, 0)}
        and checkpoint_defect_parity
        == {0: 0, 3: 1, 11: 0, 7: 0}
        and all(
            desired_selector_parity[q] == 1 ^ checkpoint_defect_parity[q]
            for q in checkpoint_word
        )
        and diagonal_signed_amplitude
        == {0: 0, 3: -1, 11: 0, 7: 0},
        "diagonal character/selector boundary changed",
    )

    # On the four seam vertices only, the desired support parity is the
    # nonzero indicator of H.  This is the canonical S3-invariant puncture
    # H\\{0}, equivalently outer OR.  It is not a character: two nonzero
    # states can add to a third nonzero state.  The q1..q5 Qempty census is
    # the sharp hostile to extending this mnemonic to an all-q rule.
    seam_nonzero_indicator = {
        q: int(
            (
                macro_truth(q, 1),
                *outer[word],
            )
            != semantic["Qempty"]
        )
        for q, word in checkpoint_word.items()
    }
    q1_to_q5_nonzero_indicator = tuple(
        int(qempty_states[q, 1] != semantic["Qempty"])
        for q in range(1, 6)
    )
    q1_to_q5_desired_parity = tuple(
        1 ^ int(q == 3) for q in range(1, 6)
    )
    require(
        seam_nonzero_indicator == desired_selector_parity
        and q1_to_q5_nonzero_indicator == (0, 0, 0, 0, 0)
        and q1_to_q5_desired_parity == (1, 1, 0, 1, 1)
        and (
            int(semantic["QB"] != semantic["Qempty"])
            ^ int(semantic["QA"] != semantic["Qempty"])
        )
        != int(semantic["QAB"] != semantic["Qempty"]),
        "punctured-H seam mnemonic or all-q hostile changed",
    )

    # The two origins have identical B=(e,u,v) states at q11 and q7.
    # Exhaust every Boolean observable on B: no function of these three
    # coordinates can distinguish the origins there.  Exactly half can
    # distinguish the two q3 states.  An origin bit is therefore the
    # minimal missing coordinate for a q11/q7 signed selector.
    ambient_index = {state: index for index, state in enumerate(ambient)}

    def boolean_value(mask, state):
        return (mask >> ambient_index[state]) & 1

    semantic_boolean_census = {"q3_distinguishing": 0, "q11_or_q7_bad": 0}
    for mask in range(1 << len(ambient)):
        if (
            boolean_value(mask, checkpoint_states[3][0])
            != boolean_value(mask, checkpoint_states[3][1])
        ):
            semantic_boolean_census["q3_distinguishing"] += 1
        if (
            boolean_value(mask, checkpoint_states[11][0])
            != boolean_value(mask, checkpoint_states[11][1])
            or boolean_value(mask, checkpoint_states[7][0])
            != boolean_value(mask, checkpoint_states[7][1])
        ):
            semantic_boolean_census["q11_or_q7_bad"] += 1
    require(
        semantic_boolean_census
        == {"q3_distinguishing": 128, "q11_or_q7_bad": 0},
        "semantic-only Boolean observable census changed",
    )

    tilted_diagonal_by_checkpoint = {
        q: tuple(
            int(chi((macro_truth(q, origin), *outer[word])) == origin)
            for origin in F2
        )
        for q, word in checkpoint_word.items()
    }
    require(
        tilted_diagonal_by_checkpoint
        == {0: (1, 0), 3: (0, 0), 11: (1, 0), 7: (1, 0)},
        "origin-tilted diagonal candidate changed",
    )

    # Projection to the two outer bits identifies H with V4.  The three
    # observed horn words are H minus its identity, not a subgroup or
    # coset.  A free H action would enlarge C169 to 676 states, but the
    # exact census supplies a section/incidence, not the missing packet
    # intertwiners needed to license that action.
    projection = {state: state[1:] for state in diagonal}
    require(
        len(set(projection.values())) == 4
        and all(state[0] == (state[1] ^ state[2]) for state in diagonal)
        and 169 * len(diagonal) == 676,
        "H-to-V4 graph identification changed",
    )

    q11_to_q7_semantic_increment = add(
        semantic["QA"], semantic["QAB"]
    )
    source_to_q11_semantic_increment = add(
        semantic["QB"], semantic["QA"]
    )
    require(
        q11_to_q7_semantic_increment == semantic["QB"]
        and source_to_q11_semantic_increment == semantic["QAB"],
        "semantic edge increments changed",
    )

    print("449-SHEET MACRO/SEMANTIC DIAGONAL AUDIT")
    print(
        f"labels=(count={len(labels)},first={labels[0]},last={labels[-1]},"
        f"residue_mod169={labels[0] % 169}); "
        f"q7_QAB_at_a+1={q7_qab_plus_one_hits}"
    )
    print(
        f"Qempty_global_counts_q0_to_q12={qempty_global_counts}; "
        f"Qempty_hits_on_449_q1_to_q5={qempty_label_hits}"
    )
    print(
        f"H=ker(e_XOR_slot7_XOR_slot8)={diagonal}=graph_of_outer_parity;"
        " H_is_V4"
    )
    print(
        f"semantic_states={semantic}; "
        "QB_XOR_QA=QAB; observed_QB_QA_QAB_is_H_minus_identity_"
        "and_is_neither_subgroup_nor_coset"
    )
    print(
        f"Qempty_chi_pairs_q1_to_q5={qempty_chi_pairs}; "
        "stepped_q3=000_completes_H; zero_q3=100_lies_in_opposite_coset"
    )
    print(
        f"checkpoint_chi_by_origin={checkpoint_chi}; "
        f"defect_parity={checkpoint_defect_parity}; "
        f"diagonal_signed_amplitude={diagonal_signed_amplitude}"
    )
    print(
        f"seam_nonzero_indicator={seam_nonzero_indicator}=desired_parity; "
        f"q1_to_q5_nonzero_hostile={q1_to_q5_nonzero_indicator}_versus_"
        f"desired={q1_to_q5_desired_parity}; "
        "H_minus_zero_is_P1(F2)_and_is_not_a_character_fibre"
    )
    print(
        f"all_256_semantic_Boolean_observables={semantic_boolean_census}; "
        f"origin_tilted_kernel_chi_XOR_o={tilted_diagonal_by_checkpoint}"
    )
    print(
        "selector_relation=diagonal_defect_parity_is_delta_q3; "
        "THM2882_positive_selector_parity_is_1_XOR_defect; "
        "diagonal_explains_q3_origin_asymmetry_but_is_origin-even_and_"
        "signed-cancelling_at_q11_and_q7"
    )
    print(
        f"edge_increments=(source_QB_to_q11_QA:{source_to_q11_semantic_increment},"
        f"q11_QA_to_q7_QAB:{q11_to_q7_semantic_increment}_with_a_to_a+1)"
    )
    print(
        "C169xV4_connection=projection_H_to_(slot7,slot8)_is_an_"
        "isomorphism_and_free_independent_H_would_have_676_states; "
        "current_449-fold_object_is_a_complete_support_section_not_a_"
        "proved_physical_V4_action"
    )
    print(
        "needed_sidecar=an_origin-coupled_affine_translate_or_external_"
        "origin_character_such_as_chi_XOR_o_at_q11/q7_plus_complete-packet_"
        "intertwiners; "
        "mere_H-membership_cannot_prevent_q11/q7_signed_cancellation"
    )
    print(
        "scope=finite-exact_sheetwise_support_and_F2_character_anatomy; "
        "no_lawful_free_H_action_on_the_physical_current_or_LRC14_conclusion"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
