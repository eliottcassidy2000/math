#!/usr/bin/env python3
"""Exact all-height two-phase probe for the LRC(14) d=2 boundary.

This is a theorem-discovery companion, not a canonized theorem dependency.
It checks the complete odd residue universe modulo 112 by two independent
representations and prints the actual four-time loneliness consequence.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations_with_replacement


LAMBDA = Fraction(1, 14)
Y_BANK = (Fraction(1, 14), Fraction(5, 56))
X_BANKS = tuple(tuple((y + j) / 2 for j in range(2)) for y in Y_BANK)
ODD_RESIDUES_112 = tuple(range(1, 112, 2))
H0 = tuple(range(1, 11)) + (12,)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def circle_norm(value):
    value %= 1
    return min(value, 1 - value)


def direct_mask(speed, bank_index):
    return tuple(
        j
        for j, x in enumerate(X_BANKS[bank_index])
        if circle_norm(speed * x) < LAMBDA
    )


def circular_residue(value, modulus):
    value %= modulus
    return min(value, modulus - value)


def modular_mask(speed, bank_index):
    # All four lift times have denominator 112 and numerators 4,60,5,61.
    numerators = ((4, 60), (5, 61))[bank_index]
    return tuple(
        j
        for j, numerator in enumerate(numerators)
        if circular_residue(speed * numerator, 112) < 8
    )


def fully_spoils_pair(alpha, beta, bank_index):
    return set(direct_mask(alpha, bank_index)) | set(
        direct_mask(beta, bank_index)
    ) == {0, 1}


def pack_safe_at_y(pack, y):
    return all(circle_norm(h * y) >= LAMBDA for h in pack)


def row_clearance(pack, exceptions, x):
    speeds = tuple(2 * h for h in pack) + tuple(exceptions)
    return min(circle_norm(v * x) for v in speeds)


def first_bank_witness(pack, exceptions):
    for bank_index, xs in enumerate(X_BANKS):
        for label, x in enumerate(xs):
            clearance = row_clearance(pack, exceptions, x)
            if clearance >= LAMBDA:
                return bank_index, label, x, clearance
    return None


def main():
    transcript = []

    def emit(line):
        transcript.append(line)
        print(line)

    for residue in ODD_RESIDUES_112:
        for bank_index in range(2):
            require(
                direct_mask(residue, bank_index)
                == modular_mask(residue, bank_index),
                f"direct/modular mask mismatch at {residue},{bank_index}",
            )
            require(
                len(direct_mask(residue, bank_index)) <= 1,
                f"one odd speed spoiled both labels at {residue},{bank_index}",
            )

    y1_label0 = tuple(
        r for r in ODD_RESIDUES_112 if direct_mask(r, 0) == (0,)
    )
    y1_label1 = tuple(
        r for r in ODD_RESIDUES_112 if direct_mask(r, 0) == (1,)
    )
    require(
        tuple(r for r in ODD_RESIDUES_112 if r % 28 in (1, 27))
        == y1_label0,
        "the +/-1 mod 28 class was not the first y1 label",
    )
    require(
        tuple(r for r in ODD_RESIDUES_112 if r % 28 in (13, 15))
        == y1_label1,
        "the +/-15 mod 28 class was not the second y1 label",
    )
    require(
        all(direct_mask(r, 1) == () for r in y1_label1),
        "a +/-15 mod 28 speed was dangerous at the second phase",
    )

    residue_pairs = tuple(combinations_with_replacement(ODD_RESIDUES_112, 2))
    cover_y1 = tuple(
        pair for pair in residue_pairs if fully_spoils_pair(*pair, 0)
    )
    cover_y2 = tuple(
        pair for pair in residue_pairs if fully_spoils_pair(*pair, 1)
    )
    cover_both = tuple(pair for pair in cover_y1 if fully_spoils_pair(*pair, 1))
    require(not cover_both, "an odd residue pair spoiled both phase banks")

    allowed_h_residues = tuple(
        h
        for h in range(56)
        if all(pack_safe_at_y((h,), y) for y in Y_BANK)
    )
    forbidden_h_residues = tuple(h for h in range(56) if h not in allowed_h_residues)
    expected_forbidden = (0, 11, 14, 22, 23, 28, 33, 34, 42, 45)
    require(
        forbidden_h_residues == expected_forbidden,
        "pack residue firewall classification changed",
    )
    require(
        all(pack_safe_at_y(H0, y) for y in Y_BANK),
        "the canonical eleven-pack misses one phase bank",
    )

    controls = ((1, 11), (1, 15), (3, 5), (9, 15), (111, 113))
    control_lines = []
    for exceptions in controls:
        witness = first_bank_witness(H0, exceptions)
        require(witness is not None, f"missing witness for {exceptions}")
        bank_index, label, x, clearance = witness
        control_lines.append(
            f"exceptions={exceptions}:bank={bank_index + 1},label={label},"
            f"x={x},clearance={clearance}"
        )

    # These are hostiles only to deleting one of the two pack-side residue
    # hypotheses. They are not Lonely Runner counterexamples.
    hostile_y2_pack = first_bank_witness((11,), (1, 15))
    hostile_y1_pack = first_bank_witness((14,), (1, 11))
    hostile_even_exceptions = first_bank_witness(H0, (22, 28))
    require(hostile_y2_pack is None, "h=11 did not break the four-time bank")
    require(hostile_y1_pack is None, "h=14 did not break the four-time bank")
    require(
        hostile_even_exceptions is None,
        "the even-exception parity hostile did not break the four-time bank",
    )

    emit("LRC14 D2 TWO-PHASE ALL-HEIGHT EXACT PROBE")
    emit("status=PROOF-COMPLETE SIGNAL + FINITE-EXACT RESIDUE AUDIT; NOT CANON")
    emit("endpoint_convention=pack_closed(>=1/14);danger_open(<1/14)")
    emit(
        "phase_bank_y=(1/14,5/56);lift_time_banks="
        "((1/28,15/28),(5/112,61/112))"
    )
    emit(
        f"odd_residue_universe_mod112={len(ODD_RESIDUES_112)};"
        f"unordered_with_repetition_pairs={len(residue_pairs)};"
        f"direct_modular_masks_match=True"
    )
    emit(
        f"pairs_spoiling_y1={len(cover_y1)};pairs_spoiling_y2={len(cover_y2)};"
        f"pairs_spoiling_both={len(cover_both)}"
    )
    emit(f"y1_label0_residues_mod112={y1_label0}")
    emit(f"y1_label1_residues_mod112={y1_label1}")
    emit("y1_label1_class=+/-15_mod28;its_y2_mask_is_empty=True")
    emit(
        f"allowed_pack_residues_mod56={allowed_h_residues};"
        f"forbidden_pack_residues_mod56={forbidden_h_residues}"
    )
    emit(f"canonical_pack_H0={H0};safe_at_both_y=True")
    for line in control_lines:
        emit(f"positive_control_{line}")
    emit(
        "hypothesis_hostile_y2_pack=H=(11,),exceptions=(1,15),"
        "all_four_bank_times_fail=True"
    )
    emit(
        "hypothesis_hostile_y1_pack=H=(14,),exceptions=(1,11),"
        "all_four_bank_times_fail=True"
    )
    emit(
        "parity_hostile=H=canonical_H0,exceptions=(22,28),"
        "all_four_bank_times_fail=True"
    )
    emit(
        "consequence=for_every_pack_H_in_allowed_mod56_and_every_two_odd_"
        "exceptions_some_displayed_lift_time_has_full_row_clearance>=1/14"
    )
    digest = sha256(("\n".join(transcript) + "\n").encode()).hexdigest()
    print(f"semantic_digest={digest}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
