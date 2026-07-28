#!/usr/bin/env python3
"""Exact companion for the THM-2835 q11 word/Bockstein horn.

It independently reconstructs the complete q-labelled outer-word graph,
isolates its 449 QB(source)->QA(target) whole-cylinder labels, and tests
every adjacent ancestry orientation at q=0,7,11.  It then separates the
resulting coefficient nerve from a lawful translation of the complete
six-factor packet and checks the rational C13 representation boundary.
"""

from __future__ import annotations

from bisect import bisect_right
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import importlib.util
from itertools import permutations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINNED = {
    COMP / "lrc14_literal_fixed_sheet_allocation_thm2806.py":
        "311d0d85500f0c65945ebe5913f09d34a16293119c942b42eeaa854fbf85f71e",
    COMP / "lrc14_169_twist_variance_opus_20260726.py":
        "c0aa9c55c3e7d38dc28b4705e58c776a979f17d2874d1b762f96b97d2364e5e9",
    COMP / "lrc14_replica_dichotomy_typed_row_opus_20260727.py":
        "6ba64a68a9fd008d2e06949b1f1cf75012f1f4e734f75f55ce0af58ae20ad7b9",
}
for path, expected in PINNED.items():
    require(
        lf_hash(path) == expected,
        f"pinned dependency changed: {path}",
    )


import lrc14_literal_fixed_sheet_allocation_thm2806 as allocation
import lrc14_169_twist_variance_opus_20260726 as atlas


def crlf_safe_present_loader():
    path = COMP / "lrc14_replica_dichotomy_typed_row_opus_20260727.py"
    spec = importlib.util.spec_from_file_location(
        "q11_word_bockstein_present", path
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


allocation.physical.target.load_present_module = crlf_safe_present_loader

P = 13
DEPTH = P**5
WORDS = ("QA", "QB", "QAB")


def word_flags(interval, offset, period, pieces, starts):
    """Return whole-cylinder flags and reject every midpoint-only hit."""
    result = bytearray(DEPTH)
    for ancestry in range(DEPTH):
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
            "midpoint and whole-cylinder word flags differ",
        )
        result[ancestry] = whole_hit
    return result


def compose(left, right):
    """Relational composition: first right, then left."""
    return {
        (source, target)
        for source, middle in right
        for middle_again, target in left
        if middle == middle_again
    }


def rational_rank(matrix):
    """Exact row rank over Q, with no optional CAS dependency."""
    work = [[Fraction(entry) for entry in row] for row in matrix]
    if not work:
        return 0
    rows = len(work)
    columns = len(work[0])
    rank = 0
    for column in range(columns):
        pivot = next(
            (
                candidate
                for candidate in range(rank, rows)
                if work[candidate][column]
            ),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        work[rank] = [entry / pivot_value for entry in work[rank]]
        for row in range(rows):
            if row == rank or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [
                left - multiple * right
                for left, right in zip(work[row], work[rank])
            ]
        rank += 1
        if rank == rows:
            break
    return rank


def modular_rank(matrix, prime):
    """Exact row rank over F_prime."""
    work = [
        [entry % prime for entry in row]
        for row in matrix
    ]
    if not work:
        return 0
    rows = len(work)
    columns = len(work[0])
    rank = 0
    for column in range(columns):
        pivot = next(
            (
                candidate
                for candidate in range(rank, rows)
                if work[candidate][column]
            ),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], prime - 2, prime)
        work[rank] = [
            entry * inverse % prime
            for entry in work[rank]
        ]
        for row in range(rows):
            if row == rank or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [
                (left - multiple * right) % prime
                for left, right in zip(work[row], work[rank])
            ]
        rank += 1
        if rank == rows:
            break
    return rank


def block_hankel(matrices):
    """Return H_(q,i),(h,j)=K_(q+h)[i,j]."""
    row_size = len(matrices[0])
    column_size = len(matrices[0][0])
    return [
        [
            matrices[(q + h) % P][source][target]
            for h in range(P)
            for target in range(column_size)
        ]
        for q in range(P)
        for source in range(row_size)
    ]


def permutation_power(permutation, exponent):
    result = tuple(range(len(permutation)))
    base = permutation
    while exponent:
        if exponent & 1:
            result = tuple(
                base[result[index]]
                for index in range(len(permutation))
            )
        base = tuple(base[base[index]] for index in range(len(base)))
        exponent //= 2
    return result


def circular_shift(interval, shift, period):
    left = (interval[0] + shift) % period
    right = (interval[1] + shift) % period
    require(left < right, "selected interval crossed the circle seam")
    return left, right


def safe_comb_contains(interval, module, speed, denominator, lo, hi):
    """Test subtract_comb containment without materializing its union."""
    require(
        module.T % (denominator * speed) == 0,
        "comb grid ceased to resolve",
    )
    unit = module.T // (denominator * speed)
    step = denominator * unit
    length = (hi - lo) * unit
    base = (lo % denominator) * unit
    left, right = interval
    k0 = (left - base) // step
    for k in range(k0 - 1, k0 + 3):
        forbidden_left = base + k * step
        forbidden_right = forbidden_left + length
        if max(left, forbidden_left) < min(right, forbidden_right):
            return False
    return True


def main():
    module, full, _details, e3, clocks, _q_pairs = (
        allocation.build_geometry()
    )
    period = full.T
    unit = period // P
    atom = allocation.ATOM_INTERVAL
    target0 = tuple(x + allocation.physical.SHIFT for x in atom)
    old = allocation.physical.relative.lift.m.core.old
    require(
        period == atlas.T_DEN == old.T
        and old.rail.DEPTH == DEPTH,
        "shared physical/ancestry scale changed",
    )

    patterns = {
        "QA": atlas.PAT_QA,
        "QB": atlas.PAT_QB,
        "QAB": atlas.PAT_QAB,
    }
    require(
        tuple(
            (patterns[word][7], patterns[word][8])
            for word in WORDS
        ) == (("in", "out"), ("out", "in"), ("in", "in")),
        "outer two-bit word atlas changed",
    )
    pieces = {
        word: tuple(old.base.build_set(pattern, old.base.ZELL))
        for word, pattern in patterns.items()
    }
    starts = {
        word: tuple(left for left, _right in pieces[word])
        for word in WORDS
    }

    cache = {}

    def target_flags(word, q):
        key = word, q
        if key not in cache:
            cache[key] = word_flags(
                atom,
                allocation.physical.SHIFT + q * unit,
                period,
                pieces[word],
                starts[word],
            )
        return cache[key]

    source = {
        word: word_flags(atom, 0, period, pieces[word], starts[word])
        for word in WORDS
    }
    source_counts = {
        word: sum(source[word])
        for word in WORDS
    }
    source_bits = {
        word: int.from_bytes(source[word], "little")
        for word in WORDS
    }
    target_bits = {
        (word, q): int.from_bytes(target_flags(word, q), "little")
        for word in WORDS
        for q in range(P)
    }
    target_counts = {
        word: tuple(sum(target_flags(word, q)) for q in range(P))
        for word in WORDS
    }
    directed = {
        (source_word, target_word): tuple(
            (
                source_bits[source_word]
                & target_bits[target_word, q]
            ).bit_count()
            for q in range(P)
        )
        for source_word in WORDS
        for target_word in WORDS
    }
    expected_directed = {
        ("QA", "QA"): (0,) * P,
        ("QA", "QB"): (0,) * P,
        ("QA", "QAB"): (0,) * P,
        ("QB", "QA"):
            (0, 0, 0, 0, 0, 0, 449, 0, 449, 449, 449, 449, 449),
        ("QB", "QB"):
            (66099, 0, 0, 0, 0, 0, 0, 65612, 0, 0, 0, 0, 0),
        ("QB", "QAB"):
            (0, 0, 0, 0, 0, 0, 0, 449, 0, 0, 0, 0, 0),
        ("QAB", "QA"):
            (0, 10786, 10785, 10785, 10783, 10780, 10779,
             0, 10329, 10329, 10328, 10328, 10328),
        ("QAB", "QB"): (0,) * P,
        ("QAB", "QAB"):
            (10786, 0, 0, 0, 0, 0, 0, 10779, 0, 0, 0, 0, 0),
    }
    require(
        source_counts == {"QA": 0, "QB": 66099, "QAB": 10786}
        and target_counts["QA"]
        == (0, 10787, 10787, 10787, 10786, 10783, 11232,
            0, 10782, 10785, 10785, 10787, 10788)
        and target_counts["QB"]
        == (66099, 0, 0, 0, 0, 0, 0, 65652, 0, 0, 0, 0, 0)
        and target_counts["QAB"]
        == (10786, 0, 0, 0, 0, 0, 0, 11232, 0, 0, 0, 0, 0)
        and directed == expected_directed,
        "complete source/target word census changed",
    )

    source_qb = source["QB"]
    labels = tuple(
        ancestry
        for ancestry in range(DEPTH)
        if source_qb[ancestry] and target_flags("QA", 11)[ancestry]
    )
    require(
        len(labels) == 449
        and labels[0] == 59306
        and labels[-1] == 311961
        and {ancestry % 169 for ancestry in labels} == {156},
        "449-sheet q11 cospan changed",
    )

    neighbor = {}
    for q in (0, 7, 11):
        for delta in (-2, -1, 0, 1, 2):
            neighbor[q, delta] = tuple(
                sum(
                    target_flags(word, q)[
                        (ancestry + delta) % DEPTH
                    ]
                    for ancestry in labels
                )
                for word in WORDS
            )
    expected_neighbor = {
        (0, -2): (0, 446, 0),
        (0, -1): (0, 446, 0),
        (0, 0): (0, 449, 0),
        (0, 1): (0, 0, 449),
        (0, 2): (0, 0, 449),
        (7, -2): (0, 446, 0),
        (7, -1): (0, 447, 0),
        (7, 0): (0, 0, 449),
        (7, 1): (0, 0, 449),
        (7, 2): (0, 0, 449),
        (11, -2): (0, 0, 0),
        (11, -1): (0, 0, 0),
        (11, 0): (449, 0, 0),
        (11, 1): (449, 0, 0),
        (11, 2): (449, 0, 0),
    }
    require(
        neighbor == expected_neighbor,
        "adjacent word/carry census changed",
    )

    # The beta-selected natural arrow is (q,h)=(11,9), hence
    # (a,11)->(a+1,7).  Every one of the 449 sheets lands in QAB and none
    # lands in QB or QA.  The other wrapped ancestry arrow h=2 lands in the
    # same QAB word at q=0.
    natural_word_states = tuple(
        tuple(
            word
            for word in WORDS
            if target_flags(word, 7)[(ancestry + 1) % DEPTH]
        )
        for ancestry in labels
    )
    require(
        set(natural_word_states) == {("QAB",)}
        and neighbor[7, 1] == (0, 0, 449)
        and neighbor[0, 1] == (0, 0, 449),
        "borrow-selected QAB completion changed",
    )
    q7_qab = target_flags("QAB", 7)
    global_qab_0 = {
        ancestry
        for ancestry in range(DEPTH)
        if source_qb[ancestry] and q7_qab[ancestry]
    }
    global_qab_plus_1 = {
        ancestry
        for ancestry in range(DEPTH)
        if source_qb[ancestry]
        and q7_qab[(ancestry + 1) % DEPTH]
    }
    global_new_plus_1 = global_qab_plus_1 - global_qab_0
    require(
        len(global_qab_0) == 449
        and len(global_qab_plus_1) == 895
        and global_qab_0 < global_qab_plus_1
        and len(global_new_plus_1) == 446
        and not any(
            target_flags("QA", 11)[ancestry]
            for ancestry in global_new_plus_1
        ),
        "global QAB carry leakage changed",
    )
    reverse_exceptions = tuple(
        ancestry
        for ancestry in labels
        if not target_flags("QB", 7)[(ancestry - 1) % DEPTH]
    )
    reverse_bad_patterns = {}
    for slot in (1, 4):
        pattern = dict(patterns["QB"])
        pattern[slot] = "in"
        reverse_bad_patterns[slot] = pattern
    reverse_bad_flags = {
        slot: word_flags(
            atom,
            allocation.physical.SHIFT + 7 * unit,
            period,
            tuple(old.base.build_set(pattern, old.base.ZELL)),
            tuple(
                left
                for left, _right in old.base.build_set(
                    pattern, old.base.ZELL
                )
            ),
        )
        for slot, pattern in reverse_bad_patterns.items()
    }
    reverse_bad_factor = {}
    for ancestry in reverse_exceptions:
        hits = tuple(
            slot
            for slot in (1, 4)
            if reverse_bad_flags[slot][(ancestry - 1) % DEPTH]
        )
        require(
            len(hits) == 1,
            "reverse exception lacks a unique low-factor defect",
        )
        reverse_bad_factor[ancestry] = (hits[0], old.base.W[hits[0]])
    require(
        reverse_exceptions == (107978, 154622)
        and reverse_bad_factor
        == {107978: (1, 14), 154622: (4, 53)},
        "reverse q7 failure mechanism changed",
    )

    # The complete q-labelled support graph is triangular.  It gives a
    # constructive nonwrapping factorization R_11=R_4 o R_7, but cannot be
    # an F13 action because the inverse composite R_2 o R_11 is empty.
    relations = {
        q: {
            (source_word, target_word)
            for source_word in WORDS
            for target_word in WORDS
            if directed[source_word, target_word][q]
        }
        for q in range(P)
    }
    expected_relations = {
        0: {("QB", "QB"), ("QAB", "QAB")},
        1: {("QAB", "QA")},
        2: {("QAB", "QA")},
        3: {("QAB", "QA")},
        4: {("QAB", "QA")},
        5: {("QAB", "QA")},
        6: {("QB", "QA"), ("QAB", "QA")},
        7: {
            ("QB", "QB"), ("QB", "QAB"), ("QAB", "QAB"),
        },
        8: {("QB", "QA"), ("QAB", "QA")},
        9: {("QB", "QA"), ("QAB", "QA")},
        10: {("QB", "QA"), ("QAB", "QA")},
        11: {("QB", "QA"), ("QAB", "QA")},
        12: {("QB", "QA"), ("QAB", "QA")},
    }
    require(
        relations == expected_relations,
        "q-labelled semantic relation graph changed",
    )
    rank = {"QB": 2, "QAB": 1, "QA": 0}
    require(
        all(
            rank[target] <= rank[source]
            for relation in relations.values()
            for source, target in relation
        )
        and compose(relations[4], relations[7]) == relations[11]
        and compose(relations[2], relations[11]) == set()
        and relations[0],
        "triangular rank/functor obstruction changed",
    )

    # No action on the three semantic words can realize this family:
    # every homomorphism C13 -> S3 is trivial.  More generally, form the
    # binary and weighted 3x3 response matrices K_q.  Their 39x39
    # block-Hankel response modules have exact rank 26.  This is both the
    # universal lower bound for a rational linear realization
    # K_q=O rho(q) B and the dimension of the canonical response-profile
    # realization.  Fourier diagonalization splits the invoice as two
    # invariant dimensions plus two copies of the 12-dimensional
    # cyclotomic/augmentation representation.
    identity = tuple(range(len(WORDS)))
    require(
        all(
            permutation_power(permutation, P) != identity
            or permutation == identity
            for permutation in permutations(range(len(WORDS)))
        ),
        "a nontrivial C13 action on three words appeared",
    )
    binary_matrices = []
    weighted_matrices = []
    for q in range(P):
        binary_matrices.append([
            [
                int((source_word, target_word) in relations[q])
                for target_word in WORDS
            ]
            for source_word in WORDS
        ])
        weighted_matrices.append([
            [
                directed[source_word, target_word][q]
                for target_word in WORDS
            ]
            for source_word in WORDS
        ])
    binary_sum = [
        [
            sum(matrix[row][column] for matrix in binary_matrices)
            for column in range(len(WORDS))
        ]
        for row in range(len(WORDS))
    ]
    weighted_sum = [
        [
            sum(matrix[row][column] for matrix in weighted_matrices)
            for column in range(len(WORDS))
        ]
        for row in range(len(WORDS))
    ]
    rank_prime = 1_000_003
    binary_hankel_rank = modular_rank(
        block_hankel(binary_matrices), rank_prime
    )
    weighted_hankel_rank = modular_rank(
        block_hankel(weighted_matrices), rank_prime
    )
    qb_binary_hankel_rank = modular_rank(block_hankel([
        [matrix[WORDS.index("QB")]]
        for matrix in binary_matrices
    ]), rank_prime)
    qb_weighted_hankel_rank = modular_rank(block_hankel([
        [matrix[WORDS.index("QB")]]
        for matrix in weighted_matrices
    ]), rank_prime)
    active_source_indices = (
        WORDS.index("QB"),
        WORDS.index("QAB"),
    )
    # The explicit regular realization has one state (residue, source)
    # for every residue and active source.  Its observation at state
    # (q,source) is exactly the corresponding response row K_q(source,*).
    # This gives the matching rational upper bounds 13 and 26.
    require(
        all(
            matrix[WORDS.index("QA")] == [0, 0, 0]
            for matrix in binary_matrices + weighted_matrices
        )
        and P == 13
        and P * len(active_source_indices) == 26,
        "regular response realization dimensions changed",
    )
    require(
        binary_sum == [[0, 0, 0], [6, 2, 1], [11, 0, 2]]
        and rational_rank(binary_sum) == 2
        and rational_rank(weighted_sum) == 2
        and binary_hankel_rank == weighted_hankel_rank == 26
        and qb_binary_hankel_rank == qb_weighted_hankel_rank == 13,
        "exact rational response-module rank changed",
    )

    # The outer word slots and the address residue are different typed
    # coordinates.  The physical carrier changes low address 8 -> 7, while
    # the q=0 semantic bank has 66099 literal QB -> QB sheets.
    n_source = (6716 + 9 * P**5) % P**6
    n_target = (n_source - 1) % P**6
    require(
        (old.base.W[7], old.base.W[8]) == (P**3, 2 * P**5)
        and (n_source % P, n_target % P) == (8, 7)
        and directed["QB", "QB"][0] == 66099,
        "factor-slot/address alias hostile changed",
    )

    # Literal translation by +9 allocation units does send J_11 to J_7,
    # but it does not preserve the complete six-factor packet.
    target11 = circular_shift(target0, 11 * unit, period)
    target7 = circular_shift(target0, 7 * unit, period)
    target11_plus9 = circular_shift(target11, 9 * unit, period)
    require(
        target11_plus9 == target7,
        "literal +9 target translation ceased to land at q7",
    )

    def section_signature(interval, s, t, clock):
        return (
            allocation.contained(interval, e3),
            allocation.contained(interval, clocks[clock]),
            safe_comb_contains(
                interval, full, full.W[1], 182,
                -14 * s - 13, -14 * s + 13,
            ),
            safe_comb_contains(
                interval, full, full.W[2], 182,
                -14 * t - 13, -14 * t + 13,
            ),
            safe_comb_contains(
                interval, full, full.C2, 182,
                14 * s - 13, 14 * s + 13,
            ),
            safe_comb_contains(
                interval, full, full.C3, 182,
                14 * t - 13, 14 * t + 13,
            ),
        )

    all_567_q7_complete = 0
    for s in allocation.COMMON_S:
        for t in allocation.COMMON_T:
            for clock in range(7):
                if all(
                    all(section_signature(interval, s, t, clock))
                    for interval in (atom, target0, target7)
                ):
                    all_567_q7_complete += 1
    joint_clock_carriers = tuple(
        clock
        for clock in range(7)
        if all(
            allocation.contained(interval, clocks[clock])
            for interval in (atom, target0, target11, target7)
        )
    )
    require(
        all_567_q7_complete == 0
        and joint_clock_carriers == (1,),
        "global semantic/clock q7 repair boundary changed",
    )
    signature_census = Counter()
    for s in allocation.COMMON_S:
        for t in range(5, 12):
            require(
                all(section_signature(target11, s, t, 1)),
                "one repairing q11 cell lost its six factors",
            )
            signature_census[
                section_signature(target7, s, t, 1)
            ] += 1
    expected_signatures = Counter({
        (False, True, True, True, True, True): 35,
        (False, True, True, False, True, True): 14,
        (False, True, False, True, True, True): 10,
        (False, True, False, False, True, True): 4,
    })
    require(
        signature_census == expected_signatures
        and not any(all(signature) for signature in signature_census),
        "literal +9 packet-loss census changed",
    )

    print("THM-2835 Q11 SEMANTIC WORD HORN AND ACTION NO-GO")
    print(
        "status=VERIFIED-EXACT COMPANION; "
        "coefficient/support theorem only; no current, row decrement, or LRC14"
    )
    print(
        f"L=(count={len(labels)},first={labels[0]},last={labels[-1]},"
        "residue_mod169=156,whole_cylinder=1)"
    )
    print(
        "word_bits_(slot7,slot8)="
        "QB=(out,in),QA=(in,out),QAB=(in,in)"
    )
    for q in (0, 7, 11):
        print(
            f"neighbor_q{q}_delta_-2_to_2="
            f"{tuple(neighbor[q, delta] for delta in (-2,-1,0,1,2))}; "
            "tuple_order=(QA,QB,QAB)"
        )
    print(
        "natural_beta_lift=(a,11)->(a+1,7); "
        "word_state=QAB_on_449/449; QA=QB=0; "
        "reverse_a-1_at_q7=(QA0,QB447,QAB0); "
        f"reverse_exceptions={reverse_exceptions}; "
        f"exception_low_factors={reverse_bad_factor}"
    )
    print(
        "q7_B_delta_QAB_counts="
        f"{tuple(neighbor[7, delta][2] for delta in (-2,-1,0,1,2))}; "
        "q7_A_delta_QB_counts="
        f"{tuple(neighbor[7, delta][1] for delta in (-2,-1,0,1,2))}; "
        "B0=Bplus1=Bplus2=449 (carry-blind support fillers)"
    )
    print(
        "distinct_446_censuses="
        f"(global_sourceQB_QAB_q7:C0={len(global_qab_0)},"
        f"Cplus1={len(global_qab_plus_1)},new={len(global_new_plus_1)},"
        "new_intersects_horn=0; "
        f"restricted_horn_reverse_minus2_QB={neighbor[7, -2][1]})"
    )
    print(
        "semantic_relations="
        f"{tuple((q, tuple(sorted(relations[q]))) for q in range(P))}"
    )
    print(
        "relation_laws=R11=R4_o_R7; "
        "inverse_failure=R2_o_R11=empty!=R0; "
        "rank=(QB2,QAB1,QA0)_nonincreasing"
    )
    print(
        f"binary_relation_sum={binary_sum}; "
        f"weighted_relation_sum_rank={rational_rank(weighted_sum)}; "
        f"response_ranks=(QB_binary={qb_binary_hankel_rank},"
        f"QB_weighted={qb_weighted_hankel_rank},"
        f"full_binary={binary_hankel_rank},"
        f"full_weighted={weighted_hankel_rank}); "
        f"lower_certificate_mod_prime={rank_prime}; "
        "matching_upper_certificate=regular_response_realization; "
        "decomposition=QB:(1+12),full:(2+2*12)"
    )
    print(
        "linearization_boundary=C13_on_3_words_is_trivial; "
        "QB_response_is_cyclic_regular_module_with_nonnegative_forward_map; "
        "full_bank_minimal_rational_state_dimension=26; "
        "missing_sidecar=typed_address_carry_to_QB_response_basepoint"
    )
    print(
        f"typed_alias_hostile=word_slot_speeds_{old.base.W[7]}_"
        f"{old.base.W[8]}_versus_address_{n_source % P}_to_"
        f"{n_target % P}; q0_QB_to_QB={directed['QB', 'QB'][0]}"
    )
    print(
        f"literal_translation=J11_plus9U_equals_J7={target7}; "
        f"J7_six_factor_signature_census={tuple(sorted(signature_census.items()))}; "
        "complete_cells=0/63; "
        f"all_567_I_J0_J7_complete={all_567_q7_complete}; "
        f"joint_I_J0_J11_J7_clock_carriers={joint_clock_carriers}; "
        "factor_order=(E3,clock,q1,q2,c2,c3)"
    )
    print(
        "CONCLUSION: the 449-sheet QB(source)->QA(q11) cospan has an "
        "exact coefficient-nerve completion at the beta-selected natural "
        "lift: every (a,11) moves with its +1 ancestry borrow to a unique "
        "QAB(q7,a+1) word.  The two-bit transition QA=(in,out) to "
        "QAB=(in,in) flips only slot8.  The full relation support even "
        "factors nonwrapping q11 as q7 then q4.  But the relation is "
        "triangular and noninvertible, while the address/future source is "
        "a groupoid; its q11 inverse composite is empty.  Moreover literal "
        "physical +9U translation sends J11 to J7 and loses E3 on all "
        "63 repairing cells (sometimes q1/q2 as well).  Hence this is a "
        "positive coefficient/word 2-simplex and a sharp semantic "
        "Bockstein candidate, not yet a lawful whole-packet action or "
        "canonical current.  The canonical QB response does carry one "
        "exact regular C13 module, and the full relation bank carries two; "
        "this settles abstract rational capacity but does not supply the "
        "typed basepoint identifying the address carry column with the "
        "physical QB response."
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
