#!/usr/bin/env python3
"""Parent-state provenance audit for THM-4065.

This no-import reconstruction uses only THM-4047's seven already-audited
physical query bits.  It rebuilds all fixed left-front tails through depth
100000, checks the complete zero census, and freezes the two absolute-phase
period-32 words used by the C transfer certificates.
"""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path


MAX_R = 100_000
CHECKS = 0
REPO_ROOT = Path(__file__).resolve().parents[1]
PARENT_SOURCE = REPO_ROOT / "04-computation/rule30_left_front_affine_monodromy_thm4047.py"
PARENT_OUTPUT = REPO_ROOT / "05-knowledge/results/rule30_left_front_affine_monodromy_thm4047.out"
PARENT_SOURCE_SHA = "d454cc5b40315c02ebf486f29e227ebc9a79c78ede18fc449a7dd32f8dc21148"
PARENT_OUTPUT_SHA = "bae500127999c350ebeff77c3145fbb4abf7b5b7292d2882eeecc6492fba75e3"
PHYSICAL_QUERIES = {
    3: (2, 0, (1, 1)),
    8: (10, 0, (1, 1)),
    29: (90, 0, (1, 1)),
    400: (3050, 1, (1, 1)),
    53208: (847962, 0, (1, 0)),
    58287: (929210, 1, (1, 0)),
    87867: (1402474, 0, (1, 1)),
}


def require(condition, label):
    global CHECKS
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")
    CHECKS += 1


class Tail:
    __slots__ = ("start", "word", "monodromy")

    def __init__(self, start, word, monodromy):
        self.start = start
        self.word = word
        self.monodromy = monodromy

    def at(self, time):
        require(time >= self.start, ("pre-tail query", time, self.start))
        return self.word[time % len(self.word)]


def rule30_step(left, right, state):
    return left ^ right ^ state ^ (right & state)


def least_period(values):
    raw_length = len(values)
    for candidate in range(1, raw_length + 1):
        if raw_length % candidate == 0 and all(
            values[index] == values[index % candidate]
            for index in range(raw_length)
        ):
            return candidate
    raise RuntimeError("least period not found")


def absolute_phase_word(start, values):
    period = least_period(values)
    answer = [None] * period
    for offset, value in enumerate(values):
        residue = (start + offset) % period
        require(
            answer[residue] is None or answer[residue] == value,
            ("phase collision", start, residue),
        )
        answer[residue] = value
    require(all(value is not None for value in answer), ("phase holes", start))
    return tuple(answer)


def packed(word):
    return sum(value << index for index, value in enumerate(word))


def reconstruct():
    tails = [Tail(0, (1,), (1, 0)), Tail(1, (1,), (0, 1))]
    used_queries = set()
    for depth in range(2, MAX_R + 1):
        left = tails[depth - 2]
        right = tails[depth - 1]
        common_start = max(left.start, right.start)
        common_period = max(len(left.word), len(right.word))
        require(
            common_period % len(left.word) == 0
            and common_period % len(right.word) == 0,
            ("dyadic period nesting", depth),
        )

        image_zero = 0
        image_one = 1
        for time in range(common_start, common_start + common_period):
            image_zero = rule30_step(left.at(time), right.at(time), image_zero)
            image_one = rule30_step(left.at(time), right.at(time), image_one)
        monodromy = (image_zero ^ image_one, image_zero)

        if monodromy[0] == 0:
            tail_start = common_start + common_period
            sample_length = common_period
            initial = monodromy[1]
        else:
            require(depth in PHYSICAL_QUERIES, ("missing phase", depth))
            query_start, initial, expected_monodromy = PHYSICAL_QUERIES[depth]
            require(
                (common_start, monodromy) == (query_start, expected_monodromy),
                ("phase query drift", depth),
            )
            used_queries.add(depth)
            tail_start = common_start
            sample_length = common_period if monodromy[1] == 0 else 2 * common_period

        values = []
        state = initial
        for time in range(tail_start, tail_start + sample_length):
            values.append(state)
            state = rule30_step(left.at(time), right.at(time), state)
        require(state == initial, ("sample wrap", depth))
        word = absolute_phase_word(tail_start, tuple(values))
        tail = Tail(tail_start, word, monodromy)
        for time in range(
            tail_start,
            tail_start + 2 * max(common_period, len(word)),
        ):
            require(
                tail.at(time + 1)
                == rule30_step(left.at(time), right.at(time), tail.at(time)),
                ("wrapped law", depth, time),
            )
        tails.append(tail)

    require(used_queries == set(PHYSICAL_QUERIES), "physical packet consumption")
    return tails


def main():
    require(
        sha256(PARENT_SOURCE.read_bytes()).hexdigest() == PARENT_SOURCE_SHA,
        "THM4047 source drift",
    )
    require(
        sha256(PARENT_OUTPUT.read_bytes()).hexdigest() == PARENT_OUTPUT_SHA,
        "THM4047 output drift",
    )
    parent_text = PARENT_OUTPUT.read_text(encoding="utf-8")
    require(
        "zero_columns=(2, 7, 28, 399, 53207, 58286, 87866)" in parent_text,
        "parent zero line",
    )
    require(
        "certificate_sha256=62c469a9af7fd4b6ed2eeed6aef384a76e1712d848cca3dd01fa157f6440f915"
        in parent_text,
        "parent semantic certificate",
    )

    tails = reconstruct()
    zeros = tuple(
        depth for depth, tail in enumerate(tails) if not any(tail.word)
    )
    require(zeros == (2, 7, 28, 399, 53207, 58286, 87866), "zero census")

    targets = []
    expected = {
        87867: (1402474, 32, 0xCF3030CF),
        99998: (1790666, 32, 0xF687520B),
        99999: (1790698, 32, 0xBE79924B),
        100000: (1790730, 32, 0x90F58380),
    }
    for depth, target in expected.items():
        tail = tails[depth]
        row = (tail.start, len(tail.word), packed(tail.word))
        require(row == target, ("target word", depth, row))
        targets.append((depth, tail.start, len(tail.word), f"{packed(tail.word):08x}"))

    left = tails[99999]
    right = tails[100000]
    for time in range(max(left.start, right.start), max(left.start, right.start) + 64):
        next_value = rule30_step(
            tails[99998].at(time), left.at(time), right.at(time)
        )
        require(next_value == right.at(time + 1), ("target recurrence", time))

    print("RULE30_TEMPORAL_CYLINDER_TRANSFER_THM4065_PARENT_AUDIT")
    print("method=no_import_query_packet_reconstruction")
    print(f"parent_source_sha256={PARENT_SOURCE_SHA}")
    print(f"parent_output_sha256={PARENT_OUTPUT_SHA}")
    print("parent_semantic_sha256=62c469a9af7fd4b6ed2eeed6aef384a76e1712d848cca3dd01fa157f6440f915")
    print(f"zero_columns={zeros}")
    print(f"target_absolute_words={tuple(targets)}")
    print(f"checks={CHECKS}")
    print("status=PASS;conditional_only_on_already_audited_THM4047_physical_query_packet")


if __name__ == "__main__":
    main()
