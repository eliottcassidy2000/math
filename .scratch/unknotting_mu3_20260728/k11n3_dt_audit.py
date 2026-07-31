#!/usr/bin/env python3
"""Independent PD -> signed DT -> KnotTheory table identification audit.

The PD convention is the one in the prompt and KnotTheory: in [a,b,c,d],
slots a,c are the incoming/outgoing under-strand ends and the four slots
are cyclically ordered.  The orientation of the over strand is propagated
from the condition that each arc has one incoming and one outgoing end.
"""

import re
from pathlib import Path


PD = [
    [4, 2, 5, 1],
    [10, 6, 11, 5],
    [8, 3, 9, 4],
    [2, 9, 3, 10],
    [11, 16, 12, 17],
    [7, 15, 8, 14],
    [15, 7, 16, 6],
    [13, 20, 14, 21],
    [17, 22, 18, 1],
    [21, 18, 22, 19],
    [19, 12, 20, 13],
]

EXPECTED_DT = [4, 8, 10, -14, 2, -16, -20, -6, -22, -12, -18]
EXPECTED_STRING = "bdeGaHJCKFI"
EXPECTED_NAME = "K11n3"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def occurrences(pd):
    occ = {}
    for crossing, row in enumerate(pd):
        require(len(row) == 4, "each PD crossing must have four slots")
        for slot, label in enumerate(row):
            occ.setdefault(label, []).append((crossing, slot))
    n = len(pd)
    require(set(occ) == set(range(1, 2 * n + 1)),
            "labels are not exactly 1,...,2n")
    require(all(len(ends) == 2 for ends in occ.values()),
            "each arc label must occur exactly twice")
    return occ


def orient_ends(pd, occ):
    # True means incoming at the crossing.
    incoming = {}
    for crossing in range(len(pd)):
        incoming[(crossing, 0)] = True
        incoming[(crossing, 2)] = False

    def other_end(label, here):
        p, q = occ[label]
        return q if p == here else p

    changed = True
    while changed:
        changed = False
        for crossing, row in enumerate(pd):
            if (crossing, 1) in incoming:
                continue
            for slot in (1, 3):
                far = other_end(row[slot], (crossing, slot))
                if far in incoming:
                    incoming[(crossing, slot)] = not incoming[far]
                    incoming[(crossing, (slot + 2) % 4)] = incoming[far]
                    changed = True
                    break
    require(len(incoming) == 4 * len(pd),
            "over-strand orientation did not propagate")
    for label, ends in occ.items():
        require(sorted(incoming[end] for end in ends) == [False, True],
                "arc %d is not once-in/once-out" % label)
    return incoming


def signed_dt(pd):
    occ = occurrences(pd)
    incoming = orient_ends(pd, occ)

    # Start by entering along arc 1.  At every encounter, leave through the
    # opposite slot and follow that arc to its other endpoint.
    current = next(end for end in occ[1] if incoming[end])
    encounters = []
    for _ in range(2 * len(pd)):
        crossing, slot = current
        encounters.append((crossing, slot, slot % 2 == 1))
        out_slot = (slot + 2) % 4
        label = pd[crossing][out_slot]
        p, q = occ[label]
        current = q if p == (crossing, out_slot) else p

    require(current == next(end for end in occ[1] if incoming[end]),
            "traversal did not close")
    require(len({(i, crossing) for i, (crossing, _, _) in
                 enumerate(encounters)}) == 2 * len(pd),
            "bad encounter traversal")

    positions = {}
    for number, (crossing, _, is_over) in enumerate(encounters, 1):
        positions.setdefault(crossing, []).append((number, is_over))
    require(all(len(pair) == 2 for pair in positions.values()),
            "a crossing was not encountered exactly twice")

    code = []
    for odd in range(1, 2 * len(pd) + 1, 2):
        crossing, _, odd_is_over = encounters[odd - 1]
        pair = positions[crossing]
        mate = next(number for number, _ in pair if number != odd)
        require(mate % 2 == 0, "DT parity condition failed")
        code.append(mate if odd_is_over else -mate)
    return code, encounters


def encode_knot_theory(dt):
    chars = []
    for value in dt:
        k = abs(value) // 2
        chars.append(chr((96 if value > 0 else 64) + k))
    return "".join(chars)


def table_name(code_string, table_path):
    text = table_path.read_text()
    body = text.split("DTStringsTo11 = {", 1)[1].split("\n}", 1)[0]
    strings = re.findall(r'"([A-Za-z]+)"', body)
    require(len(strings) == 801, "unexpected KnotTheory DT table length")
    require(strings.count(code_string) == 1, "DT string is not unique")
    index = strings.index(code_string)

    # NumberOfKnots[n,Alternating/NonAlternating] used by the package loop.
    counts = [
        (3, 1, 0), (4, 1, 0), (5, 2, 0), (6, 3, 0), (7, 7, 0),
        (8, 18, 3), (9, 41, 8), (10, 123, 42), (11, 367, 185),
    ]
    require(sum(a + non_a for _, a, non_a in counts) == len(strings),
            "crossing-count partition does not cover table")
    cursor = 0
    for crossings, alternating, nonalternating in counts:
        for kind, count in (("a", alternating), ("n", nonalternating)):
            if cursor <= index < cursor + count:
                return "K%d%s%d" % (crossings, kind, index - cursor + 1), index
            cursor += count
    raise RuntimeError("DT string index was not assigned a knot name")


def main():
    dt, encounters = signed_dt(PD)
    require(dt == EXPECTED_DT, "derived signed DT differs from expected")
    string = encode_knot_theory(dt)
    require(string == EXPECTED_STRING, "KnotTheory encoding differs")
    path = Path(__file__).with_name("DTCode4KnotsTo11.m")
    name, zero_index = table_name(string, path)
    require(name == EXPECTED_NAME, "table name differs from K11n3")

    gauss = [
        "%d:%s@X%d" % (i, "O" if over else "U", crossing + 1)
        for i, (crossing, _, over) in enumerate(encounters, 1)
    ]
    print("encounters =", " ".join(gauss))
    print("signed DT =", dt)
    print("KnotTheory string =", string)
    print("table zero-based index =", zero_index)
    print("identified knot =", name)


if __name__ == "__main__":
    main()
