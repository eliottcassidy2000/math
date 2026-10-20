"""Exact small controls for square_counts.md; standard library only.

Universe: the four-vertex tournament with forward path 0,1,2,3 and
backward chords 2->0,3->1, and its ordinal powers of lengths 0 through 4.
No pruning beyond legal directed Hamiltonian steps; independent product
words are compared with a direct DFS in each full tournament.
"""
from itertools import product
from pathlib import Path
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def power(r):
    size = 4 * r
    arcs = set()
    for i in range(size):
        for j in range(i + 1, size):
            backward = i // 4 == j // 4 and (i % 4, j % 4) in {(0, 2), (1, 3)}
            arcs.add((j, i) if backward else (i, j))
    return size, arcs


def paths(size, arcs):
    answer = []
    def visit(prefix, remaining):
        if not remaining:
            answer.append(prefix)
        for v in remaining:
            if not prefix or (prefix[-1], v) in arcs:
                visit(prefix + (v,), remaining - {v})
    visit((), set(range(size)))
    return set(answer)


def reverse_reflect(path):
    return tuple(len(path) - 1 - v for v in reversed(path))


def main():
    base = paths(*power(1))
    require(len(base) == 5, "Q must have exactly five directed Hamiltonian paths")
    rows = []
    all_paths = {}
    for r in range(5):
        size, arcs = power(r)
        direct = paths(size, arcs)
        decoded = {
            tuple(v + 4 * block for block, word in enumerate(words) for v in word)
            for words in product(base, repeat=r)
        }
        require(direct == decoded, f"product decoder failed at r={r}")
        require(all((size-1-v, size-1-u) in arcs for u, v in arcs), "anti-automorphism")
        blue = {p for p in direct if reverse_reflect(p) == p}
        require(len(direct) == 5**r, "power count")
        expected_blue = 5 ** (r//2) * (3 if r % 2 else 1)
        require(len(blue) == expected_blue, "blue formula")
        all_paths[r] = direct
        rows.append({"blocks": r, "vertices": size, "H": len(direct), "blue": len(blue)})

    x = all_paths[2]
    standard = tuple(range(8))
    require(standard in x and reverse_reflect(standard) == standard, "pointed control")
    pairs = set(product(x, repeat=2))
    diagonal = {(p, p) for p in x}
    pointed_square = {(p, q) for p, q in pairs if p != standard and q != standard}
    require(len(pairs - diagonal) == 600, "off-diagonal count")
    require(len(pointed_square) == 576, "point-deleted product count")
    require(len(pointed_square & diagonal) == 24, "retained diagonal")
    require(len(pointed_square - diagonal) == 552, "retained off-diagonal")
    same_port_blue = tuple(range(16))
    same_port_black = (0,1,2,3,4,7,5,6,8,9,10,11,12,13,14,15)
    require(same_port_blue in all_paths[4] and same_port_black in all_paths[4], "same-port legality")
    require(reverse_reflect(same_port_blue) == same_port_blue and
            reverse_reflect(same_port_black) != same_port_black, "same-port colour hostile")
    # The second half is reversed/reflected to put blue paths on a diagonal.
    decode = lambda pair: pair[0] + tuple(v+8 for v in reverse_reflect(pair[1]))
    require({decode(pair) for pair in pairs} == all_paths[4], "square path decoder")
    require(all((reverse_reflect(decode(pair)) == decode(pair)) == (pair in diagonal)
                for pair in pairs), "literal blue diagonal")
    selected_sum = sum(len(all_paths[r]) for r in (0, 2, 3, 4))
    require(selected_sum == 776, "selected disjoint union")
    result = {
        "status": "FINITE-EXACT; structural proofs in square_counts.md",
        "base_Q_paths": sorted(base), "ordinal_powers": rows,
        "square": {"all": len(pairs), "diagonal": len(diagonal), "off_diagonal": 600,
                   "point_deleted_product": 576, "retained_diagonal": 24,
                   "retained_off_diagonal": 552},
        "selected_lengths": [0, 2, 3, 4], "selected_sum": selected_sum,
        "omitted_length_one": 5,
        "same_port_hostile": {"blue": same_port_blue, "black": same_port_black},
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    Path(__file__).with_suffix(".json").write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
