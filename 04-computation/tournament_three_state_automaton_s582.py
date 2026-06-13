#!/usr/bin/env python3
"""Three-state pair automata as a bridge from event words to tournaments.

The point is intentionally small: each unordered pair reads a word over
{L, M, R}.  Its terminal state is projected to a tournament edge, while M is
resolved by a fixed Hamiltonian tie path.  This keeps the wall/tie state visible
before the usual Tournament Analysis projection erases it.
"""

from collections import Counter
from itertools import combinations, permutations, product


VERTICES = ("A", "B", "C", "D", "E")
TIE_ORDER = {v: i for i, v in enumerate(VERTICES)}

# Words are read from the lower tie-order vertex toward the higher one.
# L means the lower vertex currently owns the pair, R means the higher vertex
# owns it, and M means the pair is on the wall / middle corridor.
WORDS = {
    ("A", "B"): "LMLL",
    ("A", "C"): "RMRR",
    ("A", "D"): "MM",
    ("A", "E"): "RML",
    ("B", "C"): "LL",
    ("B", "D"): "RLR",
    ("B", "E"): "RR",
    ("C", "D"): "LML",
    ("C", "E"): "MRM",
    ("D", "E"): "RMR",
}


AUTOMATA = {
    "wall": {
        "M": {"L": "L", "M": "M", "R": "R"},
        "L": {"L": "L", "M": "M", "R": "M"},
        "R": {"L": "M", "M": "M", "R": "R"},
    },
    "last_nonmiddle": {
        "M": {"L": "L", "M": "M", "R": "R"},
        "L": {"L": "L", "M": "L", "R": "R"},
        "R": {"L": "L", "M": "R", "R": "R"},
    },
}


def run_word(table, word):
    state = "M"
    trace = [state]
    for symbol in word:
        state = table[state][symbol]
        trace.append(state)
    return state, trace


def orient(pair, terminal):
    lo, hi = pair
    if terminal == "L":
        return lo, hi, "terminal L"
    if terminal == "R":
        return hi, lo, "terminal R"
    return lo, hi, "tie path"


def build_tournament(name):
    table = AUTOMATA[name]
    matrix = {u: {v: False for v in VERTICES if v != u} for u in VERTICES}
    records = []
    for pair in combinations(VERTICES, 2):
        terminal, trace = run_word(table, WORDS[pair])
        winner, loser, reason = orient(pair, terminal)
        matrix[winner][loser] = True
        matrix[loser][winner] = False
        records.append((pair, WORDS[pair], terminal, winner, loser, reason, trace))
    return matrix, records


def score_histogram(matrix):
    scores = {u: sum(matrix[u].values()) for u in VERTICES}
    return scores, dict(sorted(Counter(scores.values()).items()))


def directed_3cycles(matrix):
    cycles = []
    for triple in combinations(VERTICES, 3):
        for perm in permutations(triple):
            a, b, c = perm
            if matrix[a][b] and matrix[b][c] and matrix[c][a]:
                cycles.append(perm)
                break
    return cycles


def strongly_connected_components(matrix):
    index = 0
    stack = []
    on_stack = set()
    indices = {}
    low = {}
    components = []

    def visit(v):
        nonlocal index
        indices[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)
        for w in VERTICES:
            if w == v or not matrix[v][w]:
                continue
            if w not in indices:
                visit(w)
                low[v] = min(low[v], low[w])
            elif w in on_stack:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            comp = []
            while True:
                w = stack.pop()
                on_stack.remove(w)
                comp.append(w)
                if w == v:
                    break
            components.append(tuple(sorted(comp, key=TIE_ORDER.get)))

    for v in VERTICES:
        if v not in indices:
            visit(v)
    return components


def hamiltonian_path_count(matrix):
    count = 0
    examples = []
    for path in permutations(VERTICES):
        if all(matrix[path[i]][path[i + 1]] for i in range(len(path) - 1)):
            count += 1
            if len(examples) < 4:
                examples.append(path)
    return count, examples


def side_to_side_transitions(table):
    total = 0
    bad = []
    for state, symbol in product(("L", "R"), ("L", "M", "R")):
        target = table[state][symbol]
        if state != target and target in ("L", "R"):
            total += 1
            bad.append((state, symbol, target))
    return total, bad


def print_tournament(name):
    matrix, records = build_tournament(name)
    scores, hist = score_histogram(matrix)
    cycles = directed_3cycles(matrix)
    sccs = strongly_connected_components(matrix)
    h_count, h_examples = hamiltonian_path_count(matrix)
    direct_flips, flip_examples = side_to_side_transitions(AUTOMATA[name])

    print(f"\nAutomaton: {name}")
    print(f"  direct side-to-side state transitions: {direct_flips} {flip_examples}")
    print("  pair terminals:")
    for pair, word, terminal, winner, loser, reason, trace in records:
        joined_trace = "".join(trace)
        print(
            f"    {pair[0]}{pair[1]} word={word:<4} trace={joined_trace:<5} "
            f"terminal={terminal} edge={winner}->{loser} ({reason})"
        )
    print(f"  scores: {scores}")
    print(f"  score histogram: {hist}")
    print(f"  directed 3-cycles: {len(cycles)} {cycles[:6]}")
    print(f"  SCCs: {sccs}")
    print(f"  Hamiltonian path count: {h_count}; examples={h_examples}")
    return matrix


def edge_flips(matrix_a, matrix_b):
    flips = []
    for u, v in combinations(VERTICES, 2):
        winner_a = u if matrix_a[u][v] else v
        winner_b = u if matrix_b[u][v] else v
        if winner_a != winner_b:
            flips.append((u, v, winner_a, winner_b))
    return flips


def main():
    print("Three-state Tournament Analysis automaton (S582)")
    print("Observable: pair event word over {L,M,R}.")
    print("Switch: terminal L gives lower->higher, terminal R gives higher->lower.")
    print("Tie Hamiltonian path: A->B->C->D->E resolves terminal M.")
    print("Interpretation: M is the wall/corridor state retained before projection.")

    wall = print_tournament("wall")
    last = print_tournament("last_nonmiddle")

    flips = edge_flips(wall, last)
    print("\nRule comparison")
    print(f"  edge flips between wall and last_nonmiddle: {len(flips)} {flips}")
    print(
        "  lesson: the wall automaton forbids direct L<->R state flips, so every "
        "orientation change must pass through M before it becomes a tournament edge."
    )

    pair_count = len(VERTICES) * (len(VERTICES) - 1) // 2
    print("\nCompression note")
    print(f"  pair automata for n={len(VERTICES)} use {pair_count} local 3-state cells.")
    print(f"  product state count is 3^{pair_count}, but one scan costs O(pair-events).")
    print(
        "  for LRC, replace A..E by components, residue states, cover arcs, or "
        "proof obligations; project only after the middle states have been audited."
    )


if __name__ == "__main__":
    main()
