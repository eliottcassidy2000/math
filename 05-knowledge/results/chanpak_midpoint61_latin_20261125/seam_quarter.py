"""All single-arc derivatives and neutral components of the stronger lower witnesses."""
from collections import deque
from math import comb
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "chanpak_e61_fibonacci_target_20261122"))
from rail_arcs import arc_changes, inflate, two_slot_kernel
from rail_plateau import canonical, transitive_blocks
from seam_hosts import decode, FOUR_SOURCE, FOUR_CASES


ROWS = tuple(sum(b << j for j, b in enumerate(row)) for row in decode(FOUR_SOURCE))


def component(initial):
    states = {canonical(initial): 0}
    todo = deque([initial])
    records = []
    while todo:
        graph = todo.popleft()
        source, ports, widths, blocks = transitive_blocks(graph)
        edges = arc_changes(source, ports, widths, two_slot_kernel(source, ports))
        deltas = []
        neutral_destinations = []
        for ((s, a), (t, b)), change in edges:
            deltas.append(change)
            if change or (s == t and abs(a - b) == 1):
                continue
            i, j = blocks[s][a], blocks[t][b]
            neighbor = list(graph)
            neighbor[i] ^= 1 << j
            neighbor[j] ^= 1 << i
            code = canonical(neighbor)
            if code not in states:
                states[code] = len(states)
                todo.append(neighbor)
            neutral_destinations.append(states[code])
        assert len(deltas) == comb(len(graph), 2) and all(d % 2 == 0 for d in deltas)
        records.append((deltas.count(2), deltas.count(0), neutral_destinations,
                        min(d for d in deltas if d > 0), max(d for d in deltas if d < 0)))
    return records


def main():
    for ports, widths in FOUR_CASES:
        print(ports, widths, component(inflate(ROWS, ports, widths)), flush=True)


if __name__ == "__main__":
    main()
