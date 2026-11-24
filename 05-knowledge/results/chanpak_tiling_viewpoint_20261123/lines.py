"""Compare two-hole reachability with literal common words, and read the colored fork."""
from itertools import combinations, permutations
from pathlib import Path

CATALOG = Path(__file__).resolve().parents[3] / "session212" / "cycle-levels"


def directed_paths(rows):
    return {word for word in permutations(range(len(rows)))
            if all(rows[a] >> b & 1 for a, b in zip(word, word[1:]))}


def change(rows, pairs, erase=False):
    rows = list(rows)
    for a, b in pairs:
        if erase:
            rows[a] &= ~(1 << b)
            rows[b] &= ~(1 << a)
        else:
            rows[a] ^= 1 << b
            rows[b] ^= 1 << a
    return tuple(rows)


def reachability(rows):
    reach = [row | 1 << i for i, row in enumerate(rows)]
    for k in range(len(rows)):
        for i in range(len(rows)):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return reach


def endpoints_comparable(rows, pairs):
    reach = reachability(rows)
    return all(reach[a] >> b & 1 or reach[b] >> a & 1 for a, b in pairs)


def census():
    for n in range(3, 8):
        edges = list(combinations(range(n), 2))
        cases = obstructions = minimal = 0
        for line in (CATALOG / f"classes{n}.rows").read_text().splitlines():
            rows = tuple(map(int, line.split()))
            words = directed_paths(rows)
            # One bit per directed Hamilton path of the full tournament.
            using = {edge: 0 for edge in edges}
            for j, word in enumerate(words):
                for x, y in zip(word, word[1:]):
                    using[tuple(sorted((x, y)))] |= 1 << j
            all_words = (1 << len(words)) - 1
            for e, f in combinations(edges, 2):
                common = (using[e] | using[f]) != all_words
                d = change(rows, (e, f), erase=True)
                assert common == endpoints_comparable(d, (e, f))
                cases += 1
                obstructions += not common
                minimal += not common and using[e] != all_words and using[f] != all_words
        print(n, "rooted at listed class representatives, choices of two marked pairs",
              cases, "no common path", obstructions, "minimal two-pair failures", minimal)


def partners(rows, path):
    chords = [(path[i], path[j]) for i in range(len(path))
              for j in range(i+2, len(path))]
    return change(rows, chords)


def blue(rows, path):
    return all(((rows[a] >> b) & 1) ==
               ((rows[path[-1-j]] >> path[-1-i]) & 1)
               for i, a in enumerate(path) for j, b in enumerate(path) if j > i+1)


def fork():
    holes = ((0, 3), (1, 4))
    start = (28, 17, 10, 18, 4)
    corners = [change(start, [e for j, e in enumerate(holes) if mask >> j & 1])
               for mask in range(4)]
    paths = [directed_paths(r) for r in corners]
    assert [len(w) for w in paths] == [13, 11, 15, 13]
    common = sorted(set.intersection(*paths))
    results = [([len(directed_paths(partners(r, p))) for r in corners],
                ''.join(str(int(blue(r, p))) for r in corners)) for p in common]
    assert results == [([13, 9, 15, 13], "0110"), ([13, 11, 5, 5], "0000"),
                       ([9, 5, 5, 3], "0001"), ([13, 11, 9, 9], "0000"),
                       ([9, 11, 9, 13], "0000"), ([5, 11, 5, 13], "0000"),
                       ([3, 5, 5, 9], "1000")]
    for word, reading in zip(common, results):
        print("".join(map(str, word)), *reading)


def obstructions():
    t = (0, 5, 9, 3)
    e, f = (1, 2), (0, 1)
    assert [len(directed_paths(change(t, edits))) for edits in ((), (e,), (f,), (e,f))] == [3,1,5,1]
    assert not directed_paths(change(t, (e,f), erase=True))
    triangle = (8, 8, 8, 16, 7)
    assert not directed_paths(triangle) and reachability(triangle) == [31]*5
    path = change((30,20,24,18,0), ((0,1),(1,4),(2,4)), erase=True)
    assert not directed_paths(path) and endpoints_comparable(path, combinations(range(5), 2))


if __name__ == "__main__":
    obstructions()
    fork()
    census()
