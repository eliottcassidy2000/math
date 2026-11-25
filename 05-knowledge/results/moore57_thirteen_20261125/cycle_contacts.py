"""Exact local contact decks on an induced cycle of length 13."""
from collections import Counter
from itertools import combinations
from ortools.sat.python import cp_model
from moore_chart import hoffman_edge


def admissible(s):
    return all(2 < j-i < 11 for i, j in combinations(s, 2))


def crossing(e, f):
    a, b = e
    c, d = f
    return a < c < b < d or c < a < d < b


def main():
    cycle = (0, 1, 2, 3, 28, 8, 7, 6, 44, 15, 16, 17, 30)
    assert all(hoffman_edge(u, v) == (j-i in (1, 12))
               for (i, u), (j, v) in combinations(enumerate(cycle), 2))
    contacts = Counter(sum(hoffman_edge(x, y) for y in cycle)
                       for x in range(50) if x not in cycle)
    assert [contacts[j] for j in range(5)] == [7, 9, 10, 8, 3]
    print("Hoffman–Singleton exterior contact counts", dict(sorted(contacts.items())))
    blocks = {k: [s for s in combinations(range(13), k) if admissible(s)]
              for k in (2, 3, 4, 5)}
    assert list(map(len, blocks.values())) == [52, 65, 13, 0]
    assert sum(crossing(e, f) for e, f in combinations(blocks[2], 2)) == 598
    for fours, expected in enumerate((15, 13, 11, 9, 5)):
        model = cp_model.CpModel()
        select = {s: model.NewBoolVar(str(s)) for k in (3, 4) for s in blocks[k]}
        for e in blocks[2]:
            model.Add(sum(x for s, x in select.items() if set(e) <= set(s)) <= 1)
        model.Add(sum(select[s] for s in blocks[4]) == fours)
        if fours:  # All four-contact supports are translates; choose the first.
            model.Add(select[0, 3, 6, 9] == 1)
        model.Maximize(sum(select[s] for s in blocks[3]))
        solver = cp_model.CpSolver()
        solver.parameters.num_search_workers = 1
        assert solver.Solve(model) == cp_model.OPTIMAL
        witness = [s for s, var in select.items() if solver.Value(var)]
        assert all(len(set(s) & set(t)) <= 1 for s, t in combinations(witness, 2))
        assert sum(len(s) == 3 for s in witness) == expected
        print("four-contact sites", fours, "maximum three-contact sites", expected,
              "minimum two-contact sites", 52 - 6*fours - 3*expected)


if __name__ == "__main__":
    main()
