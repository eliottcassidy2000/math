"""Hamilton path permutation signs and the cyclic-deletion lower bound."""
from itertools import permutations
from pathlib import Path

BANK = Path(__file__).resolve().parents[3]/"session212"/"cycle-levels"


def words(n):
    for path in permutations(range(n)):
        inversions = sum(b < a for i,a in enumerate(path) for b in path[i+1:])
        yield path, inversions % 2


def directed(rows, word):
    return all(rows[a] >> b & 1 for a,b in zip(word, word[1:]))


def minor_bound(rows):
    n = len(rows)
    m = n//2
    total = 0
    for v in range(n):
        others = [i for i in range(n) if i != v]
        start = others[0]
        for tail in permutations(others[1:]):
            ring = (start,)+tail
            if not directed(rows, ring+(start,)):
                continue
            p = sum(rows[v] >> a & 1 for a in ring[::2])
            q = sum(rows[v] >> a & 1 for a in ring[1::2])
            total += m - abs(p-q)
    return total


def main():
    for n in (3,5,7):
        bank = [list(map(int, row.split()))
                for row in (BANK/f"classes{n}.rows").read_text().splitlines()]
        assert len(bank) == {3:2,5:12,7:456}[n]
        paths = list(words(n))
        strong_pure = []
        equality = 0
        for rows in bank:
            counts = [0,0]
            strong = False
            for path, sign in paths:
                if directed(rows,path):
                    counts[sign] += 1
                    strong |= bool(rows[path[-1]] >> path[0] & 1)
            if strong and min(counts)==0:
                strong_pure.append((sum(counts), rows))
            bound = minor_bound(rows)
            assert bound <= 2*min(counts)
            if strong and 0 < bound == 2*min(counts):
                equality += 1
        assert [h for h,_ in strong_pure] == {3:[3],5:[15],7:[]}[n]
        print(n,"strong classes with one sign",strong_pure,
              "positive sharp strong cases",equality)


if __name__ == "__main__":
    main()
