"""Consume the height-capped census to give ordinary and unit-neighbor class series."""
from collections import Counter, defaultdict
from functools import lru_cache
import sys


def read(stream):
    fibres = defaultdict(dict)
    closed = False
    for line in stream:
        fields = line.split()
        if fields[:1] == ["fibre"]:
            n, h, total, lower, upper, both = map(int, fields[1:])
            fibres[h][n] = (total, lower, upper, both)
        elif fields[-1:] == ["closed"]:
            closed = True
    if not closed or not fibres:
        raise ValueError("full --fibres output through an empty generation required")
    return fibres


def grammar(fibres):
    @lru_cache(None)
    def words(h):
        # Key is (number of nonsingleton components, total strong vertices).
        if h == 1:
            return {(0, 0): 1}
        result = Counter()
        for d in range(3, h + 1, 2):
            if h % d or d not in fibres:
                continue
            for (length, total), value in words(h // d).items():
                for n, fields in fibres[d].items():
                    result[length+1, total+n] += value * fields[0]
        return result
    return words


def main():
    if len(sys.argv) not in (2, 3):
        raise SystemExit("usage: capped_cores.py cap [odd-height] < capped-fibres")
    cap = int(sys.argv[1])
    fibres = read(sys.stdin)
    if cap != max(fibres):
        raise ValueError("this summary takes a cap equal to the largest attained strong height")
    words = grammar(fibres)
    if len(sys.argv) == 3:
        h = int(sys.argv[2])
        if not 1 <= h <= cap or h % 2 == 0:
            raise ValueError("odd target in the completed range required")
        for (k, n), coefficient in sorted(words(h).items()):
            print(f"ordinary {coefficient} * z^{n} / (1-z)^{k+1}")
        if h >= 5:
            for n, (_, lower, upper, both) in sorted(fibres[h].items()):
                print("strong core", n, "lower", lower, "upper", upper, "both", both)
    else:
        missing = [h for h in range(1,cap+1,2) if not words(h)]
        no_step = [h for h in range(5,cap-1,2)
                   if words(h) and words(h+2) and
                   not any(row[2] for row in fibres[h].values())]
        unique = []
        for h in range(1,cap+1,2):
            records = words(h)
            if not records:
                continue
            min_order = min(n for _, n in records)
            if sum(count for (_, n),count in records.items() if n == min_order) == 1:
                unique.append(h)
        print("missing ordinary heights", missing)
        print("attainable adjacent heights without a single-arc reversal", no_step)
        print("heights with exactly one minimum ordinary class", unique)


if __name__ == "__main__":
    main()
