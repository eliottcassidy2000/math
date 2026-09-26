"""Exact hostile controls for incoming balance and carry-clock claims."""
from fractions import Fraction as F
from itertools import combinations
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def cell(k, a):
    words = []
    for positions in combinations(range(k), a):
        occupied = set(positions)
        word = [int(i in occupied) for i in range(k)]
        if all(3 ** sum(word[:j]) > 2 ** j for j in range(1, k + 1)):
            words.append(word)
    best = F(0)
    for i in range(a):
        for j in range(k - a):
            count = 0
            for w in words:
                one = [t for t, b in enumerate(w) if b]
                zero = [t for t, b in enumerate(w) if not b]
                count += one[i] < zero[j]
            best = max(best, F(min(count, len(words) - count), len(words)))
    return len(words), best


def main():
    rows = []
    for k, a, expected_count, expected_delta in [(5, 4, 3, F(1, 3)), (7, 5, 7, F(3, 7)), (10, 7, 30, F(7, 15))]:
        count, delta = cell(k, a)
        require((count, delta) == (expected_count, expected_delta), "exact balance hostile")
        rows.append({"k": k, "odd": a, "extensions": count, "balance": str(delta)})
    # Exact same-band segment; an internal band visit prevents overclaiming
    # this as a counterexample specifically about consecutive band visits.
    n = 4067
    states, total_depth = [n], 0
    for _ in range(29):
        x = 3 * n + 1
        a = (x & -x).bit_length() - 1
        n = x >> a
        total_depth += a
        states.append(n)
    require(n == 2051 and total_depth == 47, "carry-sensitive same-band witness")
    require(2 ** 45 < 3 ** 29 < 2 ** 46, "carry-free floors are 45 and 46")
    require(len(set(states)) == len(states), "distinct witness")
    require(2429 in states[1:-1], "consecutive-visit scope guard")
    # At the first odd step27->41, C1=82/81 and reciprocal term1/27.
    # log(1+x)<x proves the asserted <=3log(C1) has the wrong direction.
    require(F(82, 81) - 1 == F(1, 81), "clock multiplicative carry")
    print(json.dumps({"balance_hostiles": rows, "same_band": {"start": 4067, "end": 2051,
        "odd_steps": 29, "halvings": total_depth, "intermediate_band_visit": 2429},
        "log_clock_repair": "3 log C < sum 1/m <= 4 log C; natural logarithms",
        "scope": "exact controls; no general first-return counterexample claimed"}, indent=2))


if __name__ == "__main__":
    main()
