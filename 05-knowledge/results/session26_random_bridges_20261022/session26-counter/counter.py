"""Exact Hamilton-path counter for a marked transitive-deletion presentation.
Q[u][v] == 1 iff core u -> core v. Word masks use bit u for core u -> base.
Only Python standard library. No graph canonization or minimum deletion search.
"""
from itertools import permutations
from functools import lru_cache


def validate(Q, word):
    Q = tuple(tuple(row) for row in Q)
    word = tuple(word)
    k = len(Q)
    if any(len(row) != k for row in Q):
        raise ValueError('core matrix must be square')
    if any(type(x) not in (int, bool) or x not in (0, 1) for row in Q for x in row):
        raise ValueError('core entries must be bits')
    if any(Q[i][i] for i in range(k)) or any(Q[i][j]+Q[j][i] != 1 for i in range(k) for j in range(i)):
        raise ValueError('core matrix must describe a loop-free tournament')
    if any(type(x) is not int or not 0 <= x < (1 << k) for x in word):
        raise ValueError('attachment mask outside 0..2^k-1')
    return Q, word


@lru_cache(None)
def transitions(k):
    """State digits: 0 empty, 1 occupied good exit, 2 occupied bad exit.
    Terminal slot is a separate occupancy bit. Row entries are destinations
    for each allowed slot choice, retaining multiplicity if destinations agree.
    """
    powers = [3**j for j in range(k)]
    t = 3**k
    labels = []
    for mask in range(1 << k):
        rows = []
        for state in range(2*t):
            digits = [(state//p) % 3 for p in powers]
            destinations = []
            for slot in range(k):
                if slot and digits[slot] == 0 and not (mask >> (slot-1) & 1):
                    continue
                new_digit = 2 if (mask >> slot & 1) else 1
                destinations.append(state + (new_digit-digits[slot])*powers[slot])
            if state >= t or (mask >> (k-1) & 1):
                destinations.append(state if state >= t else state+t)
            rows.append(tuple(destinations))
        labels.append(tuple(rows))
    return tuple(labels)


def count_hamilton(Q, word):
    """Count directed Hamilton vertex words, including H(empty graph)=1."""
    Q, word = validate(Q, word)
    k = len(Q)
    if k == 0:
        return 1
    table = transitions(k)
    powers = [3**j for j in range(k)]
    size = 2*3**k
    answer = 0
    for sigma in permutations(range(k)):
        counts = {0: 1}
        for mask in word:
            local_mask = sum(((mask >> u) & 1) << j for j, u in enumerate(sigma))
            rows = table[local_mask]
            new = {}
            for state, weight in counts.items():
                for target in rows[state]:
                    new[target] = new.get(target, 0) + weight
            counts = new
        for state, weight in counts.items():
            good = True
            for j in range(k):
                digit = (state//powers[j]) % 3
                if digit == 2 or (j and digit == 0 and not Q[sigma[j-1]][sigma[j]]):
                    good = False
                    break
            if good:
                answer += weight
    return answer


if __name__ == '__main__':
    import argparse, json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='JSON file containing core matrix Q and integer-mask word')
    args = parser.parse_args()
    with open(args.input) as f:
        data = json.load(f)
    print(count_hamilton(data['Q'], data['word']))
