"""Exact cyclic contacts, Ramsey contraction, core extensions and Gaussian odd girths."""
from collections import Counter
from functools import cache
from itertools import combinations
from math import gcd


def unsigned(x):
    return min(x % 13, -x % 13)


def differences(points):
    return frozenset(unsigned(x-y) for x, y in combinations(points, 2))


def possible_contacts(step):
    occupied = {unsigned(step), unsigned(2*step)}
    return [s for k in range(4) for s in combinations(range(13), k)
            if len(differences(s)) == k*(k-1)//2 and not (differences(s) & occupied)]


def small_configurations():
    all_steps = set(range(1, 7))
    quartets = {s for s in combinations(range(13), 4) if len(differences(s)) == 6}
    singer = (0, 1, 3, 9)
    assert quartets == {tuple(sorted((a*x+b) % 13 for x in singer))
                        for a in range(1,13) for b in range(13)}
    assert len(quartets) == 52
    triples = {differences(s) for s in combinations(range(13), 3)
               if len(differences(s)) == 3}
    complementary = {s for s in triples if all_steps - s in triples}
    assert len(triples) == 8
    assert complementary == {frozenset((1,3,4)), frozenset((2,5,6))}
    assert {differences(s) for s in possible_contacts(1) if len(s) == 3} == {
        frozenset((3,4,6))}
    for t, allowed in enumerate(((5,6), (5,), (6,), (), (3,4)), 2):
        forbidden = {unsigned(x) for x in (1,2,t,2*t,1+t,1-t)}
        assert all_steps - forbidden == set(allowed)
    return len(quartets), sorted(tuple(sorted(v)) for v in triples)


def gaussian_rows(n, t):
    rows = [sum(1 << j for j in ((i-1) % n, (i+1) % n, n+i)) for i in range(n)]
    rows += [sum(1 << j for j in (n+(i-t) % n, n+(i+t) % n, i)) for i in range(n)]
    return rows


def cycles(rows):
    tally = Counter()
    for start in range(len(rows)):
        def extend(last, used, length):
            if rows[last] >> start & 1 and length >= 3:
                tally[length] += 1
            choices = rows[last] & ~used & ~((1 << start)-1)
            while choices:
                bit = choices & -choices
                choices -= bit
                extend(bit.bit_length()-1, used | bit, length+1)
        extend(start, 1 << start, 1)
    return {n: count//2 for n, count in sorted(tally.items())}


def check_ramsey():
    graph = [{(i+d) % 13 for d in (1, -1, 5, -5)} for i in range(13)]
    for x in range(13):
        for y in graph[x]:
            assert not graph[x] & graph[y]
    independent = lambda s: all(y not in graph[x] for x, y in combinations(s, 2))
    assert any(independent(s) for s in combinations(range(13), 4))
    assert not any(independent(s) for s in combinations(range(13), 5))


def exterior_exact_cover():
    # The low two four-bit fields are unsigned differences inside each panel.
    u = (3,4,5,6)
    v = (1,2,4,6)
    cross = (2,3,4,6,7,9,10,11)
    masks = {}
    for p in possible_contacts(1):
        for q in possible_contacts(5):
            product = {(b-a) % 13 for a in p for b in q}
            if len(product) != len(p)*len(q) or not product.issubset(cross):
                continue
            mask = sum(1 << u.index(x) for x in differences(p))
            mask += sum(1 << (4+v.index(x)) for x in differences(q))
            mask += sum(1 << (8+cross.index(x)) for x in product)
            if mask:
                size = (len(p), len(q))
                assert masks.get(mask, size) == size
                masks[mask] = size

    @cache
    def covers(free):
        if not free:
            return {()}
        pivot = free & -free
        return {tuple(sorted(tail + (size,)))
                for mask, size in masks.items() if mask & pivot and mask & free == mask
                for tail in covers(free ^ mask)}

    profiles = covers((1 << 16)-1)
    n_shared = Counter(sum(a*b != 0 for a,b in row) for row in profiles)
    assert (len(masks), len(profiles)) == (86, 76)
    assert n_shared == {2:3, 3:11, 4:18, 5:20, 6:12, 7:8, 8:4}
    return n_shared


if __name__ == "__main__":
    print("Singer blocks and triple difference sets:", small_configurations())
    check_ramsey()
    assert cycles(gaussian_rows(5, 2)) == {5:12, 6:10, 8:15, 9:20}
    counts = cycles(gaussian_rows(13, 5))
    assert set(counts) == set(range(7,27)) and counts[7] == counts[26] == 26
    print("Gaussian 13 native cycle counts:", counts)
    print("External contact size profiles by shared orbit count:", sorted(exterior_exact_cover().items()))
    for g, wanted in ((23,[221,233]), (25,[265,269])):
        possible = [a*a+b*b for a in range((g-3)//2, 0, -1) for b in (g-2-a,)
                    if gcd(a,b)==1 and 2*(a*a+b*b)<=550]
        assert possible == wanted
        print("Odd girth", g, "and allowed periods before embedding:", possible)
