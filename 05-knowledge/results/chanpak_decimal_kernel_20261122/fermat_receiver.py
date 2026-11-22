"""Exact profiles of cyclic transitive-spine Farey receivers, and target inversion."""

from argparse import ArgumentParser
from functools import lru_cache
from itertools import permutations
from array import array
from math import gcd, isqrt

from e15_pair import A, D, EDGE, flipped, join, p, tagged_data

TARGET = (1 << 60)*((1 << 61)-1)-1


def nontags(n):
    assert n > 0 and n % 4 == 0
    numerator = (3*((1 << n)+1)-1)//5
    return tuple(bool(numerator & (1 << i)) for i in reversed(range(n)))


def literal_receiver(word):
    n = len(word)
    a, b = 1 << n, 2 << n
    rows = tuple(((a-1) ^ ((2 << i)-1)) | a for i in range(n))
    tags = sum(1 << i for i, bit in enumerate(word) if not bit) | a | b
    return rows+(b, a-1), tags


def beats(u, v):
    # Characters t,u represent tagged/untagged spine vertices at a run boundary.
    # A bilateral separator | bijects paths with ordered pairs of path slots.
    if "|" in (u, v):
        return True
    if u == "*":
        return v in "tab"
    if v == "*":
        return u == "u"
    return (u in "tu" and v == "a") or (u == "a" and v == "b") or (u == "b" and v in "tu")


def contact(v, wanted):
    return wanted is None or (v != "|" and (v in "tab") == wanted)


@lru_cache(None)
def run_machine(order, first, last):
    count = len(order)+1
    weights = tuple(3**j for j in range(count))
    states = 3**count
    transitions = []
    terminal = []
    for s in range(states):
        digits = [(s//w) % 3 for w in weights]
        accept = True
        for j, d in enumerate(digits):
            if d == 2:
                accept = False
            if d == 0:
                if j == 0:
                    accept &= contact(order[0], first)
                elif j == count-1:
                    accept &= contact(order[-1], last)
                else:
                    accept &= beats(order[j-1], order[j])
        terminal.append(accept)
        moves = []
        for symbol in ("t", "u"):
            destinations = []
            for j, (w, d) in enumerate(zip(weights, digits)):
                enter = contact(symbol, first) if j == 0 else beats(order[j-1], symbol)
                if not d and not enter:
                    continue
                exit_ok = contact(symbol, last) if j == count-1 else beats(symbol, order[j])
                destinations.append(s + w*((1 if exit_ok else 2)-d))
            moves.append(destinations)
        transitions.append(moves)
    return transitions, terminal


def count(word, markers, first=None, last=None):
    answer = 0
    for order in permutations(markers):
        transitions, terminal = run_machine(order, first, last)
        dp = [0]*len(terminal)
        dp[0] = 1
        for untagged in word:
            after = [0]*len(dp)
            for i, mass in enumerate(dp):
                if mass:
                    for dest in transitions[i][untagged]:
                        after[dest] += mass
            dp = after
        answer += sum(mass for mass, good in zip(dp, terminal) if good)
    return answer


def receiver(n):
    word = nontags(n)
    h = count(word, "ab")
    q = count(word, "ab|")
    c = count(word, "ab|", first=True)
    t = count(word, "ab|", last=False)
    e = count(word, "ab", last=True)
    r = count(word, "ab", first=False)
    k = count(word, "ab*")
    assert h == (1 << n)+1 and 5*r == 3*h-1
    data = h, q, c, t, e, r, k
    if n in (4, 8):
        assert data == tagged_data(*literal_receiver(word))
    return data


def address_height(chart, x, y):
    h, q, c, t, e, r, k = chart
    return k+q*x*y+c*y+t*x+e*(y & 1)+r*(x & 1)+h*(p(x)+p(y))


def single_runs(chart, target):
    """Includes the empty word and all 1^w 0^v positive addresses on both sides."""
    h, q, c, t, e, r, k = chart
    found = []
    if target < k:
        return found

    def lefts():
        yield 0, 0
        for bits in range(1, ((target-k)//t).bit_length()+1):
            for tail in range(bits):
                ones = (1 << (bits-tail))-1
                yield ones << tail, ((1 << tail)-1)*ones+int(tail > 0)

    for x, fx in lefts():
        rem = target-k-t*x-r*(x & 1)-h*fx
        if rem < 0:
            continue
        if not rem:
            found.append((x, 0))
        slope = c+q*x
        for tail in range((rem//slope).bit_length()):
            den = ((slope+h) << tail)-h
            num = rem-(h if tail else e)
            if num < den:
                continue
            ones, residual = divmod(num, den)
            if residual == 0 and (ones & (ones+1)) == 0:
                found.append((x, ones << tail))
    assert all(address_height(chart, x, y) == target for x, y in found)
    return found


def unrestricted_face(chart, target, right, cap):
    """Completely decode a one-sided face when its address progression fits cap."""
    h, _, c, t, e, r, k = chart
    slope, odd = (c, e) if right else (t, r)
    g = gcd(slope, h)
    mod = h//g
    residue_inverse = pow(slope//g, -1, mod)
    progressions = []
    for parity in range(2):
        rem = target-k-parity*odd
        if rem < 0 or rem % g:
            continue
        lo = (rem+slope+h-1)//(slope+h)
        hi = rem//slope
        res = (rem//g)*residue_inverse % mod
        first = lo+(res-lo) % mod
        step = mod
        if (first & 1) != parity:
            if not (step & 1):
                continue
            first += step
        if step & 1:
            step *= 2
        if first <= hi:
            progressions.append((first, hi, step, rem))
    size = sum((hi-first)//step+1 for first, hi, step, _ in progressions)
    if size > cap:
        return None
    return [y for first, hi, step, rem in progressions for y in range(first, hi+1, step)
            if slope*y+h*p(y) == rem]


def both_faces(chart, target, limit=3_300_000):
    """Try the smaller positive address first, in either orientation."""
    h, q, c, t, e, r, k = chart
    if k >= target:
        return []
    bound = isqrt((target-k)//q)
    if bound > limit:
        return None
    reads = array("I", [0])*(bound+1)
    for u in range(1, bound+1):
        half = u//2
        reads[u] = reads[half] if u & 1 else reads[half]+half+(half & 1)
    hits = []
    for swapped in (False, True):
        cx, cy, ex, ey = (c, t, e, r) if swapped else (t, c, r, e)
        for x in range(1, bound+1):
            slope = cy+q*x
            base = k+cx*x+ex*(x & 1)+h*reads[x]
            for parity in (0, 1):
                rem = target-base-parity*ey
                if rem < 0:
                    continue
                hi = rem//slope
                lo = max(x+int(swapped), (rem+slope+h-1)//(slope+h))
                if lo > hi:
                    continue
                if hi-lo < 20:
                    first = lo+((parity-lo) & 1)
                    step = 2
                else:
                    g = gcd(slope, h)
                    if rem % g:
                        continue
                    mod = h//g
                    residue = (rem//g)*pow(slope//g, -1, mod) % mod
                    first = lo+(residue-lo) % mod
                    if (first & 1) != parity:
                        if not (mod & 1):
                            continue
                        first += mod
                    step = 2*mod if mod & 1 else mod
                for y in range(first, hi+1, step):
                    if slope*y+h*p(y) == rem:
                        hits.append((y, x) if swapped else (x, y))
    assert all(address_height(chart, x, y) == target for x, y in hits)
    return hits


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--through-36", action="store_true",
                        help="also exhaust the 29,408,100 smaller positive addresses at n=36")
    args = parser.parse_args()
    pair_limit = 30_000_000 if args.through_36 else 3_300_000
    ad = [join(tagged_data(A, 12), tagged_data(mid, 35)) for mid in (D, flipped(D, *EDGE))]
    cleared = [[], []]
    large = []
    cleared_both = []
    for n in range(4, 113, 4):
        v = receiver(n)
        upper, lower = (join(a, v) for a in ad)
        assert upper[:-1] == lower[:-1] and upper[-1] == lower[-1]+2
        if n == 4:
            for x, y in ((0, 0), (3, 8), (12, 60), (0, 2), (7, 0)):
                assert (x, y) in single_runs(lower, address_height(lower, x, y))
            for x, y in ((1, 1), (11, 6), (15, 37)):
                assert (x, y) in both_faces(lower, address_height(lower, x, y))
            for right in (False, True):
                ans = address_height(lower, 0, 77) if right else address_height(lower, 77, 0)
                assert 77 in unrestricted_face(lower, ans, right, 100_000)
        pair_result = both_faces(lower, TARGET, pair_limit)
        if pair_result is not None:
            print("both nonempty", n, pair_result)
            if not pair_result:
                cleared_both.append(n)
        hits = single_runs(lower, TARGET)
        print("n", n, "lower intercept", lower[-1], "single-run hits", hits)
        assert not hits
        if lower[-1] > TARGET:
            large.append(n)
        for right in (False, True):
            result = unrestricted_face(lower, TARGET, right, cap=200_000)
            if result is not None:
                print("unrestricted", "right" if right else "left", n, result)
                if not result:
                    cleared[right].append(n)
    print("unrestricted empty positive pair", cleared_both)
    print("base above target", large)
    print("unrestricted empty left", cleared[0])
    print("unrestricted empty right", cleared[1])


if __name__ == "__main__":
    main()
