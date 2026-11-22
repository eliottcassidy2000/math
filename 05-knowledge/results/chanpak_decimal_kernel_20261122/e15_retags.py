"""Finite exponent-61 congruence test on every retagging of three fixed blocks."""

from functools import reduce
from math import gcd

from e15_pair import A, D, V, EDGE, flipped, join, tagged_data


def compatible(chart, target, branch):
    h, q, c, t, e, r, k = chart
    for a in range(2):
        for b in range(2):
            if (branch == "right" and a) or (branch == "left" and b):
                continue
            origin = k + (t+r)*a + (c+e)*b + q*a*b
            increments = [2*h]
            if branch != "right":
                increments.append(2*(t+q*b))
            if branch != "left":
                increments.append(2*(c+q*a))
            if branch == "both":
                increments.append(4*q)
            if (target-origin) % reduce(gcd, increments) == 0:
                return True
    return False


def main():
    def profiles(block):
        return [tagged_data(block, s) for s in range(1 << len(block))]

    left, mid, other, right = (profiles(block)
                               for block in (A, D, flipped(D, *EDGE), V))
    # A and V are nonempty, hence the joined e,r do not read the middle tags.
    pairs = [(tag, d1, d2) for tag, (d1, d2) in enumerate(zip(mid, other))
             if d1[:4] == d2[:4]]
    legal = {}
    marked_placements = 0
    for left_tag, before in enumerate(left):
        for middle_tag, d1, d2 in pairs:
            a1, a2 = join(before, d1), join(before, d2)
            for right_tag, after in enumerate(right):
                u, v = join(a1, after), join(a2, after)
                if abs(u[-1]-v[-1]) == 2:
                    assert u[:-1] == v[:-1]
                    legal[min(u, v)] = (left_tag, middle_tag, right_tag)
                    marked_placements += 1
    target = 2**60*(2**61-1)-1
    assert len(pairs) == 4 and marked_placements == len(legal) == 48
    print("admissible middle tags", len(pairs))
    print("retaggings with permanent unit edge", marked_placements)
    print("distinct lower charts", len(legal))
    print("stationary slope range right", min(c[2] for c in legal), max(c[2] for c in legal))
    print("stationary slope range left", min(c[3] for c in legal), max(c[3] for c in legal))
    print("two-arm maximum order", 21+(target//next(iter(legal))[1]).bit_length())
    for branch in ("right", "left", "both"):
        passed = sorted(c for c in legal if compatible(c, target, branch))
        assert len(passed) == {"right": 2, "left": 30, "both": 30}[branch]
        print(branch, "pass", len(passed))
        if branch != "both":
            slope_index = 2 if branch == "right" else 3
            print("maximum order by slope", 20+(target//min(c[slope_index] for c in passed)).bit_length())
        if branch == "right":
            for c in passed:
                print("tags", legal[c], "h,c,e,k", c[0], c[2], c[4], c[-1])


if __name__ == "__main__":
    main()
