"""Exact height inversion for two parity histories of the shortened 3x-1 map."""
from argparse import ArgumentParser
from collections import defaultdict
from dataclasses import dataclass
from math import isqrt

CHARTS = {4: (2465, 2293760, 1344270, 712944, 1740, 986, 423385),
          8: (37265, 144834560, 84839910, 57602544, 24360, 14906, 33838015)}
TARGET = (1 << 60) * ((1 << 61) - 1) - 1
CYCLES = ((1,), (5, 7, 10), (17, 25, 37, 55, 82, 41, 61, 91, 136, 68, 34))


def step(x):
    return (3*x-1)//2 if x%2 else x//2


def read(x):
    result = 0
    while x:
        x, last = divmod(x, 2)
        if last == 0:
            result += x + x%2
    return result


def height(chart, x, y):
    h, q, c, t, e, r, k = chart
    return k+q*x*y+t*x+c*y+r*(x%2)+e*(y%2)+h*(read(x)+read(y))


@dataclass(frozen=True)
class Branch:
    # Inputs B+A*z, outputs D+C*z for all nonnegative z.
    A: int
    B: int
    C: int
    D: int


def children(branch):
    A, B, C, D = branch.A, branch.B, branch.C, branch.D
    for bit in range(2):
        value = D+C*bit
        yield Branch(2*A, B+A*bit, C*(3 if value%2 else 1), step(value))


def controlled(x):
    """Choose the signed odd step which makes 3*x +/- 1 divisible by 4."""
    if x%2==0:
        return x//2
    return (3*x+(1 if x%4==1 else -1))//4


def controlled_children(branch):
    A, B, C, D = branch.A, branch.B, branch.C, branch.D
    even_bit = D%2
    yield Branch(2*A, B+A*even_bit, C, (D+C*even_bit)//2)
    for bits in (1-even_bit, 3-even_bit):
        yield Branch(4*A,B+A*bits,3*C,controlled(D+C*bits))


def chart_join(chart, target, left, right, offset):
    """Return (envelope survival, post-congruence survival, exact pair)."""
    # The ternary offsets are canonical, 0 <= D < C.
    # The caller has matched D modulo the smaller power of three.
    # Right output is w+offset. Start at the first common permitted w.
    lower = max(left.D, right.D-offset)
    period = max(left.C, right.C)
    residue = left.D if left.C >= right.C else right.D-offset
    w = lower + (residue-lower)%period
    x = left.B + left.A*((w-left.D)//left.C)
    y = right.B + right.A*((w+offset-right.D)//right.C)
    a, b = left.A*(period//left.C), right.A*(period//right.C)
    h, q, c, t, e, r, k = chart
    rem = target - (k+q*x*y+t*x+c*y)
    if rem < 0:
        return 0, 0, None
    slope = t*a+c*b+q*(x*b+y*a)
    quad = q*a*b
    z = (isqrt(slope*slope+4*quad*rem)-slope)//(2*quad)
    x, y = x+a*z, y+b*z
    gap = rem-z*(slope+quad*z)-r*(x%2)-e*(y%2)
    if min(x, y) <= 0 or not 0 <= gap <= h*(x+y):
        return 0, 0, None
    if gap%(2*h):
        return 1, 0, None
    return 1, 1, (x,y) if height(chart,x,y)==target else None


def matched_pairs(left, right, offset):
    """Enumerate exactly the compatible ordered ternary cylinders."""
    def half(small, large, strict, shift):
        buckets = defaultdict(lambda: defaultdict(list))
        for b in small:
            buckets[b.C][b.D].append(b)
        levels = sorted(buckets.items())
        for b in large:
            for c, group in levels:
                if c > b.C or (strict and c == b.C):
                    break
                for a in group.get((b.D+shift)%c, ()):
                    yield a, b
    yield from half(left, right, False, -offset)
    for r, l in half(right, left, True, offset):
        yield l, r


def inventory(pairs, chart, target, offset=0):
    joined = band = congruent = 0
    hits = set()
    for left, right in pairs:
        joined += 1
        inside, residue, hit = chart_join(chart, target, left, right, offset)
        band += inside
        congruent += residue
        if hit is not None:
            hits.add(hit)
    return joined, band, congruent, sorted(hits)


def free_pairs(depth, offset=0, expand=children):
    # Any meeting by depth D can be extended in lockstep until at least
    # one clock reads D. Deduplicate the depth (D,D) on the second pass.
    levels = [[Branch(1,0,1,0)]]
    for _ in range(depth):
        levels.append([b for parent in levels[-1] for b in expand(parent)])
    combined = [b for level in levels for b in level]
    shorter = combined[:-len(levels[-1])]
    if offset:
        yield from matched_pairs(combined, combined, offset)
    else:
        yield from matched_pairs(levels[-1], combined, 0)
        yield from matched_pairs(shorter, levels[-1], 0)


def cycle_branches(bound):
    out = set()
    for cycle in CYCLES:
        for start in cycle:
            depth = odds = 0
            current = start
            # Beyond the input cutoff, only the exact start phase can
            # remain. Every terminal phase then recurs within one period.
            last = max(bound, start).bit_length() + len(cycle)-1
            while depth <= last:
                A, C = 1 << depth, 3**odds
                quotient, B = divmod(start, A)
                out.add(Branch(A,B,C,current-C*quotient))
                depth += 1
                odds += current%2
                current = step(current)
    return sorted(out, key=lambda v:(v.C,v.D,v.A,v.B))


def cycle_pairs(chart, target, offset=0):
    _, _, c, t, _, _, k = chart
    return matched_pairs(cycle_branches((target-k)//t),
                         cycle_branches((target-k)//c), offset)


def controls():
    for cyc in CYCLES:
        assert tuple(map(step,cyc)) == cyc[1:]+cyc[:1]
    phases = sum(CYCLES, ())
    assert len({x%243 for x in phases}) == len(phases) == 15
    assert 1%81 == 82%81
    for x in (3,11,2):
        y=x
        for _ in range(6): y=step(y)
        assert y==1
    differences = [{abs(x-y) for x in CYCLES[i] for y in CYCLES[j]}
                   for i,j in ((0,1),(0,2),(1,2))]
    first, second, third = differences
    assert not first & (second | third)
    assert sorted(second & third) == [24,36,54,81]
    for n in range(2,300):
        assert 4*(controlled(n)-1) <= 3*(n-1)
    for parent in [Branch(1,0,1,0),*controlled_children(Branch(1,0,1,0))]:
        expanded = list(controlled_children(parent))
        for u in range(75):
            x = parent.A*u+parent.B
            selected = [v for v in expanded if x%v.A==v.B]
            assert len(selected)==1
            v = selected[0]
            assert controlled(parent.C*u+parent.D)==v.C*(x//v.A)+v.D
    for chart in CHARTS.values():
        # Both controlled branches have the common successor 5.
        x,y=10,7
        assert controlled(x)==controlled(y)==5
        assert (x,y) in inventory(free_pairs(1,0,controlled_children),chart,height(chart,x,y))[-1]
        # 29 -> 43 -> 64, whereas 128 -> 64. Both histories end there.
        assert step(step(29))==step(128)==64
        for x,y in ((29,128),(128,29)):
            target=height(chart,x,y)
            assert (x,y) in inventory(free_pairs(2),chart,target)[-1]
        # 37 and 61 are different starts on the long cycle. Also test
        # a lifted word of the phase 1 against a physical phase.
        for x,y in ((37,61),(61,37),(85,127),(127,85)):
            target=height(chart,x,y)
            assert (x,y) in inventory(cycle_pairs(chart,target),chart,target)[-1]
        for gap in (-24,-4,4,24):
            x = 1057
            y = step(step(x))+gap
            target = height(chart,x,y)
            assert (x,y) in inventory(free_pairs(2,gap),chart,target,gap)[-1]
            x = 127
            y = x+gap
            target = height(chart,x,y)
            assert (x,y) in inventory(cycle_pairs(chart,target,gap),chart,target,gap)[-1]
    print("cycle periods, phase separation, offsets, unequal-depth and height controls passed")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--depth",type=int,help="run all independent histories through this depth")
    parser.add_argument("--controlled",type=int,help="run contracting signed histories through this depth")
    parser.add_argument("--cycles",action="store_true",help="run cycle-shadow histories at every depth")
    parser.add_argument("--offsets",type=int,nargs="+",default=[0],help="right terminal value minus left")
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    if args.check:
        controls()
    for n, chart in CHARTS.items():
        if args.controlled is not None:
            if args.controlled<0: parser.error("depth must be nonnegative")
            print("controlled",args.controlled,n,*inventory(free_pairs(args.controlled,0,controlled_children),chart,TARGET),flush=True)
        for gap in args.offsets:
            if args.cycles:
                print("cycles",n,gap,*inventory(cycle_pairs(chart,TARGET,gap),chart,TARGET,gap),flush=True)
            if args.depth is not None:
                if args.depth<0: parser.error("depth must be nonnegative")
                print("free",args.depth,n,gap,*inventory(free_pairs(args.depth,gap),chart,TARGET,gap),flush=True)
