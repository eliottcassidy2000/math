"""Evaluate the three grade-three divisor fibers above the support {3,7,11}."""
from collections import Counter
from itertools import accumulate, product
from math import comb, gcd
from pathlib import Path
import sys

PACKET = Path(__file__).resolve().parent
sys.path.insert(0, str(PACKET.parent / "chanpak_gap_star_ports_20261122"))
from replay import subpaths


def strong_blocks():
    """For the required heights, only orders 3 through 8 can contribute."""
    data = PACKET.parents[2] / "session212" / "cycle-levels"
    wanted = {3, 11, 49, 77}
    blocks = {height: [] for height in wanted}
    for n in range(3, 9):
        for line in (data / f"classes{n}.rows").read_text().splitlines():
            rows = tuple(map(int, line.split()))
            scores = sorted(row.bit_count() for row in rows)
            # A proper tournament subset with the minimum possible total score
            # is a sink set. Any such cut makes a lowest-score prefix tight.
            if any(total == comb(k, 2) for k, total in
                   enumerate(accumulate(scores[:-1]), start=1)):
                continue
            height = subpaths(rows)[-1]
            if height in wanted:
                blocks[height].append(rows)
    return blocks


def ordered_words(target, counts):
    result = Counter()
    # Repetitions are intentional: equal-height positions can carry distinct classes.
    for heights in product(counts, repeat=3):
        if heights[0] * heights[1] * heights[2] != target:
            continue
        for a, b, c in product(*(counts[h].items() for h in heights)):
            result[a[0]+b[0]+c[0]] += a[1]*b[1]*c[1]
    return result


def seam_profile(rows):
    full = (1<<len(rows))-1
    paths = subpaths(rows)
    paths[0] = 1
    q = sum(v*paths[full^s] for s,v in enumerate(paths))
    def endpoints(arcs):
        return [sum(v*paths[full^s] for s,v in enumerate(subpaths(arcs,1<<u)))
                for u in range(len(rows))]
    reverse = tuple(full ^ row ^ (1<<i) for i,row in enumerate(rows))
    return q, endpoints(rows), endpoints(reverse)


def exit_census(target, blocks):
    profiles = {row:seam_profile(row) for cores in blocks.values() for row in cores}
    class_count = exits = movable = 0
    least = None
    for heights in product(blocks, repeat=3):
        if heights[0] * heights[1] * heights[2] != target:
            continue
        for a,b,c in product(*(blocks[h] for h in heights)):
            class_count += 1
            n_exits = 0
            for start in profiles[a][1]:
                for end in profiles[c][2]:
                    gain = start*profiles[b][0]*end
                    if gcd(gain,target) == 1:
                        n_exits += 1
                        new_height = target+gain
                        least = new_height if least is None else min(least,new_height)
            exits += n_exits
            movable += bool(n_exits)
    print("exits",target,"cores",class_count,"with a coprime neighbor",movable,
          "coprime arcs",exits,"least coprime neighbor",least)


def cyclic_checks():
    from replay import flipped
    certificates = {
        3: (3,1,2,1,17),
        7: (175,117,142,157,97081),
        11: (93027,71589,80046,88537,31575675105),
    }
    for prime, expected in certificates.items():
        rows = tuple(sum(1<<((v+j)%prime) for j in range(1,(prime+1)//2))
                     for v in range(prime))
        paths = subpaths(rows)
        paths[0] = 1
        rooted = subpaths(rows,1)
        full = (1<<prime)-1
        even, odd = (sum(paths[full^subset]*count for subset,count in enumerate(rooted)
                         if subset.bit_count()%2 == parity) for parity in (0,1))
        single = subpaths(flipped(rows,0,1))[-1]
        across = paths[-1]**2 + 4*even*odd
        got = (paths[-1],even,odd,single,across)
        assert got == expected
        assert paths[-1]%prime == 0 and paths[-1]%(prime*prime) != 0
        assert single%prime and across%prime
        if prime == 3:
            second = full << prime
            two_blocks = tuple(row|second for row in rows) + tuple(row<<prime for row in rows)
            assert subpaths(flipped(two_blocks,0,prime))[-1] == across
        print("cyclic",prime,got)


if __name__ == "__main__":
    blocks = strong_blocks()
    counts = {h:Counter(map(len,cores)) for h,cores in blocks.items()}
    for h, orders in sorted(counts.items()):
        print("strong", h, sorted(orders.items()))
    for target in (63*11, 147*11, 21*11**2):
        print("grade 3", target, sorted(ordered_words(target, counts).items()))
        exit_census(target,blocks)
    cyclic_checks()
