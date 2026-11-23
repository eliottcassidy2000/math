"""Classify fixed exterior forks after two vertices expand independently."""
from collections import Counter, defaultdict
from functools import cache
from itertools import combinations
from math import comb, gcd
from sys import stdin

from physical_navigation import direct, flipped


@cache
def coefficients(rows, a, b):
    """Count exterior skeletons by their numbers of nonempty visits to each module."""
    n = len(rows)
    others = [i for i in range(n) if i not in (a, b)]
    pair = (a,b) if rows[a] & (1<<b) else (b,a)
    gap_words = ((), (a,), (b,), pair)
    # Each entry ends immediately after an exterior vertex; -1 is the empty prefix.
    stages = {(0,-1): Counter({(0,0): 1})}
    for step in range(len(others)+1):
        following = defaultdict(Counter)
        for (used,previous), ways in stages.items():
            nexts = [v for v in others if not used & (1<<v)] if step<len(others) else [-1]
            for v in nexts:
                for word in gap_words:
                    segment = (() if previous==-1 else (previous,)) + word + (() if v==-1 else (v,))
                    if any(not rows[i] & (1<<j) for i,j in zip(segment,segment[1:])):
                        continue
                    state = (used | (1<<v),v) if v>=0 else (used,-1)
                    dest = following[state]
                    for (i,j),count in ways.items():
                        dest[i+(a in word),j+(b in word)] += count
        stages = following
    return {key:value for counts in stages.values() for key,value in counts.items() if min(key)>0}


def progression(profile, rate_a, rate_b):
    """Exponential coefficients for widths (1 + rate_a*t, 1 + rate_b*t)."""
    if not rate_a:
        return profile
    out = Counter()
    for (i,j),value in profile.items():
        for x in range(1,i+1):
            for y in range(1,j+1):
                sign = (-1)**(i+j-x-y)
                out[x**rate_a*y**rate_b] += value*sign*comb(i,x)*comb(j,y)*x*y
    return {key:v for key,v in out.items() if v}


def matches(center, neighbors, rates):
    c = progression(center,*rates)
    if any(v%3 for v in c.values()):
        return False
    third = {key:v//3 for key,v in c.items()}
    low,high = c.copy(),c.copy()
    unit = 1 if rates[0] else (1,1)
    low[unit] = low.get(unit,0)-2
    high[unit] = high.get(unit,0)+2
    targets = [dict((k,v) for k,v in t.items() if v) for t in (third,low,high)]
    return all(any(progression(profile,*rates)==t for profile in group)
               for group,t in zip(neighbors,targets))


def main(stream):
    # At orders <=10, each expanded port occurs in at most nine skeleton slots.
    # For the remaining ratios of positive widths, the exponential keys cannot collide.
    rates = [(0,0)]+[(a,b) for a in range(1,4) for b in range(1,4) if gcd(a,b)==1]
    tally = Counter()
    bases = pairs = growing = 0
    for line in stream:
        if not line.startswith("base "):
            continue
        mid, *rows = map(int, line.split()[1:])
        rows = tuple(rows)
        assert len(rows)<=10
        bases += 1
        wanted = (mid//3,mid-2,mid+2)
        candidates = defaultdict(list)
        for edge in combinations(range(len(rows)),2):
            neighbor = flipped(rows,*edge)
            h = direct(neighbor)
            if h in wanted:
                candidates[h].append((edge,neighbor))
        for a,b in combinations(range(len(rows)),2):
            groups = [[r for e,r in candidates[h] if a not in e and b not in e] for h in wanted]
            if not all(groups):
                continue
            pairs += 1
            c = coefficients(rows,a,b)
            assert c[1,1]==mid
            if len(c)==1:
                continue
            growing += 1
            gs = [[coefficients(r,a,b) for r in group] for group in groups]
            for rates_ab in rates:
                if matches(c,gs,rates_ab):
                    tally[len(rows),rates_ab]+=1
                    print("solution",len(rows),a,b,rates_ab,c,rows)
        coefficients.cache_clear()
    print("base centers",bases,"eligible port pairs",pairs,"growing",growing,"successful",sorted(tally.items()))


def positive_control():
    # Adjoin a fresh ordinal source to the certified nine-vertex fork. Growing that
    # source and the original slot independently must preserve all three reversals.
    center = (494,508,504,368,481,448,128,264,64,511)
    c = coefficients(center,0,9)
    groups = [[coefficients(flipped(center,*e),0,9)] for e in ((6,7),(3,5),(7,8))]
    assert c=={(1,1):87,(2,1):63}
    assert all(matches(c,groups,r) for r in ((0,0),(1,1),(3,2)))


if __name__=="__main__":
    positive_control()
    main(stdin)
