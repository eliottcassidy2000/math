#!/usr/bin/env python3
"""Exact fixed-core, arbitrary-speed attachment compiler; standard library only."""
from fractions import Fraction as F
from math import ceil, floor, gcd, lcm
from itertools import combinations
from pathlib import Path
import json

def need(b, message):
    if not b:
        raise RuntimeError(message)

def dist(x):
    return min(x % 1, (-x) % 1)

def merge(intervals):
    out = []
    for a, b in sorted(intervals):
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(b, out[-1][1]))
        else:
            out.append((a, b))
    return out

def safe_intervals(core, radius):
    """Reference core builder: intersect all closed combs in the original clock."""
    radius = F(radius)
    need(0 < radius <= F(1, 2), 'radius')
    intervals = [(F(0), F(1))]
    for v in sorted(set(core)):
        need(isinstance(v, int) and v > 0, 'positive integer speeds')
        nxt = []
        for left, right in intervals:
            for k in range(ceil(v*left-(1-radius)), floor(v*right-radius)+1):
                a = max(left, (k+radius)/v)
                b = min(right, (k+1-radius)/v)
                if a <= b:
                    nxt.append((a, b))
        intervals = merge(nxt)
    return intervals

def forward_gap(phase, radius):
    if radius <= phase <= 1-radius:
        return F(0)
    return radius-phase if phase < radius else 1+radius-phase

def attachment(intervals, radius, speed):
    """Return a native rational witness, or None. No loop depends on speed."""
    need(isinstance(speed, int) and speed > 0, 'positive attached speed')
    for left, right in intervals:
        gap = forward_gap((speed*left) % 1, radius)
        if speed*(right-left) >= gap:
            return left+gap/speed
    return None

def residue_threshold(intervals, radius, residue):
    """None denotes infinity. Caller uses residue modulo the returned period."""
    candidates = []
    for left, right in intervals:
        gap = forward_gap((residue*left) % 1, radius)
        if gap == 0:
            return F(0)
        if right > left:
            candidates.append(gap/(right-left))
    return min(candidates) if candidates else None

def period(intervals):
    return lcm(*(left.denominator for left, _ in intervals))

def finite_rejections(intervals, radius):
    """All rejected positive speeds; valid only if a positive-width component exists."""
    width = max((right-left for left,right in intervals), default=F(0))
    need(width > 0, 'positive width required for a finite global conductor')
    # Every phase gap is <= 2r. Equality is accepted.
    cutoff = ceil(2*radius/width)
    rejected = [v for v in range(1,cutoff) if attachment(intervals,radius,v) is None]
    return cutoff,rejected

def direct_pair_maximum(speeds):
    """Independent small-input maximum: all zero, cusp and pair-crossing clocks."""
    speeds = tuple(sorted(set(speeds)))
    denominators = {2*v for v in speeds}
    denominators |= {a+b for a,b in combinations(speeds, 2)}
    denominators |= {b-a for a,b in combinations(speeds, 2)}
    points = {F(k,d) for d in denominators for k in range(d+1)}
    return max(min(dist(v*t) for v in speeds) for t in points)

def run():
    queries = 0
    cores = 0
    infinite_banks = []
    for core in combinations(range(1, 13), 3):
        if gcd(*core) != 1:
            continue
        cores += 1
        cutoff,rejected=finite_rejections(safe_intervals(core,F(1,5)),F(1,5))
        need(not rejected, ('infinite-speed extension bank',core,rejected))
        infinite_banks.append(dict(core=core,tail_starts=cutoff))
        for radius in (F(1, 5), F(1, 4), F(1, 3), F(1, 2)):
            intervals = safe_intervals(core, radius)
            Q = period(intervals)
            for v in range(1, 65):
                t = attachment(intervals, radius, v)
                reference = safe_intervals((*core, v), radius)
                need((t is not None) == bool(reference), ('literal comb test', core, radius, v))
                threshold = residue_threshold(intervals, radius, v % Q)
                need((t is not None) == (threshold is not None and v >= threshold), 'residue threshold')
                if t is not None:
                    need(min(dist(s*t) for s in (*core, v)) >= radius, 'native witness')
                queries += 1
    # Exact equality, safe endpoints, and point components, all in the same clock.
    need(attachment([(F(0), F(1, 8))], F(1, 4), 2) == F(1, 8), 'closed equality')
    point = safe_intervals((1, 2), F(1, 3))
    need(all(a == b for a,b in point), 'isolated maxima')
    need(attachment(point, F(1, 3), 3) is None, 'point obstruction')
    need(attachment(point, F(1, 3), 4) is not None, 'point survives')
    # Local-germ and independent-arrival errors are not admitted by this compiler.
    need(direct_pair_maximum((1, 6, 7)) == F(4, 13), 'germ hostile maximum')
    need(attachment(safe_intervals((1,6), F(3,8)), F(3,8), 7) is None, 'germ hostile')
    interval12 = safe_intervals((1, 2), F(1, 4))
    need(attachment(interval12, F(1,4), 3) is not None, 'individual attachment 3')
    need(attachment(interval12, F(1,4), 12) is not None, 'individual attachment 12')
    need(attachment(safe_intervals((1,2,3), F(1,4)), F(1,4), 12) is None, 'joint arrivals hostile')
    need(direct_pair_maximum((1,2,3,12)) == F(3,13), 'joint hostile maximum')
    # The prior packet's first good period is 27, beyond its deliberately finite q scan.
    core = tuple(26*j for j in range(1, 13));radius = F(1,14)
    intervals = safe_intervals(core, radius)
    rows = []
    for v in (339, 10**100+39, 10**1000+123):
        t = attachment(intervals, radius, v)
        need(t is not None and min(dist(s*t) for s in (*core,v)) >= radius, 'large exact input')
        rows.append(dict(speed_digits=len(str(v)), speed_suffix=str(v)[-20:], witness=str(t),
                         minimum_distance=str(min(dist(s*t) for s in (*core,v)))))
    # In a fixed residue, coverage is upward closed even though general v coverage need not be.
    c = (1,2);rad=F(1,3);ints=safe_intervals(c,rad)
    need(attachment(ints,rad,4) is not None and attachment(ints,rad,6) is None,
         'unrestricted monotonicity must not be asserted')
    out = dict(status='FINITE-EXACT controls for proved parametric attachment compiler',
               primitive_triple_cores=cores, literal_queries=queries,
               radii=['1/5','1/4','1/3','1/2'], attached_speeds='1..64',
               large_core=list(core), radius=str(radius), components=len(intervals),
               phase_period=period(intervals), large_inputs=rows,
               all_speed_core_banks=infinite_banks,
               hostile_germ='(1,6)+7: 4/13, not 3/8',
               hostile_independent_arrivals='(1,2)+3+12: 3/13, not 1/4')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('large_inputs','all_speed_core_banks')}, indent=2))
    print('Verified all positive attached speeds for',len(infinite_banks),'fixed cores at radius 1/5')

if __name__ == '__main__':
    run()
