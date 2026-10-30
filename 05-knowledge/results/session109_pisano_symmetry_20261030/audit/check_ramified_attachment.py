#!/usr/bin/env python3
"""Independent, stdlib-only exact audit. Checks use exceptions, not assert."""
import json
import sys
from collections import Counter
from pathlib import Path


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def census(k):
    n = 5**k
    # Integral-basis coordinates x=a+b*phi. This implementation does not
    # import any producer modules or use their orbit or valuation functions.
    points = [(a, b) for a in range(n) for b in range(n)]
    def g(x):
        a, b = x
        return ((2*a+3*b) % n, (3*a+5*b) % n)
    def h(x):
        a, b = x
        return ((a+b) % n, (a+2*b) % n)
    def neg(x):
        return ((-x[0]) % n, (-x[1]) % n)
    # Build actual images pi^s R by repeated multiplication, avoiding a
    # valuation formula whose errors could mirror the claimed theorem.
    ideals = [set(points)]
    for _ in range(2*k):
        ideals.append({((-a+2*b) % n, (2*a+b) % n) for a,b in ideals[-1]})
    require(ideals[-1] == {(0, 0)}, 'last ideal')
    depth = {}
    for s, ideal in enumerate(ideals[:-1]):
        require(len(ideal) == 5**(2*k-s), 'ideal cardinality')
        for x in ideal - ideals[s+1]:
            depth[x] = 2*k-s
    orbit_id = {}
    orbits = []
    for x in points[1:]:
        if x in orbit_id:
            continue
        orbit = []
        y = x
        while y not in orbit_id:
            orbit_id[y] = len(orbits)
            orbit.append(y)
            y = g(y)
        require(y == x, 'permutation closes at initial point')
        r = depth[x]
        require(all(depth[z] == r for z in orbit), 'depth is constant')
        require(len(orbit) == 5**(r//2), 'predicted orbit length')
        orbits.append(orbit)
    antipode = {}
    for i, orbit in enumerate(orbits):
        j = orbit_id[neg(orbit[0])]
        require(j != i, 'antipodal collision')
        require(all(orbit_id[neg(x)] == j for x in orbit), 'antipodal orbit map')
        require(all(orbit_id[h(x)] == j for x in orbit), 'Q2 reverses orbit pair')
        antipode[i] = j
    shells = Counter(depth.values())
    for r in range(1, 2*k+1):
        require(shells[r] == 4*5**(r-1), 'shell size')
    require(len(orbits) == 2*(5**k-1), 'orbit count')
    return dict(k=k, n=n, orbits=orbits, orbit_id=orbit_id,
                antipode=antipode, shells=dict(sorted(shells.items())))


def attach(lower, upper):
    n = lower['n']
    downstairs = lower['orbit_id']
    above = {i: [] for i in range(len(lower['orbits']))}
    kernel = []
    for i, orbit in enumerate(upper['orbits']):
        images = {(a % n, b % n) for a,b in orbit}
        if images == {(0,0)}:
            kernel.append(i)
        else:
            j = downstairs[next(iter(images))]
            require(images == set(lower['orbits'][j]), 'full orbit image')
            above[j].append(i)
            require(len(orbit) == 5*len(images), 'degree five orbit map')
    require(all(len(v) == 5 for v in above.values()), 'five orbits above each')
    require(len(kernel) == 8, 'eight kernel orbits')
    require(all(upper['antipode'][i] in kernel for i in kernel), 'four kernel pairs')
    B = len(lower['orbits'])//2
    require(len(upper['orbits'])//2 == 5*B+4, 'attachment recurrence')
    # Downstairs orientation picks the lexicographically earlier orbit in
    # each antipodal pair. Lift it coherently, then flip exactly one pair.
    signs = {}
    for j, lifts in above.items():
        for i in lifts:
            signs[i] = j < lower['antipode'][j]
    for i in kernel:
        signs[i] = i < upper['antipode'][i]
    require(all(signs[i] != signs[j] for i,j in upper['antipode'].items()), 'valid signs')
    j = next(j for j in above if j < lower['antipode'][j])
    changed = above[j][0]
    signs[changed] = not signs[changed]
    opposite = upper['antipode'][changed]
    signs[opposite] = not signs[opposite]
    require(all(signs[i] != signs[t] for i,t in upper['antipode'].items()), 'incoherent is valid tournament')
    d = lower['orbits'][j][0]
    witnesses = []
    for i in above[j][:2]:
        x = next(x for x in upper['orbits'][i] if (x[0] % n, x[1] % n) == d)
        witnesses.append({'difference': x, 'sign': signs[i]})
    require(witnesses[0]['sign'] != witnesses[1]['sign'], 'incoherence witnessed')
    # Enumerate all 32 local patterns, verifying exactly two can descend,
    # and exactly one matches a specified downstairs orientation.
    patterns = [[bool(mask & (1 << bit)) for bit in range(5)] for mask in range(32)]
    coherent = [p for p in patterns if len(set(p)) == 1]
    require(len(coherent) == 2 and sum(all(p) for p in coherent) == 1, 'local attachment signs')
    # All sixteen kernel assignments, independently checking antipodality.
    kernel_pairs = [(i,upper['antipode'][i]) for i in kernel if i < upper['antipode'][i]]
    for mask in range(16):
        ks = {}
        for b,(i,t) in enumerate(kernel_pairs):
            ks[i] = bool(mask & (1 << b))
            ks[t] = not ks[i]
        require(len(ks) == 8 and all(ks[i] != ks[t] for i,t in kernel_pairs), 'kernel assignment')
    return {'lower_k': lower['k'], 'upper_k': upper['k'],
            'orbits_above_each_nonzero_orbit': 5, 'kernel_antipodal_pairs': 4,
            'fixed_downstairs_compatible_lifts': 16,
            'relative_bits_discarded_by_coherence': 4*B,
            'local_patterns': 32, 'coherent_local_patterns': 2,
            'incoherent_counterexample': {'downstairs_difference': d, 'lifts': witnesses}}


def main():
    levels = [census(k) for k in range(1,4)]
    report = {'status': 'PASS', 'implementation': 'independent integral-coordinate permutation and ideal-image census',
              'checks_use_assert': False,
              'levels': [{'k': d['k'], 'points': d['n']**2, 'nonzero_orbits': len(d['orbits']),
                          'sign_bits': len(d['orbits'])//2, 'shell_sizes': d['shells'],
                          'Q2_reverses_every_orbit_pair': True} for d in levels],
              'attachments': [attach(a,b) for a,b in zip(levels,levels[1:])]}
    rendered = json.dumps(report, indent=2) + '\n'
    if len(sys.argv)>1:
        Path(sys.argv[1]).write_text(rendered)
    else:
        print(rendered, end='')

if __name__ == '__main__':
    main()
