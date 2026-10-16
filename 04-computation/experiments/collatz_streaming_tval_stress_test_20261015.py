#!/usr/bin/env python3
"""Stress test of the on-the-fly parity routine tval of collatz_census_extension_2p40_20261015.py (8-step jumps,
single uint64 steps above 2^64/26, exact two-limb continuation above the uint64 guard) against Python-integer
trajectories.  A table of parity bits below 2^22 is built; then tval(P, n, 2^22, A, Cc, jguard, guard) is compared with
the true parity of tau(n) for (i) 4,000 random n < 2^60 with the jump guard and the uint64 guard both forced low
(2^23), so that almost every orbit runs through the single-step and two-limb paths, (ii) 4,000 random n < 2^60 with
the production guards, and (iii) path-record starters whose T-orbits exceed 2^64, their neighbours and doubles.
Session opus-2026-10-15-S37.  Usage: python collatz_streaming_tval_stress_test_20261015.py out.json
"""
import sys, os, json, random
import importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('c40', os.path.join(HERE, 'collatz_census_extension_2p40_20261015.py'))
c40 = importlib.util.module_from_spec(spec)
_argv = sys.argv
sys.argv = ['x', '30', '26']
spec.loader.exec_module(c40)
sys.argv = _argv
OUT = sys.argv[1] if len(sys.argv) > 1 else 'tval_stress.json'

A, Cc = c40.jump_tables()
P, _ = c40.ext.parity_bits(22)
low = 1 << 22
rng = random.Random(7)
res = {}
cases = [('forced_low_guards', np.uint64(1 << 23), np.uint64(1 << 23)), ('production_guards', c40.JGUARD, c40.GUARD)]
for name, jg, g in cases:
    bad = 0
    for _ in range(4000):
        n = rng.randrange(low, 1 << 60)
        if c40.tval(P, n, low, A, Cc, jg, g) != 1 - 2 * c40.ext.tau_parity(n):
            bad += 1
    res[name] = dict(tested=4000, mismatches=bad)
    print(f'{name}: 4000 random n < 2^60, {bad} mismatches', flush=True)
starters = [23035537407, 45871962271, 51739336447, 59152641055, 59436135663, 70141259775, 77566362559, 110243094271,
            204430613247, 231913730799, 272025660543, 446559217279, 567839862631, 871673828443]
bad = 0
tested = 0
for n0 in starters:
    for n in (n0, n0 + 2, 2 * n0):
        for jg, g in ((c40.JGUARD, c40.GUARD), (np.uint64(1 << 23), np.uint64(1 << 23))):
            tested += 1
            if c40.tval(P, n, low, A, Cc, jg, g) != 1 - 2 * c40.ext.tau_parity(n):
                bad += 1
res['path_record_starters'] = dict(tested=tested, mismatches=bad)
print(f'path-record starters (orbits above 2^64), neighbours and doubles, both guard settings: {tested} tested, {bad} mismatches',
      flush=True)
json.dump(res, open(OUT, 'w'), indent=1)
