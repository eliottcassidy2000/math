#!/usr/bin/env python3
"""Halting opportunities of the eleven rule systems: at every configuration where a halting rule's branch is decided,
record the deciding quantity, and look for residue classes that are never visited (candidate invariants) (S56).

For each machine a function `opp(p)` returns (key, value) when p is a configuration whose next rule is chosen among
branches that include Halt, with `value` the quantity whose exact value decides halting; or None.
Usage: python halting_opportunities_20261025.py [macro_steps]
"""
import sys, json, os
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
from deepmind_cryptids_20261025 import step, START, HALT, MACHINES


def opp(nm, p):
    if nm == 'A1':
        a, b, c, d = p
        if c == 0 and a == 0 and d >= 2 and b < d - 1:
            return ('d-b (halt at 2)', d - b)
    if nm == 'A2':
        a, b, w = p
        if a == 0 and b % 2 == 0 and w <= 3:
            return ('w at even-b end (halt at 0)', w)
        if a == 0 and b % 2 == 1:
            return ('w at odd-b, a=0 (halt iff next even state has w=0, i.e. w=2)', min(w, 9))
    if nm == 'A3':
        a, c, f = p
        K, M = a // 2, c // 2
        if K >= 2 and M >= 2:
            return ('K-M (halt at 0)', K - M)
        if K == 1:
            return ('M when a=2 (halt at 1,2)', min(M, 9))
    if nm == 'A4':
        a, b = p
        if b <= 2 * a + 2:
            return ('b at end of cascade (halt at 2)', b if b <= 9 else 9)
    if nm == 'A5':
        a, z, w = p
        if a == 0 and z % 2 == 0 and w <= 3:
            return ('w at even-z end (halt at 0)', w)
    if nm == 'A6':
        u, v, w = p
        if u >= 1 and w == 1:
            return ('u when w=1 (halt at 1)', min(u, 9))
    if nm == 'A7':
        a, k, c = p
        if k == 1 and c % 2 == 1:
            return ('c when k=1, c odd (halt at 1)', min(c, 11))
    if nm == 'A8':
        a, b = p
        if b % 2 == 0 and a >= 1:
            k = (b - 6) // 2
            return ('k-(a-1) (halt at 4,5)', k - (a - 1))
        if b % 2 == 1:
            k = (b - 1) // 2
            return ('a-k (halt at 0,1)', a - k)
    if nm == 'A9':
        a, b, c = p
        if b >= 1 and c >= 1:
            return ('b-c at loop entry (halt at 0, and at 1 via (.,1,0))', b - c)
        if b == 0:
            return ('c when b=0 (halt at 0)', min(c, 9))
    if nm == 'A10':
        a, b, c = p
        if b == 0 and c >= 1:
            return ('c when b=0 (halt at 3)', min(c, 12))
    if nm == 'A11':
        a, b, c, d = p
        if c == 0 and d == 6:
            b, d = b + 3, 3
        if a == 0 and c in (0, 1) and d == 3:
            return ('b parity at a=0, c in {0,1}, d=3 (halt: c=0 & b odd, or c=1 & b even)', (c, b % 2))
    return None


def residues(vals, mods=(2, 3, 4, 5, 6, 8, 12)):
    out = {}
    for m in mods:
        out[m] = sorted(Counter(v % m for v in vals).items())
    return out


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    report = {}
    for nm in MACHINES:
        p = START[nm]
        recs = defaultdict(list)
        for i in range(N):
            if p == HALT:
                break
            o = opp(nm, p)
            if o is not None:
                recs[o[0]].append(o[1])
            p = step(nm, p)
        rep = {}
        for key, vals in recs.items():
            small = Counter(v for v in vals if isinstance(v, tuple) or abs(v) <= 12)
            entry = {'count': len(vals), 'small_values': sorted(small.items(), key=lambda t: str(t[0]))}
            if vals and not isinstance(vals[0], tuple):
                entry['residues'] = residues([v for v in vals])
            rep[key] = entry
        report[nm] = rep
        print(f'== {nm}')
        for key, e in rep.items():
            print(f'   {key}: n={e["count"]}, small={e["small_values"][:14]}')
            if 'residues' in e:
                print(f'      residues: ' + '; '.join(f'mod {m}: {r}' for m, r in list(e['residues'].items())[:5]))
    json.dump(report, open(os.path.join(HERE, 'halting_opportunities.json'), 'w'), indent=1, default=str)
