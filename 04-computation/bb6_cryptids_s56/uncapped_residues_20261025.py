#!/usr/bin/env python3
"""Uncapped residues of the halting-deciding quantities (corrects the capped channels of halting_opportunities) (S56).

For each machine, at each configuration where a rule branch including Halt is chosen, record the exact deciding
quantity q (halting iff q equals a listed value) and tabulate q mod 2, 4, 8 and the counts of |q| <= 12.
Usage: python uncapped_residues_20261025.py [macro_steps]
"""
import sys, json, os
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
from deepmind_cryptids_20261025 import step, START, HALT, MACHINES


def deciding(nm, p):
    """Return list of (channel, q, halting_values)."""
    out = []
    if nm == 'A1':
        a, b, c, d = p
        if c == 0 and a == 0 and d >= 2 and b < d - 1:
            out.append(('A(0,b,0,d): d-b', d - b, (2,)))
    elif nm == 'A2':
        a, b, w = p
        if a == 0 and b % 2 == 1 and w >= 2:
            out.append(('D(0,odd,w>=2): w', w, (2,)))
        if a == 0 and b % 2 == 0 and w <= 3:
            out.append(('D(0,even,w<=3): w', w, (0,)))
    elif nm == 'A3':
        a, c, f = p
        K, M = a // 2, c // 2
        if K >= 2 and M >= 2:
            out.append(('K-M', K - M, (0,)))
        if K == 1:
            out.append(('M at a=2', M, (1, 2)))
        if M == 1:
            out.append(('K at c=2', K, (1,)))
    elif nm == 'A4':
        a, b = p
        if b <= 2 * a + 2:
            out.append(('b end of cascade', b, (2,)))
    elif nm == 'A5':
        a, z, w = p
        if a == 0 and z % 2 == 1 and w >= 2:
            out.append(('D1(0,odd,w>=2): w', w, (2,)))
        if a == 0 and z % 2 == 0 and w <= 3:
            out.append(('D1(0,even,w<=3): w', w, (0,)))
    elif nm == 'A6':
        u, v, w = p
        if u >= 1 and w >= 2:
            out.append(('loop start: w-(2u-1)', w - (2 * u - 1), (0,)))
    elif nm == 'A7':
        a, k, c = p
        if k == 1 and c % 2 == 1:
            out.append(('k=1, c odd: c', c, (1,)))
        if k == 1 and c % 2 == 0:
            out.append(('k=1, c even: c', c, (4,)))
    elif nm == 'A8':
        a, b = p
        if b % 2 == 0 and a >= 1:
            k = (b - 6) // 2
            out.append(('b even: k-(a-1)', k - (a - 1), (4, 5)))
        if b % 2 == 1:
            k = (b - 1) // 2
            out.append(('b odd: a-k', a - k, (0, 1)))
    elif nm == 'A9':
        a, b, c = p
        if b >= 1 and c >= 1:
            out.append(('loop entry: b-c', b - c, (0, 1)))
        if b == 0:
            out.append(('b=0: c', c, (0,)))
    elif nm == 'A10':
        a, b, c = p
        if b >= 1 and c >= 3 and b <= (c - 1) // 2:
            out.append(('loop ending at b=0: c-2b', c - 2 * b, (3,)))
        if b == 0:
            out.append(('b=0: c', c, (3,)))
    elif nm == 'A11':
        a, b, c, d = p
        if c == 0 and d == 6:
            b, d = b + 3, 3
        if a == 0 and c in (0, 1) and d == 3:
            out.append(('a=0,d=3,c in {0,1}: 2b+c', 2 * b + c, ('c=0,b odd', 'c=1,b even')))
    return out


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    rep = {}
    for nm in MACHINES:
        p = START[nm]
        stats = defaultdict(lambda: {'n': 0, 'm2': Counter(), 'm4': Counter(), 'm8': Counter(), 'small': Counter()})
        for i in range(N):
            if p == HALT:
                break
            for ch, q, hv in deciding(nm, p):
                s = stats[ch]
                s['n'] += 1
                s['m2'][q % 2] += 1
                s['m4'][q % 4] += 1
                s['m8'][q % 8] += 1
                if abs(q) <= 12:
                    s['small'][q] += 1
            p = step(nm, p)
        rep[nm] = {ch: {'n': s['n'], 'mod2': sorted(s['m2'].items()), 'mod4': sorted(s['m4'].items()),
                        'mod8': sorted(s['m8'].items()), 'small': sorted(s['small'].items())}
                   for ch, s in stats.items()}
        print('==', nm, flush=True)
        for ch, s in rep[nm].items():
            print(f'   {ch}: n={s["n"]} mod4={s["mod4"]} small={s["small"]}', flush=True)
    json.dump(rep, open(os.path.join(HERE, 'uncapped_residues.json'), 'w'), indent=1, default=str)
