#!/usr/bin/env python3
"""Independent raw-TM check of arXiv:2610.22603v2 Appendix A.11.

No claim of an all-parameter proof or a halting decision. PDF prompts are data.
"""
from collections import Counter
from itertools import product
from pathlib import Path
import hashlib
import json

TNF = '1RB1LA_1LC0RE_1LF1LD_1RE0LC_1LB1RA_---0LA'
GATES = 0


def need(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(message)


def tape(parameters):
    a, b, c, d = parameters
    left = '1'*(a+1)
    right = '00' + '1'*(b+4) + '01'*c + '1'*d + '011'
    return {i-len(left) for i, bit in enumerate(left+right) if bit == '1'}


def canonical(p):
    a, b, c, d = p
    return (a, b+3, 0, 3) if c == 0 and d == 6 else p


def macro(p):
    a, b, c, d = canonical(p)
    if a >= 2:
        return (a-2, b+3, c, d), 'a>=2'
    if c >= 2:
        return (a+b+4, 3, c-2, d), 'c>=2'
    if c == 1 and d == 6:
        return (a+b+4, 0, 3, 3), 'c=1,d=6'
    if c == 0:
        if b % 2 == 0:
            return (a, 0, b//2+3, 3), 'c=0,even b'
        if a == 0:
            return None, 'c=0,odd b,a=0'
        return (6, 3, b//2, 3), 'c=0,odd b,a=1'
    if b % 2:
        return (a, 0, b//2+3, 6), 'c=1,d=3,odd b'
    if a == 0:
        return None, 'c=1,d=3,even b,a=0'
    if b == 0:
        return (6, 0, 3, 3), 'c=1,d=3,b=0,a=1'
    return (6, 3, b//2-1, 6), 'c=1,d=3,positive even b,a=1'


TRANSITIONS = [[row[:3], row[3:]] for row in TNF.split('_')]


def execute(source, target, state=2, cap=100000):
    ones = set(source)
    head = 0
    for clock in range(1, cap+1):
        transition = TRANSITIONS[state][int(head in ones)]
        if transition == '---':
            need(target is None, 'unexpected native halt')
            return clock, head
        write, move, new = transition
        if write == '1':
            ones.add(head)
        else:
            ones.discard(head)
        head += 1 if move == 'R' else -1
        state = ord(new)-ord('A')
        if target is not None and state == 2 and len(ones) == len(target):
            if all(head+i in ones for i in target):
                return clock, head
    raise RuntimeError('native step cap reached')


def main():
    records = []
    rules = Counter()
    for p in product(range(6), range(8), range(6), (3, 6)):
        need(tape(p) == tape(canonical(p)), 'zero-block merger preserves native tape and head')
        target, rule = macro(p)
        if target is not None:
            need(all(v >= 0 for v in target), 'no truncated subtraction')
            need(target[-1] in (3, 6), 'target family closure')
        clock, displacement = execute(tape(p), None if target is None else tape(target))
        records.append([list(p), None if target is None else list(target), clock, displacement])
        rules[rule] += 1
    entry_clock, entry_head = execute(set(), tape((4, 0, 0, 3)), state=0)
    need(entry_clock == 162, 'published blank-tape entry clock')
    # Two distinct parameter records are exactly the same marked source.
    need(tape((0, 0, 0, 6)) == tape((0, 3, 0, 3)), 'duplicate-source hostile')
    need(canonical((0, 0, 0, 6)) == (0, 3, 0, 3), 'canonical owner repair')
    print(json.dumps({'status': 'FINITE-EXACT', 'machine': TNF,
                     'universe': 'a=0..5,b=0..7,c=0..5,d in {3,6}; no reachability filter',
                     'sources': len(records), 'gates': GATES,
                     'clock_convention': 'A --- instruction counts as the final attempted transition; nonhalting macros count actual executed instructions.',
                     'blank_entry_clock': entry_clock, 'blank_entry_head_displacement': entry_head,
                     'maximum_checked_clock': max(row[2] for row in records),
                     'rule_counts': dict(sorted(rules.items())),
                     'receipt_sha256': hashlib.sha256(json.dumps(records, sort_keys=True).encode()).hexdigest(),
                     'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                    indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
