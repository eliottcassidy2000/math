#!/usr/bin/env python3
"""Recheck the self-contained inherited bank inventory; not a core census."""
import json
from pathlib import Path
from verify import need, profile


def main():
    obj = json.loads(Path(__file__).with_name('inventory135.json').read_text())
    records = obj['records']
    need(len(records) == 135, 'bank size')
    need([z['offset'] for z in records] == list(range(-115, 154, 2)), 'offset band')
    eligible = 0
    maximum_period = 0
    for z in records:
        p = profile(z['rows'], z['mask'])
        need(p == z['profile'], 'physical inventory interface')
        R = p['R']
        need(R > 0 and R & (R-1) == 0 and p['b'] == z['offset'], 'dyadic seed')
        r = R.bit_length()-1
        v = 0
        m = p['B']
        while m % 2 == 0:
            v += 1
            m //= 2
        need((v, m) == (z['two_adic_B'], z['odd_part_B']), 'B factorization')
        yes = r >= v and m > 1
        need(yes == z['infinitely_many_source_returns'], 'return eligibility')
        if yes:
            eligible += 1
            d = 1
            value = 2 % m
            while value != 1:
                value = 2*value % m
                d += 1
            need(d == z['return_period'] and d >= 2, 'least return period')
            need(((1 << (r+d))-R) % p['B'] == 0, 'first return exact congruence')
            maximum_period = max(maximum_period, d)
        else:
            need(z['return_period'] is None, 'no false period')
    need(eligible == obj['eligible_entries'] == 134, '134 eligible seeds')
    print(json.dumps(dict(status='PASS: finite exact 135-entry bank inventory',
                          entries=135, eligible=eligible, maximum_period=maximum_period,
                          excluded_offsets=[z['offset'] for z in records
                                            if not z['infinitely_many_source_returns']]),
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
