#!/usr/bin/env python3
"""Equal complete assembly signatures do not determine a later named flip."""
import hashlib
import json
from pathlib import Path

from audit_component import signature as permutation_signature
from component_inverse import attach, direct, require, strong


def main():
    inputs = ((66, 68, 121, 99, 107, 67, 0),
              (66, 68, 121, 67, 107, 75, 0))
    ports = 84
    common = (17, 448, 259, 217, 112, 17, 9)
    expected = ((27, 576, 349, 223, 114, 27, 9),
                (17, 448, 269, 167, 84, 17, 5))
    records = []
    for rows, out_sigma in zip(inputs, expected):
        require(direct(rows, ports) == common, 'input subset signature')
        require(permutation_signature(rows, ports) == common, 'input permutation signature')
        require(rows[4] >> 5 & 1 and not rows[5] >> 4 & 1, 'wrong action orientation')
        after = list(rows)
        after[4] ^= 1 << 5
        after[5] ^= 1 << 4
        after = tuple(after)
        require(direct(after, ports) == out_sigma, 'output subset signature')
        require(permutation_signature(after, ports) == out_sigma, 'output permutation signature')
        clocks = []
        for clock in (1, 2):
            pair = []
            for graph, sig in ((rows, common), (after, out_sigma)):
                full = attach(graph, ports, clock)
                require(strong(full), 'full graph not strong')
                observed = direct(full, 1 << (len(full) - 1))
                require(observed[0] == sig[0] - sig[4] + (1 << clock) * sig[2], 'full count')
                require(observed[2] - observed[0] == sig[4], 'physical defect')
                pair.append({'H': observed[0], 'F_root': observed[2], 'defect': sig[4]})
            clocks.append({'clock': clock, 'before_after': pair})
        records.append({'before_rows': rows, 'after_rows': after, 'ports': ports,
                        'before_sigma': common, 'after_sigma': out_sigma, 'clocks': clocks})
    receipt = {'scope': 'FINITE-EXACT explicit hostile; proves seven-coordinate quotient does not determine named flip transport',
               'action': [4, 5], 'records': records,
               'all_clock_output_H_difference': '80*2^clock-20',
               'output_defect_difference': 30}
    content = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    Path(__file__).with_suffix('.json').write_text(content)
    print(json.dumps({'PASS': True, 'sha256': hashlib.sha256(content.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
