#!/usr/bin/env python3
"""Two literal fixed-chart collisions; no target-coverage search."""
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import json
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SPEC = spec_from_file_location('local_carry_check', HERE/'check.py')
carry = module_from_spec(SPEC)
SPEC.loader.exec_module(carry)
corelib = carry.corelib
INCOMING_COMMIT = 'bf945ba947'
INCOMING_PATH = '05-knowledge/results/session255_modular_collatz_20261111/cross_frontiers/PROOF.md'


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    core, cut_mask, left = (2, 4, 1), 1, 1
    cf = corelib.coefficients(core)
    cut = cf['cuts'][cut_mask]
    need(cf['h'] == 3 and cf['q']*left+cut['s'] == 15,
         'fixed chart coefficients')
    records = []
    for address, expected_g in ((37, 15), (39, 5)):
        rows = corelib.graph(core, left, address, cut_mask)
        state = carry.direct_state(1, address)
        formula = 3*state[0]+15*address+10+(address&1)
        forward = corelib.hp(rows)
        predecessor = carry.direct_h(rows)
        need(corelib.tournament(rows) and len(rows) == 11, 'actual named tournament')
        need(state[0] == expected_g, 'g value')
        need(formula == corelib.value(cf, cut, left, address)
             == forward == predecessor == 611, 'three exact count paths')
        records.append({'address': address, 'binary': bin(address)[2:],
                        'g': state[0], 'state': state, 'H': formula,
                        'rows': rows, 'vertex_count': len(rows)})
    difference = tuple(b-a for a,b in zip(records[0]['state'], records[1]['state']))
    need(difference == (-10, 2, 0, 0), 'mixed difference')
    need(any(x < 0 for x in difference) and any(x > 0 for x in difference),
         'neither ordering comparable')
    need(3*difference[0]+15*difference[1] == 0, 'exact signed cancellation')
    suffixes = []
    for digit, expected in ((0, (1279, 1315)), (1, (1181, 1211))):
        outputs = []
        addresses = []
        matrix = carry.digit_matrix(1, digit)
        for record in records:
            address = 2*record['address']+digit
            state = carry.apply(matrix, record['state'])
            need(state == carry.direct_state(1, address), 'literal digit transition')
            count = corelib.value(cf, cut, left, address)
            need(count == 3*state[0]+15*state[1]+10+state[3], 'suffix count readout')
            addresses.append(address)
            outputs.append(count)
        need(tuple(outputs) == expected, 'suffix separates the equal counts')
        suffixes.append({'digit': digit, 'addresses': addresses, 'counts': outputs})
    incoming = subprocess.check_output(
        ['git', 'show', f'{INCOMING_COMMIT}:{INCOMING_PATH}'], cwd=ROOT)
    dependencies = [HERE/'check.py', carry.HELPER]
    print(json.dumps({
        'status': 'PASS',
        'scope': 'Two literal eleven-vertex graphs; four exact suffix readouts',
        'named_constructor': {'core_rows': core, 'outmask': cut_mask, 'left_address': left,
                              'left_names': [0], 'core_names': [1, 2, 3],
                              'right_names': list(range(4, 10)), 'apex_name': 10},
        'graphs': records, 'carry_difference': difference, 'suffixes': suffixes,
        'incoming_source': {'commit': INCOMING_COMMIT, 'path': INCOMING_PATH,
                            'sha256': sha256(incoming).hexdigest()},
        'source_sha256': {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
                          for path in dependencies},
        'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'limits': ['No inference of rooted O/Z/J-frame incomparability',
                   'No impossibility claim for richer signed compilers',
                   'No target coverage conclusion']
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
