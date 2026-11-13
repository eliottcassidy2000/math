#!/usr/bin/env python3
"""Verify the closed research packet; --replay reruns all eight exact checks.

Standard library only. All gates survive optimized Python. Historical source
pins are checked against their Git revisions, not against evolving canon.
"""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

PACKET = Path(__file__).resolve().parent
REPO = PACKET.parents[2]
RUNS = [
    ('pairings/check.py', 'pairings/receipt.json', 'pairings/receipt-optimized.json'),
    ('metagraph/replay.py', 'metagraph/normal.json', 'metagraph/optimized.json'),
    ('sporadic/check.py', 'sporadic/check.json', 'sporadic/check.optimized.json'),
    ('audit/check_small.py', 'audit/small.normal.json', 'audit/small.optimized.json'),
    ('audit_sporadic/check.py', 'audit_sporadic/receipt.json', 'audit_sporadic/receipt-optimized.json'),
    ('sporadic/check_other_lanes.py', 'sporadic/review_other_lanes.json', 'sporadic/review_other_lanes.optimized.json'),
    ('metagraph/pairing_lift_check.py', 'metagraph/pairing_lift.normal.json', 'metagraph/pairing_lift.optimized.json'),
    ('audit/slot_profiles.py', 'audit/slot_profiles.normal.json', 'audit/slot_profiles.optimized.json'),
]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(data):
    return sha256(data).hexdigest()


def read_json(name):
    return json.loads((PACKET / name).read_text())


def check_hashes(mapping, prefix=''):
    for name, expected in mapping.items():
        need(digest((PACKET / prefix / name).read_bytes()) == expected,
             'stale source binding: ' + prefix + '/' + name)


def check_packet():
    manifest = read_json('MANIFEST.json')['files']
    paths = list(PACKET.rglob('*'))
    need(not any(p.is_symlink() for p in paths), 'symlink in closed packet')
    actual = {p.relative_to(PACKET).as_posix() for p in paths
              if p.is_file() and p.name != 'MANIFEST.json'}
    need(actual == set(manifest), ('closed file set changed', actual ^ set(manifest)))
    for name, row in manifest.items():
        data = (PACKET / name).read_bytes()
        need(len(data) == row['bytes'] and digest(data) == row['sha256'],
             'manifest mismatch: ' + name)

    for name, row in read_json('PROMOTIONS.json')['files'].items():
        text = (PACKET / name).read_text()
        need(digest(text.encode()) == row['after_sha256'], 'changed promoted proof: ' + name)
        for edit in reversed(row['edits']):
            need(text.count(edit['new']) == edit['count'], 'ambiguous promotion inverse: ' + name)
            text = text.replace(edit['new'], edit['old'])
        need(digest(text.encode()) == row['before_sha256'], 'review snapshot lost: ' + name)

    for code, normal, optimized in RUNS:
        need((PACKET / normal).read_bytes() == (PACKET / optimized).read_bytes(),
             'normal/optimized receipt mismatch: ' + code)
        json.loads((PACKET / normal).read_text())
    for code, receipt in [('pairings/check.py', 'pairings/receipt.json'),
                          ('metagraph/replay.py', 'metagraph/normal.json'),
                          ('metagraph/pairing_lift_check.py', 'metagraph/pairing_lift.normal.json')]:
        need(read_json(receipt)['source_sha256'] == digest((PACKET / code).read_bytes()),
             'receipt source mismatch: ' + code)
    audit = read_json('audit_sporadic/receipt.json')
    check_hashes(audit['source_hashes'], 'sporadic')
    need(audit['audit_code_sha256'] == digest((PACKET / 'audit_sporadic/check.py').read_bytes()),
         'sporadic audit code binding')
    check_hashes(read_json('sporadic/review_other_lanes.json')['source_sha256'])
    need(read_json('metagraph/pairing_lift.normal.json')['root_receipt_sha256'] ==
         digest((PACKET / 'audit/small.normal.json').read_bytes()), 'lift input receipt')
    slot = read_json('audit/slot_profiles.normal.json')
    need(slot['root_census_sha256'] == digest((PACKET / 'audit/small.normal.json').read_bytes()),
         'slot input census')
    need(slot['checker_sha256'] == digest((PACKET / 'audit/slot_profiles.py').read_bytes()),
         'slot source binding')

    for lane in ['metagraph', 'sporadic']:
        for line in (PACKET / lane / 'SHA256SUMS').read_text().splitlines():
            expected, name = line.split('  ', 1)
            need(digest((PACKET / lane / name).read_bytes()) == expected,
                 'final lane manifest mismatch: ' + lane + '/' + name)
    proof_hash = digest((PACKET / 'PAIRING_WIRING.md').read_bytes())
    need(proof_hash in (PACKET / 'metagraph/PAIRING_LIFT_REVIEW.md').read_text(),
         'final pairing proof is not bound by its review')
    rows = next(r for r in read_json('audit/small.normal.json')['rows'] if r['n'] == 5)
    displayed = re.findall(r'^\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|(yes|no)\|$',
                           (PACKET / 'PAIRING_WIRING.md').read_text(), re.M)
    metadata = [(str(r['code']), str(r['H']), str(r['aut']), str(r['degree']),
                 'yes' if r['self_converse'] else 'no') for r in rows['node_metadata']]
    need(displayed == metadata, 'canonical-code table does not match its stated gauge')
    svg = ET.fromstring((PACKET / 'pairing-obstruction.svg').read_bytes())
    figure_meta = json.loads(svg.find('{http://www.w3.org/2000/svg}metadata').text)
    check_hashes(figure_meta['source_sha256'])

    pins = read_json('SOURCE_PINS.json')['files']
    env = dict(os.environ, GIT_NO_LAZY_FETCH='1', GIT_TERMINAL_PROMPT='0')
    for name, row in pins.items():
        data = subprocess.check_output(['git', 'show', row['git_commit'] + ':' + name],
                                       cwd=REPO, env=env)
        need(len(data) == row['bytes'] and digest(data) == row['sha256'],
             'historical source mismatch: ' + name)
    for name, value in read_json('pairings/source-pins.json')['files'].items():
        need(pins[name]['sha256'] == value, 'pairings source pin: ' + name)
    for name, row in read_json('incoming-sources.json')['sources'].items():
        need(pins[name]['sha256'] == row['sha256'], 'incoming source pin: ' + name)
    for name, row in slot['source_pins'].items():
        need(pins[name]['sha256'] == row['sha256'] and pins[name]['git_commit'] == row['commit'],
             'late slot source pin: ' + name)
    for file in PACKET.rglob('*.md'):
        for link in re.findall(r'\]\(([^\n)]+)\)', file.read_text()):
            if '://' in link or link.startswith('#'):
                continue
            target = (file.parent / link.split('#')[0]).resolve()
            if target.is_relative_to(PACKET):
                need(target.exists(), 'broken packet link: ' + str(target))
            else:
                need(target.is_relative_to(REPO), 'link escapes repository: ' + link)
                need(target.relative_to(REPO).as_posix() in pins,
                     'unbound historical reference: ' + link)
    print('PASS: closed manifest, eight receipt pairs, promotions, source/review/figure bindings.', flush=True)


def replay():
    for code, normal, optimized in RUNS:
        for mode, receipt in [(False, normal), (True, optimized)]:
            command = [sys.executable, '-I', '-B'] + (['-O'] if mode else [])
            output = subprocess.check_output(command + [str(PACKET / code)], cwd=REPO)
            need(output == (PACKET / receipt).read_bytes(),
                 'replay differs: ' + code + (' -O' if mode else ''))
        print('PASS: normal and optimized replay ' + code, flush=True)
    with tempfile.TemporaryDirectory(prefix='joint-pairing-figure-') as tmp:
        target = Path(tmp)
        for name in ['make_pairing_obstruction.py', 'PAIRING_WIRING.md', 'audit/small.normal.json']:
            (target / name).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(PACKET / name, target / name)
        subprocess.check_output([sys.executable, '-I', '-B', str(target / 'make_pairing_obstruction.py')])
        need((target / 'pairing-obstruction.svg').read_bytes() ==
             (PACKET / 'pairing-obstruction.svg').read_bytes(), 'figure replay differs')
    print('PASS: exact SVG regeneration.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    check_packet()
    if args.replay:
        replay()
