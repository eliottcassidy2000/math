#!/usr/bin/env python3
"""Check native manifests, exact reviewed inputs, and inert pinned source snapshots."""
from pathlib import Path, PurePosixPath
import hashlib
import json

ROOT = Path(__file__).resolve().parent
checks = 0


def need(ok, why):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(why)


def path(base, name):
    p = PurePosixPath(name)
    need(not p.is_absolute() and '..' not in p.parts and name == p.as_posix(),
         'unsafe native path: ' + name)
    return base.joinpath(*p.parts)


def read(base, name):
    return json.loads(path(base, name).read_text())


def check_row(base, row):
    b = path(base, row['path']).read_bytes()
    need(hashlib.sha256(b).hexdigest() == row['sha256'], 'native hash: ' + str(base / row['path']))
    if 'bytes' in row:
        need(len(b) == row['bytes'], 'native size: ' + str(base / row['path']))
    return b


def sums(base):
    for line in (base / 'SHA256SUMS').read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        check_row(base, {'path': name.lstrip('*'), 'sha256': digest})


def main():
    power = ROOT / 'session69-power'
    aliquot = ROOT / 'session69-aliquot'
    arithmetic = ROOT / 'session69-arithmetic-audit'
    tournament = ROOT / 'session69-tournament'
    audit = ROOT / 'session69-tournament-audit'
    root = ROOT / 'session69-root'
    addendum = ROOT / 'session69-arithmetic-corollary-audit'
    snapshots = ROOT / 'source-snapshots'
    for row in read(power, 'MANIFEST.json')['files']:
        check_row(power, row)
    for base in (aliquot, arithmetic, addendum):
        sums(base)
    for row in read(arithmetic, 'INPUTS.json')['files']:
        check_row({'power': power, 'aliquot': aliquot}[row['lane']], row)
    for row in read(tournament, 'FILE_MANIFEST.json')['files']:
        check_row(tournament, row)
    tm = read(audit, 'MANIFEST.json')
    for name, row in tm['files'].items():
        check_row(audit, dict(row, path=name))
    for name, row in tm['reviewed_producer_files'].items():
        check_row(tournament, dict(row, path=name))
    incoming = read(audit, 'INCOMING_SCOPE_MANIFEST.json')
    check_row(audit, dict(incoming['base_manifest'], path='MANIFEST.json'))
    for name, row in incoming['files'].items():
        b = check_row(audit, dict(row, path=name))
        if name in incoming['source_git_blobs']:
            need(hashlib.sha1(b'blob ' + str(len(b)).encode() + b'\0' + b).hexdigest()
                 == incoming['source_git_blobs'][name], 'incoming source Git blob')
    for mode in ('normal', 'optimized'):
        need((audit / ('producer.replay.' + mode + '.json')).read_bytes()
             == (tournament / 'receipt.json').read_bytes(), 'audited producer replay')
    for row in read(root, 'MANIFEST.json')['files']:
        check_row(root, row)
    review = read(addendum, 'ROOT_REVIEW.json')
    need(review['status'] == 'ACCEPT', 'root review not accepted')
    for row in review['root_payload_files']:
        check_row(root, row)
    check_row(root, {'path': 'MANIFEST.json', 'sha256': review['root_manifest_sha256']})
    for stem, native in (('root', 'rare_six'), ('firewalls', 'firewalls')):
        for mode in ('normal', 'optimized'):
            need((addendum / (stem + '.' + mode + '.json')).read_bytes()
                 == (root / (native + '.json')).read_bytes(), 'root replay receipt')
    sm = read(snapshots, 'SOURCE_SNAPSHOTS.json')
    inherited = read(tournament, 'SOURCE_MANIFEST.json')
    need(sm['input_commit'] == inherited['input_commit'], 'source snapshot commit')
    need([r['path'] for r in sm['files']] == [r['path'] for r in inherited['inputs']],
         'source snapshot universe')
    for row, native in zip(sm['files'], inherited['inputs']):
        need(all(row[k] == native[k] for k in ('path', 'bytes', 'sha256')), 'source pin differs')
        b = check_row(snapshots, row)
        need(hashlib.sha1(b'blob ' + str(len(b)).encode() + b'\0' + b).hexdigest()
             == row['git_blob'], 'source Git blob')
    print(json.dumps({'status': 'PASS', 'active_checks': checks,
                      'research_pin': inherited['input_commit'],
                      'native_manifests_and_reviewed_inputs': 'MATCH',
                      'inert_source_snapshots': len(sm['files']),
                      'snapshot_programs_executed': False}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
