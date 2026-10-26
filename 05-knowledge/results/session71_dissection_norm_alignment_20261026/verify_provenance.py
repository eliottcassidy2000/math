#!/usr/bin/env python3
"""Actively check native hash lists, reviewed bytes, and inert source Git blobs."""
import hashlib
import json
from pathlib import Path, PurePosixPath
import re

ROOT = Path(__file__).resolve().parent
CHECKS = 0

def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)

def path(name):
    p = PurePosixPath(name)
    need(not p.is_absolute() and '..' not in p.parts and '\\' not in name and name == p.as_posix() and name not in ('', '.'),
         'unsafe provenance path')
    dest = ROOT.joinpath(*p.parts)
    for q in (dest, *dest.parents):
        if q == ROOT:
            break
        need(not q.is_symlink(), 'provenance symlink')
    return dest

def check_bytes(name, digest, size=None):
    data = path(name).read_bytes()
    need(hashlib.sha256(data).hexdigest() == digest, 'provenance SHA256 mismatch: ' + name)
    if size is not None:
        need(len(data) == size, 'provenance size mismatch: ' + name)
    return data

def main():
    rules = json.loads(path('PROVENANCE_RULES.json').read_text())
    need(rules['schema'] == 'session71-provenance-rules-v1', 'provenance schema')
    source_rows = 0
    sum_manifests = [{'lane': lane, 'file': 'SHA256SUMS'} for lane in rules.get('sha256sum_lanes', [])]
    sum_manifests += rules.get('sha256sum_manifests', [])
    for record in sum_manifests:
        lane = record['lane']
        for line in path(lane + '/' + record['file']).read_text().splitlines():
            digest, file = line.split(maxsplit=1)
            check_bytes(lane + '/' + file.lstrip('*'), digest)
    for source in rules.get('source_ledgers', []):
        lane = source['lane']
        ledger = json.loads(path(lane + '/' + source['file']).read_text())
        rows = ledger
        if source.get('rows_key'):
            rows = rows[source['rows_key']]
        for row in rows:
            if source.get('require_row_key') and source['require_row_key'] not in row:
                need('url' in row, 'only external URL metadata may be omitted from snapshot verification')
                continue
            data = check_bytes(lane + '/' + row[source.get('snapshot_key', 'snapshot')], row['sha256'], row.get('bytes'))
            pin = ledger[source['ledger_pin_key']] if source.get('ledger_pin_key') else row[source.get('pin_key', 'pin')]
            need(re.fullmatch('[0-9a-f]{40}', pin) is not None, 'nonexact source pin')
            need(pin in rules.get('allowed_source_pins', []), 'unapproved source pin')
            blob = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
            need(blob == row[source.get('blob_key', 'git_blob')], 'source Git blob mismatch')
            source_rows += 1
    for record in rules.get('native_manifests', []):
        ledger = json.loads(path(record['lane'] + '/' + record['file']).read_text())
        rows = ledger[record.get('rows_key', 'files')]
        if record.get('mapping'):
            rows = [{'path': key, 'sha256': value} for key, value in rows.items()]
        for row in rows:
            check_bytes(record['lane'] + '/' + row[record.get('path_key', 'path')], row['sha256'], row.get('bytes'))
    for binding in rules.get('reviewed_bindings', []):
        data = check_bytes(binding['path'], binding['sha256'], binding.get('bytes'))
        if binding.get('equal_to'):
            need(data == path(binding['equal_to']).read_bytes(), 'review replay mismatch')
    print(json.dumps({'status': 'PASS', 'active_checks': CHECKS, 'inert_source_snapshots': source_rows,
                      'native_hash_lists': len(sum_manifests),
                      'reviewed_bindings': len(rules.get('reviewed_bindings', [])),
                      'snapshot_programs_executed': False}, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
