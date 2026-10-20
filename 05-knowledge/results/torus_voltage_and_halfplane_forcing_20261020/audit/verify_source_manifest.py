#!/usr/bin/env python3
"""Verify saved primary-source bytes without running or importing those sources.

Usage: python3 verify_source_manifest.py MANIFEST SNAPSHOT_ROOT
SNAPSHOT_ROOT contains files under the repository paths named by MANIFEST.
The publication intentionally excludes downloaded copies of primary sources.
Obtain the exact files at the pinned commit before reproducing this check.
"""
import hashlib
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
manifest = json.loads(Path(sys.argv[1]).read_text())
snapshots = Path(sys.argv[2])
checks = []
for record in manifest['files']:
    data = (snapshots/record['path']).read_bytes()
    lines = data.decode('utf-8').splitlines()
    blob = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    targets = [{'line': t['line'], 'matches': lines[t['line']-1] == t['text']}
               for t in record.get('targets', [])]
    checks.append({'path': record['path'], 'expected_blob_sha': record['sha'],
                   'recomputed_blob_sha': blob, 'blob_matches': blob == record['sha'],
                   'lines': len(lines), 'line_count_matches': len(lines) == record['lines'],
                   'targets': targets})
passed = all(r['blob_matches'] and r['line_count_matches'] and
             all(t['matches'] for t in r['targets']) for r in checks)
result = {'status': 'PASS_SOURCE_BYTES_ONLY' if passed else 'FAIL',
          'repository': manifest['repository'], 'ref': manifest['ref'],
          'files': checks, 'lean_on_path': bool(shutil.which('lean')),
          'lake_on_path': bool(shutil.which('lake')), 'build_performed': False,
          'comparator_executed': False, 'transitive_axiom_audit_performed': False,
          'trusted_code_review_performed': False,
          'scope': 'Blob hashes, physical line counts, and recorded declaration locations only; no proof verification.'}
text = json.dumps(result, indent=2)+'\n'
(HERE/'source_verification.json').write_text(text)
print(text, end='')
if not passed:
    raise SystemExit(1)
