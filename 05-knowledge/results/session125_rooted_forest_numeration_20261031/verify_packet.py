#!/usr/bin/env python3
"""Strict self-contained packet integrity check; no removable asserts."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parent

def check(value, message):
    if not value:
        raise RuntimeError(message)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

whitelist = json.loads((root/'PUBLICATION_FILES.json').read_text())
check(whitelist == sorted(set(whitelist)), 'whitelist not sorted and unique')
check(all(not Path(p).is_absolute() and '..' not in Path(p).parts for p in whitelist),
      'unsafe whitelist path')
actual = sorted(str(p.relative_to(root)) for p in root.rglob('*') if p.is_file())
check(actual == whitelist, 'missing or unexpected publication file')
manifest = json.loads((root/'PUBLICATION_MANIFEST.json').read_text())
expected_content = set(whitelist)-{'PUBLICATION_MANIFEST.json','SHA256SUMS'}
check(set(manifest['files']) == expected_content, 'manifest coverage')
for name, item in manifest['files'].items():
    check((root/name).stat().st_size == item['bytes'], 'byte count: '+name)
    check(sha(root/name) == item['sha256'], 'content hash: '+name)
sums = {}
for line in (root/'SHA256SUMS').read_text().splitlines():
    digest, name = line.split('  ',1)
    check(name not in sums, 'duplicate checksum entry')
    sums[name] = digest
check(set(sums) == set(whitelist)-{'SHA256SUMS'}, 'checksum coverage')
for name,digest in sums.items():
    check(sha(root/name) == digest, 'checksum mismatch: '+name)
for stem in ('verification','audit/independent'):
    a=root/(stem+'.normal.json');b=root/(stem+'.optimized.json')
    check(a.read_bytes() == b.read_bytes(), 'normal/-O mismatch: '+stem)
producer=json.loads((root/'verification.normal.json').read_text())
independent=json.loads((root/'audit/independent.normal.json').read_text())
check(producer['status']=='PASS' and producer['checks']==1245858, 'producer result')
check(independent['status']=='PASS' and independent['checks']==129144, 'independent result')
for name,digest in producer['source_sha256'].items():
    check(sha(root/name)==digest,'producer source binding: '+name)
for name,digest in independent['input_sha256'].items():
    check(sha(root/name)==digest,'independent input binding: '+name)
print(json.dumps({'status':'PASS','files':len(whitelist),
                  'producer_checks':producer['checks'],
                  'independent_checks':independent['checks']},sort_keys=True))
