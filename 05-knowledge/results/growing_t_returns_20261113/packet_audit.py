"""Packet integrity: links, immutable inputs, and normal/optimized receipts."""
from pathlib import Path
from urllib.parse import unquote
import gzip
import hashlib
import json
import re
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_json(path):
    data = path.read_bytes()
    if path.suffix == '.gz':
        data = gzip.decompress(data)
    return json.loads(data)


pins = {}


def check_pin(path, sha):
    target = ROOT / path
    need(target.exists(), ('missing pinned file', path))
    data = target.read_bytes()
    need(sha in (digest(data), digest(data.replace(b'\r\n', b'\n'))), ('stale pin', path))
    pins[(path, sha)] = True


links = 0
incoming_links = []
for path in HERE.rglob('*.md'):
    prose = path.read_text(encoding='utf-8')
    need('\ufffd' not in prose, ('replacement character', str(path)))
    for bad in ('\u00e2\u20ac', '\u00c3\u00a2', '\u00c2\u00b7'):
        need(bad not in prose, ('mojibake', str(path)))
    for target in re.findall(r'\[[^\]\n]*\]\(([^()\n]+)\)', prose):
        if target.startswith(('https:', 'http:', 'mailto:', '#')):
            continue
        name = unquote(target.split('#')[0]).strip('<>')
        if not name or ('.' not in name and '/' not in name):
            continue
        destination = (path.parent / name).resolve()
        if not destination.exists():
            relative = str(destination.relative_to(ROOT)).replace('\\', '/')
            status = subprocess.run(['git', 'cat-file', '-e', 'origin/main:' + relative], cwd=ROOT,
                                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            need(status.returncode == 0, ('broken link', str(path), target))
            incoming_links.append(relative)
        links += 1

pairs = []
for optimized in sorted(HERE.rglob('*optimized.json*')):
    normal = optimized.with_name(optimized.name.replace('.optimized', ''))
    need(normal.exists(), ('missing normal receipt', str(optimized)))
    need(read_json(normal) == read_json(optimized), ('normal/optimized mismatch', str(normal)))
    pairs.append(str(normal.relative_to(HERE)).replace('\\', '/'))

for item in read_json(HERE / 'geometry/source_pins.json'):
    check_pin(item['path'], item['sha256'])
for item in read_json(HERE / 'inverse/SOURCE_PINS.json')['files']:
    check_pin(item['path'], item['sha256'])
for item in read_json(HERE / 'geometry/star_kernel_audit.json')['source_pins']:
    check_pin(item['path'], item['sha256'])
for item in read_json(HERE / 'geometry/THREE_STEP_AUDIT.json')['source_pins']:
    check_pin(item['path'], item['sha256'])
for relative in ('containment/receipt.json.gz', 'star_kernel/receipt.json'):
    for path, sha in read_json(HERE / relative)['source_pins'].items():
        check_pin(path, sha)
for path, sha in read_json(HERE / 'repairs/source_manifest.json')['source_sha256'].items():
    check_pin(path, sha)

incoming = read_json(HERE / 'geometry/incoming_clocks.json')['git_source']
data = subprocess.check_output(['git', 'show', incoming['commit'] + ':' + incoming['path']], cwd=ROOT)
need(digest(data) == incoming['sha256'], 'immutable incoming mathematical source')

report = dict(status='PASS', local_markdown_links=links,
              incoming_links_pending_checkout=sorted(set(incoming_links)),
              normal_optimized_pairs=pairs, distinct_files_and_hashes_checked=len(pins),
              immutable_incoming_commit=incoming['commit'],
              scope='Integrity and receipt agreement; independent mathematical proof reviews remain in AUDIT.md. Inherited complete isomorphism censuses are pinned, not regenerated.')
(HERE / 'packet_audit.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
print(json.dumps(report, indent=2))
