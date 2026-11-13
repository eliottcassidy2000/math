"""Integrity, declared-source and receipt agreement checks for this packet."""
from pathlib import Path
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def digest(data):
    return hashlib.sha256(data).hexdigest()

def read(path):
    return json.loads(path.read_text(encoding='utf-8'))

links = 0
for path in HERE.rglob('*.md'):
    content = path.read_text(encoding='utf-8')
    need('\ufffd' not in content and '\u00e2\u20ac' not in content, ('invalid text', str(path)))
    for target in re.findall(r'\]\(([^)]+)\)', content):
        if target.startswith(('http:', 'https:', '#')):
            continue
        dest = (path.parent / target.split('#')[0]).resolve()
        need(dest.exists(), ('broken link', str(path), target))
        links += 1

pairs = []
pins = set()
for path in HERE.rglob('*.optimized.json'):
    normal = path.with_name(path.name.replace('.optimized', ''))
    need(read(normal) == read(path), ('receipt mismatch', str(path)))
    pairs.append(str(normal.relative_to(HERE)).replace('\\','/'))
    sources = read(normal).get('source_pins', [])
    if isinstance(sources, dict):
        sources = [dict(path=p, sha256=h) for p,h in sources.items()]
    for source in sources:
        data = (ROOT / source['path']).read_bytes()
        expected = source.get('sha256', source.get('sha256_lf'))
        need(expected in (digest(data),digest(data.replace(b'\r\n',b'\n'))), ('stale pin',source['path']))
        pins.add((source['path'],expected))

data = (HERE/'TREE_BLUEPRINT.md').read_bytes()
need(read(HERE/'tree_blueprint.json')['proof_sha256'] in
     (digest(data),digest(data.replace(b'\r\n',b'\n'))), 'tree proof receipt stale')
report = dict(status='PASS',local_markdown_links=links,normal_optimized_pairs=sorted(pairs),
              distinct_source_pins=len(pins),scope='Integrity checks; independent proof reviews and exact computation universes are stated in AUDIT.md.')
(HERE/'packet_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report,indent=2))
