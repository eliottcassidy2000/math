"""Pin inherited sources at explicitly read commits and hash the finished packet."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

p=Path(__file__).resolve().parent
repo=p.parents[2]
base='83b601158076669edc8729ce461f150da1c3877e'
incoming='1b86e4eadcd6bab19019a1226059f2f1fbf37e13'
paths=set()
for f in p.glob('*.md'):
    for link in re.findall(r'\]\(([^)]+)\)',f.read_text(encoding='utf-8')):
        if '://' in link: continue
        target=link.split('#')[0]
        if not target: continue
        q=(f.parent/target).resolve()
        if p not in q.parents and q!=p:
            paths.add(q.relative_to(repo).as_posix())
pins=[]
for path in sorted(paths):
    commit=incoming if ('session85_arithmetic_codes_free_bits' in path or 'free_bits_signed_incidence_20261028' in path or 'THM-4839-' in path) else base
    out=subprocess.run(['git','show',commit+':'+path],cwd=repo,capture_output=True,check=True).stdout
    pins.append({'path':path,'commit':commit,
                 'sha256_lf':hashlib.sha256(out.replace(b'\r\n',b'\n')).hexdigest()})
(p/'source_pins.json').write_bytes((json.dumps({'base_commit':base,'incoming_read_commit':incoming,
    'scope':'Inherited sources; exact sections and transfer boundaries in ATLAS.md, not full independent audits.',
    'sources':pins},indent=2,sort_keys=True)+'\n').encode('utf-8'))
manifest={}
for f in sorted(p.iterdir()):
    if not f.is_file() or f.name in ('SHA256SUMS.json','verification.json'): continue
    raw=f.read_bytes()
    if f.suffix!='.png': raw=raw.replace(b'\r\n',b'\n')
    manifest[f.name]=hashlib.sha256(raw).hexdigest()
(p/'SHA256SUMS.json').write_bytes((json.dumps(manifest,indent=2,sort_keys=True)+'\n').encode('utf-8'))
print(json.dumps({'pinned_sources':len(pins),'packet_files':len(manifest)}))
