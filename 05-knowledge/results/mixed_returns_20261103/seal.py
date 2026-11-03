"""Check artifact links and seal normalized files."""
from pathlib import Path
import hashlib,json,re
HERE=Path(__file__).resolve().parent
for f in HERE.glob('*.md'):
    for link in re.findall(r'\]\(([^)]+)\)',f.read_text(encoding='utf8')):
        if '://' not in link and not (f.parent/link.split('#')[0]).exists():raise RuntimeError((f.name,link))
hashes={f.name:hashlib.sha256(f.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for f in sorted(HERE.iterdir()) if f.is_file() and f.name!='manifest.json'}
(HERE/'manifest.json').write_text(json.dumps({'status':'VERIFIED_CHECKPOINT','normalization':'CRLF to LF','local_links':'PASS','sha256_lf':hashes},indent=2)+'\n',encoding='utf8')
print('Sealed',len(hashes),'files')
