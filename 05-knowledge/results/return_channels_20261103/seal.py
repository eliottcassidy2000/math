"""Validate local links and seal the current research artifact files."""
from pathlib import Path
import hashlib,json,re
HERE=Path(__file__).resolve().parent
bad=[]
for f in HERE.glob('*.md'):
    for link in re.findall(r'\]\(([^)]+)\)',f.read_text(encoding='utf8')):
        if '://' in link:continue
        target=link.split('#')[0]
        if target and not (f.parent/target).exists():bad.append((f.name,target))
if bad:raise RuntimeError(bad)
files={f.name:hashlib.sha256(f.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
       for f in sorted(HERE.iterdir()) if f.is_file() and f.name!='manifest.json' and not f.name.startswith('.')}
(HERE/'manifest.json').write_text(json.dumps({'status':'VERIFIED_CHECKPOINT','normalization':'CRLF to LF',
    'sha256_lf':files,'local_links':'PASS'},indent=2)+'\n',encoding='utf8')
print('Sealed',len(files),'files; all local links resolve')
