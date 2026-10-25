#!/usr/bin/env python3
"""Verify frozen sources, local links, receipts, and the packet manifest."""
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

BASE=Path(__file__).resolve().parent
ROOT=Path(subprocess.check_output(['git','rev-parse','--show-toplevel'],cwd=BASE,text=True).strip())


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    pins=json.loads((BASE/'SOURCES.json').read_text())
    for row in pins['repository_sources']:
        data=subprocess.check_output(['git','show',row['commit']+':'+row['path']],cwd=ROOT)
        if digest(data)!=row['sha256']:
            raise RuntimeError('source pin changed: '+row['path'])
    for name,script in [('normal','check.py'),('smith_normal','audit_smith.py'),('twists_normal','check_twists.py'),('ports_normal','check_ports.py')]:
        normal=BASE/(name+'.json')
        optimized=BASE/((name.replace('normal','optimized'))+'.json')
        if normal.read_bytes()!=optimized.read_bytes():
            raise RuntimeError('optimization-mode mismatch: '+name)
        receipt=json.loads(normal.read_text())
        if receipt['script_sha256']!=digest((BASE/script).read_bytes()):
            raise RuntimeError('receipt source mismatch: '+script)
    link_count=0
    for p in BASE.glob('*.md'):
        for dest in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if '://' in dest or dest.startswith('#'):
                continue
            dest=dest.split('#')[0]
            path=(p.parent/dest).resolve()
            if not path.exists():
                rel=str(path.relative_to(ROOT))
                if subprocess.run(['git','cat-file','-e','HEAD:'+rel],cwd=ROOT,capture_output=True).returncode:
                    raise RuntimeError('missing link '+dest)
            link_count+=1
    actual={str(p.relative_to(BASE)):digest(p.read_bytes()) for p in sorted(BASE.rglob('*'))
            if p.is_file() and p.name not in ('MANIFEST.json','verification.json') and '__pycache__' not in p.parts}
    if '--seal' in sys.argv:
        (BASE/'MANIFEST.json').write_text(json.dumps(actual,indent=2,sort_keys=True)+'\n')
    expected=json.loads((BASE/'MANIFEST.json').read_text())
    if actual!=expected:
        raise RuntimeError('manifest differs; reseal only after intentional changes and replays')
    print(json.dumps({'status':'PASS','files':len(actual),'source_pins':len(pins['repository_sources']),
                      'local_links':link_count,'matching_normal_optimized_receipts':4},indent=2,sort_keys=True))


if __name__=='__main__':
    main()
