"""Check local links and source receipts, then hash this packet's checks."""
from pathlib import Path
import hashlib
import json
import re

HERE=Path(__file__).resolve().parent


def main():
    links=0
    for doc in HERE.glob('*.md'):
        for raw in re.findall(r'\]\(([^)]+)\)',doc.read_text(encoding='utf-8')):
            target=raw.split('#')[0]
            if not target or '://' in target:continue
            if target=='manifest.json':continue
            if not (doc.parent/target).resolve().exists():raise RuntimeError(f'missing link {doc.name}: {target}')
            links+=1
    receipts=[]
    for path in HERE.glob('*.json'):
        if path.name=='manifest.json':continue
        data=json.loads(path.read_text(encoding='utf-8'))
        expected=data.get('source_sha256')
        if expected:
            script=path.with_suffix('.py')
            digest=hashlib.sha256(script.read_bytes()).hexdigest()
            if digest!=expected:raise RuntimeError(f'source hash mismatch {path}')
            receipts.append(path.name)
    files={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
           for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='manifest.json'}
    out=dict(status='VERIFIED local references and supplied source hashes',
             local_links=links,source_receipts=receipts,files=files)
    (HERE/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(local_links=links,source_receipts=receipts,hashed_files=len(files)),indent=2))


if __name__=='__main__':main()
