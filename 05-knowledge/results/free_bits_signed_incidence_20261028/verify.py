"""Check immutable source pins and packet receipts; optionally replay."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def need(p,s):
    if not p: raise RuntimeError(s)

def digest(data): return sha256(data.replace(b'\r\n',b'\n')).hexdigest()
def read(name): return json.loads((HERE/name).read_text(encoding='utf-8'))

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--replay',action='store_true'); args=parser.parse_args()
    names=[]
    for line in (HERE/'SHA256SUMS').read_text().splitlines():
        expected,name=line.split('  ',1)
        need(digest((HERE/name).read_bytes())==expected,'hash '+name); names.append(name)
    need(set(names)=={p.name for p in HERE.iterdir() if p.is_file() and p.name!='SHA256SUMS'},'inventory')
    for item in read('SOURCES.json')['repository_sources']:
        content=subprocess.check_output(['git','cat-file','blob',item['git_blob']],cwd=ROOT)
        need(digest(content)==item['sha256_lf'],'source '+item['path'])
    for p in HERE.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if '://' not in link and not link.startswith('#'):
                need((p.parent/link.split('#')[0]).resolve().exists(),'link '+link)
    for prefix in ('','audit-'):
        need(read(prefix+'normal.json')==read(prefix+'optimized.json'),'normal/O '+prefix)
    if args.replay:
        for script,receipt in [('check.py','normal.json'),('audit.py','audit-normal.json')]:
            for flags in ([],['-O']):
                result=subprocess.run([sys.executable,'-B',*flags,str(HERE/script)],cwd=ROOT,check=True,capture_output=True,timeout=60)
                need(json.loads(result.stdout)==read(receipt),'replay '+script)
    print(json.dumps(dict(sealed_files=len(names),source_pins='PASS',links='PASS',receipts='PASS',replay=args.replay),sort_keys=True))

if __name__=='__main__': main()
