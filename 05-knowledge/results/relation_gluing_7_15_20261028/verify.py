"""Verify receipts, pinned sources, local links, and optionally replay both implementations."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys

sys.stdout.reconfigure(encoding='utf-8', newline='\n')
p=Path(__file__).resolve().parent
repo=p.parents[2]
args=argparse.ArgumentParser()
args.add_argument('--replay',action='store_true')
args=args.parse_args()
def need(b,msg):
    if not b: raise RuntimeError(msg)
def read(path):
    path=path.resolve()
    if os.name=='nt': path=Path('\\\\?\\'+str(path))
    return path.read_bytes()
def basis(name,b):
    return b if name.endswith('.png') else b.replace(b'\r\n',b'\n')
manifest=json.loads(read(p/'SHA256SUMS.json'))
for name,digest in manifest.items():
    need(hashlib.sha256(basis(name,read(p/name))).hexdigest()==digest,'file hash: '+name)
for first,second in [('normal.json','optimized.json'),('audit-normal.json','audit-optimized.json')]:
    need(read(p/first)==read(p/second),'ordinary/optimized byte equality: '+first)
need(json.loads(read(p/'audit-normal.json'))['status']=='PASS','separate audit status')
pins=json.loads(read(p/'source_pins.json'))
for pin in pins['sources']:
    out=subprocess.run(['git','show',pin['commit']+':'+pin['path']],cwd=repo,capture_output=True,check=True).stdout
    need(hashlib.sha256(out.replace(b'\r\n',b'\n')).hexdigest()==pin['sha256_lf'],'source pin: '+pin['path'])
links=0
for f in p.glob('*.md'):
    for link in re.findall(r'\]\(([^)]+)\)',read(f).decode('utf-8')):
        if '://' in link: continue
        target=link.split('#')[0]
        if target:
            read(f.parent/target)
            links+=1
if args.replay:
    for script,normal,optimized in [('check.py','normal.json','optimized.json'),
                                    ('audit.py','audit-normal.json','audit-optimized.json')]:
        for flag,receipt in [([],normal),(['-O'],optimized)]:
            run=subprocess.run([sys.executable,'-B']+flag+[str(p/script)],cwd=p,capture_output=True,check=True)
            need(run.stdout==read(p/receipt),'replay: '+script+' '+receipt)
print(json.dumps({'status':'PASS','manifest_files':len(manifest),'source_pins':len(pins['sources']),
                  'local_links':links,'receipt_pairs_equal':2,'replay':args.replay},indent=2,sort_keys=True))
