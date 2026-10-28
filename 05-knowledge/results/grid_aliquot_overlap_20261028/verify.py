"""Verify frozen evidence; replay never overwrites expected outputs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[2]
sys.stdout.reconfigure(encoding='utf-8',newline='\n')

def need(ok,message):
    if not ok: raise RuntimeError(message)

def native(path):
    value=str(path.resolve())
    return Path('\\\\?\\'+value) if os.name=='nt' and not value.startswith('\\\\?\\') else Path(value)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--replay',action='store_true')
    parser.add_argument('--check-sources',action='store_true')
    args=parser.parse_args()
    names=[]
    for line in (ROOT/'SHA256SUMS').read_text(encoding='utf-8').splitlines():
        digest,name=line.split('  ',1)
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,'hash mismatch '+name)
        names.append(name)
    actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
            if p.is_file() and p.name!='SHA256SUMS' and '__pycache__' not in p.parts}
    need(len(names)==len(set(names)) and set(names)==actual,'manifest coverage')
    for md in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^\s)]+)\)',md.read_text(encoding='utf-8')):
            if '://' in link or link.startswith('#'): continue
            need(native(md.parent/link.split('#')[0]).exists(),'broken local link '+link)
    print('Frozen manifest and local links PASS',flush=True)
    if args.check_sources:
        pins=json.loads((ROOT/'source-pins.json').read_text(encoding='utf-8'))
        for name,digest in pins['committed_sources'].items():
            blob=subprocess.run(['git','show','HEAD:'+name],cwd=REPO,capture_output=True,check=True).stdout
            need(hashlib.sha256(blob).hexdigest()==digest,'committed source drift '+name)
        print('Inherited committed source pins PASS',flush=True)
    if args.replay:
        expected=(ROOT/'normal.json').read_bytes()
        need(expected==(ROOT/'optimized.json').read_bytes(),'frozen mode mismatch')
        for mode in ([],['-O']):
            result=subprocess.run([sys.executable,'-B']+mode+[str(ROOT/'check.py')],
                                  cwd=REPO,capture_output=True,check=True)
            need(result.stdout==expected,'replay mismatch '+str(mode))
        print('Normal and optimized exact replay PASS',flush=True)

if __name__=='__main__': main()
