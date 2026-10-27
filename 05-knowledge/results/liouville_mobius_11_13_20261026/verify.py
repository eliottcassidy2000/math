"""Verify frozen bytes before replay; never update expected receipts."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

sys.stdout.reconfigure(encoding='utf-8',newline='\n')
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[2]

def need(test,message):
    if not test: raise RuntimeError(message)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--replay',action='store_true')
    parser.add_argument('--check-sources',action='store_true')
    args=parser.parse_args()
    entries=[]
    for line in (ROOT/'SHA256SUMS').read_text(encoding='utf-8').splitlines():
        digest,name=line.split('  ',1); path=ROOT/name
        need(path.is_file(),'missing '+name)
        need(hashlib.sha256(path.read_bytes()).hexdigest()==digest,'hash mismatch '+name)
        entries.append(name)
    actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
            if p.is_file() and '__pycache__' not in p.parts and p.name!='SHA256SUMS'}
    need(len(entries)==len(set(entries)) and set(entries)==actual,'manifest coverage')
    for md in ROOT.glob('*.md'):
        for target in re.findall(r'\]\(([^\s)]+)\)',md.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'): continue
            need((md.parent/target.split('#')[0]).resolve().exists(),'broken relative link '+target)
    print('Manifest and local links PASS: '+str(len(entries))+' files',flush=True)
    if args.check_sources:
        pins=json.loads((ROOT/'source-pins.json').read_text(encoding='utf-8'))
        for name,digest in pins['sources'].items():
            data=(REPO/name).read_bytes().replace(b'\r\n',b'\n')
            need(hashlib.sha256(data).hexdigest()==digest,'inherited source drift '+name)
        print('Inherited source pins PASS',flush=True)
    if args.replay:
        for script,stem in [('check.py',''),('audit.py','audit-')]:
            expected=(ROOT/(stem+'normal.json')).read_bytes()
            need(expected==(ROOT/(stem+'optimized.json')).read_bytes(),'frozen optimization pair differs')
            for optimized in (False,True):
                cmd=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(ROOT/script)]
                result=subprocess.run(cmd,cwd=REPO,capture_output=True,check=True)
                need(result.stdout==expected,'replay differs '+script+(' -O' if optimized else ''))
            print(script+': normal and optimized PASS',flush=True)

if __name__=='__main__': main()
