#!/usr/bin/env python3
"""Replay scoped controls and verify the frozen mathematical source surface."""
from hashlib import sha256
from pathlib import Path
import argparse
import json
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SOURCE_HEAD='ab5098e5b8e3'
JOBS=('mutation/check','arithmetic/check','returns/check',
      'matrix/check','arithmetic/audit_matrix')

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def digest(p):return sha256(p.read_bytes()).hexdigest()
def dump(p,obj):p.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--rerun',action='store_true')
    parser.add_argument('--refresh',action='store_true')
    args=parser.parse_args()
    sources=set()
    for job in JOBS:
        script=HERE/(job+'.py')
        for optimized in (False,True):
            receipt=HERE/(job+('.optimized' if optimized else '')+'.json')
            if args.rerun or args.refresh:
                cmd=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(script)]
                out=subprocess.run(cmd,cwd=ROOT,capture_output=True,check=True)
                obj=json.loads(out.stdout)
                if args.refresh:dump(receipt,obj)
                else:need(obj==json.loads(receipt.read_text()),('fresh replay changed',job,optimized))
        normal=HERE/(job+'.json');optimized=HERE/(job+'.optimized.json')
        need(normal.read_bytes()==optimized.read_bytes(),('optimized mismatch',job))
        obj=json.loads(normal.read_text())
        need(str(obj['status']).startswith(('PASS','FINITE-EXACT')),('receipt status',job))
        need(obj['script_sha256']==digest(script),('script changed',job))
        for name,expected in obj.get('source_sha256',{}).items():
            path=(ROOT/name).resolve()
            need(digest(path)==expected,('receipt source changed',name))
            if not path.is_relative_to(HERE):sources.add(str(path.relative_to(ROOT)))
    links=0
    for md in HERE.rglob('*.md'):
        for link in re.findall(r'\]\(([^\s)]+)\)',md.read_text()):
            if '://' in link or link.startswith('#'):continue
            path=(md.parent/link.split('#',1)[0]).resolve()
            need(path.is_file(),('broken link',str(md.relative_to(HERE)),link))
            links+=1
            if path.is_relative_to(ROOT) and not path.is_relative_to(HERE):
                sources.add(str(path.relative_to(ROOT)))
    pins_path=HERE/'SOURCE_PINS.json'
    if args.refresh:
        dump(pins_path,{'source_HEAD':SOURCE_HEAD,
             'source_sha256':{p:digest(ROOT/p) for p in sorted(sources)},
             'primary_source':{'url':'https://arxiv.org/abs/math/0104151',
                 'title':'Fomin–Zelevinsky, Cluster algebras I: Foundations',
                 'checked':'2026-11-12','use':'Attribution of ordinary matrix mutation'}})
    pins=json.loads(pins_path.read_text())
    for name,expected in pins['source_sha256'].items():
        path=ROOT/name
        if path.is_file() and digest(path)==expected:continue
        raw=subprocess.run(['git','show',pins['source_HEAD']+':'+name],
                           cwd=ROOT,capture_output=True,check=True).stdout
        need(sha256(raw).hexdigest()==expected,('historical pin changed',name))
    manifest={str(p.relative_to(HERE)):digest(p) for p in sorted(HERE.rglob('*'))
              if p.is_file() and p.name!='MANIFEST.json' and '__pycache__' not in p.parts}
    if args.refresh:dump(HERE/'MANIFEST.json',manifest)
    need(json.loads((HERE/'MANIFEST.json').read_text())==manifest,'artifact manifest changed')
    print(json.dumps({'status':'PASS','checker_pairs':len(JOBS),'source_pins':len(pins['source_sha256']),
                      'local_links':links,'artifacts':len(manifest),
                      'fresh_replay':args.rerun or args.refresh},indent=2,sort_keys=True))

if __name__=='__main__':main()
