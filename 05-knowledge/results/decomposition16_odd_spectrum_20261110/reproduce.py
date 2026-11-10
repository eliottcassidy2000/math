#!/usr/bin/env python3
"""Validate frozen evidence; --rerun also executes all seven checker pairs."""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
JOBS=('induction/check','attachment/check','history/check',
      'history/attachment_audit','induction/audit_joint','joint_check','core_audit')

def require(x,msg):
    if not x:raise RuntimeError(msg)

def digest(p):return sha256(p.read_bytes()).hexdigest()

def resolve_pin(s):
    candidates={p.resolve() for p in (ROOT/s,HERE/s) if p.is_file()}
    require(len(candidates)==1,('unresolved or ambiguous source',s))
    return next(iter(candidates))

def files():
    return sorted(p for p in HERE.rglob('*') if p.is_file()
                  and p.name!='manifest.json' and '__pycache__' not in p.parts)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--rerun',action='store_true')
    parser.add_argument('--write-manifest',action='store_true')
    args=parser.parse_args()
    pins=set()
    for stem in JOBS:
        script=HERE/(stem+'.py')
        normal=HERE/(stem+'.json');optimized=HERE/(stem+'.optimized.json')
        require(normal.read_bytes()==optimized.read_bytes(),('optimization mismatch',stem))
        data=json.loads(normal.read_text())
        require(data['status']=='PASS',('receipt failed',stem))
        require(data['script_sha256']==digest(script),('script changed',stem))
        for key in ('source_sha256','reviewed_sha256'):
            for name,expected in data.get(key,{}).items():
                path=resolve_pin(name)
                require(digest(path)==expected,('pinned source changed',stem,name))
                pins.add(str(path.relative_to(ROOT)))
        if 'proof_sha256' in data:
            require(digest(script.parent/'PROOF.md')==data['proof_sha256'],('proof changed',stem))
        if args.rerun:
            for opt in (False,True):
                command=[sys.executable,'-B']+(['-O'] if opt else [])+[str(script)]
                result=subprocess.run(command,cwd=ROOT,capture_output=True,check=True)
                require(result.stdout==normal.read_bytes(),('fresh replay differs',stem,opt))
    links=0
    for md in HERE.rglob('*.md'):
        for link in re.findall(r'\]\(([^\s)]+)\)',md.read_text()):
            if '://' in link or link.startswith('#'):continue
            target=link.split('#',1)[0]
            require((md.parent/target).is_file(),('broken local link',str(md.relative_to(HERE)),link))
            links+=1
    actual={str(p.relative_to(HERE)):digest(p) for p in files()}
    manifest=HERE/'manifest.json'
    if args.write_manifest:
        manifest.write_text(json.dumps(actual,indent=2,sort_keys=True)+'\n')
    else:
        require(json.loads(manifest.read_text())==actual,'manifest differs')
    print(json.dumps({'status':'PASS','normal_optimized_pairs':len(JOBS),
        'distinct_pinned_inputs':len(pins),'local_links':links,'artifacts':len(actual),
        'mode':'fresh mathematical reruns' if args.rerun else 'saved evidence validation'},
        indent=2,sort_keys=True))

if __name__=='__main__':main()
