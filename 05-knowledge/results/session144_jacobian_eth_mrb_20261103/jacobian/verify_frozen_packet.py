#!/usr/bin/env python3
"""Verify an explicitly frozen packet and replay its pinned exact checks.

No Python assert gates. Run normally or under -O: the child interpreters keep
that mode. A file edit or a different test output is a hard failure.
"""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path

def need(ok,message):
    if not ok:raise RuntimeError(message)

def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,default=Path(__file__).with_name('MANIFEST.json'));p.add_argument('--hashes-only',action='store_true');a=p.parse_args()
    manifest=a.manifest.resolve();root=manifest.parent;m=json.loads(manifest.read_text());seen=set()
    for row in m['files']:
        rel=Path(row['path']);need(not rel.is_absolute() and '..' not in rel.parts,'unsafe manifest path')
        need(row['path'] not in seen,'duplicate manifest path');seen.add(row['path']);data=(root/rel).read_bytes()
        need(len(data)==row['bytes'],'size mismatch '+row['path']);need(hashlib.sha256(data).hexdigest()==row['sha256'],'hash mismatch '+row['path'])
    checks=[]
    if not a.hashes_only:
        for row in m.get('replay',[]):
            need(row['script'] in seen,'unpinned replay script')
            cmd=[sys.executable,'-B']+(['-O'] if sys.flags.optimize else [])+[str(root/row['script'])]+row.get('args',[])
            result=subprocess.run(cmd,cwd=root,capture_output=True,check=False)
            need(result.returncode==0,'replay failed '+row['script']+': '+result.stderr.decode(errors='replace'))
            sha=hashlib.sha256(result.stdout).hexdigest();need(sha==row['stdout_sha256'],'replay output mismatch '+row['script'])
            checks.append({'script':row['script'],'stdout_sha256':sha})
    print(json.dumps({'status':'PASS','scope':m['scope'],'files_verified':len(seen),'replay':checks,'assert_gates_used':False},sort_keys=True,indent=2))
if __name__=='__main__':main()
