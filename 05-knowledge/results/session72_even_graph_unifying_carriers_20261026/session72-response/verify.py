#!/usr/bin/env python3
"""Verify the fixed local artifact set, deterministic replay, and hostile gates."""
import hashlib, json, pathlib, subprocess, sys
ROOT=pathlib.Path(__file__).resolve().parent

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def run(args):
    return subprocess.run([sys.executable,'-I','-B']+args,cwd=ROOT,capture_output=True,text=True)

def main():
    manifest=json.loads((ROOT/'FILE_MANIFEST.json').read_text())
    for name,h in manifest.items():
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,'artifact hash mismatch: '+name)
    expected=(ROOT/'normal.json').read_text()
    need(expected==(ROOT/'optimized.json').read_text(),'frozen normal/optimized receipts differ')
    replays=[]
    for mode in [[],['-O']]:
        r=run(mode+['check.py'])
        need(r.returncode==0,'positive replay failed: '+r.stderr)
        need(r.stdout==expected,'positive replay differs from frozen receipt')
        replays.append({'optimized':bool(mode),'returncode':r.returncode,'matches_frozen':True})
    hostile=[]
    for mutation in ['rank','trace','cone','return']:
        for mode in [[],['-O']]:
            r=run(mode+['check.py','--max-n','5','--corrupt',mutation])
            need(r.returncode!=0 and 'intentional corruption:' in r.stderr,'hostile gate failed: '+mutation)
            hostile.append({'mutation':mutation,'optimized':bool(mode),'returncode':r.returncode,'expected_exception':True})
    print(json.dumps({'manifest_files_verified':len(manifest),'positive_replays':replays,'hostile_gates':hostile},sort_keys=True,indent=2))
if __name__=='__main__':main()
