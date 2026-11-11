#!/usr/bin/env python3
"""Replay only Session257's four finite new-claim engines."""
import argparse, hashlib, json, pathlib, subprocess, sys
p=argparse.ArgumentParser();p.add_argument('--optimized',action='store_true');a=p.parse_args();a.optimized = a.optimized or not __debug__
root=pathlib.Path(__file__).resolve().parent
jobs=[('chambers/check.py','chambers/results.json',False),('permutohedron/reproduce.py','permutohedron/certificate.json',False),('gram/check.py','gram/results.json',False),('filtration/check.py','filtration/exact.json',True)]
records=[]
for script,receipt,stdout_json in jobs:
    path=root/script; out=root/receipt
    before=out.read_bytes() if out.exists() else None
    cmd=[sys.executable]+(['-O'] if a.optimized else [])+['-B',path.name]
    done=subprocess.run(cmd,cwd=path.parent,capture_output=True,text=True)
    if done.returncode:raise RuntimeError(script+' failed:\n'+done.stderr+'\n'+done.stdout)
    if stdout_json:
        parsed=json.loads(done.stdout)
        if before is not None and json.loads(before)!=parsed:raise RuntimeError(script+' differs from saved receipt')
        if before is None:out.write_text(json.dumps(parsed,indent=2)+'\n')
    after=out.read_bytes()
    if before is not None and before!=after:raise RuntimeError(script+' changed its deterministic receipt')
    records.append({'script':script,'receipt':receipt,'status':'PASS','source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'receipt_sha256':hashlib.sha256(after).hexdigest()})
record={'scope':'Four Session257 new-claim engines; historical inputs are inherited','optimized':a.optimized,'checks':'Explicit runtime checks remain active with -O','engines':records}
(root/('REPLAY_OPTIMIZED.json' if a.optimized else 'REPLAY.json')).write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
