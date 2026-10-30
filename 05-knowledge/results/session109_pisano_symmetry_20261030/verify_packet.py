#!/usr/bin/env python3
"""Isolated, optimization-safe replay of the complete Session109 packet."""
import ast,hashlib,json,pathlib,shutil,subprocess,sys,tempfile
ROOT=pathlib.Path(__file__).resolve().parent
CASES=[('affine/check.py','affine/results.json'),('clocks/check_clocks.py','clocks/results.json'),('clocks/check_metallic.py','clocks/metallic_results.json'),('tournaments/check.py','tournaments/check.json'),('audit/check_ramified_attachment.py','audit/ramified_attachment_audit.json'),('audit_clocks/check_independent.py','audit_clocks/results.json'),('synthesis/check.py','synthesis/results.json')]
def check(ok,msg):
    if not ok:raise RuntimeError(msg)
manifest=json.loads((ROOT/'manifest.json').read_text())
for f,digest in manifest.items():
    check(hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==digest,'manifest '+f)
for script,_ in CASES:
    check(not any(isinstance(node,ast.Assert) for node in ast.walk(ast.parse((ROOT/script).read_text()))),'optimized-away assertions '+script)
rows=[]
for script,expected in CASES:
    ref=json.loads((ROOT/expected).read_text()); outputs=[]
    for opt in (False,True):
        with tempfile.TemporaryDirectory(prefix='pisano-replay-') as temp:
            dest=pathlib.Path(temp)/'packet'
            dest.mkdir()
            for f in manifest:
                target=dest/f;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/f,target)
            command=[sys.executable]+(['-O'] if opt else [])+[str(dest/script)]
            proc=subprocess.run(command,cwd=dest,capture_output=True,text=True,timeout=180)
            check(proc.returncode==0,script+' '+proc.stderr)
            output=json.loads(proc.stdout)
            check(output==ref,'recorded output differs '+script)
            outputs.append(output)
    check(outputs[0]==outputs[1],'optimization mismatch '+script)
    rows.append({'script':script,'recorded':expected,'normal':'PASS','optimized':'PASS','output_sha256':hashlib.sha256(json.dumps(ref,sort_keys=True,separators=(',',':')).encode()).hexdigest()})
report={'status':'PASS','isolated_replays':2*len(rows),'manifest_entries':len(manifest),'checks_remain_active_under_optimization':True,'cases':rows,'scope':'Exact finite tests plus independent proof reviews; no CI claim'}
print(json.dumps(report,indent=2))
