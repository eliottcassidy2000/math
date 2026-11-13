"""Replay every packet verifier normally and with optimization; compare mathematical receipts."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess,sys,json,hashlib,os
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
# file=None means a JSON receipt is emitted on stdout.
SPECS=[
 ('module/module_certificates.py','module/verification.json','module/verification.optimized.json'),
 ('module/escape_probes.py','module/escapes.json','module/escapes.optimized.json'),
 ('curves/verify_curves.py','curves/verification.json','curves/verification.json'),
 ('curves/module_audit.py','curves/module_audit.json','curves/module_audit.json'),
 ('curves/verify_881_profiles.py','curves/881_profiles.json','curves/881_profiles.json'),
 ('anchor/check.py',None,'anchor/receipt'),
 ('anchor/audit_loops.py',None,'anchor/audit_loops'),
 ('loops/verify.py',None,'loops/receipt'),
 ('loops/audit_anchor.py',None,'loops/audit_anchor')]

def run(spec):
    script,normal,optimized=spec;receipts=[]
    for opt in (False,True):
        args=[sys.executable]+(['-O'] if opt else [])+['-B',str(HERE/script)]
        env=os.environ.copy();env['PYTHONIOENCODING']='utf-8'
        result=subprocess.run(args,cwd=ROOT,env=env,text=True,encoding='utf-8',capture_output=True)
        if result.returncode:raise RuntimeError(script+' failed: '+result.stderr+' '+result.stdout)
        if normal is None:
            data=json.loads(result.stdout)
            name=optimized+('.optimized.json' if opt else '.json')
            (HERE/name).write_text(json.dumps(data,sort_keys=True,indent=2)+'\n',encoding='utf-8',newline='\n')
        else:data=json.loads((HERE/(optimized if opt else normal)).read_text(encoding='utf-8'))
        receipts.append(data)
    if receipts[0]!=receipts[1]:raise RuntimeError(script+' normal/optimized mismatch')
    print('PASS',script,flush=True)
    return dict(script=script,normal_optimized_equal=True,script_text_lf_sha256=hashlib.sha256((HERE/script).read_bytes().replace(b'\r\n',b'\n')).hexdigest(),mathematical_receipt_sha256=hashlib.sha256(json.dumps(receipts[0],sort_keys=True,separators=(',',':')).encode()).hexdigest())

if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(run,SPECS))
    result=dict(status='PASS',engines=len(rows),executions=2*len(rows),python=sys.version.split()[0],checks=rows)
    (HERE/'REPLAY.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS',len(rows),'engines;',2*len(rows),'executions')
