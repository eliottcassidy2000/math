#!/usr/bin/env python3
"""Replay accepted checkpoint-one kernels from a fresh relocated directory.
Use: python3 replay.py [--output report.json]. All controls survive Python -O.
"""
from pathlib import Path
import argparse, hashlib, json, os, shutil, subprocess, sys, tempfile

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);args=ap.parse_args()
    src=Path(__file__).resolve().parent
    def check(ok,label):
        if not ok:raise RuntimeError(label)
    scripts=['stability/verify_n748.py','independent_threshold/audit_n748.py','independent/audit_bridges.py','riemann_five_eleven/verify_history.py']
    runs=[]
    with tempfile.TemporaryDirectory(prefix='session139-checkpoint1-relocated-') as name:
        dst=Path(name)
        for f in scripts:
            (dst/f).parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src/f,dst/f)
        mutant=(src/scripts[0]).read_text().replace('h=eta-C*eta*eta','h=F(0)')
        check(mutant!=(src/scripts[0]).read_text(),'mutation applied')
        (dst/'mutated_producer.py').write_text(mutant)
        outputs={}
        def run(label,script,flags=(),expected=0,needle=None,opt=False):
            command=[sys.executable]+(['-O'] if opt else [])+[script,*flags]
            proc=subprocess.run(command,cwd=dst,text=True,capture_output=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},timeout=300)
            check(proc.returncode==expected,f'{label}: exit {proc.returncode}: {proc.stderr[-1500:]}')
            if needle:check(needle in proc.stdout+proc.stderr,f'{label}: missing expected rejection')
            runs.append({'label':label,'mode':'optimized' if opt else 'normal','exit_code':proc.returncode,'stdout_sha256':hashlib.sha256(proc.stdout.encode()).hexdigest(),'expected_rejection':expected!=0})
            return proc.stdout
        for opt in (False,True):
            for f in scripts:outputs[f,opt]=run(f,f,opt=opt)
            run('threshold747-negative',scripts[1],('--threshold','747'),1,'shared_budget_strict_contradiction',opt)
            run('producer-h-zero-negative','mutated_producer.py',expected=1,needle='shared_budget_contradiction',opt=opt)
            for mutation in ('hamiltonian','paley'):
                run('bridge-'+mutation+'-negative',scripts[2],('--mutate',mutation),1,'RuntimeError',opt)
            for mutation in ('tetranacci','paley','incidence'):
                run('history-'+mutation+'-negative',scripts[3],('--mutate',mutation),1,'RuntimeError',opt)
        for f in scripts:check(outputs[f,False]==outputs[f,True],f+' output identity')
        check(outputs[scripts[0],False]==(src/'stability/verification_n748.json').read_text(),'producer frozen output')
        check(outputs[scripts[1],False]==(src/'independent_threshold/normal.json').read_text(),'independent frozen output')
        check(outputs[scripts[2],False]==(src/'independent/audit_bridges.out').read_text(),'bridge frozen output')
        check(outputs[scripts[3],False]==(src/'riemann_five_eleven/history.normal.json').read_text(),'history exact output')
        report={'status':'PASS','fresh_relocation':True,'producer_imported_by_independent_audit':False,'normal_optimized_identical':scripts,'scientific_executions':len(runs),'positive_runs':8,'expected_rejections':14,'independent_endpoint_named_checks':len(json.loads(outputs[scripts[1],False])['active_checks']),'source_pins':{f:hashlib.sha256((src/f).read_bytes()).hexdigest() for f in scripts},'runs':runs}
    body=json.dumps(report,sort_keys=True,indent=2)+'\n'
    if args.output:args.output.write_text(body)
    print(body,end='')
if __name__=='__main__':main()
