#!/usr/bin/env python3
"""Replay producer and independent checks, then repeat from a relocated package."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parent

def require(ok,message):
    if not ok:raise RuntimeError(message)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--skip-relocation',action='store_true',help='Internal: do not recursively repeat relocation')
    args=parser.parse_args();work=args.work_dir.resolve();work.mkdir(parents=True,exist_ok=True)
    require(ROOT!=work and ROOT not in work.parents,'Work directory must be outside the package')
    receipts=work/'receipts';receipts.mkdir(exist_ok=True)
    started=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());records=[]
    def run(name,command):
        t=time.monotonic()
        result=subprocess.run([str(x) for x in command],cwd=work,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True)
        (receipts/(name+'.out')).write_text(result.stdout)
        if result.stderr:(receipts/(name+'.stderr')).write_text(result.stderr)
        require(result.returncode==0,(name,result.returncode,result.stderr[-2000:]))
        records.append({'check':name,'returncode':0,'seconds':round(time.monotonic()-t,3)})
        print(name,'PASS',flush=True)
        return result.stdout
    def pair(name,script,extra=()):
        a=run(name+'.normal',[sys.executable,script,*extra])
        b=run(name+'.optimized',[sys.executable,'-O',script,*extra])
        require(a==b,('Normal/-O output differs',name))
    for stem in ('check','check_macro_inverse'):pair('producer.'+stem,ROOT/(stem+'.py'))
    for mode in ('normal','optimized'):
        flags=[] if mode=='normal' else ['-O']
        run('producer.fibers.'+mode,[sys.executable,*flags,ROOT/'fibers.py','--output',work/('fibers.'+mode+'.json')])
        require((work/('fibers.'+mode+'.json')).read_bytes()==(ROOT/'fiber_boundaries.json').read_bytes(),'Regenerated fiber data differ')
    require((receipts/'producer.fibers.normal.out').read_bytes()==(receipts/'producer.fibers.optimized.out').read_bytes(),'Fiber outputs differ')
    pair('independent',ROOT/'audit/independent_check.py')
    pair('structure',ROOT/'audit/structure_check.py')
    # Current mathematical audit outputs must match the original independent
    # audit exactly, apart from the intentionally changed producer source hash.
    for stem,oldname in (('independent','original_independent.json'),('structure','original_structure.json')):
        old=json.loads((ROOT/'provenance'/oldname).read_text())
        new=json.loads((receipts/(stem+'.normal.out')).read_text())
        old.pop('source_sha256',None);new.pop('source_sha256',None)
        require(old==new,('Independent mathematical audit output changed',stem))
    graphs=sorted((ROOT/'audit').glob('adjacency_*.txt'))
    require(len(graphs)==4,'Require all four actual frontier matrices')
    for opt in ('O0','O2'):
        run('compile.frontier.'+opt,['g++','-'+opt,'-std=c++17','-Wall','-Wextra','-pedantic',ROOT/'audit/adjacency_dp.cpp','-o',work/('frontier.'+opt)])
        run('frontier.'+opt,[work/('frontier.'+opt),*graphs])
    require((receipts/'frontier.O0.out').read_bytes()==(receipts/'frontier.O2.out').read_bytes(),'Frontier O0/O2 differs')
    pair('frontier_mutation',ROOT/'audit/mutation.py',[work/'frontier.O2'])
    relocation=None
    if not args.skip_relocation:
        relocated=work/'relocated-package'
        require(not relocated.exists(),'Use a fresh work directory for relocation')
        shutil.copytree(ROOT,relocated,ignore=shutil.ignore_patterns('payloads','__pycache__'))
        relocated_work=work/'relocated-run'
        run('relocation',[sys.executable,relocated/'run_checks.py','--work-dir',relocated_work,'--skip-relocation'])
        other=relocated_work/'receipts'
        common=[p.name for p in receipts.glob('*.out') if p.name!='relocation.out']
        for name in common:require((receipts/name).read_bytes()==(other/name).read_bytes(),('Relocation output differs',name))
        relocated_verification=json.loads((other/'verification.json').read_text())
        require(relocated_verification['status']=='PASS','Relocated replay did not pass')
        relocation={'status':'PASS','all_stdout_receipts_identical':len(common),'full_normal_optimized_replay':True,'cwd_outside_package':True}
    fingerprint={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(ROOT.rglob('*')) if p.is_file() and p.suffix in ('.py','.cpp')}
    inputs={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('actual_seed_roots.json','fiber_boundaries.json','frontier_repairs.jsonl')}
    result={'status':'PASS','started_utc':started,'finished_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'python':platform.python_version(),'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],'checks':records,'source_sha256':fingerprint,'input_sha256':inputs,'relocation':relocation,'scope':'Exact actual-seed language identities, fibers and finite inverse completeness. No universal positive-answer or infinite tournament spectrum theorem.'}
    (receipts/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print('ALL CHECKS PASS',flush=True)
if __name__=='__main__':main()
