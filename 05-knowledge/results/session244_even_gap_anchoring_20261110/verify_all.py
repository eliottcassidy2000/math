#!/usr/bin/env python3
"""Replay the sealed Session244 packet. Full order-12 regeneration is optional.

Default: all theorem checks in normal/-O, strict native controls, every stored
forest witness recounted. --full-census additionally regenerates all 12!
histories and compares every CSV, which takes roughly half an hour.
"""
import argparse
import json
from pathlib import Path
import resource
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
FOREST = ROOT / 'logic/forest_census/release'
JOBS = {
 'selector/verify.py':'selector/RESULT.json',
 'selector/verify_round2.py':'selector/ROUND2.json',
 'selector/verify_round3.py':'selector/ROUND3.json',
 'selector/verify_round4.py':'selector/ROUND4.json',
 'composition/verify.py':'composition/RESULTS.json',
 'composition/verify_round2.py':'composition/ROUND2_RESULTS.json',
 'composition/verify_round3.py':'composition/ROUND3_RESULTS.json',
 'composition/verify_structure.py':'composition/STRUCTURE_RESULTS.json',
 'composition/verify_interfaces.py':'composition/INTERFACES_CHECK.json',
 'logic/verify.py':'logic/CHECK.json',
 'logic/verify_round2.py':'logic/ROUND2_CHECK.json',
 'logic/verify_round3.py':'logic/ROUND3_CHECK.json',
 'dyadic/verify.py':'dyadic/CHECK.json',
 'audit/independent_verify.py':'audit/INDEPENDENT.json',
 'audit/verify_cores.py':'audit/CORES.json',
 'audit/final_independent_verify.py':'audit/FINAL_INDEPENDENT.json',
 'audit/verify_hostile217.py':'audit/HOSTILE217.json',
 'audit/verify_interfaces_independent.py':'audit/INTERFACES_INDEPENDENT.json',
 'audit/verify_exception97309.py':'audit/EXCEPTION97309.json',
 'audit/verify_normalization.py':'audit/NORMALIZATION.json',
 'logic/forest_census/release/verify_depth_exception.py':'logic/forest_census/release/DEPTH_EXCEPTION_CHECK.json',
}
FAIL = ['composition/verify.py','composition/verify_interfaces.py','dyadic/verify.py',
        'audit/final_independent_verify.py','audit/verify_hostile217.py',
        'audit/verify_interfaces_independent.py','audit/verify_exception97309.py','audit/verify_normalization.py']
FLAGS = ['-std=c++17','-O3','-Wall','-Wextra','-Werror','-pedantic']


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def run(command, good=True):
    p = subprocess.run(list(map(str,command)),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    need((p.returncode == 0) == good, ('command outcome',list(map(str,command)),p.returncode,p.stderr.decode(errors='replace')[-3000:]))
    return p.stdout


def python(script, optimized=False, args=(), good=True):
    return run([sys.executable,'-I','-B']+(['-O'] if optimized else [])+[ROOT/script,*args],good)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--full-census',action='store_true')
    ap.add_argument('--without-census',action='store_true',help='explicitly partial development check')
    args=ap.parse_args()
    need(not(args.full_census and args.without_census),'incompatible modes')
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    result={'status':'PARTIAL' if args.without_census else 'PASS','python_checks':[],'native_checks':[],'full_census_regenerated':False}
    for script,receipt in JOBS.items():
        expected=json.loads((ROOT/receipt).read_text())
        for opt in (False,True):
            need(json.loads(python(script,opt))==expected,('stored receipt mismatch',script,opt))
        if script in FAIL:
            for opt in (False,True):
                python(script,opt,['--fail'],good=False)
        result['python_checks'].append({'script':script,'normal':'PASS','optimized':'PASS','hostile':'PASS' if script in FAIL else 'not applicable'})
    with tempfile.TemporaryDirectory(prefix='session244-replay-') as td:
        td=Path(td)
        for source,receipt,name,kind in [('face/verify.cpp','face/CHECK.json','face','json'),('dyadic/probe.cpp','dyadic/PROBE.txt','dyadic','text')]:
            outputs=[]
            for mode in ('normal','ndebug'):
                binary=td/(name+'-'+mode)
                run(['g++',*FLAGS,*(['-DNDEBUG'] if mode=='ndebug' else []),ROOT/source,'-o',binary])
                outputs.append(run([binary]))
                if name=='face':run([binary,'--fail'],good=False)
            need(outputs[0]==outputs[1],('native mode mismatch',name))
            expected=(ROOT/receipt).read_bytes()
            need((json.loads(outputs[0])==json.loads(expected)) if kind=='json' else outputs[0]==expected,('native receipt mismatch',name))
            result['native_checks'].append({'source':source,'normal_and_ndebug':'PASS'})
        if not args.without_census:
            for opt in (False,True):
                checked=json.loads(python('logic/forest_census/release/audit.py',opt))
                need(checked==json.loads((FOREST/'AUDIT.json').read_text()),'forest Python audit receipt')
            control=json.loads(python('logic/forest_census/release/verify_controls.py'))
            need(control==json.loads((FOREST/'CONTROLS.json').read_text()),'forest strict controls receipt')
            binary=td/'all-witnesses'
            run(['g++',*FLAGS,'-DNDEBUG',FOREST/'audit_all_witnesses.cpp','-o',binary])
            witness_output=run([binary,'12',FOREST])
            need(witness_output==(FOREST/'ALL_WITNESSES.jsonl').read_bytes(),'all forest witness receipt')
            small=json.loads(python('audit/verify_small_forests.py',args=[FOREST]))
            need(small==json.loads((ROOT/'audit/SMALL_FORESTS.json').read_text()),'small independent permutations')
            result['native_checks'].append({'source':'logic/forest_census/release/census.cpp','order_8_regeneration_and_hostile_controls':'PASS','order_12_stored_witnesses':'PASS','order_9_normal_form_obstruction':'PASS'})
            shallow=td/'shallow'
            run(['g++',*FLAGS,'-DNDEBUG',FOREST/'shallow_census.cpp','-o',shallow])
            for normalized in ('0','1'):
                rec=json.loads(run([shallow,'8',FOREST,normalized]))
                need(rec['status']=='PASS' and rec['distinct_shallow_values']==192 and not rec['full_census_values_absent_shallow'],'small shallow replay')
            stored=json.loads((FOREST/'NORMALIZED_SHALLOW_AUDIT.json').read_text())
            need(stored['status']=='PASS' and stored['n']==12 and stored['histories']==23670073 and stored['distinct_shallow_values']==53548 and stored['full_census_values_absent_shallow']==[97309], 'stored normalized shallow receipt')
            result['native_checks'].append({'source':'logic/forest_census/release/shallow_census.cpp','order_8_both_modes':'PASS','order_12_receipt':'checked; full negative replay available with --full-census'})
            if args.full_census:
                binary=td/'census';out=td/'regenerated';out.mkdir()
                run(['g++',*FLAGS,'-DNDEBUG',FOREST/'census.cpp','-o',binary])
                output=run([binary,'12',out])
                need(output==(FOREST/'SUMMARY.jsonl').read_bytes(),'full census summary')
                for n in range(1,13):
                    name=f'n{n}_witnesses.csv'
                    need((out/name).read_bytes()==(FOREST/name).read_bytes(),('full census CSV',n))
                need(json.loads(run([shallow,'12',out,'1']))==stored,'full independent normalized shallow replay')
                result['full_census_regenerated']=True
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
