#!/usr/bin/env python3
"""Replay the independent Session134 audit from any cwd, normal and -O."""
import hashlib,json,os,pathlib,subprocess,sys,tempfile
ROOT=pathlib.Path(__file__).resolve().parent
AUDITS=(
 ('audit_two_apex.py','two_apex'),
 ('audit_two_apex_codes.py','two_apex_codes'),
 ('audit_migration.py','migration'),
 ('audit_migration_compiler.py','migration_compiler'),
 ('audit_zero_tail_migration.py','zero_tail_migration'),
 ('audit_migration_closure.py','migration_closure'),
 ('audit_boundaries.py','boundaries'),
 ('audit_core_budget.py','core_budget'),
 ('audit_spine_bound.py','spine_bound'),
 ('audit_nonnesting_family.py','nonnesting_family'),
 ('audit_radix_macros.py','radix_macros'),
)
def need(c,msg):
    if not c:raise RuntimeError(msg)
def main():
    manifest=json.loads((ROOT/'MANIFEST.json').read_text());checked=0
    for name,expected in manifest['sha256'].items():
        got=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        need(got==expected,'manifest mismatch: '+name);checked+=1
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1');records=[]
    with tempfile.TemporaryDirectory(prefix='session134-independent-replay-') as d:
        for script,stem in AUDITS:
            outputs=[]
            for flag,suffix in (([], 'normal'),(['-O'],'optimized')):
                result=subprocess.run([sys.executable,*flag,str(ROOT/script)],cwd=d,env=env,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                frozen=(ROOT/(stem+'.'+suffix+'.json')).read_bytes()
                need(result.stdout==frozen,'receipt mismatch: '+stem+'.'+suffix);outputs.append(result.stdout)
            need(outputs[0]==outputs[1],'normal/-O disagreement: '+stem)
            records.append({'audit':script,'receipt_sha256':hashlib.sha256(outputs[0]).hexdigest(),'normal_optimized_match':True})
    print(json.dumps({'status':'PASS','manifest_payload_files':checked,'python_consumers':records,'foreign_cwd_replay':True,'producer_imports':False,'global_completeness':'OPEN'},indent=2,sort_keys=True))
if __name__=='__main__':main()
