#!/usr/bin/env python3
"""Verify frozen bytes and replay CP2 in a temporary copy, normally or under -O."""
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def need(ok,message):
    if not ok:raise RuntimeError(message)

def normalized_precision(data):
    data=json.loads(data);data.pop('optimized_python',None)
    for row in data['certificates']:row.pop('seconds',None)
    return data

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--hashes-only',action='store_true');args=parser.parse_args()
    manifest=json.loads((ROOT/'MANIFEST.json').read_text());seen=set()
    for row in manifest['files']:
        path=Path(row['path']);need(not path.is_absolute() and '..' not in path.parts,'unsafe manifest path')
        need(row['path'] not in seen,'duplicate manifest path');seen.add(row['path']);data=(ROOT/path).read_bytes()
        need(len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'],'frozen bytes changed '+row['path'])
    records=[]
    if not args.hashes_only:
        with tempfile.TemporaryDirectory(prefix='s144-cp2-replay-') as temp:
            dest=Path(temp)
            for rel in seen:
                p=dest/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/rel,p)
            def run(script,argv=()):
                path=dest/script;cmd=[sys.executable,'-B']+(['-O'] if sys.flags.optimize else [])+[str(path),*argv]
                r=subprocess.run(cmd,cwd=path.parent,capture_output=True)
                need(r.returncode==0,'replay failed '+script+': '+r.stderr.decode(errors='replace'))
                return r.stdout
            pairs=[
                ('jacobian_followup/verify_padded_boundary.py','jacobian_followup/padded_boundary.json'),
                ('jacobian_followup/verify_equivariant_gauge.py','jacobian_followup/equivariant_gauge.json'),
                ('jacobian_followup/verify_coordinate_fibers.py','jacobian_followup/coordinate_fibers.json'),
                ('jacobian_followup/independent/verify_quartic_boundary.py','jacobian_followup/independent/quartic_checks.json'),
                ('eth/verify_eth_interfaces.py','eth/verification.json'),
                ('mrb_independent/check_mrb_boundary.py','mrb_independent/mrb_boundary.normal.json'),
                ('mrb/check.py','mrb/check.out'),
                ('mrb/test_rejections.py','mrb/rejections.out'),
            ]
            for script,expected in pairs:
                out=run(script);need(out==(ROOT/expected).read_bytes(),'output mismatch '+script)
                records.append({'script':script,'output_sha256':hashlib.sha256(out).hexdigest()})
            need((dest/'mrb/results.json').read_bytes()==(ROOT/'mrb/results.normal.json').read_bytes(),'MRB mathematical output changed')
            run('mrb/compare_certificates.py');need((dest/'mrb/certificate_comparison.json').read_bytes()==(ROOT/'mrb/certificate_comparison.json').read_bytes(),'MRB certificate comparison changed')
            precision_path=dest/'mrb/audit/replayed_precision.json'
            run('mrb/audit/exact_precision_audit.py',['--precisions','1','8','16','32','--output',str(precision_path)])
            norm=normalized_precision(precision_path.read_text())
            need(norm==normalized_precision((ROOT/'mrb/audit/exact_precision_results.json').read_text()),'rational precision certificate changed')
            need(norm==json.loads((ROOT/'mrb_independent/precision.normalized.json').read_text()),'independent rational precision certificate changed')
            digest=hashlib.sha256(json.dumps(norm,sort_keys=True,separators=(',',':')).encode()).hexdigest()
            records.append({'script':'mrb/audit/exact_precision_audit.py','precisions':[1,8,16,32],'normalized_sha256':digest,'metadata_ignored':['optimized_python','seconds']})
    print(json.dumps({'status':'PASS','files_verified':len(seen),'assert_gates_used':False,'replays':records,'frozen_files_modified':False,'scope':manifest['scope']},indent=2,sort_keys=True))
if __name__=='__main__':main()
