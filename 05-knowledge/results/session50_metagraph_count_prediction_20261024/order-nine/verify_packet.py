#!/usr/bin/env python3
"""Manifest integrity and manifest-only full standalone relocation replay."""
import argparse, ast, hashlib, json, shutil, subprocess, sys, tempfile
from pathlib import Path
CHECKS=0
def require(c,m):
    global CHECKS
    CHECKS+=1
    if not c:raise RuntimeError(m)
def run(command,cwd):
    result=subprocess.run(command,cwd=cwd,check=True,text=True,capture_output=True)
    return result.stdout.strip()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--full',action='store_true');ap.add_argument('--output',default=None);args=ap.parse_args()
    here=Path(__file__).resolve().parent;manifest=json.loads((here/'manifest.json').read_text())
    for item in manifest['files']:
        path=here/item['path'];data=path.read_bytes();require(len(data)==item['bytes'] and hashlib.sha256(data).hexdigest()==item['sha256'],'manifest member '+item['path'])
    normal=json.loads((here/'normal.json').read_text());optimized=json.loads((here/'optimized.json').read_text());compiled=json.loads((here/'upper-compiled.json').read_text());regular=json.loads((here/'regular-normal.json').read_text())
    require(normal==optimized,'normal and -O exact equality')
    require(regular==json.loads((here/'regular-optimized.json').read_text()),'regular normal and -O exact equality')
    for key in ('pointed_extension_states','upper_refinement_trace','stable_root_cell_histogram','per_parent'):require(normal[key]==compiled[key],'independent compiled '+key)
    require(normal['regular_sector']=={'pointed_states':regular['pointed_regular_states'],'trace':regular['regular_trace'],'classes':regular['regular_classes'],'error':regular['regular_error']},'separate regular trace implementation')
    require(normal['upper_refinement_trace']-normal['total_error_correction']==9*normal['predicted_class_count'],'final exact trace arithmetic')
    # No calibration using target answer, target upper trace, trace correction,
    # or regular-sector output is permitted in executable prediction source.
    forbidden={normal['predicted_class_count'],normal['upper_refinement_trace'],normal['exact_unpointing_trace'],normal['total_error_correction'],regular['regular_error'],regular['regular_trace'],regular['regular_classes'],regular['pointed_regular_states'],normal['pointed_extension_states']}
    tree=ast.parse((here/'predictor.py').read_text())
    constants={n.value for n in ast.walk(tree) if isinstance(n,ast.Constant) and type(n.value)is int}
    require(not forbidden&constants,'no derived target or correction constants in predictor')
    require(normal['child_canonicalizations']==0 and normal['general_child_reroots']==0 and compiled['child_deletions']==0,'target operation boundary receipts')
    replay=None
    if args.full:
        with tempfile.TemporaryDirectory(prefix='order-nine-manifest-only-') as temp:
            root=Path(temp)
            for item in manifest['files']:
                dst=root/item['path'];dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(here/item['path'],dst)
            shutil.copyfile(here/'manifest.json',root/'manifest.json')
            stdout=[]
            stdout.append(run([sys.executable,'prepare_inputs.py'],root))
            require((root/'parents-eight.json').read_bytes()==(here/'parents-eight.json').read_bytes(),'relocated complete parent atlas regenerated exactly')
            require((root/'parents-eight.txt').read_bytes()==(here/'parents-eight.txt').read_bytes(),'relocated compiled parent atlas regenerated exactly')
            stdout.append(run([sys.executable,'predictor.py','--output','replayed-normal.json'],root))
            stdout.append(run([sys.executable,'-O','predictor.py','--output','replayed-optimized.json'],root))
            require(json.loads((root/'replayed-normal.json').read_text())==normal,'relocated normal prediction exact')
            require(json.loads((root/'replayed-optimized.json').read_text())==normal,'relocated -O prediction exact')
            stdout.append(run([sys.executable,'regular_trace.py','--output','replayed-regular.json'],root))
            require(json.loads((root/'replayed-regular.json').read_text())==regular,'relocated regular trace exact')
            stdout.append(run(['g++','-O3','-std=c++17','-Wall','-Wextra','-pedantic','upper_trace.cpp','-o','upper_trace'],root))
            stdout.append(run(['./upper_trace','parents-eight.txt','replayed-compiled.json'],root))
            require(json.loads((root/'replayed-compiled.json').read_text())==compiled,'relocated compiled control exact')
            replay={'status':'PASS','copied_only_manifest_members':True,'stdout':stdout}
    result={'status':'PASS','checks':CHECKS,'manifest_members':len(manifest['files']),'manifest_sha256':hashlib.sha256((here/'manifest.json').read_bytes()).hexdigest(),'full_relocation':replay}
    if args.output:Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
