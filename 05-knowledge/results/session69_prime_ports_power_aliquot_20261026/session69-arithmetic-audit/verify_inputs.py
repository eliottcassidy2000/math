#!/usr/bin/env python3
"""Read-only input, receipt, negative-control, and complete-fiber cross-check."""
import argparse, hashlib, json, os, subprocess, sys
from pathlib import Path

def need(ok, why):
    if not ok:raise RuntimeError(why)

ap=argparse.ArgumentParser()
ap.add_argument('--power',required=True,type=Path)
ap.add_argument('--aliquot',required=True,type=Path)
args=ap.parse_args()
here=Path(__file__).resolve().parent
manifest=json.loads((args.power/'MANIFEST.json').read_text())
for row in manifest['files']:
    b=(args.power/row['path']).read_bytes()
    need(len(b)==row['bytes'] and hashlib.sha256(b).hexdigest()==row['sha256'],'power manifest '+row['path'])
for line in (args.aliquot/'SHA256SUMS').read_text().splitlines():
    digest,name=line.split(maxsplit=1);name=name.lstrip('*')
    need(hashlib.sha256((args.aliquot/name).read_bytes()).hexdigest()==digest,'aliquot manifest '+name)
for lane,source in [('power',args.power/'results.json'),('aliquot',args.aliquot/'normal.json')]:
    need((here/(lane+'.normal.json')).read_bytes()==source.read_bytes(),lane+' supplied receipt')
    need((here/(lane+'.normal.json')).read_bytes()==(here/(lane+'.optimized.json')).read_bytes(),lane+' optimized receipt')
need((here/'independent.normal.json').read_bytes()==(here/'independent.optimized.json').read_bytes(),'independent optimized receipt')
for lane,error in [('aliquot','negative-control target6 completeness'),('independent','negative-control canonical first-six')]:
    for mode in ('normal','optimized'):
        need(error in (here/(lane+'.negative.'+mode+'.stderr')).read_text(),lane+' '+mode+' intended negative gate')
        need(not (here/(lane+'.negative.'+mode+'.stdout')).read_bytes(),lane+' '+mode+' no success receipt')
independent=json.loads((here/'independent.normal.json').read_text())
code='''import hashlib,json,sys
sys.path.insert(0,sys.argv[1])
from aliquot_fibers import fiber
x={str(m):fiber(m) for m in range(2,1001)}
print(hashlib.sha256(json.dumps(x,sort_keys=True,separators=(",",":")).encode()).hexdigest())
'''
env=dict(os.environ);env['PYTHONDONTWRITEBYTECODE']='1'
proc=subprocess.run([sys.executable,'-c',code,str(args.aliquot/'source')],env=env,text=True,capture_output=True,check=True)
need(proc.stdout.strip()==independent['complete_aliquot_fibers_sha256'],'full 2..1000 independently computed fibers differ')
inputs=[]
for label,root in [('power',args.power),('aliquot',args.aliquot)]:
    names=['REPORT.md','check.py']+(['MANIFEST.json','results.json','results.optimized.json'] if label=='power' else ['SOURCE_PINS.json','SHA256SUMS','source/aliquot_fibers.py','normal.json','optimized.json'])
    for name in names:
        b=(root/name).read_bytes();inputs.append({'lane':label,'path':name,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
(here/'INPUTS.json').write_text(json.dumps({'source_commit':manifest['source_commit'],'files':inputs},indent=2)+'\n')
print(json.dumps({'status':'PASS','power_manifest_entries':len(manifest['files']),
 'aliquot_manifest_entries':len((args.aliquot/'SHA256SUMS').read_text().splitlines()),
 'original_receipts_match':True,'normal_optimized_receipts_match':True,
 'negative_control_gates_match':True,'full_aliquot_target_range':[2,1000],
 'independent_complete_fibers_match_inherited_algorithm':True,
 'complete_fibers_sha256':proc.stdout.strip(),'independent_active_checks':independent['checks']},indent=2)+'\n')
