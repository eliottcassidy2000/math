#!/usr/bin/env python3
"""Bind the integral-order addendum to its independent audit and source packet."""
import ast,hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent

def need(ok,message):
 if not ok:raise RuntimeError(message)

def check(path,row):
 f=root/path;need(f.is_file() and not f.is_symlink(),'missing audited file '+path);data=f.read_bytes()
 need(len(data)==row['bytes'],'audited byte length '+path)
 need(hashlib.sha256(data).hexdigest()==row['sha256'],'audited SHA256 '+path)
 return data
pins=json.loads((root/'independent_integral_order/reviewed_source_pins.json').read_text())
for row in pins['files']:check(pins['producer_bundle']+'/'+row['path'],row)
manifests=0;files=0
for lane in ('cyclic_integral_order','independent_integral_order'):
 m=json.loads((root/lane/'manifest.json').read_text())
 for row in m['files']:check(lane+'/'+row['path'],row);files+=1
 manifests+=1
source=json.loads((root/'cyclic_integral_order/source_pins.json').read_text())
sha=hashlib.sha256((root/'tournaments/manifest.json').read_bytes()).hexdigest()
need(source['original_packet_manifest_sha256']==pins['original_tournaments_manifest_sha256']==sha,'original rational-source manifest identity')
f20=json.loads((root/'independent_f20_embedding/reviewed_source_pins.json').read_text())
for row in f20['files']:check(f20['producer_bundle']+'/'+row['path'],row)
for row in f20['unchanged_dependencies']:
 need(hashlib.sha256((root/row['packet']/'manifest.json').read_bytes()).hexdigest()==row['manifest_sha256'],'F20 unchanged dependency '+row['packet'])
scripts=('f20_embedding/check_corollary.py','f20_embedding/test_corollary.py','independent_f20_embedding/check_f20.py','cyclic_integral_order/check_integral_order.py','cyclic_integral_order/test_negative_controls.py','independent_integral_order/check_integral_order.py','independent_integral_order/audit_producer_mutations.py')
for script in scripts:need(not any(isinstance(n,ast.Assert) for n in ast.walk(ast.parse((root/script).read_text()))),'optimization-sensitive assert '+script)
print(json.dumps({'status':'PASS','audited_integral_addendum_files':len(pins['files']),'audited_f20_addendum_files':len(f20['files']),'lane_manifests':manifests,'lane_files':files,'optimization_safe_scripts':len(scripts),'inherited_tournament_manifest_sha256':sha,'scope':'Exact proof/checker binding to additive integral-order and F20 audits; existing science untouched, external subgroup/character facts remain cited.'},indent=2,sort_keys=True))
