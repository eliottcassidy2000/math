#!/usr/bin/env python3
"""Verify audit, frozen native packet and original dependency hashes."""
import argparse,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--packet',type=Path,required=True);a=p.parse_args();here=Path(__file__).resolve().parent
count=0
def check(path,entry):
 global count
 data=path.read_bytes()
 if len(data)!=entry['bytes']or hashlib.sha256(data).hexdigest()!=entry['sha256']:raise RuntimeError('manifest mismatch: '+str(path))
 count+=1
manifest=json.loads((here/'MANIFEST.json').read_text())
for scope,root in [('audit',here),('producer',a.packet),('repository',a.repo)]:
 for entry in manifest[scope]:check(root/entry['path'],entry)
for entry in json.loads((here/'PUBLISH_MANIFEST.json').read_text())['files']:check(here/entry['path'],entry)
for entry in json.loads((a.packet/'MANIFEST.json').read_text())['files']:check(a.packet/entry['path'],entry)
for n,o in [('normal.json','optimized.json'),('dyadic-replay.normal.json','dyadic-replay.optimized.json'),('producer.normal.json','producer.optimized.json'),('cost-boundary.normal.json','cost-boundary.optimized.json'),('cost-buckets.normal.json','cost-buckets.optimized.json'),('cycles.normal.json','cycles.optimized.json')]:
 if(here/n).read_bytes()!=(here/o).read_bytes():raise RuntimeError('normal/optimized mismatch: '+n)
if(here/'producer.normal.json').read_bytes()!=(a.packet/'normal.json').read_bytes():raise RuntimeError('producer replay differs from frozen output')
report=json.loads((here/'normal.json').read_text());old=json.loads((a.repo/'05-knowledge/results/session45_injective_multisource_moves_20261023/audit-dyadic/normal.json').read_text());replay=json.loads((here/'dyadic-replay.normal.json').read_text())
if any(old[k]!=replay[k]for k in old if k!='source_hashes'):raise RuntimeError('original dyadic audit mathematical replay mismatch')
ledger=json.loads((here/'proper_balanced_separators.json').read_text())
for r in report['circuits']:
 selected=ledger['circuits'][str(r['slots']//2)];sha=hashlib.sha256(json.dumps(selected,separators=(',',':')).encode()).hexdigest()
 if len(selected)!=r['proper_balanced_rejected']or sha!=r['separator_ledger_sha256']:raise RuntimeError('separator ledger mismatch')
census=json.loads((here/'cost-boundary.normal.json').read_text());buckets=json.loads((here/'cost-buckets.normal.json').read_text())
if hashlib.sha256((here/'cost-boundary.normal.json').read_bytes()).hexdigest()!=buckets['source_sha256']['census']:raise RuntimeError('cost bucket census pin')
for c,b in zip(census['results'],buckets['results']):
 if c['cost']!=b['cost']or c['complete_relation_set_sha256']!=b['complete_relation_set_sha256']:raise RuntimeError('relation set crosscheck')
case=report['circuits'][0];partition=sorted([sorted(case['plus']),sorted(case['minus'])])
if census['results'][1]['canonical_primitive_partitions']!=[partition]:raise RuntimeError('unique relation native realization')
print(json.dumps({'status':'PASS','verified_file_entries':count,'native_checks':report['checks'],'filtration_replay_checks':replay['checks']},sort_keys=True))
