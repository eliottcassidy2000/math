#!/usr/bin/env python3
from pathlib import Path
import subprocess,json,hashlib
BASE=Path(__file__).resolve().parent
REPO=Path('/tmp/session142-source')
COMMIT='9d5bae74f05a685d523df0c5559d71c3cd53cd6a'
root='05-knowledge/results/'
files={
 'square_sum_spectrum_overnight_20261031c':['D_profiles.tsv','DD_profiles.tsv']+[f'DT{k}D_profiles.tsv' for k in range(1,6)],
 'square_sum_spectrum_overnight_20261031d':['README.md','VERIFICATION.json','MANIFEST.json','bank_query.py','source_inverse.py','source_inverse.json','large_point_queries.jsonl','large_point_adjacency_audit.json','audit_large_points.py','capped_scan.cpp','audit_capped_scan.cpp'],
 'session124_million_path_spectrum_20261031':['README.md','construct_witness.py','constructors/seed_roots.json','constructors/block_target_witnesses.jsonl'],
}
manifest=[]
for package,names in files.items():
 for name in names:
  path=root+package+'/'+name
  data=subprocess.check_output(['git','-C',str(REPO),'show',COMMIT+':'+path],timeout=30)
  dest=BASE/'sources'/package/name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
  oid=subprocess.check_output(['git','-C',str(REPO),'rev-parse',COMMIT+':'+path],text=True).strip()
  manifest.append({'repository_path':path,'commit':COMMIT,'blob':oid,'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data),'local_snapshot':str(dest.relative_to(BASE))})
  print(path,len(data),flush=True)
(BASE/'sources.json').write_text(json.dumps(manifest,indent=2)+'\n')
