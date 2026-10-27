#!/usr/bin/env python3
"""Producer-only explicit provenance refresh; no traversal or repository writes."""
from pathlib import Path
import hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parent
REPO=Path(sys.argv[1]) if len(sys.argv)>1 else Path('/workspace/scratch/07d5532f2fce/session73-repo')
BASE='92e7f12f3b95840cff92c5e6cf40b2ae8952cae4'
INCOMING='64a7985e7054821c5958dc85f6dd489c844111f2'
files=[]
for group,names in [
('05-knowledge/results/weighted_renewal_moving_modulus_20261026',['README.md','AUDIT.md','verify.py','SOURCE_PINS.json']),
('05-knowledge/results/clock_even_guarded_composition_20261026',['README.md','THEORY.md','CLOCKS-MACHINES.md','FOURIER.md','check.py','SOURCES.json'])]:
    for name in names:
        path=group+'/'+name
        files.append((BASE,path,'sources/'+path.replace('/','_')))
for name in ['S69.md','S70.md','S71.md','S72.md']:
    path='00-navigation/'+name
    files.append((BASE,path,'sources/'+path.replace('/','_')))
files.extend([
(BASE,'05-knowledge/results/native_return_necklaces_20261026/README.md','sources/native_return-README.md'),
(BASE,'05-knowledge/results/session72_even_graph_unifying_carriers_20261026/session72-response/REPORT.md','sources/session72-response-REPORT.md'),
(BASE,'05-knowledge/results/collatz_source_refinement_floor_20261005.md','sources/collatz_source_refinement_floor_20261005.md'),
(INCOMING,'05-knowledge/hypotheses/HYP-9419-unifying-candidates-three-layer-composition-rank-laws-dissection-obstructions-pattern-clocks.md','sources/HYP-9419-INCOMING.md'),
(INCOMING,'05-knowledge/results/three_quantities_fortyfour_five_dudeney_fermat_clock_arrows_20261026.md','sources/INCOMING-S62-note.md'),
(INCOMING,'01-canon/theorems/THM-4809-incidence-cycle-law-degree-bound-on-restored-guards-and-clock-pattern-circulants.md','sources/THM-4809-INCOMING.md')])
result=[]
for commit,path,local in files:
    data=subprocess.check_output(['git','-C',str(REPO),'show',commit+':'+path])
    p=ROOT/local
    if p.read_bytes()!=data:
        raise RuntimeError('Source snapshot mismatch: '+local)
    blob=subprocess.check_output(['git','-C',str(REPO),'rev-parse',commit+':'+path],text=True).strip()
    result.append({'commit':commit,'path':path,'blob':blob,'snapshot':local,
                   'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
(ROOT/'SOURCE_PINS.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':'PASS','source_snapshots':len(result)},sort_keys=True))
