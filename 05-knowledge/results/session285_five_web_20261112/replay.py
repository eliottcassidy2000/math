#!/usr/bin/env python3
"""Replay exact finite controls; no claim this alone proves all-order statements."""
import gzip,hashlib,json,pathlib,platform,subprocess,sys
import sympy
P=pathlib.Path(__file__).resolve().parent
pairs=[('stress/verify.py','stress/receipt.json'),
       ('bridges/verify_periods.py','bridges/periods.json'),
       ('tournaments/verify.py','tournaments/replay.json'),
       ('tournaments/scan_orbits.py','tournaments/orbit_scan.json'),
       ('holonomy/verify.py','holonomy/RESULTS.json'),
       ('history/check_mechanisms.py','history/mechanisms.normal.json'),
       ('audit/independent.py','audit/independent.json'),
       ('audit/loop_rank.py','audit/loop_rank.json'),
       ('audit/tournament_counts.py','audit/tournament_counts.json'),
       ('audit/inventory_scope.py','audit/inventory_scope.json')]
rows=[]
for script,receipt in pairs:
    outputs=[]
    for opts in [[],['-O']]:
        r=subprocess.run([sys.executable,*opts,str(P/script)],capture_output=True,check=True)
        outputs.append(r.stdout)
    saved=(P/receipt).read_bytes()
    if outputs[0]!=outputs[1] or outputs[0]!=saved:
        raise ValueError('receipt mismatch '+script)
    rows.append({'script':script,'receipt':receipt,'sha256':hashlib.sha256(saved).hexdigest(),'normal_optimized':'identical'})
for name,meta in json.loads((P/'history/catalog_hashes.json').read_text()).items():
    compressed=(P/'history'/(name+'.gz')).read_bytes();raw=gzip.decompress(compressed)
    if hashlib.sha256(compressed).hexdigest()!=meta['gzip_sha256']:raise ValueError(name+' gzip hash')
    if hashlib.sha256(raw).hexdigest()!=meta['raw_sha256']:raise ValueError(name+' raw hash')
    if len(raw)!=meta['raw_bytes']:raise ValueError(name+' size')
result={'status':'PASS','python':platform.python_version(),'sympy':sympy.__version__,'pairs':rows,'catalog_hashes':'PASS'}
(P/'REPLAY.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps(result,indent=2,sort_keys=True))
