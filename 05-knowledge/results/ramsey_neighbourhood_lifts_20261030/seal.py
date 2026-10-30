"""LF-normalized artifact hashes and immutable source pins; integrity only."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
RUNTIME='05-knowledge/results/regular_parent_symmetry_20261029/verify.py'
BANK='05-knowledge/results/ramsey_saturated_ports_20261030/verification.json'
PINS={
 '16a4cdb898c7':[
   RUNTIME,BANK,
   '05-knowledge/results/ramsey_circuit_interfaces_20261030/RAMSEY.md',
   '05-knowledge/results/ramsey_saturated_ports_20261030/RECURRENCES.md',
   '05-knowledge/results/odd_descent_local_relaxations_20261030/README.md',
   '05-knowledge/results/session118_shared_gate_deletions_20261030/README.md',
   '05-knowledge/results/session119_odd_transport_classifications_20261030/README.md',
   '01-canon/theorems/THM-455-trans-subtournament-doubling-and-tower.md'],
 '1df1292225':[
   '05-knowledge/results/session120_major_frontier_interfaces_20261030/README.md',
   '05-knowledge/results/session120_major_frontier_interfaces_20261030/ramsey/SPECTRAL_COMPILER.md']}

def digest(data):return hashlib.sha256(data.replace(b'\r\n',b'\n')).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT)
def snapshot():
    artifacts={p.name:digest(p.read_bytes()) for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='manifest.json'}
    pins=[]
    for rev,paths in PINS.items():
        commit=git('rev-parse',rev).decode().strip()
        for path in paths:
            blob=git('show',commit+':'+path)
            if path in (RUNTIME,BANK) and digest(blob)!=digest((ROOT/path).read_bytes()):raise RuntimeError('Runtime dependency changed: '+path)
            pins.append({'commit':commit,'path':path,'lf_sha256':digest(blob)})
    return {'scope':'Integrity only; same-agent proof review; no formalization','artifacts':artifacts,'source_pins':pins,'live_runtime_inputs':[RUNTIME,BANK]}
def main():
    p=argparse.ArgumentParser();p.add_argument('--record',action='store_true');args=p.parse_args()
    data=snapshot();path=HERE/'manifest.json'
    if args.record:path.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    elif data!=json.loads(path.read_text()):raise RuntimeError('Manifest mismatch')
    print(json.dumps({'status':'PASS','artifacts':len(data['artifacts']),'pins':len(data['source_pins'])}))
if __name__=='__main__':main()
