"""LF-normalized integrity seal and immutable source pins; not a proof checker."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
ARTIFACTS=('.gitattributes','README.md','RAMSEY.md','CIRCUITS.md','RECOVERY.md',
           'AUDIT.md','check.py','audit.py','receipt.json','explore.py','seal.py')
RUNTIME='05-knowledge/results/regular_parent_symmetry_20261029/verify.py'
PINS={
    'b435f833bc':(
        RUNTIME,
        '05-knowledge/results/tournament_open_problem_routes_20261029/PROOFS.md',
        '05-knowledge/results/tournament_open_problem_routes_20261029/ATLAS.md',
        '01-canon/theorems/THM-455-trans-subtournament-doubling-and-tower.md',
        '01-canon/theorems/THM-3489-rule30-packed-restart-and-pointed-pascal-face.md',
        '05-knowledge/results/collatz_catalan_ramsey_20260927.md',
        '05-knowledge/results/session112_counting_interfaces_20261030/counting/README.md',
        '05-knowledge/results/scaled_count_fft_carriers_20261030/SOURCE_AUDIT.md',
        '05-knowledge/results/scaled_count_fft_carriers_20261030/source_receipt.json',
        '05-knowledge/results/cayley_dickson_integral_20261030/PROOFS.md',
        '05-knowledge/results/session116_dimension_owner_towers_20261030/README.md'),
    '25626f0220':(
        '05-knowledge/results/ramsey_saturated_ports_20261030/PROOFS.md',
        '05-knowledge/results/ramsey_saturated_ports_20261030/RECURRENCES.md',
        '05-knowledge/results/return_kernel_core_transplants_20261030/README.md')}

def digest(blob):return hashlib.sha256(blob.replace(b'\r\n',b'\n')).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT)

def snapshot():
    pins=[]
    for revision,paths in PINS.items():
        commit=git('rev-parse',revision).decode().strip()
        for path in paths:
            blob=git('show',commit+':'+path)
            if path==RUNTIME and digest(blob)!=digest((ROOT/path).read_bytes()):
                raise RuntimeError('Changed inherited runtime: '+path)
            pins.append({'commit':commit,'path':path,'lf_sha256':digest(blob)})
    return {'status':'Integrity only; same-agent proof review',
            'artifacts':{name:digest((HERE/name).read_bytes()) for name in ARTIFACTS},
            'source_pins':pins,'live_runtime_inputs':[RUNTIME],
            'external_rules_revision':'2cd3a17f3f296b7f0ac248a63a9b54a0fabdf851',
            'external_pdf_sha256':'e57f24481d2474a6f2c342d86312572f024710835ca583beb31e2bb1f9ceae40'}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true')
    args=parser.parse_args();current=snapshot();path=HERE/'manifest.json'
    if args.record:path.write_text(json.dumps(current,indent=2,sort_keys=True)+'\n',encoding='utf8',newline='\n')
    elif current!=json.loads(path.read_text(encoding='utf8')):raise RuntimeError('Manifest mismatch')
    print(json.dumps({'match':True,'artifacts':len(ARTIFACTS),'source_pins':len(current['source_pins'])}))

if __name__=='__main__':main()
