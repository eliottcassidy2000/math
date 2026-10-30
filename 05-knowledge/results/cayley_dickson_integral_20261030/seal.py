"""LF-normalized integrity seal and immutable source pins; not a proof checker."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
ARTIFACTS=('README.md','PROOFS.md','RECOVERY.md','AUDIT.md','check.py','receipt.json','explore.py','seal.py')
RUNTIME='05-knowledge/results/regular_parent_symmetry_20261029/verify.py'
PINS={
    'f128982d82':(
        RUNTIME,
        '05-knowledge/results/session112_counting_interfaces_20261030/norms/NORM_SIGNS.md',
        '05-knowledge/results/session112_counting_interfaces_20261030/norms/NORM_DEFECT.md',
        '05-knowledge/results/collatz_mod6_20260922_paley_fano_octonion_design.md',
        '05-knowledge/results/session115_odd_eighth_classifications_20261030/quadratic/ODD_EIGHTH_PHASES.md',
        '05-knowledge/results/prime_dilation_response_20261030/PROOFS.md',
        '05-knowledge/results/cd_tower_verify_s20dn.out',
        '05-knowledge/results/lrc_sixteen_sedenion_proof_search_s390.out',
        '01-canon/theorems/THM-869-milgram-corona-sedenion.md',
        '07-reflections/cd-tower-architecture.md',
        '07-reflections/two-towers-cayley-dickson-and-hyperoperation-meet-at-apex7-but-the-octonion-structure-is-only-numerical.md'),
    'dcf45496b4':(
        '05-knowledge/results/session115_odd_eighth_classifications_20261030/collatz/GENERAL_RESET.md',
        '05-knowledge/results/session115_odd_eighth_classifications_20261030/SYNTHESIS.md'),
    'effbfd26fd':(
        '05-knowledge/results/odd_eighth_owner_families_20261030/README.md',
        '05-knowledge/results/entry_cosets_resonant_echoes_20261030/README.md',
        '01-canon/theorems/THM-724-single-killer-covering-min-rigidity.md')}

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
            'source_pins':pins,'live_runtime_inputs':[RUNTIME]}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true')
    args=parser.parse_args();current=snapshot();path=HERE/'manifest.json'
    if args.record:path.write_text(json.dumps(current,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    elif current!=json.loads(path.read_text(encoding='utf-8')):raise RuntimeError('Manifest mismatch')
    print(json.dumps({'match':True,'artifacts':len(ARTIFACTS),'source_pins':len(current['source_pins'])}))

if __name__=='__main__':main()
