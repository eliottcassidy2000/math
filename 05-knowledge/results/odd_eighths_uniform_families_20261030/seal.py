"""LF-normalized artifact and immutable inheritance seal; not a proof checker."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
ARTIFACTS=('README.md','LRC.md','COLLATZ.md','PHASE_BRIDGE.md','RECOVERY.md','AUDIT.md',
           'check.py','receipt.json','seal.py')
RUNTIMES=(
    '05-knowledge/results/session114_continuation_symmetry_20261030/lrc/check_resonant_tail.py',
    '05-knowledge/results/regular_parent_symmetry_20261029/verify.py')
PINS={
    'b1c353a48c': RUNTIMES+(
        '05-knowledge/results/session114_continuation_symmetry_20261030/lrc/GENERAL_CORE.md',
        '05-knowledge/results/session114_continuation_symmetry_20261030/lrc/RESONANT_TAIL.md',
        '05-knowledge/results/session114_continuation_symmetry_20261030/collatz/PERIODIC_CONTENT_RATES.md',
        '05-knowledge/results/phase_budget_shared_orbits_20261030/COLLATZ.md',
        '05-knowledge/results/source_address_two_clocks_20261004.md',
        '05-knowledge/results/affine_tree_recurrences_20261028/REFINEMENT.md',
        '05-knowledge/results/chessboard_weave_20261006.md',
        '07-reflections/reduced-fraction-harmonic-k4-fibonacci-cassini-codex-20260816.md',
        '05-knowledge/results/analogy_reduction_atlas_20261025/QR_TOURNAMENT.md',
        '01-canon/theorems/THM-3012-central-binomial-quarter-series-and-the-spherical-euclidean-wall.md',
        '05-knowledge/results/amm12592_procgen_20260923_uniform_frontier.md',
        '05-knowledge/results/prime_dilation_response_20261030/PROOFS.md'),
    '893f84da21':(
        '05-knowledge/results/session115_odd_eighth_classifications_20261030/unitary/WORD_UNITARY_CLASSIFICATION.md',
        '05-knowledge/results/session115_odd_eighth_classifications_20261030/lrc/HARMONIC_COMPLETION.md',
        '05-knowledge/results/session115_odd_eighth_classifications_20261030/quadratic/ODD_EIGHTH_PHASES.md',
        '05-knowledge/results/session115_odd_eighth_classifications_20261030/switching/ODD_EIGHTH_GROUPS.md')}

def digest(blob):return hashlib.sha256(blob.replace(b'\r\n',b'\n')).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT)

def snapshot():
    pins=[]
    for revision,paths in PINS.items():
        commit=git('rev-parse',revision).decode().strip()
        for path in paths:
            blob=git('show',commit+':'+path)
            if path in RUNTIMES and digest(blob)!=digest((ROOT/path).read_bytes()):
                raise RuntimeError('Changed inherited runtime: '+path)
            pins.append({'commit':commit,'path':path,'lf_sha256':digest(blob)})
    return {'status':'Integrity only; same-agent proof review',
            'artifacts':{name:digest((HERE/name).read_bytes()) for name in ARTIFACTS},
            'source_pins':pins,'live_runtime_inputs':list(RUNTIMES)}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true')
    args=parser.parse_args();current=snapshot();path=HERE/'manifest.json'
    if args.record:path.write_bytes((json.dumps(current,indent=2,sort_keys=True)+'\n').encode())
    elif current!=json.loads(path.read_text(encoding='utf-8')):raise RuntimeError('Manifest mismatch')
    print(json.dumps({'match':True,'artifacts':len(ARTIFACTS),
                      'source_pins':len(current['source_pins']),'live_runtime_inputs':len(RUNTIMES)}))

if __name__=='__main__':main()
