"""Seal exact packet bytes and inherited Git source blobs; verify without mutation."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BASE='16a4cdb898c71d54dcb2c6cf58784c8e13692d22'
BASE_PATHS=[
 '05-knowledge/results/paley_multicoloring_and_carry_20261021/session23b-palettes/PALEY_MULTICOLORING_AND_CARRY.md',
 '01-canon/theorems/THM-1850-directed-wowii-domination-transitivity-inequality-and-three-refutations.md',
 '05-knowledge/results/ramsey_circuit_interfaces_20261030/RAMSEY.md',
 '05-knowledge/results/ramsey_circuit_interfaces_20261030/README.md',
 '05-knowledge/results/terminal_fees_phase_covers_20261030/LRC.md',
 '05-knowledge/results/odd_eighths_uniform_families_20261030/LRC.md',
 '05-knowledge/results/odd_descent_local_relaxations_20261030/PROOFS.md',
 '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
 '01-canon/theorems/THM-658-chi-c-GW-below-14-linearization-fails.md',
 '05-knowledge/results/clock_even_guarded_composition_20261026/THEORY.md',
 '05-knowledge/results/session119_odd_transport_classifications_20261030/RESEARCH_SYNTHESIS.md',
]
SOURCES=dict.fromkeys(BASE_PATHS,BASE)
SOURCES.update(dict.fromkeys([
 '05-knowledge/results/precision_lifts_ramsey_deficits_20261030/LRC.md',
 '05-knowledge/results/session120_major_frontier_interfaces_20261030/lrc/LRC.md',
], 'fd5a8148a4'))
SOURCES.update(dict.fromkeys([
 '05-knowledge/results/ramsey_neighbourhood_lifts_20261030/README.md',
 '05-knowledge/results/ramsey_neighbourhood_lifts_20261030/GRAPH_CAPACITY.md',
], '4eb3c4555d'))
REPAIRED='01-canon/theorems/THM-658-chi-c-GW-below-14-linearization-fails.md'
EXECUTED='05-knowledge/results/regular_parent_symmetry_20261029/verify.py'
RECEIPTS=['normality_receipt','palette_receipt','lrc_receipt','lrc_families',
          'capacity_receipt','capacity_audit','corruption_receipt','all_pair_receipt']

def git(*args):
    return subprocess.check_output(['git',*args],cwd=ROOT)

def sha(data):
    return hashlib.sha256(data).hexdigest()

def canonical(path):
    return path.read_bytes().replace(b'\r\n',b'\n')

def build():
    artifacts={}
    for path in sorted(HERE.rglob('*')):
        if not path.is_file() or path.name=='manifest.json' or '__pycache__' in path.parts:
            continue
        data=path.read_bytes()
        if b'\r' in data: raise RuntimeError('Non-LF artifact: '+path.name)
        artifacts[path.relative_to(HERE).as_posix()]=sha(data)
    sources={path:{'revision':git('rev-parse',rev).decode().strip(),
                   'sha256':sha(git('show',rev+':'+path))} for path,rev in SOURCES.items()}
    return {'sources':sources,'artifacts':artifacts,
            'repaired_current_sources':{REPAIRED:sha(canonical(ROOT/REPAIRED))},
            'status':'PROVED scoped with exact certificate / FINITE-EXACT / OPEN general targets',
            'review':'Same-agent review; separate exact algorithms; no formal build.',
            'runtime':'Python standard library; NumPy/SciPy only for optional discovery.',
            'basis':'Exact Git source blobs; exact LF packet bytes; explicit repaired-source bytes.'}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--write',action='store_true')
    ap.add_argument('--staged',action='store_true'); args=ap.parse_args()
    result=build(); manifest=HERE/'manifest.json'
    if args.write:
        manifest.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    elif json.loads(manifest.read_text())!=result:
        raise RuntimeError('Manifest differs')
    if sha(canonical(ROOT/EXECUTED))!=result['sources'][EXECUTED]['sha256']:
        raise RuntimeError('Executed inherited helper differs from its pin')
    for name in RECEIPTS:
        if (HERE/(name+'.json')).read_bytes()!=(HERE/(name+'.optimized.json')).read_bytes():
            raise RuntimeError('Normal/optimized mismatch: '+name)
    if args.staged:
        for name,digest in result['artifacts'].items():
            rel=(HERE/name).relative_to(ROOT).as_posix()
            if sha(git('show',':'+rel))!=digest: raise RuntimeError('Staged mismatch: '+name)
        if git('show',':'+manifest.relative_to(ROOT).as_posix())!=manifest.read_bytes():
            raise RuntimeError('Staged manifest differs')
        if sha(git('show',':'+REPAIRED))!=result['repaired_current_sources'][REPAIRED]:
            raise RuntimeError('Staged repaired source differs')
    print(json.dumps({'sources':len(SOURCES),'artifacts':len(result['artifacts']),
                      'normal_optimized_pairs':len(RECEIPTS),'staged':args.staged,'write':args.write}))

if __name__=='__main__': main()
