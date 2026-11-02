"""Hash the exact text artifacts and source dependency without line-ending drift."""
import argparse
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
INPUTS=[
    '05-knowledge/results/session142_tournament_odd_spectrum_20261102/jacobian_bridge/PROOFS.md',
    '05-knowledge/results/eleven_wedge_honeycomb_20261102/NATIVE_26.md',
    '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
    '05-knowledge/results/merged_metagraph_realizability_20261020/REALIZABILITY.txt',
    '05-knowledge/results/collatz_carry_interfaces_20261004.md',
    '05-knowledge/results/twenty_curvature_ports_20261102/THREE_WINDOW_CARRIER.md',
    '01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md',
    '01-canon/theorems/THM-589-metagraph-H-variance-closed-form-nud-succession.md',
    '01-canon/theorems/THM-781-hamiltonian-path-inverse-metagraph-fibre.md',
    '01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md',
]

def digest(path):
    return hashlib.sha256(path.read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest()

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    a=HERE/'verification.json';b=HERE/'verification_optimized.json'
    if a.read_bytes()!=b.read_bytes():raise RuntimeError('normal and optimized receipts differ')
    files=[p for p in HERE.iterdir() if p.is_file() and p.name!='manifest.json']
    files+=[ROOT/p for p in INPUTS]
    record={'hash_convention':'SHA256 of UTF-8 text normalized to LF',
            'files':{p.relative_to(ROOT).as_posix():digest(p) for p in sorted(files)}}
    target=HERE/'manifest.json'
    if args.check:
        if json.loads(target.read_text(encoding='utf-8'))!=record:raise RuntimeError('manifest mismatch')
    else:target.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'files':len(files),'receipts_identical':True,'manifest_checked':args.check}))

if __name__=='__main__':main()
