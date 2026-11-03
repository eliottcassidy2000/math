"""Pin exact text artifacts and the mathematical/code sources they use."""
import argparse
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
INPUTS=[
 '05-knowledge/results/boolean_realization_71_17_20261102/scout.py',
 '05-knowledge/results/boolean_realization_71_17_20261102/THEORY.md',
 '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
 '05-knowledge/results/eleven_wedge_honeycomb_20261102/extension26.py',
 '05-knowledge/results/eleven_wedge_honeycomb_20261102/NATIVE_26.md',
 '05-knowledge/results/session142_tournament_odd_spectrum_20261102/minimal_gap/PROOF.md',
 '05-knowledge/results/prime_inverse_ports_20261102/MOVES.md',
 '05-knowledge/results/session55_parker_fields_integer_lifts_20261025/padic-lifts/PROOF.md',
 '05-knowledge/results/session57_square_sum_magic_paths_20261025/biclique-bridge/THEOREM.md',
 '05-knowledge/results/collatz_mod6_20260922_paley_fano_octonion_design.md',
]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    if (HERE/'verification.json').read_bytes()!=(HERE/'verification_optimized.json').read_bytes():
        raise RuntimeError('normal and optimized outputs differ')
    files=[p for p in HERE.iterdir() if p.is_file() and p.name!='manifest.json']+[ROOT/p for p in INPUTS]
    result={'convention':'SHA256 UTF8 normalized LF',
            'files':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in sorted(files)}}
    target=HERE/'manifest.json'
    if args.check:
        if json.loads(target.read_text(encoding='utf-8'))!=result:raise RuntimeError('source or artifact hash changed')
    else:target.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'files':len(files),'normal_optimized_identical':True,'checked':args.check}))

if __name__=='__main__':main()
