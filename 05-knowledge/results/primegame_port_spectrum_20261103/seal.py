"""Pin this packet and the explicitly inherited mathematical sources."""
from pathlib import Path
import argparse,hashlib,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
INPUTS=[
 '01-canon/theorems/THM-001-redei.md',
 '05-knowledge/results/zero_heavy_ports_20261020g/profile_growth.md',
 '01-canon/theorems/THM-1900-the-insertion-response-calculus.md',
 '01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md',
 '05-knowledge/results/session142_tournament_odd_spectrum_20261102/frontier/proof.md',
 '05-knowledge/results/session142_tournament_odd_spectrum_20261102/constructions/BINARY_NO_GO.md',
 '05-knowledge/results/two_edge_transport_20261102/ESCAPES_AND_CARRIERS.md',
 '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
 '05-knowledge/results/groups_of_five_primegame_merged_metagraph_20261031.md',
 '05-knowledge/results/session129_primegame_tournament_partition_20261031/partition/README.md',
 '05-knowledge/results/stability_priority_moduli_20261102/FRACTRAN_PRECISION.md',
 '05-knowledge/results/primegame_path_orbits_20261102/THEORY.md',
 '05-knowledge/results/primegame_path_orbits_20261102/CONNECTIONS.md',
 '05-knowledge/results/metagraph_primegame_20261102/TWO_BIRTH_THEOREM.md',
 '05-knowledge/results/metagraph_primegame_20261102/INDUCED_OBSTRUCTIONS.md',
 '05-knowledge/results/session143_magic_square_cyclic_bridges_20261103/cyclic_integral_order/PROOF.md',
]
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    a=(HERE/'verification.json').read_bytes();b=(HERE/'verification_optimized.json').read_bytes()
    if a!=b:raise RuntimeError('Normal and optimized verifier outputs differ')
    paths=[p for p in HERE.iterdir() if p.is_file() and p.name!='manifest.json']+[ROOT/p for p in INPUTS]
    result={'convention':'SHA256 of UTF8 text, leading BOM removed, CRLF normalized to LF','files':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_text(encoding='utf-8-sig').replace('\r\n','\n').encode()).hexdigest() for p in sorted(paths)}}
    target=HERE/'manifest.json'
    if args.check:
        if json.loads(target.read_text())!=result:raise RuntimeError('Source or artifact hash mismatch')
    else:target.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'files':len(paths),'normal_optimized_identical':True,'checked':args.check}))
if __name__=='__main__':main()
