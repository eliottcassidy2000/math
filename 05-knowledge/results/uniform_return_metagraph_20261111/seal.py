"""Create or verify the packet's LF-normalized artifact and dependency manifest."""
from pathlib import Path
import argparse,hashlib,json

BASE=Path(__file__).resolve().parent
REPO=BASE.parents[2]
PINS=[
 '01-canon/theorems/THM-4167-tournament-exposure-capacity-deletion-support-moment-and-parity-holonomy.md',
 '05-knowledge/results/two_entry_funnels_20261103/READY.md',
 '05-knowledge/results/ready_return_pentagon_20261103/ATLAS.md',
 '05-knowledge/results/ready_return_pentagon_20261103/OBSTRUCTION.md',
 '05-knowledge/results/session181_multiport_defect_20261104/synthesis/NOTE.md',
 '05-knowledge/results/session176_tournament_information_20261104/spectrum/PROOF.md',
 '05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py',
 '05-knowledge/results/session252_fivefold_bridges_20261111/algebra/PROOF.md',
 '05-knowledge/results/chanpak_five_core_bridges_20261110/cover/positive_raw_signatures.json',
 '05-knowledge/results/chanpak_five_core_bridges_20261110/algebra/LOW_RATIO_FRONTIER.md',
 '05-knowledge/results/infinite_birth_families_20261111/BIRTH_THEOREMS.md',
 '05-knowledge/results/h13_paired_descent_20261111/dynamics/PROOF.md',
 '05-knowledge/results/h13_paired_descent_20261111/modular_cycles/ADDITIVITY.md',
 '05-knowledge/results/session254_four_frontiers_20261111/chan_pak/PROOF.md',
 '05-knowledge/results/h13_paired_descent_20261111/square_graceful/RESIDUE_PHASE.md',
 '05-knowledge/results/h13_paired_descent_20261111/dynamics/INCOMING_REPAIR.md',
 '05-knowledge/results/session255_modular_collatz_20261111/finite_infinite/PROOF.md',
]


def sha(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def collect():
    artifacts={p.relative_to(BASE).as_posix():sha(p) for p in sorted(BASE.rglob('*'))
               if p.is_file() and p.name!='MANIFEST.json' and '__pycache__' not in p.parts}
    pins={name:sha(REPO/name) for name in PINS}
    return dict(hash_basis='LF-normalized bytes',artifacts=artifacts,inherited_source_pins=pins)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');args=ap.parse_args()
    report=collect();target=BASE/'MANIFEST.json'
    if args.check:
        old=json.loads(target.read_text(encoding='utf-8'))
        if old!=report:raise RuntimeError('Artifact or inherited dependency hash mismatch')
    else:target.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('PASS:',len(report['artifacts']),'packet artifacts and',len(report['inherited_source_pins']),'inherited pins')


if __name__=='__main__':main()
