"""Explicit maintainer utility: pin dependencies or seal this reviewed packet."""
from hashlib import sha256
from pathlib import Path
import argparse,json

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
DEPENDENCIES=[
 '04-computation/a000568_exact_burnside.py',
 '04-computation/a000568_enum.py',
 '04-computation/a000568_gmp_enum_v4.c',
 '01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md',
 '01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md',
 '05-knowledge/results/tournament_future_boundary_20261020/verify.py',
 '05-knowledge/results/tournament_future_boundary_20261020/receipt.json',
 '05-knowledge/results/tournament_future_boundary_20261020/a000568-reference.txt',
 '05-knowledge/results/degree_local_tournament_enumeration_20261021/session22-enumeration/local_scale.py',
 '05-knowledge/results/unprotected_ports_cancellation_20261025/counting/README.md',
 '05-knowledge/results/strong_atom_continuation_vitali_20261026/core.py',
 '05-knowledge/results/unnamed_carrier_dissection_20261026/README.md',
 '05-knowledge/results/session72_even_graph_unifying_carriers_20261026/session72-cross-domain/native-paley/README.md',
 '05-knowledge/results/session79_tournament_enumeration_20261027/ROOT_SWITCH.md',
 '05-knowledge/results/relation_tower_frobenius_ports_20261028/check_root_switch.py',
 '05-knowledge/results/relation_tower_frobenius_ports_20261028/root_switch.json',
 '05-knowledge/results/relation_tower_frobenius_ports_20261028/tournament_bank.json.gz',
 '05-knowledge/results/frobenius_lca_metagraph_20261028/fields.py',
 '05-knowledge/results/frobenius_lca_metagraph_20261028/THEORY.md',
]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--sources',action='store_true');args=parser.parse_args()
    if args.sources:
        data={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in DEPENDENCIES};target=HERE/'sources.json'
    else:
        data={p.name:sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.iterdir())
              if p.is_file() and p.name!='manifest.json'};target=HERE/'manifest.json'
    target.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(target.name,len(data))

if __name__=='__main__':main()
