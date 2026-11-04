"""No inherited receipts are modified. Complete declared local menus only."""
import argparse
import json
import subprocess
import tempfile
from pathlib import Path

HERE=Path(__file__).resolve().parent
OLD=HERE.parent.parent/'zero_ready_strong_switch_20261103'/'literal'

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--depth',type=int,default=2)
    p.add_argument('--force-collar',action='store_true')
    a=p.parse_args()
    seed=json.loads((OLD/'realized_generation2_seed4_depth3_all_prefixes.json').read_text())[1]
    names=seed['suffix_names']
    forced=sorted([names.index(22),names.index(32)]) if a.force_collar else [-1,-1]
    data=[len(names),names.index(seed['target_inport']),len(seed['prefix_names']),a.depth,*forced,seed['target_H'],*seed['suffix_rows']]
    with tempfile.TemporaryDirectory(prefix='neutral-collar-') as tmp:
        exe=Path(tmp)/'scan.exe'
        subprocess.run(['g++','-std=c++20','-O3',str(HERE/'search.cpp'),'-o',str(exe)],check=True)
        receipt=json.loads(subprocess.check_output([str(exe)],input=' '.join(map(str,data)),text=True))
    receipt.update(status='FINITE-EXACT',source_names=seed['target_names'],suffix_names=names,
                   source_rows=seed['target_rows'],suffix_rows=seed['suffix_rows'],source_inport=seed['target_inport'],
                   required_property='strong suffix with at least one zero-ready vertex, identical full H, no prefix shrink',
                   source='zero_ready_strong_switch_20261103/literal/realized_generation2_seed4_depth3_all_prefixes.json index1')
    outfile=HERE/f"depth{a.depth}{'_forced_collar' if a.force_collar else ''}.json"
    outfile.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ('source_rows','suffix_rows','source_names','suffix_names')}))

if __name__=='__main__':main()
