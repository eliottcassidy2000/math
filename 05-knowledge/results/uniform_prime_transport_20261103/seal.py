"""Pin inherited bytes and seal this packet; does not publish or send messages."""
import hashlib
import json
from pathlib import Path
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SOURCES=[
 ('ce57cb0d03','01-canon/theorems/THM-002-ocf.md'),
 ('ce57cb0d03','01-canon/theorems/THM-155-disjoint-3cycle-identity.md'),
 ('ce57cb0d03','01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md'),
 ('ce57cb0d03','01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md'),
 ('ce57cb0d03','01-canon/theorems/THM-4810-h-spectrum-atlas-strong-fibres-exact-low-spectra-complement-parity-invariant-score-monotonicity-paley-budgets.md'),
 ('ce57cb0d03','05-knowledge/results/prime_inverse_ports_20261102/SYMMETRY.md'),
 ('ce57cb0d03','05-knowledge/results/prime_inverse_ports_20261102/REGULAR_BALANCE.md'),
 ('ce57cb0d03','05-knowledge/results/tournament_composition_classification_20261026/README.md'),
 ('ce57cb0d03','05-knowledge/results/fano_tournament_move_quotients_20261021/session23b-fano/README.md'),
 ('ce57cb0d03','05-knowledge/results/primegame_path_orbits_20261102/THEORY.md'),
 ('ce57cb0d03','05-knowledge/results/metagraph_primegame_20261102/INDUCED_OBSTRUCTIONS.md'),
 ('ce57cb0d03','05-knowledge/results/session142_tournament_odd_spectrum_20261102/minimal_gap/PROOF.md'),
 ('5e1cedf452','05-knowledge/results/primegame_port_spectrum_20261103/FOREST_RESPONSE_RECURRENCE.md'),
 ('b687cd28ca','05-knowledge/results/session144_jacobian_eth_mrb_20261103/jacobian/CORRECTION_LINEAGE.md')]


def git(*args):
    return subprocess.check_output(['git',*args],cwd=ROOT)


def main():
    pins=[]
    for ref,path in SOURCES:
        data=git('show',ref+':'+path)
        pins.append({'commit':git('rev-parse',ref).decode().strip(),'path':path,
                     'blob':git('rev-parse',ref+':'+path).decode().strip(),
                     'sha256':hashlib.sha256(data).hexdigest()})
    (HERE/'SOURCE_PINS.json').write_text(json.dumps({'sources':pins,
        'primary_external':{'author':'Arthur H. Busch','year':2006,'result':'Corollary1',
        'url':'https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf',
        'status':'CITED; primary source opened and checked in this session'}},indent=2)+'\n',encoding='utf-8')
    files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.iterdir())
           if p.is_file() and p.name!='MANIFEST.json'}
    (HERE/'MANIFEST.json').write_text(json.dumps({'status':'same-agent proof audit with independent algorithms',
                                                 'files_sha256':files},indent=2)+'\n',encoding='utf-8')
    print('sealed',len(files),'files and',len(pins),'inherited sources')


if __name__=='__main__':main()
