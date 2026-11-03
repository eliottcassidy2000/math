"""Freeze source identities and local proof/evidence bytes; no publishing."""
from pathlib import Path
import subprocess,json,hashlib
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
BASE='6a47c3bac4'
PATHS=[
 '01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md',
 '01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md',
 '01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md',
 '01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md',
 '01-canon/theorems/THM-001-redei.md',
 '01-canon/theorems/THM-1862-order-join-reduction-principle.md',
 '05-knowledge/hypotheses/HYP-3053-tournament-diagonal-layer-flow.md',
 '05-knowledge/results/session146_carrier_primegame_20261103/machines/README.md',
 '05-knowledge/results/session146_carrier_primegame_20261103/prime290/README.md',
 '05-knowledge/results/session146_carrier_primegame_20261103/prime290/WHEEL_GENERAL.md',
 '05-knowledge/results/chart_clique_norm_20261103/CLIQUE_CARRIER.md',
 '05-knowledge/results/period3_graph_square_bridge_20261019.md',
 '05-knowledge/results/quadratic_cycle_phase_20261015b.md']
SOURCES=[(BASE,p) for p in PATHS]+[
 ('917c974679','05-knowledge/results/quadratic_2916_odd_pair_triangles_n4_20261102.md'),
 ('f139791ef3','05-knowledge/results/session146_carrier_primegame_20261103/carrier9/PROOF.md'),
 ('f139791ef3','05-knowledge/results/session146_carrier_primegame_20261103/carrier10_upper/PROOF.md')]
def git(*a):return subprocess.check_output(['git',*a],cwd=ROOT)
def main():
    pins=[]
    for ref,p in SOURCES:
        data=git('show',ref+':'+p)
        pins.append({'commit':git('rev-parse',ref).decode().strip(),'path':p,
                     'blob':git('rev-parse',ref+':'+p).decode().strip(),
                     'sha256':hashlib.sha256(data).hexdigest()})
    data={'repo_sources':pins,'primary_records_opened':[
        {'url':'https://arxiv.org/abs/0808.1408','scope':'CITED endpoint-prime existence; original Dirichlet translation'},
        {'url':'https://arxiv.org/abs/math/9512217','scope':'context only; general classification is conditional in its abstract'}]}
    (HERE/'SOURCE_PINS.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='MANIFEST.json'}
    (HERE/'MANIFEST.json').write_text(json.dumps({'scope':'same-agent audit, scoped proofs and exact controls','sha256':files},indent=2)+'\n',encoding='utf-8')
    print('sealed',len(files),'files and',len(pins),'repo sources')
if __name__=='__main__':main()
