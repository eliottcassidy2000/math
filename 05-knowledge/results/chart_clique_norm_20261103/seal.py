"""Pin inherited sources and seal this packet without publishing or messaging."""
from pathlib import Path
import hashlib,json,subprocess
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
BASE='aaf2c12555'
PATHS=[
 '05-knowledge/results/session143_frobenius_magic_followup_20261103/producer/PROOF.md',
 '05-knowledge/results/session143_frobenius_magic_followup_20261103/independent/PROOF.md',
 '05-knowledge/results/preprint_connections_solfeggio_digits_20261017.md',
 '05-knowledge/results/carmichael_successor_blocks_20261014b.md',
 '05-knowledge/results/session129_primegame_tournament_partition_20261031/partition/README.md',
 '05-knowledge/results/primegame_path_orbits_20261102/CONNECTIONS.md',
 '05-knowledge/results/primegame_port_spectrum_20261103/FOREST_RESPONSE_RECURRENCE.md',
 '05-knowledge/results/uniform_prime_transport_20261103/CLASSIFICATION.md',
 '05-knowledge/results/ramsey_cascade_barrier_and_port_rigidity_20261030/README.md',
 '05-knowledge/results/coloring_clock_session_20261020/FUNCTIONAL_GRAPHS_AND_REGULARITY.md',
 '05-knowledge/results/session144_jacobian_eth_mrb_20261103/jacobian_followup/independent/QUARTIC_BOUNDARY_AUDIT.md',
 '01-canon/theorems/THM-3441-weighted-quartic-jelonek-components-and-inertia-parity.md',
 '01-canon/theorems/THM-2176-gordian-continuation-profile-and-interaction-cocycle.md',
 '01-canon/theorems/THM-2242-tournament-complement-transport-and-knot-kernel-green-rigidity.md']
SOURCES=[(BASE,p) for p in PATHS]+[
 ('2a836157af','05-knowledge/results/session145_endpoint_fence_ramsey_20261103/carrier/POOLED_PARITY_CARRIER.md'),
 ('d40bde5797','05-knowledge/results/jacobian_eth_prime_wait_20261102/WEIGHTED-CLASSIFICATION.md'),
 ('14823946db','05-knowledge/results/session146_carrier_primegame_20261103/carrier7/PROOF.md'),
 ('84154edbc1','05-knowledge/results/session146_carrier_primegame_20261103/carrier8/PROOF.md'),
 ('84154edbc1','05-knowledge/results/session146_carrier_primegame_20261103/machines/README.md'),
 ('4c6f0ca135','05-knowledge/results/session145_endpoint_fence_ramsey_20261103/portfolio/fences/PROOF.md'),
 ('4c6f0ca135','05-knowledge/results/session145_endpoint_fence_ramsey_20261103/portfolio/fences/theta_addendum/PROOF.md'),
 ('4c6f0ca135','05-knowledge/results/session145_endpoint_fence_ramsey_20261103/portfolio/magic_frequency/PROOF.md'),
 ('4c6f0ca135','05-knowledge/results/session145_endpoint_fence_ramsey_20261103/portfolio/ramsey_knots/PROOFS.md'),
 ('4c6f0ca135','05-knowledge/results/session145_endpoint_fence_ramsey_20261103/portfolio/ramsey_knots/SOURCES.md')]
EXTERNAL=[
 {'title':'Maximum Fence Area','url':'https://erich-friedman.github.io/packing/fence/','scope':'problem statement and known constructions, not optimum proofs'},
 {'title':'Maximum Product of Areas','url':'https://erich-friedman.github.io/packing/product/','scope':'known n12 construction and its credit'},
 {'title':'Gullerud and Mbirika, An Euler Phi Function for the Eisenstein Integers and Some Applications','url':'https://math.colgate.edu/~integers/u20/u20.pdf','scope':'Euclidean domain and prime classification, section2'},
 {'title':'Razborov, Flag Algebras','url':'https://www.ias.edu/sites/default/files/math/csdm/05-06/arazborov_flag_algebras.pdf','scope':'classical pattern-count framework; not claimed to supply our exact finite carrier theorem'},
 {'title':'Radziszowski, Small Ramsey Numbers, DS1.18, April24 2026','url':'https://www.cs.rit.edu/~spr/ElJC/ejcram18.pdf','scope':'page57 survey baseline57 and91; later86 bound inherited from pinned Session145 audit'}]
def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT)
def main():
    pins=[]
    for ref,path in SOURCES:
        data=git('show',ref+':'+path)
        pins.append({'commit':git('rev-parse',ref).decode().strip(),'path':path,
                     'blob':git('rev-parse',ref+':'+path).decode().strip(),
                     'sha256':hashlib.sha256(data).hexdigest()})
    (HERE/'SOURCE_PINS.json').write_text(json.dumps({'repo_sources':pins,'primary_sources_opened':EXTERNAL},indent=2)+'\n',encoding='utf-8')
    files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='MANIFEST.json'}
    (HERE/'MANIFEST.json').write_text(json.dumps({'status':'same-agent audit with independent algorithms','sha256':files},indent=2)+'\n',encoding='utf-8')
    print('sealed',len(files),'files;',len(pins),'repo source pins')
if __name__=='__main__':main()
