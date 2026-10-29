"""Explicit maintainer operation. The verifier never invokes this file."""
from pathlib import Path
from hashlib import sha256
import json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def digest(path):
    return sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def main():
    runtime=json.loads((HERE/'inherited.json').read_text(encoding='utf-8'))['runtime_sources']
    texts=[
        'moving_selector_atoms_20261028/THEORY.md',
        'session90_even_graphs_odd_descent_20261029/carriers/TRANSPORT_AND_HOLONOMY.md',
        'session92_joint_order_carriers_20261029/algebra/JOINT_CLOSURE_ALGEBRA.md',
        'session92_joint_order_carriers_20261029/observers/FOLDED_OBSERVER_COSTS.md',
        'session91_packing_golden_carriers_20261029/operators/GUARDED_FLIP_LEAKAGE.md',
        'lucas_mirror_gap_spaces_20261028/LUCAS_MIRRORS.md',
        'lucas_mirror_gap_spaces_20261028/GAP_TOPOLOGY.md',
        'golden_prime_clocks_20261003.md',
        'endpoint_pythagorean_fibers_20261019.md',
        'square_sum_minus_bounded_section_contacts_20261019.md',
        'berggren_ap_square_nonregular_20261019.md',
        'idoneal_information_channels_20261021/IDONEAL_INFORMATION_CHANNELS.md',
        'session29_framed_profiles_20261022/session29-dehn/REPORT.md',
        'halfperiod_twin_classes_cm_log_law_20261018.md',
        'two_grids_ring_or_ideal_class_13_37_20261027.md',
        'session27_repunits_lattices_20261022/session27-recurrence/RECIPROCAL_POWER_THEOREM.md']
    names=sorted(set(runtime)|{'05-knowledge/results/'+p for p in texts})
    pins={name:digest(ROOT/name) for name in names}
    (HERE/'sources.json').write_text(json.dumps(pins,indent=2)+'\n',encoding='utf-8')
    excluded={'manifest.json','verification.json','verification.optimized.json'}
    manifest={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.is_file() and p.name not in excluded}
    (HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'pinned_sources':len(pins),'sealed_artifacts':len(manifest)}))

if __name__=='__main__':main()
