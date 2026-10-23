"""Independent orbit truth against the source-only canonicalizer."""
from itertools import product
from pathlib import Path
import importlib.util
import json
import random
import sys

HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


five = load('five_independent_controls', HERE/'check.py')
canonical = load('seven_canonical_under_audit', HERE.parent/'seven/canonical.py')


def main():
    codes = {}
    representatives = {}
    rng = random.Random(523031)
    blocks_checked = 0
    for degree in range(1, 5):
        for B in five.regular_blocks(5, degree):
            truth = min(five.shift(B, a, b) for a, b in product(range(5), repeat=2))
            A = five.tournament({(0, 1): B})
            key, order, metadata = canonical.canonicalize(A, [list(range(5)), list(range(5, 10))])
            if key in codes and codes[key] != truth:
                raise ValueError('canonical code merges nonisomorphic five blocks')
            codes[key] = truth
            representatives.setdefault(truth, A)
            if metadata[0]['canonicalizer_count'] != len(five.stabilizer(B)):
                raise ValueError('canonicalizer fibre differs from independent stabilizer')
            blocks_checked += 1
    if len(codes) != 192 or len(representatives) != 192:
        raise ValueError('canonicalizer split one genuine block orbit')

    port_controls = 0
    for B, A in representatives.items():
        cells = [list(range(5)), list(range(5, 10))]
        path = tuple(five.hamiltonian_path(A))
        permutation = list(range(10))
        rng.shuffle(permutation)
        target, target_cells = canonical.relabel(A, cells, permutation)
        target_path = tuple(permutation[v] for v in path)
        key, order = canonical.canonical_port(A, cells, path)
        other_key, other_order = canonical.canonical_port(target, target_cells, target_path)
        if key != other_key:
            raise ValueError('source/port code not invariant')
        mapping = dict(zip(order, other_order))
        if tuple(mapping[v] for v in path) != target_path:
            raise ValueError('canonical path map changed shared owners')
        K = five.complement_along_path(A, path)
        L = five.complement_along_path(target, target_path)
        if any(K[i][j] != L[mapping[i]][mapping[j]] for i in range(10) for j in range(10)):
            raise ValueError('literal port complement did not transport')
        port_controls += 1

    rigid_matching_reps = sorted(B for B in representatives
                                  if B[0].bit_count() == 1 and len(five.stabilizer(B)) == 1)
    if len(rigid_matching_reps) != 4:
        raise ValueError('wrong local rigid matching alphabet')
    tree_keys = set()
    fibre_sizes = []
    for left, right in product(rigid_matching_reps, repeat=2):
        keys = set()
        for phase in range(5):
            B = {(0, 1): left, (1, 2): five.shift(right, phase, 0), (0, 2): (0,)*5}
            A = five.tournament(B)
            cells = [list(range(5*i, 5*i+5)) for i in range(3)]
            key, order, metadata = canonical.canonicalize(A, cells)
            keys.add(key)
            if metadata[0]['canonicalizer_count'] != 1:
                raise ValueError('two rigid attachments lost actual rigidity')
        if len(keys) != 5:
            raise ValueError('shared middle-cell phase did not give five inequivalent trees')
        fibre_sizes.append(len(keys))
        tree_keys.update(keys)
    if len(tree_keys) != 80:
        raise ValueError('different fixed edge signatures collided')
    result = {'status': 'FINITE-EXACT independent audit',
              'all_regular_five_blocks': blocks_checked,
              'independent_double_orbits_and_canonical_codes': len(codes),
              'actual_port_and_relabelling_controls': port_controls,
              'two_rigid_edge_tree_signatures': len(fibre_sizes),
              'joint_classes_per_tree_signature': fibre_sizes,
              'joint_tree_classes_total': len(tree_keys)}
    destination = HERE/'audit_canonical.json'
    if '--write' in sys.argv:
        destination.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    elif json.loads(destination.read_text()) != result:
        raise ValueError('audit receipt changed')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
