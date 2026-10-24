"""Literal branching obstruction to a single global Hamiltonian observer."""
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent


def main():
    supply = json.loads((RESULTS/'compatible_charged_helpers_20261023/supply/selected_supply.json').read_text())
    host = json.loads((RESULTS/'guarded_component_completion_20261024/native/xor_tree_receipt.json').read_text())
    by = {r['source']:r for r in supply['helpers']}
    paths = {}
    for s, option in ((368,'before'), (776,'after'), (3139,'after')):
        local = by[s]['selected'][option]['path']
        frame = host['source_frames'][str(s)]
        paths[s] = [frame[v] for v in local]
    if not paths[368][1] == paths[776][4] == paths[3139][7] == 4:
        raise RuntimeError('unexpected external owner boundary')
    required = {(a,b) for p in paths.values() for a,b in zip(p,p[1:])}
    if not {(4,6),(4,11)} <= required or 6 == 11:
        raise RuntimeError('missing literal successor conflict')
    reflected = {(b,a) for a,b in required}
    if not {(6,4),(11,4)} <= reflected:
        raise RuntimeError('converse must retain the predecessor conflict')
    print(json.dumps(dict(status='PASS: fixed host cannot admit a single global observer',
                          fixed_paths=paths, required_successors_of_4=sorted(b for a,b in required if a==4),
                          scope='global observer preserving all local K-partner flags; local controller remains valid')))


if __name__ == '__main__':
    main()
