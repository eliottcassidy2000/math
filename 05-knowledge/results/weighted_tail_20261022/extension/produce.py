"""Produce at most sixteen literal native edges per inherited pending owner.

Discovery is separate from check.py, which never searches an unknown suffix.
The consumer binds every row to the already authenticated ordinary prefix.
"""
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
spec = importlib.util.spec_from_file_location('extension_native_producer', ROOT/'session18_native_annulus/verify_annulus.py')
native = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = native
spec.loader.exec_module(native)


def produce():
    original = json.loads((ROOT/'future_tail_blocks_20261021b/receipt.json').read_bytes())
    V = set(json.loads((ROOT/'session18_native_annulus/frontier.json').read_bytes())['frontier'])
    cases = []
    for row in original['cases']:
        owners = []
        for pending in row['pending']:
            u = pending[0]
            path = []
            for _ in range(16):
                if u in V:
                    break
                target, clock = native.edge(u)
                path.append([u, target, clock])
                u = target
            owners.append({'owner': pending[0], 'path': path})
        cases.append({'N': row['N'], 'D': row['D'], 'owners': owners})
    return {'kind': 'ordinary native extension with immutable original owners', 'depth': 16, 'cases': cases}


if __name__ == '__main__':
    result = produce()
    (HERE/'paths.json').write_bytes((json.dumps(result, sort_keys=True, separators=(',', ':'))+'\n').encode())
    print(json.dumps({'cases': [[row['N'], len(row['owners']), sum(len(v['path']) for v in row['owners'])] for row in result['cases']]}))
