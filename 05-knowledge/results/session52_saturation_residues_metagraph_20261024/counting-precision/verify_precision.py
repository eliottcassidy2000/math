#!/usr/bin/env python3
"""Check pinned patch application in a temporary tree and exact eigenfunction controls."""
import argparse, hashlib, json, shutil, subprocess, tempfile
from itertools import combinations
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
checks = 0

def require(ok, why):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(why)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def cyclic_count(n, mask):
    scores = [0] * n
    for bit, (i, j) in enumerate(combinations(range(n), 2)):
        scores[i if (mask >> bit) & 1 else j] += 1
    return comb(n, 3) - sum(s * (s - 1) // 2 for s in scores)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--repo', type=Path, default=HERE.parent / 'math')
    p.add_argument('--receipt', type=Path, default=HERE / 'precision-receipt.json')
    args = p.parse_args()
    pins = json.loads((HERE / 'source-pins.json').read_text())['sources']
    with tempfile.TemporaryDirectory(prefix='counting-precision-') as td:
        temp = Path(td)
        for pin in pins:
            source = args.repo / pin['path']
            require(sha(source) == pin['source_sha256'], ('changed source', str(source)))
            target = temp / pin['path']
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
        # Source tree remains untouched: patch runs only in temp.
        r = subprocess.run(['patch', '--batch', '--fuzz=0', '-p1', '-i', str(HERE/'precision.patch')], cwd=temp, capture_output=True, text=True)
        require(r.returncode == 0, ('patch failed', r.stdout, r.stderr))
        for pin in pins:
            require(sha(temp / pin['path']) == pin['result_sha256'], ('unexpected patched output', pin['path']))
            require(sha(args.repo / pin['path']) == pin['source_sha256'], ('source changed', pin['path']))
    rows = []
    for n in range(3, 7):
        d = comb(n, 2)
        counts = [cyclic_count(n, mask) for mask in range(1 << d)]
        require(4 * sum(counts) == comb(n, 3) * (1 << d), ('uniform mean', n))
        fs = [4*c-comb(n,3) for c in counts]
        require(any(fs), ('zero mode', n))
        for mask, f in enumerate(fs):
            require(sum(fs[mask ^ (1 << bit)] for bit in range(d)) == (d-4)*f, ('eigenfunction', n, mask))
        rows.append({'n':n,'labelled_tournaments':1 << d,'centered_numerator':'4*c(T)-binom(n,3)','adjacency_eigenvalue':d-4})
    # Two-dimensional order-three class operator: raw cyclicity is not an eigenfunction.
    a = ((2, 1), (3, 0))
    c = (0, 1)
    f = (-1, 3)
    require(tuple(sum(a[i][j]*c[j] for j in range(2)) for i in range(2)) != tuple(-x for x in c), 'raw cyclicity hostile')
    require(tuple(sum(a[i][j]*f[j] for j in range(2)) for i in range(2)) == tuple(-x for x in f), 'centered n3 control')
    out={'status':'PASS','checks':checks,'source_pins':pins,'patch_sha256':sha(HERE/'precision.patch'),
         'scope':'Temporary-tree exact patch application and all labelled eigenfunction checks through n=6; no shared source edits.',
         'eigenfunction_controls':rows}
    args.receipt.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'checks':checks,'status':'PASS','receipt':str(args.receipt)},sort_keys=True))

if __name__=='__main__': main()
