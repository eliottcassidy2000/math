"""Compare the actual producer polynomial with independent actual-child DP.

No producer direct-DP function is used. The C++ single-graph mode checks every
cut by direct DP before emitting its weighted-cut response for comparison.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import subprocess
import sys

HERE = Path(__file__).resolve().parent

def need(condition, message):
    if not condition:
        raise ValueError(message)

def decode(code, n):
    graph = [0] * n
    k = 0
    for u in range(n):
        for v in range(u + 1, n):
            winner, loser = (u, v) if code & (1 << k) else (v, u)
            graph[winner] |= 1 << loser
            k += 1
    return tuple(graph)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path, default=HERE / 'census.normal')
    parser.add_argument('--producer-dir', type=Path, default=HERE / 'producer_snapshot')
    args = parser.parse_args()
    sys.path.insert(0, str(args.producer_dir.resolve()))
    spec = importlib.util.spec_from_file_location('producer_polynomial', args.producer_dir / 'cut_polynomial.py')
    producer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(producer)
    rng = random.Random(12470531)
    summaries = []
    cases = 0
    cuts = 0
    for n in range(1, 11):
        bits = n * (n - 1) // 2
        # Exhaustive labeled parents through n=4; 20 seeded random parents per
        # larger order, plus the two opposite transitive encodings.
        codes = range(1 << bits) if n <= 4 else sorted({0, (1 << bits) - 1} | {rng.getrandbits(bits) for _ in range(20)})
        checked = 0
        for code in codes:
            raw = subprocess.check_output([str(args.binary.resolve()), '--graph-code', str(n), str(code)], text=True)
            independent = json.loads(raw)
            actual = producer.cut_counts(decode(code, n))
            need(independent['cuts'] == actual, f'polynomial mismatch at n={n}, code={code}')
            need(independent['direct_cut_checks'] == 1 << n, 'missing direct checks')
            checked += 1
            cases += 1
            cuts += 1 << n
        summaries.append({'n': n, 'parents': checked, 'all_cut_comparisons': checked * (1 << n)})
    print(json.dumps({'method': 'producer_polynomial_vs_independent_actual_child_DP',
                      'parents': cases, 'all_cut_comparisons': cuts, 'orders': summaries,
                      'producer_sha256': {name: hashlib.sha256((args.producer_dir / name).read_bytes()).hexdigest()
                                          for name in ['cut_polynomial.py', 'full_binary.py']}}, sort_keys=True))

if __name__ == '__main__':
    main()
