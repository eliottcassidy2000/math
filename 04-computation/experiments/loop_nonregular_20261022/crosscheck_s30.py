"""Run the Session30 decorated-core predictor (weight grading w = e + leaves; 5 weight-one and 26 weight-two cores;
05-knowledge/results/session30_nonregular_loops_20261022/session30-decorated-engine/predict.py, imported unchanged) on
this lane's complete class lists and samples, with exact H from ar_tool.  Session opus-2026-10-22-S54.
Writes one row per tournament: code H y y1 y2 A1 A2 (y = n (1 - kappa/rho)).
Usage: python crosscheck_s30.py n FILE OUT [START COUNT]"""
import sys, os, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE = os.path.normpath(os.path.join(HERE, '..', '..', '..', '05-knowledge', 'results', 'session30_nonregular_loops_20261022',
                                       'session30-decorated-engine'))
sys.path.insert(0, ENGINE)
from predict import predict          # Session30 engine (unchanged)

n = int(sys.argv[1]); path = sys.argv[2]; out = sys.argv[3]
start = int(sys.argv[4]) if len(sys.argv) > 4 else 0
count = int(sys.argv[5]) if len(sys.argv) > 5 else 10 ** 9
iu = np.triu_indices(n, 1)
with open(out, 'w') as fh:
    fh.write('code H y y1 y2 A1 A2\n')
    for i, line in enumerate(open(path)):
        if i < start or i >= start + count or not line.strip():
            continue
        code, H = line.split()[0], int(line.split()[1])
        A = np.zeros((n, n), dtype=int); b = np.array([1 if ch == '1' else 0 for ch in code])
        A[iu[0], iu[1]] = b; A[iu[1], iu[0]] = 1 - b
        p = predict(A - A.T)
        kap = H * 2 ** (n - 1) / math.factorial(n)
        y = n * (1 - kap / p['rho'])
        fh.write(f"{code} {H} {y:.12g} {p['y1']:.12g} {p['y2']:.12g} {p['A1']:.12g} {p['A2']:.12g}\n")
        fh.flush()
