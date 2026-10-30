#!/usr/bin/env python3
"""Independent integer, chunked NumPy checker; does not run the relaxation."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def audit(E, path):
    from math import gcd
    q = 5**E
    N, D = q-1, 8*5**(E-1)
    g = gcd(N, D)
    N, D = N//g, D//g
    count = q*(q-q//5)
    need(path.stat().st_size == 4*count, 'complete int32 certificate')
    H = np.memmap(path, dtype='<i4', mode='r')
    maximum = int(H.max())
    need(int(H.min()) >= 0 and int(H[0]) == maximum, 'identity dominates')
    cap = 4*E
    need(N*(cap+1) >= maximum+D*E, 'every omitted positive letter')
    checked, tight, min_slack = 0, 0, 10**15
    # Decode the independently specified compact row order. All calculations
    # are exact int64, with products bounded by 4*q*q in these receipts.
    for start in range(0, count, 262144):
        index = np.arange(start, min(count, start+262144), dtype=np.int64)
        row, b = index//q, index % q
        m = row+1+row//4
        for a in range(1, cap+1):
            inv = pow(2**a, -1, q)
            mm, bb = (3*m*inv) % q, ((3*b+1)*inv) % q
            target = (mm-mm//5-1)*q+bb
            R = np.zeros(len(index), dtype=np.int64)
            for e in range(1, E+1):
                R += ((mm-1) % (5**e) == 0) & (bb % (5**e) == 0)
            slack = H[target].astype(np.int64)-H[index]-D*R+N*a
            need(bool(np.all(slack >= 0)), 'full edge inequality')
            min_slack = min(min_slack, int(slack.min()))
            tight += int(np.count_nonzero(slack == 0))
            checked += len(index)
    # Independent scalar equality word, using unnormalized coefficients.
    t = 2*5**(E-1)
    P, Q, B, reward = 1, 1, 0, 0
    for _ in range(t):
        P, Q, B = 3*P % q, 2*Q % q, (3*B+Q) % q
        reward += sum((P-Q) % (5**e) == 0 and B % (5**e) == 0
                      for e in range(1, E+1))
    need(D*reward == N*t, 'sharp constant-one word')
    sha = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1<<20), b''):
            sha.update(chunk)
    return dict(E=E, q=q, states=count, numerator=N, denominator=D,
                maximum=maximum, minimum=int(H.min()), nonzero=int(np.count_nonzero(H)),
                checked_inequalities=checked, tight_inequalities=tight,
                minimum_slack=min_slack, tail_slack=N*(cap+1)-maximum-D*E,
                equality_word_length=t, equality_reward=reward,
                sha256=sha.hexdigest(), arithmetic='exact int64; no floating point')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('E', type=int)
    parser.add_argument('certificate', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = json.dumps(audit(args.E, args.certificate), indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(report)
    print(report, end='')
