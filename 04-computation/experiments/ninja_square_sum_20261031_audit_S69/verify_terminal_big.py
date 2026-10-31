"""S69 audit: independent verifier for terminal-form square-sum cycles printed by ss_terminal_big.exe.

For every line 'N=.. T=.. attempts=.. nodes=.. path=v1,...,vN' it checks, with exact integer arithmetic:
  the path is a permutation of 1..N; it starts 1,3; it ends at T = 8 (N even) or 15 (N odd); every adjacent sum is a
  perfect square; the closing sum T + 1 is a perfect square (so the cycle passes 8-1-3 or 15-1-3).
It also requires that the set of N with a verified path is exactly the requested range, and writes a compact JSON
receipt: the SHA-256 of the whole certificate file (14.9 MB for 63..2620, kept out of the repository; it is regenerated
byte-for-byte by the deterministic ss_terminal_big.exe) and summary statistics.

usage: python verify_terminal_big.py PATHFILE LO HI RECEIPT.json
"""
import hashlib
import json
import sys
from math import isqrt


def is_sq(x):
    r = isqrt(x)
    return r * r == x


def main():
    src, lo, hi, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    seen = {}
    max_nodes = (0, None)
    attempts_hist = {}
    with open(src) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            head, _, p = line.partition(' path=')
            if not p:
                raise SystemExit('line without path (EXHAUSTED/UNDECIDED): ' + line[:80])
            kv = dict(tok.split('=', 1) for tok in head.split())
            n, t = int(kv['N']), int(kv['T'])
            P = [int(x) for x in p.split(',')]
            if t != (8 if n % 2 == 0 else 15):
                raise SystemExit(f'N={n}: wrong target {t}')
            if len(P) != n or sorted(P) != list(range(1, n + 1)):
                raise SystemExit(f'N={n}: not a permutation of 1..N')
            if P[0] != 1 or P[1] != 3 or P[-1] != t:
                raise SystemExit(f'N={n}: wrong prefix or end')
            for a, b in zip(P, P[1:]):
                if not is_sq(a + b):
                    raise SystemExit(f'N={n}: non-square edge {a}+{b}')
            if not is_sq(P[-1] + P[0]):
                raise SystemExit(f'N={n}: closing edge not square')
            if n in seen:
                raise SystemExit(f'N={n}: duplicate')
            seen[n] = {'sha256': hashlib.sha256(p.encode()).hexdigest(), 'nodes': int(kv['nodes']),
                       'attempts': int(kv['attempts'])}
            attempts_hist[kv['attempts']] = attempts_hist.get(kv['attempts'], 0) + 1
            if int(kv['nodes']) > max_nodes[0]:
                max_nodes = (int(kv['nodes']), n)
    if sorted(seen) != list(range(lo, hi + 1)):
        missing = sorted(set(range(lo, hi + 1)) - set(seen))
        raise SystemExit(f'coverage failure, missing {missing[:20]}')
    with open(src, 'rb') as fh:
        file_digest = hashlib.sha256(fh.read()).hexdigest()
    receipt = {'status': 'passed', 'range': [lo, hi], 'verified_cycles': len(seen),
               'odd_N_verified': sum(1 for n in seen if n % 2), 'even_N_verified': sum(1 for n in seen if n % 2 == 0),
               'max_nodes_last_attempt': {'nodes': max_nodes[0], 'N': max_nodes[1]},
               'attempts_histogram': dict(sorted(attempts_hist.items(), key=lambda kv: int(kv[0]))),
               'certificate_file_sha256': file_digest,
               'certificate_generator': f'ss_terminal_big.exe {lo} {hi}  (default budget 200000, 40 attempts)'}
    with open(out, 'w') as fh:
        json.dump(receipt, fh, indent=1)
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
