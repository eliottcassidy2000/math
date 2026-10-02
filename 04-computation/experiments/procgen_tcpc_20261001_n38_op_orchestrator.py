#!/usr/bin/env python3
"""procgen_tcpc_20261001_n38_op_orchestrator.py -- the orchestrator's second-engine confirmation of THM-4532 at r = 19
(N = 38): every arc of D_19 lies on an odd number of Hamiltonian paths.

The tcpc lane ran karp4 on D_19. This driver runs the lane's other inclusion-exclusion engine, karp3
(procgen_tcpc_20261001_karp3.c), on the CONVERSE D_19^op, so engine and input both differ. Reversing every
Hamiltonian path gives c_{T^op}(v -> u) = c_T(u -> v), so T is all-odd iff T^op is all-odd. This is the same two-engine
standard THM-4532 used at N = 34 (karp3 on T, karp4 on T^op).

D_r is built from the definition used in the orchestrator's audit (procgen_tcpc_20261001_orchestrator_check.py):
  sheets a_j = 2j, b_j = 2j + 1 (j in Z/r), so translation by 2 is an automorphism (karp3 requires tau = +2 invariance);
  a_j -> a_k iff k - j in {1..m};  b_j -> b_k iff j - k in {1..m};  m = (r - 1)/2;
  a_j -> b_k iff k - j in Y, else b_k -> a_j;  Y = {+-1, .., +-floor(m/2)}, plus 0 when m is odd.
One representative per tau-orbit of arcs (37 orbits at N = 38). The work is split into PARTS processes (karp3's
"part total" arguments); the parities of the parts XOR together.

usage: procgen_tcpc_20261001_n38_op_orchestrator.py [r] [parts]     (default r = 19, parts = 4; about 45 CPU-minutes
       per part at r = 19 on the orchestrator's machine)
"""
import os
import re
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def build_input(r):
    m = (r - 1) // 2
    S = set(range(1, m + 1))
    Y = set()
    for d in range(1, m // 2 + 1):
        Y |= {d % r, (-d) % r}
    if m % 2 == 1:
        Y |= {0}
    N = 2 * r
    A = [[0] * N for _ in range(N)]
    for j in range(r):
        for k in range(r):
            if j != k:
                if (k - j) % r in S:
                    A[2 * j][2 * k] = 1
                if (j - k) % r in S:
                    A[2 * j + 1][2 * k + 1] = 1
            if (k - j) % r in Y:
                A[2 * j][2 * k + 1] = 1
            else:
                A[2 * k + 1][2 * j] = 1
    assert all(A[u][v] + A[v][u] == 1 for u in range(N) for v in range(N) if u != v)
    AT = [[A[j][i] for j in range(N)] for i in range(N)]          # the converse D_r^op
    reps, seen = [], set()
    for u in range(N):
        for v in range(N):
            if u != v and AT[u][v]:
                key = (u % 2, v % 2, (v // 2 - u // 2) % r)
                if key not in seen:
                    seen.add(key)
                    reps.append((u, v))
    assert len(reps) == (N * (N - 1) // 2) // r
    return N, '%d\n' % N + '\n'.join(' '.join(map(str, row)) for row in AT) + '\n%d\n' % len(reps) + \
        '\n'.join('%d %d' % e for e in reps) + '\n', len(reps)


def main():
    r = int(sys.argv[1]) if len(sys.argv) > 1 else 19
    parts = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    N, inp, K = build_input(r)
    t0 = time.time()
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, 'karp3')
        subprocess.run(['cc', '-O3', '-march=native', '-o', exe, os.path.join(HERE, 'procgen_tcpc_20261001_karp3.c')], check=True)
        procs = [subprocess.Popen(['nice', exe, str(p), str(parts)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
                 for p in range(parts)]
        outs = [pr.communicate(inp)[0] for pr in procs]
    combine(outs, r, N, K, time.time() - t0)


def combine(outs, r, N, K, secs=None):
    H = 0
    C = {}
    reps = set()
    for p, t in enumerate(outs):
        print(f'---- part {p} of {len(outs)} ----')
        print(t.strip())
        H ^= int(re.search(r'^H (\d)', t, re.M).group(1))
        for u, v, b in re.findall(r'^C (\d+) (\d+) (\d)', t, re.M):
            C[(int(u), int(v))] = C.get((int(u), int(v)), 0) ^ int(b)
        reps.add(int(re.search(r'REPS (\d+)', t).group(1)))
    even = [k for k, x in C.items() if x == 0]
    print(f'==== D_{r}^op (N = {N}) by karp3 in {len(outs)} parts: H mod 2 = {H}; {len(C)} of {K} arc orbits combined; '
          f'odd: {sum(C.values())}, even: {even}; canonical subsets (total, as reported by every part): {sorted(reps)}')
    if secs is not None:
        print(f'elapsed {secs:.0f} s')
    print('ALL ARC ORBITS ODD' if H == 1 and len(C) == K and not even else 'NOT ALL ODD')


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--combine':
        # combine existing part outputs: --combine r file0 file1 ...
        r = int(sys.argv[2])
        N, _, K = build_input(r)
        combine([open(f).read() for f in sys.argv[3:]], r, N, K)
    else:
        main()
