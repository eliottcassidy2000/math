"""Auditor: smallest zeros of Phi_T on larger lists, with independent forest numbers (forest_audit.exe), run in
parallel chunks.  Usage:
  python zeros_big.py file n [--complete] [--limit K] [--tag name]     (codes in the first column of `file`)
  python zeros_big.py random n count seed [--tag name]                 (my own uniform random tournaments)
--complete: replace every almost-regular tournament by its regular completion (order n + 1).
Writes <tag>.rows (code H t0 kind w-root N_2 ...) and prints summary statistics."""
import os, sys, math, random, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tourn

HERE = os.path.dirname(os.path.abspath(__file__))


def run_chunk(args):
    codes, n = args
    return tourn.forest_numbers(codes, n)


def main():
    a = sys.argv[1:]
    tag = None
    if '--tag' in a:
        i = a.index('--tag'); tag = a[i + 1]; del a[i:i + 2]
    limit = None
    if '--limit' in a:
        i = a.index('--limit'); limit = int(a[i + 1]); del a[i:i + 2]
    complete = '--complete' in a
    if complete:
        a.remove('--complete')
    Hfile = None
    if a[0] == 'random':
        n, count, seed = int(a[1]), int(a[2]), int(a[3])
        rng = random.Random(seed)
        E = n * (n - 1) // 2
        codes = [''.join('1' if rng.getrandbits(1) else '0' for _ in range(E)) for _ in range(count)]
        tag = tag or f'random{n}_s{seed}'
    else:
        fn, n = a[0], int(a[1])
        lines = [l.split() for l in open(fn) if l.strip()]
        if limit:
            lines = lines[:limit]
        codes = [l[0] for l in lines]
        Hfile = [int(l[1]) for l in lines] if len(lines[0]) > 1 else None
        if complete:
            codes = [tourn.to_code(tourn.complete_regular(tourn.from_code(c))) for c in codes]
            n += 1
            Hfile = None
        tag = tag or os.path.basename(fn).split('.')[0] + ('_completed' if complete else '')
    nch = 12
    chunks = [(codes[i::nch], n) for i in range(nch)]
    with ThreadPoolExecutor(nch) as ex:
        res_chunks = list(ex.map(run_chunk, chunks))
    res = [None] * len(codes)
    for i in range(nch):
        for j, r in enumerate(res_chunks[i]):
            res[i + j * nch] = r
    rows = []
    with open(os.path.join(HERE, tag + '.rows'), 'w') as f:
        for k, (code, (H, c, _)) in enumerate(zip(codes, res)):
            N = tourn.N_from_c(c, n)
            t0, wr, kind = tourn.smallest_zero(N, n)
            rows.append((code, H, N, t0, kind, wr))
            f.write(f'{code} {H} {t0:.6f} {kind} {wr} ' + ' '.join(map(str, N)) + '\n')
    t = np.array([r[3] for r in rows])
    kinds = [r[4] for r in rows]
    print(f'{tag}: n={n}, {len(rows)} tournaments')
    if Hfile:
        print(f'  H agrees with the file for {sum(r[1] == h for r, h in zip(rows, Hfile))}/{len(rows)}')
    print(f'  t0: min {t.min():.4f}  5% {np.quantile(t, 0.05):.4f}  median {np.median(t):.4f}  max {t.max():.4f}')
    for R in (1.0, 1.1, 1.2):
        print(f'  zero in |z| <= {R}: {int((t <= R).sum())} ({(t <= R).mean():.4%})')
    print(f'  smallest zero kind: ' + ', '.join(f'{k} {kinds.count(k)} ({kinds.count(k) / len(kinds):.1%})'
                                            for k in ('neg', 'pos', 'complex', 'none')))
    i = int(np.argmin(t))
    print(f'  minimiser {rows[i][0]} H={rows[i][1]} N={rows[i][2]} w-root {rows[i][5]}')


if __name__ == '__main__':
    main()
