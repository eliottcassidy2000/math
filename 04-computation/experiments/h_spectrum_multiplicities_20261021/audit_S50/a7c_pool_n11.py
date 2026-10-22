"""Audit S50 - a7c: pool the n = 11 (and n = 12) mod-8 masses of all available runs.
Session numbers are copied from h_residue_sampling.out and h_residue_sampling_run1_crashed_at_n12.out (data);
audit numbers are parsed from a7_residue_replication.out and a7b_residue_replication.out."""
import os, re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
runs = {11: [('session run 1', 2_000_000, (1.0010, 0.9998, 1.0012, 0.9980)),
             ('session run 2', 2_000_000, (0.9982, 1.0003, 1.0010, 1.0005))],
        12: [('session run 2', 500_000, (1.0032, 0.9990, 1.0002, 0.9977))]}
for fn in ('a7_residue_replication.out', 'a7b_residue_replication.out'):
    txt = open(os.path.join(HERE, fn)).read()
    for blk in re.finditer(r'n = (\d+) \((\d+) samples, seed (\d+)\).*?mod-8 masses \(1,3,5,7\): 1: ([\d.]+), 3: ([\d.]+), 5: ([\d.]+), 7: ([\d.]+)', txt, re.S):
        n, N = int(blk.group(1)), int(blk.group(2))
        if n in runs:
            runs[n].append((f'audit {fn[:3]} seed {blk.group(3)}', N, tuple(float(blk.group(i)) for i in range(4, 8))))
for n, rs in runs.items():
    tot = sum(N for _, N, _ in rs)
    pooled = np.sum([np.array(m) * N for _, N, m in rs], axis=0) / tot
    se = 4 * np.sqrt(0.25 * 0.75 / tot)
    sed = 4 * np.sqrt(0.5 / tot)
    print(f'n = {n}: runs')
    for name, N, m in rs:
        print(f'    {name:28s} N = {N:>9,d}  masses 1,3,5,7 = {m}   leader: residue {(1, 3, 5, 7)[int(np.argmax(m))]}')
    print(f'    pooled (N = {tot:,d}): ' + ', '.join(f'{r}: {v:.4f}' for r, v in zip((1, 3, 5, 7), pooled)) +
          f'; s.e. of a mass {se:.4f}, of a difference {sed:.4f}; residue 5 - residue 1 = {pooled[2] - pooled[0]:+.4f} '
          f'({(pooled[2] - pooled[0]) / sed:+.1f} s.e.); residue 1 excess {pooled[0] - 1:+.4f} ({(pooled[0] - 1) / se:+.1f} s.e.)')
