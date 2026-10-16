#!/usr/bin/env python3
"""Concatenate the signed block moments of consecutive resumed dense-census runs into one moments file and write the
census JSON (N, S, H, t_N) used by collatz_pole_certificate_general_N_20261015.py and its relatives.

Each input is an OUT_moments.npz of collatz_census_dense_resumable_20261016.py (keys moments, block0, N, B); the runs
must be consecutive (the first starts at block 0, each next one at the previous N).  The output moments cover
[0, N_last) with B = 2^17.  Checks: contiguity, and the total sum of t (+ t(N)) against the final run's S.
Session opus-2026-10-16-S38.  Usage: python collatz_combine_moments_20261016.py out_moments.npz out_census.json
       final_run.json run1_moments.npz [run2_moments.npz ...]
"""
import sys, json
import numpy as np

OUTM, OUTJ, FINALJ = sys.argv[1:4]
parts = sys.argv[4:]
rows = []
expected_row0 = 0
N = None
for p in parts:
    with np.load(p) as Z:
        m = np.array(Z['moments'])
        b0 = int(Z['block0'])
        B = int(Z['B'])
        N = int(Z['N'])
    row0 = b0 * ((1 << 20) // B)
    if row0 != expected_row0:
        raise SystemExit(f'{p}: starts at row {row0}, expected {expected_row0}')
    rows.append(m)
    expected_row0 = row0 + len(m)
v = np.concatenate(rows, axis=0)
if len(v) != N // B:
    raise SystemExit(f'row count {len(v)} != N/B = {N // B}')
J = json.load(open(FINALJ))
tN = int(J['t_N'])
S = int(v[:, 0].sum()) + tN
if S != int(J['S']):
    raise SystemExit(f'moment total {S} != S(N) = {J["S"]}')
np.savez_compressed(OUTM, moments=v, N=np.int64(N), B=np.int64(B))
json.dump(dict(N=N, S=int(J['S']), H=int(J['H']), t_N=tN), open(OUTJ, 'w'), indent=1)
print(f'wrote {OUTM} ({len(v)} rows, N = 2^{N.bit_length() - 1}) and {OUTJ}; S(N) = {S}')
